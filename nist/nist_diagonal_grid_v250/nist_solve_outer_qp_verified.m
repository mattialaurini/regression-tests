function result = nist_solve_outer_qp_verified(G,cMinus,cPlus,options,purpose)
%NIST_SOLVE_OUTER_QP_VERIFIED Recover the SAME paper QP when necessary.
% Ordinary, consistent OPTIMAL solves keep the original numerical policy.
% A numerical primal objective is not assumed to be a feasible upper bound.
if nargin<5, purpose="lower"; end
purpose=string(purpose);
assert(isscalar(purpose) && any(purpose==["lower","screening"]),'Invalid QP purpose.');
[m,r]=size(G); A=G+cMinus; B=G+cPlus;
qp=struct('Q',blkdiag(sparse(r,r),speye(m)),'f',zeros(r+m,1), ...
    'objcon',0,'Aineq',sparse([A,-eye(m);-B,eye(m)]),'bineq',zeros(2*m,1), ...
    'Aeq',sparse([ones(1,r),zeros(1,m)]),'beq',1, ...
    'lb',[zeros(r,1);-inf(m,1)],'ub',inf(r+m,1));
start=tic;
first=solve_nist_qp(qp,options);
result=localAssess(first,A,B,1,1,options,false,purpose);
attempts={localAttempt(first,result,"original",1,1,struct())};
canRetry=isfield(first,'solver') && string(first.solver)=="gurobi";
if purpose=="screening"
    % Screening policy is unchanged: it needs only a feasible primal value.
    variants="scaled_divided_objective";
elseif ~first.hasSolution
    % Preserve the validated MGH10 recovery for failures without any iterate.
    variants=["scaled_divided_objective","dual_simplex_original","scaled_keep_objective"];
else
    variants=["dual_simplex_original","scaled_keep_objective","scaled_divided_objective"];
end
for variant=variants
    if ~canRetry || ~localNeedsRecovery(result,A,B,options,purpose), break, end
    remaining=double(options.TimeLimit)-toc(start);
    if remaining<=0, break, end
    scale=1; objectiveFactor=1; qCoefficient=1;
    if variant~="dual_simplex_original"
        scale=max(1,max(abs([A(:);B(:)])));
        if variant=="scaled_keep_objective"
            qCoefficient=scale*scale;
        else
            objectiveFactor=scale*scale;
        end
    end
    if ~isfinite(qCoefficient) || ~isfinite(objectiveFactor), continue, end
    model=struct('A',sparse([A/scale,-eye(m);-B/scale,eye(m);ones(1,r),zeros(1,m)]), ...
        'obj',zeros(r+m,1),'Q',blkdiag(sparse(r,r),qCoefficient*speye(m)), ...
        'rhs',[zeros(2*m,1);1],'sense',[repmat('<',2*m,1);'='], ...
        'lb',qp.lb,'ub',qp.ub,'modelsense','min');
    params=struct('OutputFlag',double(logical(options.Display)), ...
        'Threads',double(options.Threads), ...
        'FeasibilityTol',double(options.FeasibilityTolerance), ...
        'OptimalityTol',double(options.OptimalityTolerance), ...
        'BarConvTol',min(1e-8,double(options.OptimalityTolerance)), ...
        'DualReductions',0,'NumericFocus',3);
    if variant=="dual_simplex_original", params.Method=1; end
    if isfinite(remaining), params.TimeLimit=max(0.001,remaining); end
    try
        raw=gurobi(model,params);
    catch exception
        raw=struct('status','ERROR','error',getReport(exception,'extended','hyperlinks','off'));
    end
    retry=localWrap(raw,r+m);
    candidate=localAssess(retry,A,B,scale,objectiveFactor,options,true,purpose);
    attempts{end+1}=localAttempt(retry,candidate,variant,scale,objectiveFactor,params); %#ok<AGROW>
    if purpose=="screening"
        choose=~localNeedsRecovery(candidate,A,B,options,purpose) || ...
            (~result.hasSolution && candidate.hasSolution);
    else
        choose=(candidate.boundAvailable && (~result.boundAvailable || ...
            candidate.safeLowerBound>result.safeLowerBound)) || ...
            (~result.hasSolution && candidate.hasSolution && ~result.boundAvailable);
    end
    if choose, result=candidate; end
end
result.qpRetryCount=numel(attempts)-1;
result.qpAttempts=attempts;
result.originalQPStatus=string(first.status);
result.wallTime=toc(start);
end

function attempt=localAttempt(native,result,variant,scale,factor,parameters)
attempt=struct('status',string(native.status),'variant',variant, ...
    'residualScale',scale,'objectiveFactor',factor, ...
    'boundAvailable',result.boundAvailable, ...
    'boundVerified',result.boundValidation.verified, ...
    'parameters',parameters,'raw',native.raw);
end

function native=localWrap(raw,n)
native=struct('solver',"gurobi",'status',string(raw.status), ...
    'hasSolution',isfield(raw,'x') && numel(raw.x)==n, ...
    'x',[],'objective',Inf,'objectiveBound',-Inf,'raw',raw);
if native.hasSolution
    native.x=double(raw.x(:));
    if isfield(raw,'objval'), native.objective=double(raw.objval); end
end
if isfield(raw,'objbound') && isfinite(raw.objbound)
    native.objectiveBound=double(raw.objbound);
elseif native.status=="OPTIMAL" && isfinite(native.objective)
    native.objectiveBound=native.objective;
end
end

function needed=localNeedsRecovery(result,A,B,options,purpose)
if purpose=="screening"
    needed=~isfinite(result.feasiblePrimalUpperBound);
elseif ~result.boundAvailable
    needed=true;
elseif ~result.isRecovery || result.safeLowerBound>0
    needed=false;
else
    % A zero dual bound after recovery may be numerically uninformative.
    % It is enough only if a feasible primal point is also close to zero.
    [upper,~]=nist_outer_qp_primal_upper(A,B,result,options);
    needed=~isfinite(upper) || upper>double(options.BoundSafety)*(1+abs(upper));
end
end

function result=localAssess(native,A,B,scale,factor,options,isRecovery,purpose)
result=native; [m,r]=size(A);
result.nativeResidualScale=scale; result.nativeObjectiveFactor=factor;
result.isRecovery=isRecovery;
result.boundAvailable=false; result.safeLowerBound=0;
result.boundValidation=struct('reason',"unusable_solver_status",'verified',false);
result.objective=double(native.objective)*factor;
result.objectiveBound=double(native.objectiveBound)*factor;
if native.hasSolution && numel(native.x)==r+m && all(isfinite(native.x))
    result.x=double(native.x(:)); result.x(r+1:end)=result.x(r+1:end)*scale;
    result.hasSolution=all(isfinite(result.x));
else
    result.hasSolution=false; result.x=[];
end
if purpose=="screening"
    [result.feasiblePrimalUpperBound,result.weights]= ...
        nist_outer_qp_primal_upper(A,B,result,options);
    return
end
status=string(native.status);
if ~any(status==["OPTIMAL","SUBOPTIMAL"]) || ~result.hasSolution || ...
        ~isfinite(result.objective) || ~isfinite(result.objectiveBound)
    return
end
candidate=max(0,result.objectiveBound-double(options.BoundSafety)*(1+abs(result.objectiveBound)));
rawConsistent=candidate<=max(0,result.objective)+ ...
    double(options.BoundSafety)*(1+abs(result.objective));
checkDual=isRecovery || status~="OPTIMAL" || ~rawConsistent;
if candidate==0
    ok=status=="OPTIMAL";
    check=struct('reason',"existing_nonnegativity",'verified',ok);
elseif checkDual
    dualRaw=native.raw;
    if isfield(dualRaw,'pi'), dualRaw.pi=double(dualRaw.pi)*(factor/scale); end
    [ok,check]=nist_verify_qp_dual(A,B,dualRaw,candidate,1);
    check.qp_residual_scale=scale; check.qp_objective_factor=factor;
else
    % Preserve the benchmark's ordinary OPTIMAL + BoundSafety policy.
    ok=true;
    check=struct('reason',"original_optimal_solver_policy", ...
        'verified',false,'accepted_by_solver_status',true);
end
check.raw_primal_bound_consistent=rawConsistent;
if ok && checkDual && candidate>0
    % Raw objval can be infeasible. Compare only against a reconstructed
    % feasible point; never discard a checked dual bound just for raw objval.
    [upper,w]=nist_outer_qp_primal_upper(A,B,result,options);
    result.feasiblePrimalUpperBound=upper; result.weights=w;
    check.feasible_primal_upper=upper;
    if isfinite(upper) && candidate>upper
        ok=false; check.reason="inconsistent_with_feasible_primal";
        check.verified=false;
    end
end
result.boundValidation=check;
if ok, result.boundAvailable=true; result.safeLowerBound=candidate; end
end
