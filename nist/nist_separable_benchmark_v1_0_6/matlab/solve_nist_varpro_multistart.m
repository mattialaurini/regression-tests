function result = solve_nist_varpro_multistart(problem,options)
%SOLVE_NIST_VARPRO_MULTISTART Local bounded variable projection baseline.
%
% This is not a global solver. At every nonlinear iterate it solves the
% bounded linear least-squares subproblem exactly, then fmincon optimizes the
% reduced objective over the nonlinear box from several starts.
arguments
    problem (1,1) struct
    options (1,1) struct = nist_varpro_options()
end
problem = validate_nist_problem(problem);
startTimer = tic;
if exist("fmincon","file")==0
    result = localUnavailable(problem,"fmincon is not available.");
    return
end
starts = zeros(problem.nNonlinear,0);
if options.UseNistStarts
    starts(:,end+1) = min(max(problem.start1Z,problem.zLower),problem.zUpper);
    starts(:,end+1) = min(max(problem.start2Z,problem.zLower),problem.zUpper);
end
starts(:,end+1) = (problem.zLower+problem.zUpper)/2;
rng(options.Seed,"twister");
for k=1:options.RandomStarts
    starts(:,end+1) = problem.zLower + rand(problem.nNonlinear,1).*(problem.zUpper-problem.zLower); %#ok<AGROW>
end
% Remove duplicate starts after clipping.
starts = unique(round(starts.',14),'rows','stable').';

bestObjective = Inf; bestLinear=[]; bestZ=[]; bestResidual=[]; bestPrediction=[];
runs = repmat(struct("start",[],"z",[],"objective",Inf,"exitflag",NaN, ...
    "iterations",NaN,"funcCount",NaN,"elapsed",NaN,"error",""),size(starts,2),1);
numberCompleted = 0;
for k=1:size(starts,2)
    remaining = options.TimeLimit-toc(startTimer);
    if remaining<=0, break, end
    runStart = tic;
    fminOptions = optimoptions("fmincon","Algorithm","sqp", ...
        "Display",char(options.Display),"OutputFcn",@(x,o,state)localTimeStop(startTimer,options.TimeLimit,x,o,state), ...
        "MaxIterations",1000,"MaxFunctionEvaluations",5000, ...
        "OptimalityTolerance",1e-9,"StepTolerance",1e-10, ...
        "ConstraintTolerance",1e-9);
    objective = @(z)localReducedObjective(problem,z,options.QPOptions);
    try
        [z,fval,exitflag,output] = fmincon(objective,starts(:,k),[],[],[],[], ...
            problem.zLower,problem.zUpper,[],fminOptions);
        linearResult = solve_nist_linear_block_qp(problem,z,options.QPOptions);
        if linearResult.hasSolution
            [fvalCheck,residual,prediction] = evaluate_nist_objective(problem,linearResult.linear,z);
            if fvalCheck<bestObjective
                bestObjective=fvalCheck; bestLinear=linearResult.linear; bestZ=z(:);
                bestResidual=residual; bestPrediction=prediction;
            end
        end
        runs(k).start=starts(:,k); runs(k).z=z(:); runs(k).objective=fval;
        runs(k).exitflag=exitflag; runs(k).iterations=localField(output,"iterations",NaN);
        runs(k).funcCount=localField(output,"funcCount",NaN); runs(k).elapsed=toc(runStart);
    catch exception
        runs(k).start=starts(:,k); runs(k).elapsed=toc(runStart);
        runs(k).error=string(exception.message);
    end
    numberCompleted=k;
end
runs=runs(1:numberCompleted);
elapsed=toc(startTimer);
if isempty(bestLinear)
    result=localUnavailable(problem,"No local run returned a feasible point.");
    result.wall_time=elapsed; result.runs=runs; return
end
result=struct(); result.solver="VarPro"; result.status="local_solution";
result.has_solution=true; result.solution_source="bounded_variable_projection_multistart";
result.objective=bestObjective; result.primal_bound=bestObjective;
result.dual_bound=NaN; result.paper_certificate_satisfied=false;
result.linear=bestLinear; result.z=bestZ; result.residual=bestResidual;
result.prediction=bestPrediction; result.residual_original=bestResidual*problem.residualScale;
result.prediction_original=bestPrediction*problem.residualScale;
result.wall_time=elapsed; result.nodes=NaN;
result.total_nodes=NaN; result.runs=runs; result.number_of_starts=size(starts,2);
result.number_completed=numberCompleted; result.problem_id=string(problem.id);
result.instance_id=string(problem.instance_id); result.n_linear=problem.nLinear;
result.n_nonlinear=problem.nNonlinear; result.n_residual=problem.nResidual;
result.domain_profile=string(problem.domainProfile); result.objective_scaled=bestObjective;
result.rss_original=bestObjective*problem.objectiveScale; result.dual_bound_original=NaN;
result.certified_rss=problem.certifiedRSS;
result.reference_relative_error=abs(result.rss_original-problem.certifiedRSS)/ ...
    max(problem.certifiedRSS,realmin); result.objective_scale=problem.objectiveScale;
end

function value=localReducedObjective(problem,z,qpOptions)
try
    inner=solve_nist_linear_block_qp(problem,z,qpOptions);
    if inner.hasSolution && isfinite(inner.objective), value=inner.objective; else, value=realmax('double')/1e100; end
catch
    value=realmax('double')/1e100;
end
end

function stop=localTimeStop(startTimer,timeLimit,~,~,~)
stop=toc(startTimer)>=timeLimit;
end

function result=localUnavailable(problem,message)
result=struct("solver","VarPro","status","unavailable","has_solution",false, ...
    "solution_source","none","objective",NaN,"primal_bound",NaN, ...
    "dual_bound",NaN,"paper_certificate_satisfied",false,"linear",[],"z",[], ...
    "residual",[],"prediction",[],"wall_time",0,"nodes",NaN,"total_nodes",NaN, ...
    "problem_id",string(problem.id),"instance_id",string(problem.instance_id), ...
    "n_linear",problem.nLinear,"n_nonlinear",problem.nNonlinear, ...
    "n_residual",problem.nResidual,"domain_profile",string(problem.domainProfile), ...
    "rss_original",NaN,"dual_bound_original",NaN,"certified_rss",problem.certifiedRSS, ...
    "reference_relative_error",NaN,"message",string(message));
end
function value=localField(s,name,defaultValue)
if isstruct(s)&&isfield(s,name), value=s.(name); else, value=defaultValue; end
end
