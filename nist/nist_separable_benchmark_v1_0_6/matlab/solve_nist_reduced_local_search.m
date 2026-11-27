function result = solve_nist_reduced_local_search(problem,region,startZ,options)
%SOLVE_NIST_REDUCED_LOCAL_SEARCH Budgeted primal heuristic for NIST BnB.
%
% Optimizes only the nonlinear block z with fmincon. For every trial z,
% the bounded linear least-squares block is solved exactly. This routine is
% a primal heuristic only: it returns feasible points and never contributes
% a lower bound or a global certificate.
arguments
    problem (1,1) struct
    region (1,1) struct
    startZ (:,1) double
    options (1,1) struct
end

result=struct("available",false,"hasSolution",false,"objective",Inf, ...
    "linear",[],"z",[],"residual",[],"prediction",[], ...
    "wallTime",0,"qpTime",0,"functionEvaluations",0, ...
    "iterations",0,"exitflag",NaN,"message","");

if exist("fmincon","file")==0
    result.message="fmincon is not available.";
    return
end

problem=validate_nist_problem(problem);
zLower=region.zLower(:); zUpper=region.zUpper(:);
z0=min(max(startZ(:),zLower),zUpper);
callTimer=tic;
qpTime=0; objectiveCalls=0;

fminOptions=optimoptions("fmincon","Algorithm","sqp", ...
    "Display","off", ...
    "MaxIterations",options.MaxIterations, ...
    "MaxFunctionEvaluations",options.MaxFunctionEvaluations, ...
    "OptimalityTolerance",options.OptimalityTolerance, ...
    "StepTolerance",options.StepTolerance, ...
    "ConstraintTolerance",options.ConstraintTolerance, ...
    "OutputFcn",@localStopForTime);

try
    [z,~,exitflag,output]=fmincon(@localReducedObjective,z0,[],[],[],[], ...
        zLower,zUpper,[],fminOptions);
    finalQP=solve_nist_linear_block_qp(problem,z,localQPOptions());
    qpTime=qpTime+localField(finalQP,"wallTime",0);
    if finalQP.hasSolution
        [objective,residual,prediction]=evaluate_nist_objective(problem,finalQP.linear,z);
        result.hasSolution=isfinite(objective);
        result.objective=objective;
        result.linear=finalQP.linear;
        result.z=z(:);
        result.residual=residual;
        result.prediction=prediction;
    end
    result.available=true;
    result.exitflag=exitflag;
    result.iterations=localField(output,"iterations",0);
    result.functionEvaluations=max(objectiveCalls,localField(output,"funcCount",objectiveCalls));
catch exception
    result.available=true;
    result.message=string(exception.message);
end
result.wallTime=toc(callTimer);
result.qpTime=qpTime;

    function value=localReducedObjective(zTrial)
        objectiveCalls=objectiveCalls+1;
        if toc(callTimer)>=options.TimeLimit
            value=realmax("double")/1e100;
            return
        end
        qpStart=tic;
        try
            inner=solve_nist_linear_block_qp(problem,zTrial,localQPOptions());
            qpTime=qpTime+toc(qpStart);
            if inner.hasSolution && isfinite(inner.objective)
                value=inner.objective;
            else
                value=realmax("double")/1e100;
            end
        catch
            qpTime=qpTime+toc(qpStart);
            value=realmax("double")/1e100;
        end
    end

    function stop=localStopForTime(~,~,~)
        stop=toc(callTimer)>=options.TimeLimit;
    end

    function qpOptions=localQPOptions()
        qpOptions=options.QPOptions;
        remaining=max(0.001,options.TimeLimit-toc(callTimer));
        if isfield(qpOptions,"TimeLimit") && isfinite(qpOptions.TimeLimit)
            remaining=min(remaining,double(qpOptions.TimeLimit));
        end
        qpOptions.TimeLimit=remaining;
    end
end

function value=localField(s,name,defaultValue)
if isstruct(s) && isfield(s,name), value=s.(name); else, value=defaultValue; end
end
