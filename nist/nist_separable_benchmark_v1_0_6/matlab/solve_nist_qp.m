function result = solve_nist_qp(problem,options)
%SOLVE_NIST_QP Solve a convex continuous QP with Gurobi or quadprog.
arguments
    problem (1,1) struct
    options (1,1) struct = nist_qp_options()
end
solver = lower(string(options.Solver));
if solver == "auto"
    if exist("gurobi","file")~=0
        solver = "gurobi";
    elseif exist("quadprog","file")~=0
        solver = "quadprog";
    else
        error("Neither Gurobi nor quadprog is available.");
    end
end
startTime = tic;
switch solver
    case "gurobi"
        model = struct();
        model.A = sparse([problem.Aineq;problem.Aeq]);
        model.obj = double(problem.f(:));
        model.Q = sparse(problem.Q);
        model.objcon = double(localField(problem,"objcon",0));
        model.rhs = double([problem.bineq;problem.beq]);
        model.sense = [repmat('<',size(problem.Aineq,1),1); ...
            repmat('=',size(problem.Aeq,1),1)];
        model.lb = double(problem.lb(:));
        model.ub = double(problem.ub(:));
        model.modelsense = 'min';
        params = struct();
        params.OutputFlag = double(logical(options.Display));
        params.Threads = double(options.Threads);
        params.FeasibilityTol = double(options.FeasibilityTolerance);
        params.OptimalityTol = double(options.OptimalityTolerance);
        params.BarConvTol = min(1e-8,double(options.OptimalityTolerance));
        if isfinite(options.TimeLimit)
            params.TimeLimit = max(0.001,double(options.TimeLimit));
        end
        raw = gurobi(model,params);
        result.status = string(raw.status);
        result.hasSolution = isfield(raw,"x") && ~isempty(raw.x);
        result.x = [];
        result.objective = Inf;
        result.objectiveBound = -Inf;
        if result.hasSolution
            result.x = double(raw.x(:));
            if isfield(raw,"objval")
                result.objective = double(raw.objval);
            else
                result.objective = result.x.'*problem.Q*result.x + ...
                    problem.f(:).'*result.x + localField(problem,"objcon",0);
            end
        end
        if isfield(raw,"objbound") && isfinite(raw.objbound)
            result.objectiveBound = double(raw.objbound);
        elseif result.hasSolution && result.status=="OPTIMAL"
            result.objectiveBound = result.objective;
        end
        result.raw = raw;
    case "quadprog"
        displayValue = "none";
        if options.Display, displayValue = "iter"; end
        qpOptions = optimoptions("quadprog","Display",displayValue);
        try, qpOptions.ConstraintTolerance = options.FeasibilityTolerance; catch, end
        try, qpOptions.OptimalityTolerance = options.OptimalityTolerance; catch, end
        % quadprog minimizes 0.5*x'*H*x+f'*x.
        [x,fval,exitflag,output] = quadprog(2*problem.Q,problem.f(:), ...
            problem.Aineq,problem.bineq,problem.Aeq,problem.beq, ...
            problem.lb,problem.ub,[],qpOptions);
        fval = fval+localField(problem,"objcon",0);
        result.status = "QUADPROG_"+string(exitflag);
        result.hasSolution = ~isempty(x) && isfinite(fval);
        result.x = double(x(:));
        result.objective = double(fval);
        result.objectiveBound = -Inf;
        if exitflag>0, result.objectiveBound = result.objective; end
        result.raw = output;
    otherwise
        error("Unknown QP solver '%s'.",solver);
end
result.solver = solver;
result.wallTime = toc(startTime);
end

function value = localField(s,name,defaultValue)
if isfield(s,name), value = s.(name); else, value = defaultValue; end
end
