function result = solve_nist_scip(problem,options)
%SOLVE_NIST_SCIP Solve one bounded NIST problem with SCIP/PySCIPOpt.
arguments
    problem (1,1) struct
    options (1,1) struct = nist_scip_options()
end
problem = validate_nist_problem(problem);
defaults = nist_scip_options();
fields = fieldnames(defaults);
for k = 1:numel(fields)
    if ~isfield(options,fields{k}), options.(fields{k}) = defaults.(fields{k}); end
end
startTimer = tic;
if ~isfield(problem,"initial") || isempty(problem.initial)
    initial = solve_nist_linear_block_qp(problem, ...
        (problem.zLower+problem.zUpper)/2,nist_qp_options());
    if ~initial.hasSolution, error("Could not construct the common initial point."); end
    problem.initial = struct("linear",initial.linear, ...
        "z",(problem.zLower+problem.zUpper)/2, ...
        "objective",initial.objective,"residual",initial.residual);
end
[initialObjective,initialResidual,initialPrediction] = evaluate_nist_objective( ...
    problem,problem.initial.linear,problem.initial.z);
relative = localOption(options,"paper_epsilon_relative",0.01);
absolute = localOption(options,"paper_epsilon_absolute",1e-10);
if 0 >= (1-relative)*initialObjective-absolute
    result = localInitialCertificate(problem,options,initialObjective, ...
        initialResidual,initialPrediction);
    result.matlab_wall_time = toc(startTimer);
    result.wall_time = result.matlab_wall_time;
    return
end
reloadBridge = logical(localOption(options,"reload_python_module",false));
[module,bridgeInfo] = import_nist_scip_bridge(Reload=reloadBridge);
pythonOptions = options;
for name = ["save_artifacts","output_directory","reload_python_module"]
    if isfield(pythonOptions,name), pythonOptions = rmfield(pythonOptions,name); end
end
pythonOptions.model_name = char("nist_"+lower(string(problem.id)));
outputDirectory = string(localOption(options,"output_directory",""));
saveArtifacts = logical(localOption(options,"save_artifacts",false));
if saveArtifacts && strlength(outputDirectory)>0
    if ~isfolder(outputDirectory), mkdir(outputDirectory); end
    pythonOptions.log_file = char(fullfile(outputDirectory,"scip.log"));
    pythonOptions.write_problem = char(fullfile(outputDirectory,"instance.cip"));
end
rawJson = module.solve_nist_json( ...
    char(problem.modelId),double(problem.predictors),double(problem.response), ...
    double(problem.residualScale),double(problem.linearLower), ...
    double(problem.linearUpper),double(problem.zLower),double(problem.zUpper), ...
    double(problem.initial.linear),double(problem.initial.z),jsonencode(pythonOptions));
raw = jsondecode(string(rawJson));
if ~isfield(raw,"bridge_api_version") || string(raw.bridge_api_version) ~= "1.0.0"
    error("Unexpected NIST SCIP bridge result.");
end
result = raw;
result.solver = "SCIP";
result.wall_time = toc(startTimer);
result.matlab_wall_time = result.wall_time;
result.bridge_info = bridgeInfo;
result.problem_id = string(problem.id);
result.instance_id = string(problem.instance_id);
result.n_linear = problem.nLinear;
result.n_nonlinear = problem.nNonlinear;
result.n_residual = problem.nResidual;
result.domain_profile = string(problem.domainProfile);
result.linear = double(result.linear(:));
result.z = double(result.z(:));
result.residual = double(result.residual(:));
result.prediction = double(result.prediction(:));
result.residual_original = result.residual*problem.residualScale;
result.prediction_original = result.prediction*problem.residualScale;
[checkObjective,checkResidual,checkPrediction] = evaluate_nist_objective( ...
    problem,result.linear,result.z);
objectiveDifference = abs(checkObjective-result.objective);
residualDifference = max(abs(checkResidual-result.residual));
predictionDifference = max(abs(checkPrediction-result.prediction));
linearViolation = max([0;problem.linearLower-result.linear; ...
    result.linear-problem.linearUpper]);
zViolation = max([0;problem.zLower-result.z;result.z-problem.zUpper]);
feasibilityTolerance = double(localOption(options,"feasibility_tolerance",1e-6));
geometryTolerance = max(1e-8,10*feasibilityTolerance);
validationPassed = objectiveDifference <= 1e-7*(1+abs(checkObjective)) && ...
    residualDifference <= 1e-7*(1+max(abs(checkResidual))) && ...
    predictionDifference <= 1e-7*(1+max(abs(checkPrediction))) && ...
    linearViolation <= geometryTolerance && zViolation <= geometryTolerance;
result.matlab_validation = struct( ...
    "objective_recomputed",checkObjective, ...
    "objective_abs_difference",objectiveDifference, ...
    "max_abs_residual_difference",residualDifference, ...
    "max_abs_prediction_difference",predictionDifference, ...
    "maximum_linear_bound_violation",linearViolation, ...
    "maximum_z_bound_violation",zViolation, ...
    "feasibility_tolerance",feasibilityTolerance, ...
    "geometry_validation_tolerance",geometryTolerance, ...
    "passed",validationPassed);
if ~validationPassed
    error("NISTSCIP:ValidationFailure", ...
        "SCIP result failed independent MATLAB validation.");
end
result.objective_scaled = result.objective;
result.dual_bound_scaled = result.dual_bound;
result.rss_original = result.objective*problem.objectiveScale;
result.dual_bound_original = result.dual_bound*problem.objectiveScale;
result.certified_rss = problem.certifiedRSS;
result.reference_relative_error = abs(result.rss_original-problem.certifiedRSS) / ...
    max(problem.certifiedRSS,realmin);
result.objective_scale = problem.objectiveScale;
if saveArtifacts && strlength(outputDirectory)>0
    save(fullfile(outputDirectory,"problem.mat"),"problem","-v7.3");
    save(fullfile(outputDirectory,"result.mat"),"result","-v7.3");
end
end

function value = localOption(options,name,defaultValue)
if isfield(options,name), value=options.(name); else, value=defaultValue; end
end

function result = localInitialCertificate(problem,options,objective,residual,prediction)
result = struct();
result.solver = "SCIP"; result.status = "paper_gap_initial";
result.scip_status = "not_started"; result.has_solution = true;
result.solution_source = "external_initial";
result.objective = objective; result.primal_bound = objective; result.dual_bound = 0;
result.paper_certificate_satisfied = true;
result.paper_certificate_margin = -((1-options.paper_epsilon_relative)*objective ...
    - options.paper_epsilon_absolute);
result.linear = double(problem.initial.linear(:));
result.z = double(problem.initial.z(:));
result.residual = residual; result.prediction = prediction;
result.residual_original = residual*problem.residualScale;
result.prediction_original = prediction*problem.residualScale;
result.nodes = 0; result.total_nodes = 0; result.solve_time = 0;
result.python_total_time = 0;
result.problem_id = string(problem.id); result.instance_id = string(problem.instance_id);
result.n_linear = problem.nLinear; result.n_nonlinear = problem.nNonlinear;
result.n_residual = problem.nResidual; result.domain_profile = string(problem.domainProfile);
result.objective_scaled = objective; result.dual_bound_scaled = 0;
result.rss_original = objective*problem.objectiveScale;
result.dual_bound_original = 0; result.certified_rss = problem.certifiedRSS;
result.reference_relative_error = abs(result.rss_original-problem.certifiedRSS)/ ...
    max(problem.certifiedRSS,realmin);
result.objective_scale = problem.objectiveScale;
result.matlab_validation = struct("objective_recomputed",objective, ...
    "objective_abs_difference",0,"max_abs_residual_difference",0, ...
    "max_abs_prediction_difference",0,"maximum_linear_bound_violation",0, ...
    "maximum_z_bound_violation",0,"passed",true);
end
