function row = nist_result_to_row(result,timeLimit,errorMessage)
%NIST_RESULT_TO_ROW Convert one solver result into a flat table row.
arguments
    result (1,1) struct
    timeLimit (1,1) double
    errorMessage (1,1) string = ""
end
solver=string(localField(result,"solver",""));
status=string(localField(result,"status","unknown"));
certified=logical(localField(result,"paper_certificate_satisfied",false));
objective=double(localField(result,"objective",NaN));
dual=double(localField(result,"dual_bound",NaN));
rss=double(localField(result,"rss_original",NaN));
dualRss=double(localField(result,"dual_bound_original",NaN));
reference=double(localField(result,"certified_rss",NaN));
referenceError=double(localField(result,"reference_relative_error",NaN));
wall=double(localField(result,"wall_time",NaN));
if ~isfinite(wall), wall=NaN; end
par2=wall;
if ~certified || ~isfinite(wall), par2=2*timeLimit; end
nodes=double(localField(result,"nodes",NaN));
stats=localField(result,"stats",struct());
m1=double(localField(stats,"m1Evaluations",NaN));
m2=double(localField(stats,"m2Evaluations",NaN));
qpTime=double(localField(stats,"qpTime",NaN));
solutionSource=string(localField(result,"solution_source",""));
problemId=string(localField(result,"problem_id",localField(result,"instance_id","")));
instanceId=string(localField(result,"instance_id",problemId));
row=table(instanceId,problemId,solver,status,certified,objective,dual,rss,dualRss, ...
    reference,referenceError,wall,timeLimit,par2,nodes,m1,m2,qpTime, ...
    double(localField(result,"n_linear",NaN)), ...
    double(localField(result,"n_nonlinear",NaN)), ...
    double(localField(result,"n_residual",NaN)), ...
    string(localField(result,"domain_profile","study")),solutionSource,errorMessage, ...
    'VariableNames',["instance_id","problem_id","solver","status","certified", ...
    "objective_scaled","dual_bound_scaled","rss","dual_bound_rss", ...
    "nist_certified_rss","reference_relative_error","wall_time","time_limit", ...
    "par2","nodes","m1_evaluations","m2_evaluations","qp_time", ...
    "n_linear","n_nonlinear","n_residual","domain_profile", ...
    "solution_source","error_message"]);
end
function value=localField(s,name,defaultValue)
if isstruct(s)&&isfield(s,name), value=s.(name); else, value=defaultValue; end
end
