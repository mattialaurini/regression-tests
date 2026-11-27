function problem = load_nist_problem(problemId,args)
%LOAD_NIST_PROBLEM Load one separable NIST StRD benchmark instance.
arguments
    problemId (1,1) string
    args.DataDirectory (1,1) string = fullfile( ...
        fileparts(fileparts(mfilename("fullpath"))),"data","raw")
    args.AutoDownload (1,1) logical = true
    args.DomainProfile (1,1) string {mustBeMember(args.DomainProfile,["study"])} = "study"
    args.QPOptions (1,1) struct = nist_qp_options()
end
catalog = nist_problem_catalog();
index = find(strcmpi([catalog.id],problemId),1);
if isempty(index)
    error("Unknown NIST problem '%s'.",problemId);
end
meta = catalog(index);
filename = fullfile(args.DataDirectory,meta.dataFile);
if ~isfile(filename) && args.AutoDownload
    download_nist_data(DataDirectory=args.DataDirectory,ProblemIds=meta.id);
end
if ~isfile(filename)
    error("NISTBenchmark:MissingData", ...
        "Missing %s. Run download_nist_data first.",filename);
end
[yRaw,X] = parse_nist_dat(filename,meta.nObservations,meta.nPredictors);
if meta.modelId == "nelson"
    response = log(yRaw);
else
    response = yRaw;
end
b = meta.studyBounds;
problem = meta;
problem.responseRaw = yRaw;
problem.responseOriginal = response;
problem.objectiveScale = max(double(meta.certifiedRSS),1e-8);
problem.residualScale = sqrt(problem.objectiveScale);
problem.response = response/problem.residualScale;
problem.predictors = X;
problem.linearLower = b.linearLower;
problem.linearUpper = b.linearUpper;
problem.zLower = b.zLower;
problem.zUpper = b.zUpper;
problem.linearVertices = nist_box_vertices(problem.linearLower,problem.linearUpper);
problem.nLinear = numel(problem.linearLower);
problem.nNonlinear = numel(problem.zLower);
problem.nResidual = numel(problem.response);
problem.nVertices = size(problem.linearVertices,2);
problem.domainProfile = args.DomainProfile;
problem.instance_id = lower(meta.id)+"_"+args.DomainProfile;
problem.certifiedLinear = meta.certifiedParameters(meta.linearIndices);
problem.certifiedZ = meta.certifiedParameters(meta.nonlinearIndices);
problem.certifiedObjective = meta.certifiedRSS;
problem.certifiedObjectiveScaled = meta.certifiedRSS/problem.objectiveScale;
if any(problem.certifiedLinear<problem.linearLower) || ...
        any(problem.certifiedLinear>problem.linearUpper) || ...
        any(problem.certifiedZ<problem.zLower) || any(problem.certifiedZ>problem.zUpper)
    error("The study domain for %s does not contain the NIST certified parameters.",problem.id);
end
problem.start1Linear = meta.start1(meta.linearIndices);
problem.start1Z = meta.start1(meta.nonlinearIndices);
problem.start2Linear = meta.start2(meta.linearIndices);
problem.start2Z = meta.start2(meta.nonlinearIndices);
problem = validate_nist_problem(problem);

z0 = min(max(problem.start2Z,problem.zLower),problem.zUpper);
initial = solve_nist_linear_block_qp(problem,z0,args.QPOptions);
if ~initial.hasSolution
    z0 = (problem.zLower+problem.zUpper)/2;
    initial = solve_nist_linear_block_qp(problem,z0,args.QPOptions);
end
if ~initial.hasSolution
    error("Could not construct a feasible initial point for %s.",problem.id);
end
problem.initial = struct("linear",initial.linear,"z",z0, ...
    "objective",initial.objective,"residual",initial.residual);
end
