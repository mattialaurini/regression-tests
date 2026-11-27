function report = download_nist_data(args)
%DOWNLOAD_NIST_DATA Download official NIST StRD ASCII files.
arguments
    args.DataDirectory (1,1) string = fullfile( ...
        fileparts(fileparts(mfilename("fullpath"))),"data","raw")
    args.ProblemIds string = strings(0,1)
    args.Overwrite (1,1) logical = false
end
catalog = nist_problem_catalog();
if ~isempty(args.ProblemIds)
    requested = upper(args.ProblemIds(:));
    keep = arrayfun(@(p) any(upper(p.id)==requested),catalog);
    catalog = catalog(keep);
    if isempty(catalog)
        error("NISTBenchmark:UnknownProblem", ...
            "None of the requested problem identifiers are in the catalog.");
    end
end
if ~isfolder(args.DataDirectory), mkdir(args.DataDirectory); end
rows = cell(numel(catalog),1);
options = weboptions("Timeout",60,"ContentType","text");
for k = 1:numel(catalog)
    problem = catalog(k);
    destination = fullfile(args.DataDirectory,problem.dataFile);
    status = "existing";
    message = "";
    if args.Overwrite || ~isfile(destination)
        status = "downloaded";
        urls = [problem.dataUrl, ...
            replace(problem.dataUrl,"https://www.itl.nist.gov","https://itl.nist.gov"), ...
            replace(problem.dataUrl,"https://","http://")];
        errors = strings(0,1);
        succeeded = false;
        for candidate = unique(urls,'stable')
            try
                websave(destination,candidate,options);
                succeeded = true;
                break
            catch downloadError
                errors(end+1,1) = string(downloadError.message); %#ok<AGROW>
                if isfile(destination), delete(destination); end
            end
        end
        if ~succeeded
            status = "failed";
            message = strjoin(errors," | ");
        end
    end
    bytes = NaN;
    if isfile(destination)
        info = dir(destination); bytes = info.bytes;
    end
    rows{k} = table(problem.id,status,string(destination),bytes,message, ...
        'VariableNames',["problem_id","status","file","bytes","message"]);
end
report = vertcat(rows{:});
disp(report);
failed = report.status=="failed";
if any(failed)
    warning("NISTBenchmark:DownloadFailure", ...
        ["Some NIST files could not be downloaded. The official URLs are " + ...
         "listed in docs/SOURCES.md; copy the files manually into %s."], ...
        args.DataDirectory);
end
end
