function [response,predictors] = parse_nist_dat(filename,nObservations,nPredictors)
%PARSE_NIST_DAT Parse the data block of an official NIST StRD .dat file.
arguments
    filename (1,1) string
    nObservations (1,1) double {mustBeInteger,mustBePositive}
    nPredictors (1,1) double {mustBeInteger,mustBePositive}
end
if ~isfile(filename)
    error("NISTBenchmark:MissingData","NIST data file not found: %s",filename);
end
lines = splitlines(string(fileread(filename)));
pattern = "[-+]?(?:(?:\d+\.?\d*)|(?:\.\d+))(?:[Ee][-+]?\d+)?";
rows = zeros(0,nPredictors+1);
for line = lines.'
    tokens = regexp(char(line),pattern,'match');
    if numel(tokens) ~= nPredictors+1
        continue
    end
    values = str2double(tokens);
    if all(isfinite(values))
        rows(end+1,:) = values; %#ok<AGROW>
    end
end
if size(rows,1) < nObservations
    error("NISTBenchmark:ParseFailure", ...
        "Found %d numeric data rows in %s; expected at least %d.", ...
        size(rows,1),filename,nObservations);
end
rows = rows(end-nObservations+1:end,:);
response = rows(:,1);
predictors = rows(:,2:end);
end
