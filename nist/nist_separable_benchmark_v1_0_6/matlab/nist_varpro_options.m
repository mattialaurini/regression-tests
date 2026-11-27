function options = nist_varpro_options(varargin)
%NIST_VARPRO_OPTIONS Options for the bounded variable-projection baseline.
p = inputParser;
addParameter(p,"TimeLimit",120);
addParameter(p,"RandomStarts",8);
addParameter(p,"Seed",1);
addParameter(p,"Display","off");
addParameter(p,"UseNistStarts",true);
addParameter(p,"QPOptions",nist_qp_options());
parse(p,varargin{:});
options = struct();
options.TimeLimit = double(p.Results.TimeLimit);
options.RandomStarts = double(p.Results.RandomStarts);
options.Seed = double(p.Results.Seed);
options.Display = string(p.Results.Display);
options.UseNistStarts = logical(p.Results.UseNistStarts);
options.QPOptions = p.Results.QPOptions;
end
