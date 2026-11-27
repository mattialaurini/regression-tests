function opts = parse_opts(defaults, varargin)
%PARSE_OPTS Name/value option parser.
%   OPTS = PARSE_OPTS(DEFAULTS, 'Name1', Value1, ...) overrides the fields of
%   the struct DEFAULTS. In MATLAB, calls written as f(Name=Value) arrive here
%   as name/value pairs, so both calling styles work. Unknown names are errors,
%   so a misspelt option can no longer be silently ignored.
opts = defaults;
if mod(numel(varargin), 2) ~= 0
    error('parse_opts:pairs', 'Options must be given as name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = varargin{k};
    if ~ischar(name)
        name = char(name);
    end
    if ~isfield(defaults, name)
        error('parse_opts:unknown', 'Unknown option "%s".', name);
    end
    opts.(name) = varargin{k+1};
end
end
