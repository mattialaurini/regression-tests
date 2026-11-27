function cfg = paper_cfg(cfg, defaults)
%PAPER_CFG Fill the missing fields of cfg with defaults, and reject unknown fields.
%   A field of cfg that is not among the defaults is almost always a leftover of another
%   script (cfg stays in the workspace after PAPER_TESTS_LIN or PAPER_TESTS_ML): it
%   would silently change the output file, the time limit, etc. In that case an error
%   asks to clear cfg.
if isempty(cfg), cfg = struct(); end
f = fieldnames(cfg);
extra = setdiff(f, [fieldnames(defaults); {'resumed_from'}]);
if ~isempty(extra)
    error('paper_cfg:unknown', ['cfg contains fields that this script does not use: %s.\n', ...
        'It is probably a leftover of another script: run "clear cfg" and define only the ', ...
        'fields you want to change.'], strjoin(extra', ', '));
end
f = fieldnames(defaults);
for i = 1:numel(f)
    if ~isfield(cfg, f{i}), cfg.(f{i}) = defaults.(f{i}); end
end
end
