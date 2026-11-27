function paper_status(batch, what)
%PAPER_STATUS Header line for one run of a batch: position, elapsed time and an
%   estimate of the remaining time (mean time per completed run).
el = toc(batch.t0);
if batch.done > 0
    eta = el / batch.done * (batch.n - batch.done); etas = hms(eta);
else
    etas = 'n/a';
end
fprintf('\n==== run %d/%d | elapsed %s | remaining ~%s | %s\n', batch.done + 1, batch.n, hms(el), etas, what);
end

function s = hms(t)
s = sprintf('%02d:%02d:%02d', floor(t / 3600), floor(mod(t, 3600) / 60), floor(mod(t, 60)));
end
