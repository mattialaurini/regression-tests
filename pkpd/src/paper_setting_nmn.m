function T = paper_setting_nmn(settings)
%PAPER_SETTING_NMN (N, M, N^) of each row of a settings list (N = M = N^ without screening).
T = zeros(size(settings, 1), 3);
for s = 1:size(settings, 1)
    o = struct('screening', false, 'Nnoscreen', 3, 'Ngrid', 3, 'q0', 2, 'Nhat', 3);
    a = settings{s, 2};
    for i = 1:2:numel(a), o.(a{i}) = a{i+1}; end
    if o.screening, T(s, :) = [o.Ngrid, 1 + o.q0 * (o.Ngrid - 1), o.Nhat];
    else, T(s, :) = o.Nnoscreen * [1 1 1]; end
end
end
