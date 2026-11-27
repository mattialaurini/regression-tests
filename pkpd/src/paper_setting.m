function row = paper_setting(N, M, Nh)
%PAPER_SETTING Settings row {label, options} for IDENTIFY_ARX from (N, M, N^).
%   N = M = N^: without screening (valid bound m1 on N points per coordinate).
%   Otherwise: screening m2(R, N, M), valid bound once M <= N^; M must satisfy
%   M = 1 + q (N - 1) with q integer (grid nesting) and N <= N^ <= M.
if N == M && M == Nh
    row = {sprintf('NS N=%d', N), {'screening', false, 'Nnoscreen', N}};
    return
end
q = (M - 1) / (N - 1);
if N < 2 || q ~= round(q) || q < 1
    error('paper_setting:nesting', 'M = %d is not of the form 1 + q(N-1) for N = %d.', M, N);
end
if Nh < N || Nh > M
    error('paper_setting:Nhat', 'N^ = %d must satisfy N <= N^ <= M.', Nh);
end
row = {sprintf('S N=%d M=%d Nh=%d', N, M, Nh), {'screening', true, 'Ngrid', N, 'q0', q, 'Nhat', Nh}};
end
