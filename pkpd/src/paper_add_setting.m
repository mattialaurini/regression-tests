function [settings, idx] = paper_add_setting(settings, N, M, Nh)
%PAPER_ADD_SETTING Append the setting (N, M, N^) unless it is already in the list;
%   idx is its row.
row = paper_setting(N, M, Nh);
idx = find(strcmp(settings(:, 1), row{1}), 1);
if isempty(idx), settings(end+1, :) = row; idx = size(settings, 1); end
end
