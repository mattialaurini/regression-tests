function test_pkpd_plot_errors()
%TEST_PKPD_PLOT_ERRORS Plotting must not truncate small positive errors.
root=string(fileparts(mfilename('fullpath')));
oldPath=path; cleanup=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(fullfile(root,'pkpd','src'),'-begin');
original=[0;1e-8;1e-5;0.2;NaN];
[shown,level]=pkpd_error_plot_values(original);
assert(level<1e-8 && shown(1)==level);
assert(isequaln(shown(2:end),original(2:end)));
assert(original(1)==0,'Plot preparation changed the source data.');
second=pkpd_error_plot_values([0;0.5],level);
assert(second(1)==level && second(2)==0.5);
[shown,level]=pkpd_error_plot_values(zeros(3,1));
assert(level>0 && all(shown==level));
fprintf('PK/PD plot checks passed: only exact zeros are relocated.\n');
end
