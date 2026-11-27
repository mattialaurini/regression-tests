function test_portable_inputs()
%TEST_PORTABLE_INPUTS Check frozen data, scaling, boxes and common initial points.
state=prepare_nist_reproduction();
oldPath=path; cleanup=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(fullfile(state.baseRoot,'matlab'));
for k=1:8
    p=state.problems{k};
    assert(p.objectiveScale==max(p.certifiedRSS,1e-8));
    assert(norm(p.response*p.residualScale-p.responseOriginal,Inf)<1e-10*max(1,norm(p.responseOriginal,Inf)));
    [value,residual]=evaluate_nist_objective(p,p.initial.linear,p.initial.z);
    assert(abs(value-p.initial.objective)<1e-10*max(1,value));
    assert(norm(residual-p.initial.residual,Inf)<1e-10*max(1,norm(residual,Inf)));
    assert(all(p.initial.linear>=p.linearLower)&all(p.initial.linear<=p.linearUpper));
    assert(all(p.initial.z>=p.zLower)&all(p.initial.z<=p.zUpper));
    assert(state.baseOptions{k}.EnableLocalSearch&&state.baseOptions{k}.LocalSearchAtRoot);
end
fprintf('Eight portable NIST inputs and shared initial points verified. No solver called.\n');
end
