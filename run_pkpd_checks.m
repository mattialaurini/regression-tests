function run_pkpd_checks()
%RUN_PKPD_CHECKS Check saved fits, then derivative/containment routines.
root=string(fileparts(mfilename('fullpath')));
oldPath=path; oldRng=rng;
cleanup=onCleanup(@()restore(oldPath,oldRng)); %#ok<NASGU>
addpath(fullfile(root,'pkpd','src'),'-begin');
test_pkpd_plot_errors();
stored=load(fullfile(root,'pkpd','reference','patient13_fits.mat'),'R'); R=stored.R;
assert(numel(R.fits)==6&&R.cfg.sigmas==1&&R.cfg.id==13);
sim=R.data{1};
for m=1:6
    fit=R.fits{1,m};
    err=sqrt(mean((fit.yhat-sim.b_noiseless).^2));
    assert(abs(err-fit.rmse_true)<1e-10&&fit.converged);
    if m<=4
        c=max(fit.parameters.Q*fit.parameters.theta,0);
        y=sim.E0-fit.Emax*c.^fit.gamma./(1+c.^fit.gamma);
    else
        y=fit.E0-fit.Emax*ml_psi(fit.parameters.z,sim.u,fit.parameters.model);
    end
    assert(norm(y-fit.yhat,Inf)<1e-8,'Saved fit does not reproduce its output.');
end
% The original numerical regression scripts print their individual diagnostics.
tests_fd
assert(worst<1e-8,'ARX containment check failed.');
skip_solver=true;
tests_ml
assert(nviol==0&&worst<1e-8,'Double-pole derivative/containment check failed.');
tests_ml_e0
assert(worst<1e-8,'Estimated-baseline containment check failed.');
tests_ml_oc
assert(nviol==0&&worst<1e-8,'One-compartment derivative/containment check failed.');
tests_ml_third
assert(nv==0,'Third-order interpolation check failed.');
fprintf('PK/PD saved-output and containment checks passed.\n');
end
function restore(previousPath,previousRng)
path(previousPath); rng(previousRng);
end
