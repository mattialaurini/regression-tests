function R = paper_patient13_six_models(varargin)
%PAPER_PATIENT13_SIX_MODELS Fit six models from scratch; default sigma=1.
% Copy beside IDENTIFY_ARX and IDENTIFY_ML, then run:
%   R = paper_patient13_six_models;
% Re-running reuses completed fits in the dedicated checkpoint.
% Options: sigmas, outdir, timeout_arx, timeout_ml, solver, redo_unconverged.
o = parse_opts(struct('sigmas', 1, 'outdir', 'patient13_comparison_sigma1', ...
    'timeout_arx', 1800, 'timeout_ml', 3600, 'solver', 'auto', ...
    'redo_unconverged', false), varargin{:});
validateattributes(o.sigmas, {'numeric'}, {'vector','real','finite','nonnegative','nonempty'});
o.sigmas = unique(o.sigmas(:)');
assert(exist('identify_arx', 'file') == 2 && exist('identify_ml', 'file') == 2, ...
    'Run from pkpd/src, with its functions on the MATLAB path.');
assert(exist('ss', 'file') > 0 && exist('c2d', 'file') > 0, ...
    'Control System Toolbox is required by simulate_patient.');
if strcmp(o.solver, 'auto')
    if exist('gurobi', 'file') >= 2, o.solver = 'hybrid'; else, o.solver = 'matlab'; end
end
assert(any(strcmp(o.solver, {'hybrid','gurobi','matlab'})), ...
    'solver must be auto, hybrid, gurobi or matlab.');
if strcmp(o.solver, 'matlab')
    assert(exist('linprog', 'file') > 0 && exist('quadprog', 'file') > 0, ...
        'Without Gurobi, Optimization Toolbox (linprog and quadprog) is required.');
else
    assert(exist('gurobi', 'file') >= 2, 'Gurobi is not on the MATLAB path.');
end
if ~exist(o.outdir, 'dir'), mkdir(o.outdir); end
file = fullfile(o.outdir, 'patient13_checkpoint.mat');
cfg = struct('protocol_version', 1, 'id', 13, 'sigmas', o.sigmas, ...
    'seed', 13, 'omega', 0.2, 'Reltol', 0.05, 'Abstol', 1e-3, ...
    'estimate_e0', true, 'solver', o.solver, 'arx_grid', [3 5 5]);
if exist(file, 'file') == 2
    previous = load(file, 'R'); R = previous.R;
    assert(isequaln(R.cfg, cfg), ...
        'Checkpoint has another protocol/solver. Choose another outdir.');
else
    R = struct('cfg', cfg, 'labels', {{'ARX-1', 'ARX-2', 'ARX-3', 'ARX-4', ...
        'ML-DP', 'ML-1C'}}, 'data', {cell(1,numel(cfg.sigmas))}, 'fits', {cell(numel(cfg.sigmas),6)});
end
rng0 = rng; rng_guard = onCleanup(@() rng(rng0)); %#ok<NASGU>
s0 = qp_solver(); th0 = gurobi_threads(); pr0 = gurobi_presolve();
solver_guard = onCleanup(@() restore_settings(s0, th0, pr0)); %#ok<NASGU>
qp_solver(o.solver); gurobi_threads(0); gurobi_presolve(-1);
patient_profiles
setting = paper_setting(cfg.arx_grid(1), cfg.arx_grid(2), cfg.arx_grid(3));
for p = 1:numel(cfg.sigmas)
    if isempty(R.data{p})
        R.data{p} = simulate_patient(13, cfg.sigmas(p), 'seed', cfg.seed, 'omega', cfg.omega);
        save(file, 'R');
    end
    sim = R.data{p};
    for m = 1:6
        if ~isempty(R.fits{p,m})
            if ~o.redo_unconverged || R.fits{p,m}.converged, continue; end
        end
        fprintf('\nPatient 13 | sigma=%g | %s (%d/%d)\n', cfg.sigmas(p), R.labels{m}, (p-1)*6+m, 6*numel(cfg.sigmas));
        ge = struct('y', sim.b, 'u', sim.u, 'E0', sim.E0);
        common = {'Reltol', cfg.Reltol, 'Abstol', cfg.Abstol, ...
            'pc_score', 'product', 'progress', 15};
        if m <= 4
            % Same clipping and input-driven prediction as FD_TEST_BNB.
            ge.y = min(sim.b, sim.E0);
            [gh, eh, theta, st] = identify_arx(ge, Patients(13,2:5), ...
                'model', 'linearized', 'order', m, common{:}, ...
                'timeout', o.timeout_arx, setting{2}{:});
            c_raw = st.Q * theta; c = max(c_raw, 0);
            yhat = sim.E0 - eh * c.^gh ./ (1 + c.^gh); e0hat = sim.E0;
            parameters = struct('theta', theta, 'Q', st.Q, 'c_unclipped', c_raw);
        else
            if m == 5, model = 'doublepole'; else, model = 'onecomp'; end
            % Raw (unclipped) data and estimated E0, as in PAPER_TESTS_ML.
            [gh, eh, z, st] = identify_ml(ge, 'model', model, common{:}, ...
                'estimate_e0', cfg.estimate_e0, 'timeout', o.timeout_ml);
            e0hat = st.e0; yhat = e0hat - eh * ml_psi(z, sim.u, model);
            parameters = struct('z', z, 'model', model);
        end
        F = struct('label', R.labels{m}, 'gamma', gh, 'Emax', eh, 'E0', e0hat, ...
            'yhat', yhat, 'parameters', parameters, 'stat', st, ...
            'converged', st.converged, 'fobj', st.fobj, 'lb', st.lb, ...
            'rmse_raw', sqrt(mean((yhat-sim.b).^2)), ...
            'rmse_used', sqrt(mean((yhat-ge.y).^2)), ...
            'rmse_true', sqrt(mean((yhat-sim.b_noiseless).^2)), ...
            'err_gamma', abs(gh-sim.gamma_ob)/sim.gamma_ob, ...
            'err_Emax', abs(eh-sim.Emax_ob)/sim.Emax_ob);
        R.fits{p,m} = F;
        % Save each fit before plotting; interruptions preserve earlier fits.
        save(file, 'R');
        fprintf('RMSE true = %.5g BIS; RMSE raw = %.5g BIS; certified = %d\n', ...
            F.rmse_true, F.rmse_raw, F.converged);
        if ~F.converged
            warning('patient13:timeout', '%s, sigma=%g: prescribed gap not reached; incumbent retained.', ...
                F.label, cfg.sigmas(p));
        end
    end
end
export_patient13_outputs(R, o.outdir);
end

function restore_settings(s,th,pr)
qp_solver(s); gurobi_threads(th); gurobi_presolve(pr);
end
