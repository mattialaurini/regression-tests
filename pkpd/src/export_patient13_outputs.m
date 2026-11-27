function export_patient13_outputs(R, outdir)
%EXPORT_PATIENT13_OUTPUTS Write the two manuscript figures and auditable CSVs.
nsig = numel(R.cfg.sigmas);
colors = [0.80 0.15 0.15; 0.90 0.50 0.05; 0.25 0.65 0.30; ...
          0.10 0.50 0.80; 0.60 0.25 0.70; 0.05 0.15 0.35];
styles = {'-', '--', ':', '-.', '--', '-'};
% The same styles are used in both figures and noise panels.
figure('Color','w','Units','centimeters','Position',[2 2 19 7*(nsig+1)]);
for p = 1:nsig
    sim = R.data{p}; t = (0:numel(sim.u)-1)';
    subplot(nsig+1,1,p); hold on
    plot(t,sim.b,'.','Color',[0.65 0.65 0.65],'MarkerSize',6);
    labels = [{'Noisy BIS'}, R.labels];
    for m = 1:6
        F = R.fits{p,m};
        plot(t,F.yhat,styles{m},'Color',colors(m,:),'LineWidth',1.4);
        if ~F.converged, labels{m+1} = [labels{m+1} ' *']; end
    end
    title(sprintf('Patient 13, $\\sigma=%g$',R.cfg.sigmas(p)), 'Interpreter','latex');
    ylabel('$y$ (BIS)','Interpreter','latex');
    legend(labels,'Location','best','NumColumns',2,'FontSize',8,'Box','off');
    xlim([t(1) t(end)]); grid on; box on
end
subplot(nsig+1,1,nsig+1); stairs(t,R.data{1}.u,'k-','LineWidth',1.3);
xlabel('$t$ [s]','Interpreter','latex'); ylabel('$u$ [mg/s]','Interpreter','latex');
xlim([t(1) t(end)]); grid on; box on
save_figure_pdf(fullfile(outdir,'patient13_six_models_BIS.pdf'));
savefig(gcf,fullfile(outdir,'patient13_six_models_BIS.fig'));
% Noise-free simulated output is available: use it to reveal close curves.
figure('Color','w','Units','centimeters','Position',[3 3 19 8*nsig]);
for p = 1:nsig
    subplot(nsig,1,p); hold on
    sim = R.data{p}; t = (0:numel(sim.u)-1)';
    for m = 1:6
        err_abs = abs(R.fits{p,m}.yhat-sim.b_noiseless);
        err_abs(err_abs == 0) = NaN; % exact zeros are not plotted on a log axis
        semilogy(t,err_abs,styles{m}, ...
            'Color',colors(m,:),'LineWidth',1.4);
    end
    set(gca,'YScale','log');
    title(sprintf('$\\sigma=%g$',R.cfg.sigmas(p)), 'Interpreter','latex');
    xlabel('$t$ [s]','Interpreter','latex');
    ylabel('$|\widehat y-y_{\mathrm{true}}|$ (BIS)','Interpreter','latex');
    legend(R.labels,'Location','best','NumColumns',3,'FontSize',8,'Box','off');
    xlim([t(1) t(end)]); grid on; box on
end
save_figure_pdf(fullfile(outdir,'patient13_six_models_errors.pdf'));
savefig(gcf,fullfile(outdir,'patient13_six_models_errors.fig'));
rows = cell(6*nsig,11); j = 0;
for p = 1:nsig
    sim = R.data{p}; t = (0:numel(sim.u)-1)'; Y = zeros(numel(t),6);
    for m = 1:6
        F = R.fits{p,m}; Y(:,m) = F.yhat; j=j+1;
        rows(j,:) = {R.cfg.sigmas(p), R.labels{m}, F.converged, ...
            F.rmse_raw, F.rmse_used, F.rmse_true, F.err_gamma, F.err_Emax, ...
            F.fobj, F.lb, F.stat.time};
    end
    TT = array2table([t sim.u sim.b min(sim.b,sim.E0) sim.b_noiseless Y], ...
        'VariableNames', {'time_s','u','BIS_raw','BIS_ARX_used','BIS_true', ...
        'ARX1','ARX2','ARX3','ARX4','ML_DP','ML_1C'});
    tag = [strrep(sprintf('patient13_signals_sigma%g',R.cfg.sigmas(p)),'.','p') '.csv'];
    writetable(TT,fullfile(outdir,tag));
end
metrics = cell2table(rows,'VariableNames',{'sigma','model','certified', ...
    'RMSE_raw_BIS','RMSE_used_BIS','RMSE_true_BIS','rel_error_gamma', ...
    'rel_error_Emax','objective_norm','lower_bound','search_time_s'});
writetable(metrics,fullfile(outdir,'patient13_metrics.csv')); disp(metrics);
fprintf('Figures and data saved in %s\n',outdir);
end

