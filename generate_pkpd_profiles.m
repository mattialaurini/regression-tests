function generate_pkpd_profiles(outputDirectory)
%GENERATE_PKPD_PROFILES Recreate PatientsBIS.pdf from the nominal 13 profiles.
arguments
    outputDirectory (1,1) string = ""
end
root=string(fileparts(mfilename('fullpath')));
if strlength(outputDirectory)==0, outputDirectory=fullfile(root,'outputs','paper','Figures'); end
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
oldPath=path; oldRng=rng;
cleanup=onCleanup(@()restore(oldPath,oldRng)); %#ok<NASGU>
addpath(fullfile(root,'pkpd','src'),'-begin');
figure('Color','w'); hold on
for id=1:13
    sim=simulate_patient(id,0,'omega',0,'seed',id);
    plot((0:numel(sim.u)-1)',sim.b_noiseless,'LineWidth',1);
end
xlabel('$t$ [s]','Interpreter','latex'); ylabel('BIS'); grid on; box on
legend(arrayfun(@(id)sprintf('Patient %d',id),1:13,'UniformOutput',false),'Location','bestoutside');
save_figure_pdf(char(fullfile(outputDirectory,'PatientsBIS.pdf')));
end
function restore(previousPath,previousRng)
path(previousPath); rng(previousRng);
end
