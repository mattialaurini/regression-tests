function record_pkpd_environment(output)
%RECORD_PKPD_ENVIRONMENT Record execution details before timing numerical runs.
environment=struct('matlab',version,'architecture',computer('arch'), ...
    'control_toolbox',ver('control'),'optimization_toolbox',ver('optim'), ...
    'gurobi_interface',which('gurobi'),'gurobi_threads',0, ...
    'gurobi_presolve',-1,'qp_solver','hybrid','sequential',true);
% Gurobi's automatic thread policy is retained for PK/PD; NIST uses one thread.
model=struct('A',sparse(0,1),'obj',0,'lb',0,'ub',1,'modelsense','min');
params=struct('OutputFlag',1,'Threads',0);
log=evalc('gurobi(model,params);');
token=regexp(log,'Gurobi Optimizer version ([^\r\n]+)','tokens','once');
if ~isempty(token), environment.gurobi_version=token{1}; end
folder=fullfile(output,'sessions'); if ~isfolder(folder), mkdir(folder); end
file=fullfile(folder,'environment_'+string(datetime('now','Format','yyyyMMdd_HHmmssSSS'))+'.mat');
save(file,'environment','-v7');
end
