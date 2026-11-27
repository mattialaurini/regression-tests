function record_reproduction_environment(folder,includeSCIP)
%RECORD_REPRODUCTION_ENVIRONMENT Capture versions outside measured solve calls.
arguments
    folder (1,1) string
    includeSCIP (1,1) logical = false
end
if ~isfolder(folder), mkdir(folder); end
info=struct('matlab',version,'architecture',computer('arch'), ...
    'optimization_toolbox',ver('optim'),'gurobi_interface',which('gurobi'));
model=struct('A',sparse(0,1),'obj',0,'lb',0,'ub',1,'modelsense','min');
params=struct('OutputFlag',1,'Threads',1);
captured=evalc('gurobi(model,params);');
token=regexp(captured,'Gurobi Optimizer version ([^\r\n]+)','tokens','once');
if isempty(token), info.gurobi_version="not available in optimizer log";
else, info.gurobi_version=string(token{1}); end
if ismac
    [~,info.cpu]=system('sysctl -n machdep.cpu.brand_string');
    [~,info.ram_bytes]=system('sysctl -n hw.memsize');
    [~,info.operating_system]=system('sw_vers');
elseif isunix
    [~,info.cpu]=system('lscpu');
    [~,info.memory]=system('head -n 3 /proc/meminfo');
    [~,info.operating_system]=system('uname -a');
else
    info.cpu=getenv('PROCESSOR_IDENTIFIER'); info.operating_system=computer;
end
if includeSCIP
    base=fullfile(fileparts(mfilename('fullpath')),'nist_separable_benchmark_v1_0_6');
    [bridge,~]=import_nist_scip_bridge(PackageRoot=base,Reload=true);
    info.python_scip=jsondecode(string(bridge.versions_json()));
end
file=fullfile(folder,'environment_'+string(datetime('now','Format','yyyyMMdd_HHmmss'))+'.json');
fid=fopen(file,'w'); assert(fid>=0,'Cannot write environment JSON.');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(info,PrettyPrint=true));
end
