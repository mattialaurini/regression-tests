function summary = write_nist_tolerance_outputs(results,outputDirectory)
%WRITE_NIST_TOLERANCE_OUTPUTS Scaled-objective CSV and manuscript tables.
arguments
    results table
    outputDirectory (1,1) string
end
localCSV(results,fullfile(outputDirectory,'nist_tolerance_results.csv'));
summary=table();
cases=unique(results.case_id,'stable');
for c=1:numel(cases)
    caseId=cases(c); allCase=results(results.case_id==caseId,:);
    relative=allCase.epsilon_relative(1); timeLimit=allCase.time_limit(1);
    lines=strings(0,1);
    for solver=["BnB","SCIP"]
        t=allCase(allCase.solver==solver,:);
        if isempty(t), continue, end
        ok=t.certified;
        meanSuccessfulTime=NaN;
        if any(ok), meanSuccessfulTime=mean(t.wall_time(ok),'omitnan'); end
        row=table(caseId,relative,timeLimit,solver,height(t),sum(ok),100*mean(ok), ...
            meanSuccessfulTime,mean(t.par2,'omitnan'), ...
            median(t.nodes_visited,'omitnan'), ...
            localSum(t.pruned_tolerance_only),localSum(t.pruned_at_incumbent), ...
            'VariableNames',{'case_id','epsilon_relative','time_limit','solver', ...
            'runs','certified_runs','success_percent','mean_time_successes', ...
            'mean_par2','median_nodes_visited','total_pruned_tolerance_only', ...
            'total_pruned_at_incumbent'});
        summary=[summary;row]; %#ok<AGROW>
        lines(end+1)="\begin{table}[tb]";
        lines(end+1)="\centering";
        lines(end+1)="\small";
        lines(end+1)="\setlength{\tabcolsep}{3pt}";
        lines(end+1)="\begin{tabular}{lrrrrrc}";
        lines(end+1)="\hline";
        lines(end+1)="Problem & $U$ & $L$ & Time [s] & Nodes & $N_{\rm tol}$ & Gap met \\";
        lines(end+1)="\hline";
        for i=1:height(t)
            success="No"; if t.certified(i), success="Yes"; end
            fields=[t.problem_id(i),localNumber(t.objective_scaled(i)), ...
                localNumber(t.dual_bound_scaled(i)),localDecimal(t.wall_time(i)), ...
                localInteger(t.nodes_visited(i)),localInteger(t.pruned_tolerance_only(i)),success];
            lines(end+1)=strjoin(fields," & ")+" \\";
        end
        lines(end+1)="\hline";
        lines(end+1)="\end{tabular}";
        caption=solver+" results for $\epsilon_r="+string(sprintf('%.3g',relative))+ ...
            "$ and a nominal time limit of "+string(sprintf('%.0f',timeLimit))+ ...
            "~s. Bounds $U$ and $L$ are reported directly in RSS$/s$ units. "+ ...
            "The configured gap test is $L\ge(1-\epsilon_r)U-10^{-10}$. "+ ...
            "Times are measured wall times, including solver-call overhead.";
        if solver=="BnB"
            caption=caption+" Nodes counts distinct regions visited at least once. "+ ...
                "$N_{\rm tol}$ counts existing nodes discarded with "+ ...
                "$(1-\epsilon_r)UB-10^{-10}\le LB<UB$, using the incumbent "+ ...
                "at discard time; it includes queue closure at global certification.";
        else
            caption=caption+" Nodes is the processed-node count reported by SCIP. "+ ...
                "The tolerance-only discard count is unavailable and shown as a dash.";
        end
        lines(end+1)="\caption{"+caption+"}";
        lines(end+1)="\label{tab:nist_"+caseId+"_"+lower(solver)+"}";
        lines(end+1)="\end{table}";
        lines(end+1)="";
    end
    tableFolder=fullfile(outputDirectory,'tables');
    if ~isfolder(tableFolder), mkdir(tableFolder); end
    file=fullfile(tableFolder,'nist_'+caseId+'.tex');
    temporary=file+'.tmp.tex';
    writelines(lines,temporary); movefile(temporary,file,'f');
end
localCSV(summary,fullfile(outputDirectory,'nist_tolerance_summary.csv'));
end

function localCSV(t,file)
temporary=string(file)+'.tmp.csv';
writetable(t,temporary); movefile(temporary,file,'f');
end

function total=localSum(values)
if all(isnan(values)), total=NaN; else, total=sum(values,'omitnan'); end
end

function text=localNumber(value)
if ~isfinite(value), text="--"; return, end
text=string(sprintf('%.8g',value));
if contains(lower(text),'e')
    parts=split(lower(text),'e');
    text="$"+parts(1)+"\times10^{"+string(str2double(parts(2)))+"}$";
else
    text="$"+text+"$";
end
end

function text=localDecimal(value)
if isfinite(value), text=string(sprintf('%.3f',value)); else, text="--"; end
end

function text=localInteger(value)
if isfinite(value), text=string(sprintf('%.0f',value)); else, text="--"; end
end
