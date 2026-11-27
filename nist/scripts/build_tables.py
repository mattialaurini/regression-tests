#!/usr/bin/env python3
"""Build normalized-RSS LaTeX tables from reference or fresh CSVs (stdlib only)."""
import argparse,csv,math
from pathlib import Path
IDS=['BoxBOD','MGH17','Lanczos2','Kirby2','Hahn1','Gauss1','Nelson','MGH10']
META=dict(zip(IDS,[(6,1,1),(33,3,2),(24,3,3),(151,3,2),(236,4,3),(250,3,5),(128,2,1),(16,1,2)]))
def rows(path):
 with open(path,newline='') as f:return list(csv.DictReader(f))
def num(r,k):return float(r[k])
def success(r):return str(r['certified']).lower() in ('1','true','1.0')
def nfmt(x):return f'{x:.5g}' if math.isfinite(x) else '--'
def timecell(r,budget):
 if r['status']=='timelimit' and not success(r):return f'{budget:g}'+r'$^{\dagger}$'
 t=nfmt(num(r,'wall_time'))
 return t if success(r) else t+' ('+r['status'].replace('_',r'\_')+')'
def table(headers,body,caption):
 return '\n'.join([r'\begin{table}',r'\centering',r'\caption{'+caption+'}',
  r'\begin{tabular}{'+'l'*len(headers)+'}',r'\hline',' & '.join(headers)+r' \\',r'\hline']+
  [' & '.join(row)+r' \\' if row else r'\hline' for row in body]+[r'\hline',r'\end{tabular}',r'\end{table}',''])
def build(comp,ab,out):
 out.mkdir(parents=True,exist_ok=True); summary=[]
 for er,T in [(.1,100),(.01,600),(.001,3600)]:
  body=[];selected=[r for r in comp if num(r,'epsilon_relative')==er]
  for solver in ('BnB','SCIP'):
   q=[r for r in selected if r['solver']==solver]
   if q:summary.append([str(er),solver,f'{sum(success(r) for r in q)}/{len(q)}',
      f'{sum(num(r,"wall_time") if success(r) else 2*T for r in q)/len(q):.2f}'])
  for pid in IDS:
   for solver in ('BnB','SCIP'):
    rr=[r for r in selected if r['problem_id']==pid and r['solver']==solver]
    if not rr:continue
    assert len(rr)==1; r=rr[0]
    t=timecell(r,T)
    body.append([pid,solver,nfmt(num(r,'objective_scaled')),nfmt(num(r,'dual_bound_scaled')),t,
        nfmt(num(r,'nodes_visited')),nfmt(num(r,'pruned_tolerance_only')) if solver=='BnB' else '--'])
  (out/f'comparison_{er:g}.tex').write_text(table(['Problem','Solver','$U$','$L$','$t$ [s]','Nodes','$N_{\\mathrm{tol}}$'],body,
    f'Normalized RSS bounds, relative tolerance {er:g}, budget {T} s. '+r'$\dagger$: configured limit; raw wall times remain in the CSV.'))
 (out/'comparison_summary.tex').write_text(table(['$\\epsilon_r$','Solver','Certified','Mean PAR2 [s]'],summary,'Certification and penalized computation time.'))
 for pid in IDS:
  body=[];nr,nl,nn=META[pid]
  for er,T in [(.1,100),(.01,1000),(.001,3600)]:
   for M in (3,5,9,17):
    cells=[]
    for v in ('with_m2','without_m2'):
     rr=[r for r in ab if r['problem_id']==pid and num(r,'epsilon_relative')==er and num(r,'M')==M and r['variant']==v]
     if not rr:
      blocked=24*8*nr*(2**nl)*M**nn/2**30>8
      cells.append('Memory limit' if blocked else 'Missing');continue
     assert len(rr)==1;r=rr[0]
     t=timecell(r,num(r,'time_limit'))+' s'
     nodes=nfmt(num(r,'nodes_visited'))
     if v=='with_m2':nodes+=' ('+nfmt(num(r,'nodes_cheap'))+', '+nfmt(num(r,'nodes_expensive'))+')'
     lines=[t,nodes]
     if not success(r):lines.append('gap '+nfmt(num(r,'gap_relative_percent'))+r'\%')
     cells.append(r'\shortstack{'+r'\\'.join(lines)+'}')
    body.append([str(er),str(M),f'{(2**nl)*M**nn+nr}/{2*nr+1}']+cells)
    if M!=17:body.append([])
  (out/f'ablation_{pid}.tex').write_text(table(['$\\epsilon_r$','$M$','$m_1$ vars./rows','With $m_2$','Without $m_2$'],body,
   f'{pid}: time and visited nodes (cheap, costly with screening). '+
   r'$N=2$, $\widehat N=M$. Final relative gap is shown at time limits. '+
   f'The screening QP has {(2**nl)*2**nn+nr} variables and {2*nr+1} rows.'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--comparison',default='reference/comparison_results.csv');p.add_argument('--ablation',default='reference/ablation_results.csv');p.add_argument('--output',default='reports');a=p.parse_args()
 build(rows(a.comparison),rows(a.ablation),Path(a.output));print('Tables written to',a.output)
