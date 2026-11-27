#!/usr/bin/env python3
"""Verify frozen numerical sources and published result invariants (stdlib only)."""
import csv,hashlib,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def rows(name):
 with (ROOT/'reference'/name).open(newline='') as f:return list(csv.DictReader(f))
def f(r,k):return float(r[k])
def yes(r,k='certified'):return r[k].lower() in ('true','1','1.0')
def main():
 hashes=json.loads((ROOT/'reference/numerical_source_sha256.json').read_text())
 for name,expected in hashes.items():
  assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,('Source changed',name)
 comp=rows('comparison_results.csv'); ab=rows('ablation_results.csv'); progress=rows('ablation_progress.csv')
 assert len(comp)==48 and len(ab)==180 and len(progress)==192
 assert len({(r['problem_id'],r['epsilon_relative'],r['solver']) for r in comp})==48
 assert len({r['task_id'] for r in ab})==180
 assert sum(r['execution_state']=='resource_limit' for r in progress)==12
 for r in comp+ab:
  U,L,er,ea=[f(r,k) for k in ['objective_scaled','dual_bound_scaled','epsilon_relative','epsilon_absolute']]
  assert yes(r)==(0<=L<=U and L>=(1-er)*U-ea),(r['problem_id'],r['epsilon_relative'])
  assert not r['error_message']
 for r in ab:
  assert f(r,'M')==f(r,'Nhat') and f(r,'M') in (3,5,9,17)
  assert f(r,'nodes_cheap')+f(r,'nodes_expensive')==f(r,'nodes_visited')
  assert f(r,'pruned_tolerance_only')+f(r,'pruned_at_incumbent')==f(r,'pruned_total')
  assert abs(f(r,'gap_absolute_scaled')-(f(r,'objective_scaled')-f(r,'dual_bound_scaled')))<1e-10
  if r['variant']=='with_m2':assert f(r,'N')==2
  else:assert f(r,'m2_evaluations')==0
  assert f(r,'time_limit')==f(r,'target_time_limit') or (yes(r) and f(r,'epsilon_relative')==.01 and f(r,'time_limit')==600)
 for solver,counts in [('BnB',[7,8,8]),('SCIP',[1,0,0])]:
  assert [sum(yes(r) for r in comp if r['solver']==solver and f(r,'epsilon_relative')==e) for e in [.1,.01,.001]]==counts
 print(f'PASS: {len(hashes)} frozen source files, 48 comparison results, 180 ablation results, 12 memory exclusions.')
if __name__=='__main__':main()
