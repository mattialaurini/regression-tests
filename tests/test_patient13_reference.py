"""Check the supplied illustration independently of MATLAB and solvers."""
import csv,math,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/'pkpd/reference'
def read(name):
 with (ROOT/name).open(newline='') as f:return list(csv.DictReader(f))
class Patient13(unittest.TestCase):
 def test_complete_certified_six_model_reference(self):
  rows=read('patient13_metrics.csv')
  self.assertEqual([r['model'] for r in rows],['ARX-1','ARX-2','ARX-3','ARX-4','ML-DP','ML-1C'])
  self.assertTrue(all(r['certified'].lower() in ('1','true') for r in rows))
  self.assertTrue(all(float(r['sigma'])==1 for r in rows))
 def test_input_and_reconstruction_metrics(self):
  rows=read('patient13_signals_sigma1.csv');self.assertEqual(len(rows),150)
  metrics=read('patient13_metrics.csv')
  for k,r in enumerate(rows):
   self.assertEqual(float(r['time_s']),k)
   self.assertEqual(float(r['u']),10 if k<10 else 3 if k<35 else 0)
   self.assertAlmostEqual(float(r['BIS_ARX_used']),min(float(r['BIS_raw']),93.1),places=10)
  for column,m in zip(['ARX1','ARX2','ARX3','ARX4','ML_DP','ML_1C'],metrics):
   rmse=math.sqrt(sum((float(r[column])-float(r['BIS_true']))**2 for r in rows)/len(rows))
   self.assertAlmostEqual(rmse,float(m['RMSE_true_BIS']),places=9)
 def test_paper_illustration_numbers(self):
  rows=read('patient13_metrics.csv')
  for row,target in zip(rows,[29.909,.360,.976,.951,.226,.222]):
   self.assertEqual(round(float(row['RMSE_true_BIS']),3),target)
  self.assertAlmostEqual(float(rows[4]['search_time_s']),8.237077833333334,places=9)
  self.assertAlmostEqual(float(rows[5]['search_time_s']),236.39300966666667,places=9)
if __name__=='__main__':unittest.main()
