#!/usr/bin/env python3
"""Check the distribution inventory, frozen sources and saved-result structure."""
import csv,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
 inventory=json.loads((ROOT/'PACKAGE_SHA256.json').read_text())
 for name,expected in inventory.items():
  file=ROOT/name
  assert file.is_file(),f'Missing distribution file: {name}'
  assert hashlib.sha256(file.read_bytes()).hexdigest()==expected,f'Changed file: {name}'
 mapping=json.loads((ROOT/'pkpd/source_mapping.json').read_text())
 for row in mapping:
  file=ROOT/'pkpd/src'/row['distributed']
  assert hashlib.sha256(file.read_bytes()).hexdigest()==row['distributed_sha256']
 for file in ROOT.rglob('*'):
  if not file.is_file():continue
  assert file.name not in ('.DS_Store','octave-workspace')
  assert '__MACOSX' not in file.parts and 'legacy_not_in_paper' not in file.parts
 assert (ROOT/'licenses/PriorityQueue-BSD.txt').is_file()
 assert (ROOT/'LICENSE').read_text().startswith('# PolyForm Noncommercial License 1.0.0')
 print(f'PASS: {len(inventory)} distribution files and {len(mapping)} mapped PK/PD source files.')
if __name__=='__main__':main()
