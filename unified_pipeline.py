#!/usr/bin/env python3
"""Enhanced discovery followed by auditable annotation and reporting."""
import argparse, json, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
 p=argparse.ArgumentParser(description=__doc__, epilog="Unrecognized options are forwarded to enhanced discovery. Use --existing-results for reporting without rerunning searches.")
 p.add_argument('--enhanced', action='store_true', help=argparse.SUPPRESS)
 p.add_argument('--existing-results', type=Path)
 p.add_argument('--out', required=True, type=Path)
 p.add_argument('--vmr', required=True, type=Path)
 p.add_argument('--nt-fasta', required=True, type=Path)
 p.add_argument('--protein-fasta', required=True, type=Path)
 p.add_argument('--plots', choices=['all','supported','none'], default='supported')
 p.add_argument('--pdf', action='store_true')
 p.add_argument('--mapping-breadth', type=float, default=0.8)
 p.add_argument('--mapping-depth', type=float, default=5.0)
 a,extra=p.parse_known_args()
 for f in [a.vmr,a.nt_fasta,a.protein_fasta]:
  if not f.is_file():p.error('Missing file: '+str(f))
 if not 0<=a.mapping_breadth<=1 or a.mapping_depth<0:p.error('Invalid coverage thresholds')
 if a.existing_results and extra:p.error('Discovery options cannot be used with --existing-results: '+' '.join(extra))
 if '--preflight' in extra:
  if a.existing_results:p.error('--preflight applies to discovery only')
  subprocess.run([sys.executable,str(ROOT/'enhanced_v1/pipeline.py'),'--out',str(a.out/'discovery')]+extra,check=True)
  return
 a.out.mkdir(parents=True,exist_ok=False)
 manifest={'status':'running','arguments':{k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},'discovery_arguments':extra}
 def save(): (a.out/'run_manifest.json').write_text(json.dumps(manifest,indent=2))
 save()
 try:
  source=a.existing_results
  if source is None:
   source=a.out/'discovery'
   subprocess.run([sys.executable,str(ROOT/'enhanced_v1/pipeline.py'),'--out',str(source)]+extra,check=True)
  from report import build
  result=build(source,a.out/'report',a.vmr,a.nt_fasta,a.protein_fasta,a.plots,a.pdf,a.mapping_breadth,a.mapping_depth)
  manifest.update(status='completed',validation=result);save()
 except BaseException as exc:
  manifest.update(status='failed',error=repr(exc));save();raise
if __name__=='__main__':main()
