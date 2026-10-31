"""Exact UniProt OX taxonomy join; kingdom labels are not full lineages."""
import csv,gzip,hashlib,json,re,datetime,argparse
from pathlib import Path

def digest(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def classify(lineage):
 tokens={x.strip() for x in lineage.split(';')}
 viral='Viruses' in tokens
 cellular=bool(tokens & {'Bacteria','Archaea','Eukaryota'})
 return 'viral' if viral and not cellular else 'cellular' if cellular and not viral else 'unclassified'

def load_map(path):
 rows={}
 with open(path) as f:
  rd=csv.DictReader(f,delimiter='\t')
  if not {'subject_id','lineage'}<=set(rd.fieldnames or []):raise ValueError('taxonomy map requires subject_id and lineage')
  for r in rd:
   sid=r['subject_id']
   if not sid or sid in rows:raise ValueError('Empty/duplicate taxonomy subject: '+sid)
   r['classification']=classify(r['lineage']);rows[sid]=r
 return rows

def build(fasta,species,out):
 taxa={};release='';codes={'V':'Viruses','B':'Bacteria','A':'Archaea','E':'Eukaryota','O':'Other'}
 with open(species) as f:
  for line in f:
   if line.startswith('Release:'):release=line.strip()
   m=re.match(r'^\S+\s+([AB EVO])\s+(\d+): N=(.+)$',line)
   if m:
    code,tid,name=m.groups();value=(codes[code],name)
    if tid in taxa and taxa[tid][0]!=value[0]:raise ValueError('Conflicting species taxid '+tid)
    taxa[tid]=value
 if not release or len(taxa)<1000:raise ValueError('Invalid official species table')
 counts={};seen=set();fields=['subject_id','accession','taxid','lineage','scientific_name','mapping_status','source']
 with gzip.open(fasta,'rt') if str(fasta).endswith('.gz') else open(fasta) as src,open(out,'x') as dst:
  w=csv.DictWriter(dst,fields,delimiter='\t');w.writeheader()
  for line in src:
   if not line.startswith('>'):continue
   sid=line[1:].split()[0]
   if sid in seen:raise ValueError('Duplicate FASTA ID '+sid)
   seen.add(sid);m=re.search(r'\bOX=(\d+)\b',line);tid=m.group(1) if m else ''
   kingdom,name=taxa.get(tid,('',''))
   status='matched_taxid' if kingdom else 'taxid_not_in_species_table' if tid else 'missing_OX'
   group=classify(kingdom);counts[group]=counts.get(group,0)+1
   w.writerow(dict(subject_id=sid,accession=sid.split('|')[1] if '|' in sid else sid,taxid=tid,lineage=kingdom,scientific_name=name,mapping_status=status,source='UniProt speclist; '+release))
 meta=dict(source_url='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/speclist.txt',source_release=release,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),fasta_path=str(Path(fasta).resolve()),fasta_sha256=digest(fasta),species_sha256=digest(species),map_sha256=digest(out),entries=len(seen),classification_counts=counts,taxonomic_resolution='kingdom_only',note='Current taxonomy applied to local FASTA OX IDs; not a matched historical release. Missing IDs remain unclassified.')
 Path(str(out)+'.json').write_text(json.dumps(meta,indent=2))
 assert len(load_map(out))==len(seen)
 print(json.dumps(meta,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--species',required=True);p.add_argument('--out',required=True);a=p.parse_args();build(a.fasta,a.species,a.out)
