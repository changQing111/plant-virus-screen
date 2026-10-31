#!/usr/bin/env python3
"""Conservative candidate discovery; no species/host confirmation is inferred."""
import argparse,csv,gzip,json,os,shutil,subprocess,sys,time
from pathlib import Path
from collections import defaultdict
from taxonomy import load_map, digest
from Bio import SeqIO
from Bio.Seq import Seq
BASE='/data/changq/miniforge3/'
def table(path,fields,rows):
 with open(path,'w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
def records(path):
 with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:
  yield from SeqIO.parse(f,'fasta')
def hits(path,n):
 with open(path) as f:
  for line in f:
   if line.strip() and not line.startswith('#'):
    a=line.rstrip('\n').split('\t')
    if len(a)!=n: raise ValueError(f'{path}: expected {n} columns, got {len(a)}')
    yield a
def homologs(nt,aa,ids):
 evidence=defaultdict(set)
 for path,n,kind in [(nt,17,'nt'),(aa,15,'protein')]:
  for a in hits(path,n):
   if a[0] not in ids: raise ValueError('Hit ID absent from FASTA: '+a[0])
   if int(a[2])!=ids[a[0]]: raise ValueError('Query length mismatch: '+a[0])
   ev,score=float(a[-2]),float(a[-1])
   if int(a[5]) >= (300 if kind=='nt' else 100) and ev <= (1e-10 if kind=='nt' else 1e-5) and score >= (200 if kind=='nt' else 80): evidence[a[0]].add(kind)
 return evidence
def domain_hits(path,orfmap):
 out=[]
 with open(path) as f:
  for l in f:
   if l.startswith('#') or not l.strip():continue
   a=l.split(); oid=a[0]
   if oid not in orfmap:raise ValueError('Unknown ORF '+oid)
   cov=(int(a[16])-int(a[15])+1)/int(a[5])
   if float(a[12])<=1e-5 and cov>=0.25:
    out.append(dict(query_id=orfmap[oid],orf_id=oid,model=a[3],accession=a[4],i_evalue=a[12],hmm_coverage=cov,aa_start=a[17],aa_end=a[18]))
 return out
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--query',required=True);p.add_argument('--out',required=True)
 p.add_argument('--reads1');p.add_argument('--reads2');p.add_argument('--blastn');p.add_argument('--blastx')
 p.add_argument('--threads',type=int,default=8)
 p.add_argument('--nt-db',default='/data/changq/Document/viral/viral_blastn_db/viral')
 p.add_argument('--viral-db',default='/data/changq/Document/viral/uniprot_viral.dmnd')
 p.add_argument('--competitive-db',default='/data/changq/Document/uniport_db/uniprot_sprot.dmnd')
 p.add_argument('--competitive-fasta',default='/data/changq/Document/uniport_db/uniprot_sprot.fasta.gz')
 p.add_argument('--taxonomy-map',help='TSV: subject_id, lineage (semicolon-delimited exact taxon names). Verified bundled map auto-used for default Swiss-Prot.')
 p.add_argument('--pfam',default='/data/changq/Document/Pfam_hmm_db/Pfam-A.hmm')
 p.add_argument('--no-mapping',action='store_true',help='Explicitly mark mapping not performed')
 p.add_argument('--preflight',action='store_true')
 a=p.parse_args()
 taxonomy_provenance={}
 default_map=Path(__file__).resolve().parent/'taxonomy_data/swissprot_taxonomy.tsv'
 if not a.taxonomy_map and default_map.is_file() and str(Path(a.competitive_fasta).resolve())=='/data/changq/Document/uniport_db/uniprot_sprot.fasta.gz' and str(Path(a.competitive_db).resolve())=='/data/changq/Document/uniport_db/uniprot_sprot.dmnd':
  a.taxonomy_map=str(default_map)
 if a.taxonomy_map:
  meta=Path(str(a.taxonomy_map)+'.json')
  if meta.is_file():
   taxonomy_provenance=json.loads(meta.read_text())
   if taxonomy_provenance.get('fasta_sha256')!=digest(a.competitive_fasta):raise ValueError('Taxonomy FASTA digest mismatch')
   if taxonomy_provenance.get('map_sha256')!=digest(a.taxonomy_map):raise ValueError('Taxonomy map digest mismatch')
   if taxonomy_provenance.get('database_sha256') and taxonomy_provenance['database_sha256']!=digest(a.competitive_db):raise ValueError('Taxonomy database digest mismatch')
  else:taxonomy_provenance={'source':'user_supplied_map','map_sha256':digest(a.taxonomy_map)}
 taxonomy_rows=load_map(a.taxonomy_map) if a.taxonomy_map else {}
 if a.threads<1:p.error('threads must be positive')
 if bool(a.blastn)!=bool(a.blastx):p.error('provide both existing search tables')
 if not a.no_mapping and not(a.reads1 and a.reads2):p.error('paired reads required, or explicitly --no-mapping')
 tools={'diamond':BASE+'envs/kofamscan/bin/diamond','hmmsearch':BASE+'envs/kofamscan/bin/hmmsearch','blastn':BASE+'bin/blastn','bwa':BASE+'bin/bwa','samtools':BASE+'envs/kofamscan/bin/samtools'}
 for k,v in tools.items():
  tools[k]=os.environ.get(k.upper()+'_BIN',v)
  if a.no_mapping and k in ('bwa','samtools'):continue
  if a.blastn and k=='blastn':continue
  if not os.access(tools[k],os.X_OK):raise FileNotFoundError(tools[k])
 for f in [a.query,a.pfam,a.viral_db,a.competitive_db,a.competitive_fasta,a.reads1,a.reads2,a.blastn,a.blastx,a.taxonomy_map]:
  if f and not Path(f).is_file():raise FileNotFoundError(f)
 if not a.blastn and not list(Path(a.nt_db).parent.glob(Path(a.nt_db).name+'.n*')):raise FileNotFoundError(a.nt_db)
 if a.preflight:print(json.dumps({'status':'preflight_passed','tools':tools,'arguments':vars(a)},indent=2));return
 out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False)
 manifest={'status':'running','arguments':vars(a),'tools':tools,'started':time.time(),'commands':[],'taxonomy_provenance':taxonomy_provenance}
 def save(): (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
 def run(args,stdout=None):
  args=list(map(str,args));manifest['commands'].append(args);save()
  with open(out/'commands.log','a') as log:
   log.write(repr(args)+'\n');log.flush();subprocess.run(args,check=True,stdout=stdout or log,stderr=log)
 save()
 try:
  ids={}
  for r in records(a.query):
   if r.id in ids:raise ValueError('Duplicate FASTA ID: '+r.id)
   ids[r.id]=len(r)
  if not ids:raise ValueError('Empty FASTA')
  nt=Path(a.blastn) if a.blastn else out/'blastn.tsv';aa=Path(a.blastx) if a.blastx else out/'blastx.tsv'
  if not a.blastn:
   run([tools['blastn'],'-task','blastn','-query',a.query,'-db',a.nt_db,'-evalue','1e-5','-num_threads',a.threads,'-outfmt','6 qseqid sseqid qlen slen pident length qcovs qcovhsp qcovus mismatch gapopen qstart qend sstart send evalue bitscore','-out',nt])
   run([tools['diamond'],'blastx','--very-sensitive','--query',a.query,'--db',a.viral_db,'--threads',a.threads,'--evalue','1e-5','--max-target-seqs','25','--outfmt','6','qseqid','sseqid','qlen','slen','pident','length','qcovhsp','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','--out',aa])
  evidence=homologs(nt,aa,ids)
  # Select RdRP/RNA replicase profiles by curated Pfam descriptions; retain model provenance.
  models=[];block=[]
  with open(a.pfam) as src,open(out/'rdrp.hmm','w') as dst:
   for line in src:
    block.append(line)
    if line.strip()=='//':
     desc=' '.join(x.strip() for x in block if x.startswith(('NAME ','DESC ')))
     if any(x in desc.lower() for x in ('rdrp','rna-dependent rna polymerase','rna dependent rna polymerase','rna replicase')):
      dst.writelines(block);models.append(desc)
     block=[]
  if not models:raise ValueError('No RdRP profiles selected')
  (out/'selected_models.txt').write_text('\n'.join(models)+'\n')
  orfmap={};orfrows=[]
  with open(out/'orfs.faa','w') as f:
   for r in records(a.query):
    for strand,seq in [('+',r.seq),('-',r.seq.reverse_complement())]:
     for frame in range(3):
      pep=str(seq[frame:len(seq)-(len(seq)-frame)%3].translate());start=0
      for part in pep.split('*'):
       if len(part)>=80:
        oid='orf'+str(len(orfmap)+1);orfmap[oid]=r.id
        lo=frame+3*start;hi=lo+3*len(part)
        left,right=(lo+1,hi) if strand=='+' else (len(r)-hi+1,len(r)-lo)
        orfrows.append(dict(orf_id=oid,query_id=r.id,strand=strand,frame=frame+1,nt_start=left,nt_end=right,aa_length=len(part)))
        f.write('>'+oid+'\n'+part+'\n')
       start+=len(part)+1
  table(out/'orf_coordinates.tsv',['orf_id','query_id','strand','frame','nt_start','nt_end','aa_length'],orfrows)
  dom=out/'domains.domtbl'
  if orfmap:run([tools['hmmsearch'],'--cpu',a.threads,'--domtblout',dom,'-E','1e-3','--domE','1e-3',out/'rdrp.hmm',out/'orfs.faa'])
  else:dom.write_text('')
  dh=domain_hits(dom,orfmap)
  table(out/'domain_evidence.tsv',['query_id','orf_id','model','accession','i_evalue','hmm_coverage','aa_start','aa_end'],dh)
  for d in dh:evidence[d['query_id']].add('RdRP_domain')
  candidates=set(evidence)
  with open(out/'candidates.fasta','w') as f:
   for r in records(a.query):
    if r.id in candidates:SeqIO.write(r,f,'fasta')
  (out/'candidate_ids.txt').write_text(''.join(q+'\n' for q in ids if q in candidates))
  comp=out/'competitive.tsv';best={};comp_rows=[]
  if candidates:
   run([tools['diamond'],'blastx','--very-sensitive','--query',out/'candidates.fasta','--db',a.competitive_db,'--threads',a.threads,'--evalue','1e-5','--max-target-seqs','50','--outfmt','6','qseqid','sseqid','pident','length','qstart','qend','sstart','send','evalue','bitscore','--out',comp])
   subjects={h[1] for h in hits(comp,10)};annotations={}
   for r in records(a.competitive_fasta):
    if r.id in subjects:annotations[r.id]=r.description
   for h in hits(comp,10):
    tax=taxonomy_rows.get(h[1],{})
    group=tax.get('classification','unclassified')
    row=dict(query_id=h[0],subject_id=h[1],bitscore=h[9],evalue=h[8],qstart=h[4],qend=h[5],classification=group,annotation=annotations.get(h[1],''),taxid=tax.get('taxid',''),lineage=tax.get('lineage',''),scientific_name=tax.get('scientific_name',''),mapping_status=tax.get('mapping_status','provided_lineage' if tax else 'subject_not_mapped'),taxonomy_source=tax.get('source','user_supplied_map' if tax else ''))
    comp_rows.append(row)
    if h[0] not in best or float(h[9])>float(best[h[0]]['bitscore']):best[h[0]]=row
  else:comp.write_text('')
  table(out/'competitive_annotated.tsv',['query_id','subject_id','bitscore','evalue','qstart','qend','classification','annotation','taxid','lineage','scientific_name','mapping_status','taxonomy_source'],comp_rows)
  manifest['competitive_taxonomy_counts']={k:sum(r['classification']==k for r in comp_rows) for k in ('viral','cellular','unclassified')}
  manifest['competitive_unmapped_subjects']=len({r['subject_id'] for r in comp_rows if r['mapping_status']=='subject_not_mapped'})
  if not a.no_mapping and candidates:
   # Full assembly background, no host reference required. Source assembly remains untouched.
   with open(out/'mapping_reference.fasta','w') as ref:
    SeqIO.write(records(a.query),ref,'fasta')
   run([tools['bwa'],'index',out/'mapping_reference.fasta'])
   with open(out/'mapping.sam','w') as f:run([tools['bwa'],'mem','-t',a.threads,out/'mapping_reference.fasta',a.reads1,a.reads2],f)
   run([tools['samtools'],'sort','-@',a.threads,'-o',out/'mapping.bam',out/'mapping.sam'])
   run([tools['samtools'],'index',out/'mapping.bam'])
   run([tools['samtools'],'quickcheck',out/'mapping.bam'])
   with open(out/'mapping.flagstat.txt','w') as f:run([tools['samtools'],'flagstat',out/'mapping.bam'],f)
   bed=out/'candidates.bed';bed.write_text(''.join(f'{q}\t0\t{ids[q]}\n' for q in ids if q in candidates))
   for name,mq in [('all',0),('mapq20',20)]:
    with open(out/f'depth_{name}.tsv','w') as f:run([tools['samtools'],'depth','-a','-b',bed,'-Q',mq,out/'mapping.bam'],f)
    stats={q:[0,0,0] for q in candidates}
    for h in hits(out/f'depth_{name}.tsv',3):
     if h[0] in stats:
      d=int(h[2]);s=stats[h[0]];s[0]+=d;s[1]+=d>0;s[2]+=d>=5
    table(out/f'coverage_{name}.tsv',['query_id','length','mean_depth','breadth_1x','breadth_5x'],[dict(query_id=q,length=ids[q],mean_depth=s[0]/ids[q],breadth_1x=s[1]/ids[q],breadth_5x=s[2]/ids[q]) for q,s in sorted(stats.items())])
   (out/'mapping.sam').unlink()
  rows=[dict(query_id=q,length=ids[q],retained=int(q in candidates),evidence=';'.join(sorted(evidence.get(q,[]))),competitive_best_class=best[q]['classification'] if q in best else 'no_hit',mapping_status='not_requested' if a.no_mapping else ('completed' if candidates else 'no_candidates'),host_status='undetermined',interpretation='candidate_requires_validation' if q in candidates else 'no_evidence_under_current_search') for q in ids]
  table(out/'all_contigs_audit.tsv',list(rows[0]),rows)
  table(out/'final_candidates.tsv',list(rows[0]),[r for r in rows if r['retained']])
  observed={r.id for r in records(out/'candidates.fasta')}
  assert observed==candidates and len(rows)==len(ids)
  manifest.update(status='completed',input_contigs=len(ids),retained_contigs=len(candidates),domain_hits=len(dh),mapping_performed=not a.no_mapping and bool(candidates),taxonomy_status='provided' if a.taxonomy_map else 'unclassified_requires_lineage_map',finished=time.time());save()
  print(json.dumps(manifest))
 except BaseException as e:
  manifest.update(status='failed',error=repr(e));save();raise
if __name__=='__main__':main()
