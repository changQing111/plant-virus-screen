"""Evidence prioritization, not infection diagnosis."""
import csv,json,hashlib,sys,subprocess
from pathlib import Path
from collections import Counter,defaultdict
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'legacy'))
import select_final_candidates as legacy
import summarize_final_retained_viruses as summary

def read(p):
 with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields):
 with open(p,'w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def grade(evidence,competition,mapping_ok):
 if competition=='cellular':return 'D_conflicting','cellular_competitive_hit'
 n=len(set(evidence.split(';')) & {'nt','protein','RdRP_domain'})
 if n>=2 and competition=='viral' and mapping_ok:return 'A_supported','multiple_evidence_viral_competition_mapping'
 if n>=2:return 'B_supported','multiple_evidence_incomplete_validation'
 return 'C_provisional','single_evidence_or_unresolved'
def build(src,out,vmr,ntf,aaf,plots,pdf,breadth,depth):
 src=Path(src);out.mkdir(parents=True,exist_ok=False)
 required=['final_candidates.tsv','all_contigs_audit.tsv','blastn.tsv','blastx.tsv','competitive_annotated.tsv','domain_evidence.tsv','manifest.json']
 for n in required:
  if not (src/n).is_file():raise ValueError('Missing '+n)
 manifest=json.loads((src/'manifest.json').read_text())
 if not isinstance(manifest,dict) or manifest.get('status')!='completed':raise ValueError('Invalid discovery manifest')
 rows=read(src/'final_candidates.tsv'); ids=[r['query_id'] for r in rows];idset=set(ids)
 if len(ids)!=len(idset):raise ValueError('Duplicate candidate IDs')
 audit=read(src/'all_contigs_audit.tsv')
 if {r['query_id'] for r in audit if r['retained']=='1'}!=idset:raise ValueError('Audit/candidate mismatch')
 best={};wanted=set()
 for name,n,kind,minlen,ev,score in [('blastn.tsv',17,'nt',300,1e-10,200),('blastx.tsv',15,'protein',100,1e-5,80)]:
  with open(src/name) as f:
   for line in f:
    a=line.rstrip('\n').split('\t')
    if len(a)!=n:raise ValueError('Unexpected search format')
    if a[0] not in idset:continue
    if int(a[5])<minlen or float(a[-2])>ev or float(a[-1])<score:continue
    key=(a[0],kind)
    if key not in best or float(a[-1])>float(best[key][-1]):best[key]=a
 for a in best.values():wanted.add(a[1])
 labels={}
 for path in [ntf,aaf]:
  with open(path) as f:
   for line in f:
    if line.startswith('>'):
     text=line[1:].strip();acc=text.split()[0]
     if acc in wanted:labels[acc]=text
 coverage={}
 if (src/'coverage_mapq20.tsv').exists():
  coverage={r['query_id']:r for r in read(src/'coverage_mapq20.tsv')}
  if set(coverage)!=idset:raise ValueError('Coverage/candidate mismatch')
 annotated=[]
 for r in rows:
  q=r['query_id'];x=dict(r);x['contig_id']=q;x['representative_query_id']=q
  for kind,key in [('nt','best_nt_hit'),('protein','best_protein_hit')]:
   a=best.get((q,kind));x[key]=''
   if a:
    label=labels.get(a[1],a[1]);x[key]=(a[1]+' | '+label[len(a[1]):].strip()) if kind=='nt' else (a[1]+' | '+label[len(a[1]):].strip())
    x['blastn_anno' if kind=='nt' else 'blastx_anno']=label[len(a[1]):].strip()
  x['putative_virus']=legacy.infer_putative_virus_name(x)
  if x['putative_virus']=='Unclassified virus-like candidate':x['putative_virus']+=' ['+q+']'
  x['putative_segment_or_protein']=legacy.infer_segment_or_protein(x)
  x['candidate_name']=x['putative_virus'];x['group_id']=x['putative_virus']
  c=coverage.get(q,{})
  ok=bool(c) and float(c['breadth_1x'])>=breadth and float(c['mean_depth'])>=depth
  x['evidence_grade'],x['grade_reason']=grade(r['evidence'],r['competitive_best_class'],ok)
  x['mapq20_breadth_1x']=c.get('breadth_1x','');x['mapq20_mean_depth']=c.get('mean_depth','')
  x['priority']='primary' if x['evidence_grade'].startswith(('A_','B_')) else 'secondary'
  x['host_status']='undetermined';x['classification_scope']='reference_similarity_not_species_confirmation'
  annotated.append(x)
 legacy.annotate_with_vmr(annotated,vmr)
 for x in annotated:
  x['vmr_annotation_scope']='legacy_accession_name_or_genus_mapping;reference_host_only'
  if x.get('plant_virus_status')=='non-plant virus':x['plant_virus_status']='plant_host_not_documented_by_this_mapping'
 fields=list(dict.fromkeys(k for r in annotated for k in r)) or ['contig_id','evidence_grade']
 write(out/'annotated_candidates.tsv',annotated,fields)
 summary.write_summary(out/'virus_summary_legacy.tsv',summary.summarize_rows(annotated))
 grouped=defaultdict(list)
 for x in annotated:grouped[x['putative_virus']].append(x)
 groups=[]
 for name,rr in sorted(grouped.items()):
  groups.append({'reference_related_group':name,'contig_count':len(rr),'grade_counts':json.dumps(dict(Counter(r['evidence_grade'] for r in rr)),sort_keys=True),'contig_ids':';'.join(r['contig_id'] for r in rr),'interpretation':'reference_related_group_not_confirmed_species'})
 write(out/'virus_evidence_summary.tsv',groups,['reference_related_group','contig_count','grade_counts','contig_ids','interpretation'])
 selected=[r for r in annotated if plots=='all' or (plots=='supported' and r['evidence_grade'].startswith(('A_','B_')))] if plots!='none' else []
 write(out/'plot_candidates.tsv',selected,fields)
 if selected:
  for kind,script,mod,table in [('nt','batch_visualize_references.py','visualize_blastn_reference.py','blastn.tsv'),('protein','batch_visualize_proteins.py','visualize_blastx_protein.py','blastx.tsv')]:
   cmd=[sys.executable,str(ROOT/'legacy'/script),'--final-tsv',str(out/'plot_candidates.tsv'),'--blastn' if kind=='nt' else '--blastx',str(src/table),'--module',str(ROOT/'legacy'/mod),'--outdir',str(out/(kind+'_plots'))]
   if kind=='protein':cmd+=['--protein-fasta',str(aaf)]
   if not pdf:cmd+=['--no-pdf']
   with open(out/(kind+'_plots.log'),'w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
 figures=list(out.rglob('*.html'))
 if selected and any(r.get('best_nt_hit') or r.get('best_protein_hit') for r in selected) and not figures:raise ValueError('Plotting requested but no HTML figures generated')
 counts=dict(Counter(r['evidence_grade'] for r in annotated));check=read(out/'annotated_candidates.tsv')
 assert [r['contig_id'] for r in check]==ids
 assert sum(g['contig_count'] for g in groups)==len(ids)
 validation={'status':'passed','input_candidates':len(ids),'output_candidates':len(check),'grade_counts':counts,'plot_candidates':len(selected),'html_figures':len(figures),'input_sha256':{n:sha(src/n) for n in required},'vmr_sha256':sha(vmr),'thresholds':{'mapq20_breadth_1x':breadth,'mapq20_mean_depth':depth},'grading_version':'0.1.0-heuristic'}
 (out/'validation.json').write_text(json.dumps(validation,indent=2))
 (out/'REPORT.md').write_text('# 病毒候选证据报告\n\n候选数：'+str(len(ids))+'。全部保留，无静默剔除。\n\n分级计数：'+json.dumps(counts)+'\n\nA：多类证据、病毒竞争最佳命中、回贴达到阈值；B：多类证据但验证不完整；C：单类或未解决证据；D：细胞来源竞争冲突。分级是启发式复核优先级，不是感染、物种或新病毒确认。\n\nICTV 注释复用传统 accession/name/genus 映射，宿主字段仅描述参考病毒。汇总按推断参考相关名称分组，组数不是病毒物种数。竞争证据当前只使用最高分分类，尚未完成区域级病毒/细胞分数竞争，高等级仍需人工复核。\n',encoding='utf-8')
 return validation
