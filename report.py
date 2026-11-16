"""Evidence prioritization, not infection diagnosis."""
import csv,json,hashlib,sys,subprocess
from enhanced_v1.evidence_policy import regional_competition, discovery_priority
from pathlib import Path
from collections import Counter,defaultdict
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'legacy'))
import select_final_candidates as legacy
import summarize_final_retained_viruses as summary
import plot_index
import visualize_summary as summary_figs
import visualize_composition as composition_figs

def read(p):
 with open(p,encoding='utf-8') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields):
 with open(p,'w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def grade(evidence,competition,mapping_ok):
 if competition=='cellular':return 'D_conflicting','cellular_competitive_hit'
 kinds={x.strip() for x in (evidence or '').split(';') if x.strip()}
 n=len(kinds & {'nt','protein','RdRP_domain'})
 # DIAMOND BLASTX and the RdRp HMM generally hit the same polymerase region,
 # so protein+RdRP_domain is not independent confirmation. Require a
 # nucleotide hit as well before calling a candidate A_supported.
 if n>=2 and 'nt' in kinds and competition=='viral' and mapping_ok:return 'A_supported','multiple_evidence_viral_competition_mapping_with_nt'
 if n>=2:return 'B_supported','multiple_evidence_incomplete_validation'
 return 'C_provisional','single_evidence_or_unresolved'

def evidence_values(text):return {x.strip() for x in (text or '').split(';') if x.strip()}

def validate_discovery_semantics(src,rows,audit):
 audit_by={}
 for r in audit:
  q=r.get('query_id','')
  if not q or q in audit_by:raise ValueError('Empty/duplicate audit query ID: '+q)
  audit_by[q]=r
 candidate_by={r['query_id']:r for r in rows}
 best={};observed=defaultdict(set)
 specs=[('blastn.tsv',17,'nt',300,1e-10,200),('blastx.tsv',15,'protein',100,1e-5,80)]
 for name,n,kind,minlen,ev,score in specs:
  with open(src/name,encoding='utf-8') as f:
   for line_no,line in enumerate(f,1):
    if not line.strip() or line.startswith('#'):continue
    a=line.rstrip('\n').split('\t')
    if len(a)!=n:raise ValueError('%s:%s: expected %s columns, got %s'%(name,line_no,n,len(a)))
    q=a[0]
    if q not in audit_by:raise ValueError('%s:%s: query absent from audit: %s'%(name,line_no,q))
    if int(a[2])!=int(audit_by[q]['length']):raise ValueError('%s:%s: query length mismatch: %s'%(name,line_no,q))
    if int(a[5])<minlen or float(a[-2])>ev or float(a[-1])<score:continue
    observed[q].add(kind);key=(q,kind)
    if key not in best or float(a[-1])>float(best[key][-1]):best[key]=a
 search_evidence_queries={q for q,values in observed.items() if values}
 domain_rows=read(src/'domain_evidence.tsv')
 for i,r in enumerate(domain_rows,2):
  q=r.get('query_id','')
  if q not in audit_by:raise ValueError('domain_evidence.tsv:%s: query absent from audit: %s'%(i,q))
  if float(r['i_evalue'])>1e-5 or float(r['hmm_coverage'])<0.25:raise ValueError('domain_evidence.tsv:%s: hit does not meet discovery thresholds'%i)
  observed[q].add('RdRP_domain')
 expected_ids={q for q,values in observed.items() if values}
 for q,r in audit_by.items():
  expected=observed.get(q,set());actual=evidence_values(r.get('evidence',''))
  if actual!=expected:raise ValueError('Audit evidence mismatch for %s: expected %s, got %s'%(q,';'.join(sorted(expected)),';'.join(sorted(actual))))
  if (r.get('retained')=='1')!=bool(expected):raise ValueError('Audit retained/evidence mismatch for '+q)
 if set(candidate_by)!=expected_ids:raise ValueError('Candidate set does not match recomputed evidence union')
 for q,r in candidate_by.items():
  audit_row=audit_by[q]
  differing=sorted(k for k in set(r)|set(audit_row) if r.get(k)!=audit_row.get(k))
  if differing:raise ValueError('Candidate/audit row mismatch for %s: %s'%(q,','.join(differing)))
  actual=evidence_values(r.get('evidence',''));expected=observed[q]
  if actual!=expected:raise ValueError('Candidate evidence mismatch for %s: expected %s, got %s'%(q,';'.join(sorted(expected)),';'.join(sorted(actual))))
 competitive_rows=read(src/'competitive_annotated.tsv');competitive_best={}
 for i,r in enumerate(competitive_rows,2):
  q=r.get('query_id','');classification=r.get('classification','')
  if q not in candidate_by:raise ValueError('competitive_annotated.tsv:%s: non-candidate query: %s'%(i,q))
  if classification not in {'viral','cellular','unclassified'}:raise ValueError('competitive_annotated.tsv:%s: invalid classification: %s'%(i,classification))
  if float(r['evalue'])>1e-5:raise ValueError('competitive_annotated.tsv:%s: hit does not meet discovery threshold'%i)
  score=float(r['bitscore'])
  if q not in competitive_best or score>competitive_best[q][0]:competitive_best[q]=(score,classification)
 meta=json.loads((src/'manifest.json').read_text()) if (src/'manifest.json').exists() else {}
 regional,_=regional_competition(competitive_rows) if meta.get('competition_method')=='regional_v1' else ({},[])
 for q,r in audit_by.items():
  expected=competitive_best[q][1] if q in competitive_best else 'no_hit'
  if meta.get('competition_method')=='regional_v1':expected=regional.get(q,'no_hit')
  if r.get('competitive_best_class')!=expected:raise ValueError('Audit competition mismatch for %s: expected %s, got %s'%(q,expected,r.get('competitive_best_class')))
  if q in candidate_by and candidate_by[q].get('competitive_best_class')!=expected:raise ValueError('Candidate competition mismatch for %s: expected %s, got %s'%(q,expected,candidate_by[q].get('competitive_best_class')))
 return best,{'search_evidence_queries':len(search_evidence_queries),'candidate_evidence_queries':len(expected_ids),'domain_rows':len(domain_rows),'competitive_rows':len(competitive_rows),'evidence_consistent':True,'competition_consistent':True}

def build(src,out,vmr,ntf,aaf,plots,pdf,breadth,depth):
 src=Path(src);out.mkdir(parents=True,exist_ok=False)
 required=['final_candidates.tsv','all_contigs_audit.tsv','blastn.tsv','blastx.tsv','competitive_annotated.tsv','domain_evidence.tsv','manifest.json']
 for n in required:
  if not (src/n).is_file():raise ValueError('Missing '+n)
 manifest=json.loads((src/'manifest.json').read_text(encoding='utf-8'))
 if not isinstance(manifest,dict) or manifest.get('status')!='completed':raise ValueError('Invalid discovery manifest')
 rows=read(src/'final_candidates.tsv'); ids=[r['query_id'] for r in rows];idset=set(ids)
 if len(ids)!=len(idset):raise ValueError('Duplicate candidate IDs')
 audit=read(src/'all_contigs_audit.tsv')
 if {r['query_id'] for r in audit if r['retained']=='1'}!=idset:raise ValueError('Audit/candidate mismatch')
 best,semantic_checks=validate_discovery_semantics(src,rows,audit);wanted=set()
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
  nt_ref=legacy.clean_reference_name(x.get('blastn_anno',''))
  os_name=legacy.extract_os_name(x.get('blastx_anno',''))
  protein_ref=legacy.clean_reference_name(os_name) if os_name else ''
  x['nt_reference_name']=nt_ref;x['protein_reference_name']=protein_ref
  x['reference_name_concordance']=legacy.reference_name_concordance(nt_ref,protein_ref)
  x['putative_virus_scope']='protein_reference_preferred;nt_and_protein_references_disagree' if x['reference_name_concordance']=='discordant' else 'reference_similarity_only'
  x['candidate_name']=x['putative_virus'];x['group_id']=x['putative_virus']
  c=coverage.get(q,{})
  ok=bool(c) and float(c['breadth_1x'])>=breadth and float(c['mean_depth'])>=depth
  x['evidence_grade'],x['grade_reason']=grade(r['evidence'],r['competitive_best_class'],ok)
  x['mapq20_breadth_1x']=c.get('breadth_1x','');x['mapq20_mean_depth']=c.get('mean_depth','')
  x['discovery_priority'],x['discovery_priority_reason']=discovery_priority(r['evidence'],r['competitive_best_class'],ok,manifest.get('hmm_policy')=='curated_core_v1')
  x['novelty_status']='not_assessed_requires_reference_comparison_and_phylogeny'
  x['assembly_validation_status']='not_assessed'
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
 selected=[r for r in annotated if plots=='all' or (plots=='supported' and (r['evidence_grade'].startswith(('A_','B_')) or r.get('discovery_priority','').startswith('novelty_review')))] if plots!='none' else []
 write(out/'plot_candidates.tsv',selected,fields)
 if selected:
  for kind,script,mod,table in [('nt','batch_visualize_references.py','visualize_blastn_reference.py','blastn.tsv'),('protein','batch_visualize_proteins.py','visualize_blastx_protein.py','blastx.tsv')]:
   cmd=[sys.executable,str(ROOT/'legacy'/script),'--final-tsv',str(out/'plot_candidates.tsv'),'--blastn' if kind=='nt' else '--blastx',str(src/table),'--module',str(ROOT/'legacy'/mod),'--outdir',str(out/(kind+'_plots'))]
   if kind=='protein':cmd+=['--protein-fasta',str(aaf)]
   if not pdf:cmd+=['--no-pdf']
   with open(out/(kind+'_plots.log'),'w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
  depth_tsv=src/'depth_mapq20.tsv'
  if depth_tsv.is_file():
   cmd=[sys.executable,str(ROOT/'legacy'/'batch_visualize_coverage.py'),'--final-tsv',str(out/'plot_candidates.tsv'),'--depth',str(depth_tsv),'--module',str(ROOT/'legacy'/'visualize_coverage.py'),'--outdir',str(out/'coverage_plots')]
   if not pdf:cmd+=['--no-pdf']
   with open(out/'coverage_plots.log','w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
 figures=[p for p in out.rglob('*.html') if p.name not in('index.html','report_summary.html','composition_summary.html')]
 if selected and any(r.get('best_nt_hit') or r.get('best_protein_hit') for r in selected) and not figures:raise ValueError('Plotting requested but no HTML figures generated')
 summary_figs.build_summary(out,breadth,depth,pdf)
 composition_figs.build_composition(out,pdf)
 if figures:plot_index.build_index(out)
 pdfs=[p for p in out.rglob('*.pdf') if p.name not in('report_summary.pdf','composition_summary.pdf')]
 counts=dict(Counter(r['evidence_grade'] for r in annotated));check=read(out/'annotated_candidates.tsv')
 assert [r['contig_id'] for r in check]==ids
 assert sum(g['contig_count'] for g in groups)==len(ids)
 validation={'status':'passed','input_candidates':len(ids),'output_candidates':len(check),'grade_counts':counts,'plot_candidates':len(selected),'html_figures':len(figures),'pdf_figures':len(pdfs),'input_sha256':{n:sha(src/n) for n in required},'vmr_sha256':sha(vmr),'thresholds':{'mapq20_breadth_1x':breadth,'mapq20_mean_depth':depth},'semantic_checks':semantic_checks,'grading_version':'0.3.0-separate-discovery-priority','discovery_policy':manifest.get('evidence_policy_version','legacy')}
 (out/'validation.json').write_text(json.dumps(validation,indent=2),encoding='utf-8')
 (out/'REPORT.md').write_text('# 病毒候选证据报告\n\n候选数：'+str(len(ids))+'。全部保留，无静默剔除。\n\n分级计数：'+json.dumps(counts)+'\n\nA：多类证据且包含核酸命中、病毒竞争最佳命中、回贴达到阈值；B：多类证据但缺少核酸命中或验证不完整；C：单类或未解决证据；D：细胞来源竞争冲突。分级是启发式复核优先级，不是感染、物种或新病毒确认。\n\n报告汇总图（分级分布、分组构成、覆盖总览）见 report_summary.html；A/B/C 候选分类组成图（宿主来源、基因组组成、植物病毒状态）见 composition_summary.html。图形摘要索引见 index.html（参考比对图与候选覆盖深度图）；PDF 为矢量导出，缩放不失真。覆盖图在有回贴深度数据，或报告表明确记录零深度时生成；后者显示全零曲线并在覆盖图清单标记 ok_zero_depth。\n\nICTV 注释复用传统 accession/name/genus 映射，宿主字段仅描述参考病毒。汇总按推断参考相关名称分组，组数不是病毒物种数。新发现使用区域级竞争，旧结果保留原始判定。请结合 discovery_priority 与来源 manifest 解读；新颖性、宿主和物种仍需独立验证。\n',encoding='utf-8')
 return validation
