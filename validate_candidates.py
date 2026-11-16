#!/usr/bin/env python3
"""Optional assembly diagnostics and homolog-specific phylogeny. No species calls."""
import argparse, csv, hashlib, json, shutil, subprocess
from pathlib import Path
from collections import defaultdict

def digest(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()

def write(path, rows, fields):
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)

def coverage_metrics(length, positions):
    if length <= 0: raise ValueError('Invalid sequence length')
    total=n1=n5=0; last=0; gap=0
    for pos,depth in sorted(positions.items()):
        if not 1<=pos<=length or depth<0: raise ValueError('Invalid depth coordinate/value')
        total+=depth;n1+=depth>0;n5+=depth>=5
        if depth:
            gap=max(gap,pos-last-1);last=pos
    return dict(mean_depth=total/length,breadth_1x=n1/length,breadth_5x=n5/length,
                longest_zero_run=max(gap,length-last))

def assembly(a,out):
    with open(a.discovery/'final_candidates.tsv') as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    lengths={r['query_id']:int(r['length']) for r in rows}
    depths=defaultdict(dict)
    # Recompute to avoid mixing depth conventions from historical runs.
    bam=a.discovery/'mapping.bam'
    if not bam.is_file(): raise ValueError('Assembly diagnostics require mapping.bam')
    bed=out/'candidates.bed'
    bed.write_text(''.join(f'{q}\t0\t{n}\n' for q,n in lengths.items()))
    cmd=[a.samtools,'depth','-s','-a','-Q','20','-b',str(bed),str(bam)]
    with open(out/'depth.tsv','w') as f: subprocess.run(cmd,stdout=f,check=True)
    with open(out/'depth.tsv') as f:
        for line in f:
            q,p,d=line.split()
            if q in lengths:depths[q][int(p)]=int(d)
    counts=defaultdict(lambda:defaultdict(int))
    import re
    with subprocess.Popen([a.samtools,'view','-F','2308','-q','20',str(bam)],stdout=subprocess.PIPE,text=True) as proc:
        for line in proc.stdout:
            r=line.split('\t');q=r[2]
            if q not in lengths:continue
            flag=int(r[1]);c=counts[q];c['primary_alignments']+=1
            c['paired_alignments']+=bool(flag&1)
            c['proper_pair_alignments']+=bool(flag&2)
            c['softclipped_alignments']+=any(int(n)>=10 for n in re.findall(r'(\d+)S',r[5]))
        if proc.wait():raise RuntimeError('samtools view failed')
    results=[]
    for q,n in lengths.items():
        x=dict(query_id=q,length=n,**coverage_metrics(n,depths[q]));c=counts[q]
        x['primary_alignments']=c['primary_alignments']
        x['softclip_fraction']=c['softclipped_alignments']/max(1,c['primary_alignments'])
        x['proper_pair_fraction']=c['proper_pair_alignments']/c['paired_alignments'] if c['paired_alignments'] else ''
        x['assembly_status']='review' if x['breadth_5x']<0.8 or x['longest_zero_run']>=100 or x['softclip_fraction']>0.1 else 'coverage_screen_passed'
        x['limitation']='not_junction_validation;not_infection_confirmation'
        results.append(x)
    fields=list(results[0]) if results else ['query_id','assembly_status']
    write(out/'assembly_validation.tsv',results,fields)
    return dict(bam=str(bam),depth_command=cmd,candidates=len(results))

def phylogeny(a,out):
    from Bio import SeqIO
    seqs=list(SeqIO.parse(a.proteins,'fasta'))
    ids=[s.id for s in seqs]
    if len(seqs)<4 or len(set(ids))!=len(ids):raise ValueError('Need >=4 homologous proteins with unique IDs')
    if any(not s.seq or '*' in str(s.seq) for s in seqs):raise ValueError('Empty or stop-containing protein')
    wanted=set(a.candidate_ids.read_text().split())
    if not wanted or not wanted<set(ids):raise ValueError('Candidate IDs must be present, alongside references')
    if not a.marker.strip():raise ValueError('Specify homologous marker name')
    for tool in (a.mafft,a.iqtree):
        if not shutil.which(tool):raise FileNotFoundError(tool)
    aln=out/'aligned.faa'
    with open(aln,'w') as f,open(out/'mafft.log','w') as log:
        subprocess.run([a.mafft,'--thread',str(a.threads),'--auto',str(a.proteins.resolve())],stdout=f,stderr=log,check=True)
    aligned=list(SeqIO.parse(aln,'fasta'));size=len(aligned[0])
    keep=[i for i in range(size) if sum(str(s.seq)[i] not in '-X?' for s in aligned)/len(aligned)>=a.occupancy]
    if len(keep)<80:raise ValueError('Fewer than 80 columns remain; review homolog selection')
    trimmed=out/'trimmed.faa'
    trimmed.write_text(''.join('>'+s.id+'\n'+''.join(str(s.seq)[i] for i in keep)+'\n' for s in aligned))
    cmd=[a.iqtree,'-s',str(trimmed),'-st','AA','-m','MFP','-B','1000','--alrt','1000','-T',str(a.threads),'-seed','42','--prefix',str(out/'tree')]
    with open(out/'iqtree.log','w') as log:subprocess.run(cmd,stdout=log,stderr=log,check=True)
    if not (out/'tree.treefile').is_file():raise RuntimeError('Missing inferred tree')
    return dict(marker=a.marker,protein_sha256=digest(a.proteins),candidate_ids=sorted(wanted),
                retained_columns=len(keep),occupancy=a.occupancy,command=cmd,
                interpretation='unrooted_gene_tree;user_curated_homologs;not_species_demarcation')

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='mode',required=True)
    s=sub.add_parser('assembly');s.add_argument('--discovery',required=True,type=Path);s.add_argument('--samtools',default='samtools')
    t=sub.add_parser('phylogeny');t.add_argument('--proteins',required=True,type=Path);t.add_argument('--candidate-ids',required=True,type=Path)
    t.add_argument('--marker',required=True);t.add_argument('--mafft',default='mafft');t.add_argument('--iqtree',default='iqtree2')
    t.add_argument('--threads',type=int,default=4);t.add_argument('--occupancy',type=float,default=0.5)
    for parser in (s,t):parser.add_argument('--out',required=True,type=Path)
    a=p.parse_args()
    if a.mode=='phylogeny' and (a.threads<1 or not 0<a.occupancy<=1):p.error('Invalid threads/occupancy')
    a.out=a.out.resolve();a.out.mkdir(parents=True,exist_ok=False)
    manifest={'status':'running','mode':a.mode}
    try:
        manifest.update(assembly(a,a.out) if a.mode=='assembly' else phylogeny(a,a.out));manifest['status']='completed'
    except BaseException as exc:
        manifest.update(status='failed',error=repr(exc));raise
    finally:(a.out/'manifest.json').write_text(json.dumps(manifest,indent=2))

if __name__=='__main__':main()
