import csv
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from report import grade, legacy, validate_discovery_semantics
import batch_visualize_coverage as coverage_batch


class GradingTests(unittest.TestCase):
 def test_conflict_overrides_support(self):self.assertEqual(grade('nt;protein','cellular',True)[0],'D_conflicting')
 def test_mapping_required_for_A(self):self.assertEqual(grade('nt;protein','viral',False)[0],'B_supported')
 def test_A(self):self.assertEqual(grade('nt;protein','viral',True)[0],'A_supported')
 def test_domain_only(self):self.assertEqual(grade('RdRP_domain','no_hit',False)[0],'C_provisional')
 def test_unclassified(self):self.assertEqual(grade('nt;protein','unclassified',True)[0],'B_supported')


class VmrAccessionTests(unittest.TestCase):
 def test_uniprot_fragments_are_not_accessions(self):
  text='tr|A0A8S5TBU1|A0A8S5TBU1_9CAUD 30S ribosomal protein S4'
  self.assertEqual(legacy.extract_accessions(text),[])

 def test_supported_accessions_preserve_order(self):
  text='NC_001477.1; L: AJ410618; RNA3: AB771723; NC_001477.1'
  self.assertEqual(legacy.extract_accessions(text),['NC_001477.1','NC_001477','AJ410618','AB771723'])

 def test_accession_join_uses_nucleotide_hit_only(self):
  accession={'AB123456':{'dsRNA'}}
  accession_host={'AB123456':{'plants'}}
  maps=(accession,{},{},{},accession_host,{},{},{})
  rows=[{'best_nt_hit':'','best_protein_hit':'sp|P12345| protein AB123456','putative_virus':'Unknown'}]
  with mock.patch.object(legacy,'load_vmr_genome_maps',return_value=maps):
   legacy.annotate_with_vmr(rows,'unused.xlsx')
  self.assertEqual(rows[0]['genome_composition'],'')
  self.assertEqual(rows[0]['host_source'],'')

 def test_nucleotide_accession_join_still_works(self):
  accession={'AB123456':{'dsRNA'}}
  accession_host={'AB123456':{'plants'}}
  maps=(accession,{},{},{},accession_host,{},{},{})
  rows=[{'best_nt_hit':'AB123456 | example virus','best_protein_hit':'','putative_virus':'Unknown'}]
  with mock.patch.object(legacy,'load_vmr_genome_maps',return_value=maps):
   legacy.annotate_with_vmr(rows,'unused.xlsx')
  self.assertEqual(rows[0]['genome_composition'],'dsRNA')
  self.assertEqual(rows[0]['plant_virus_status'],'plant virus')


class DiscoverySemanticTests(unittest.TestCase):
 def write_tsv(self,path,fields,rows):
  with Path(path).open('w',encoding='utf-8',newline='') as handle:
   writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t')
   writer.writeheader();writer.writerows(rows)

 def fixture(self,root,evidence='nt;protein',competition='viral',competitive_evalue='1e-20'):
  root=Path(root)
  nt=['q1','AB123456','1000','1000','90','400','40','40','40','0','0','1','400','1','400','1e-20','300']
  aa=['q1','P1','1000','300','50','120','36','0','0','1','360','1','120','1e-20','100']
  (root/'blastn.tsv').write_text('\t'.join(nt)+'\n',encoding='utf-8')
  (root/'blastx.tsv').write_text('\t'.join(aa)+'\n',encoding='utf-8')
  self.write_tsv(root/'domain_evidence.tsv',['query_id','orf_id','model','accession','i_evalue','hmm_coverage','aa_start','aa_end'],[])
  comp={'query_id':'q1','subject_id':'P1','bitscore':'100','evalue':competitive_evalue,'qstart':'1','qend':'360','classification':'viral','annotation':'','taxid':'1','lineage':'Viruses','scientific_name':'Example virus','mapping_status':'matched','taxonomy_source':'test'}
  self.write_tsv(root/'competitive_annotated.tsv',list(comp),[comp])
  row={'query_id':'q1','length':'1000','retained':'1','evidence':evidence,'competitive_best_class':competition,'mapping_status':'completed','host_status':'undetermined','interpretation':'candidate_requires_validation'}
  return [dict(row)],[dict(row)]

 def test_consistent_files_pass(self):
  with tempfile.TemporaryDirectory() as tmp:
   rows,audit=self.fixture(tmp)
   best,checks=validate_discovery_semantics(Path(tmp),rows,audit)
  self.assertEqual(set(best),{('q1','nt'),('q1','protein')})
  self.assertTrue(checks['evidence_consistent'])
  self.assertTrue(checks['competition_consistent'])

 def test_evidence_mismatch_fails(self):
  with tempfile.TemporaryDirectory() as tmp:
   rows,audit=self.fixture(tmp,evidence='nt')
   with self.assertRaisesRegex(ValueError,'evidence mismatch'):
    validate_discovery_semantics(Path(tmp),rows,audit)

 def test_competition_mismatch_fails(self):
  with tempfile.TemporaryDirectory() as tmp:
   rows,audit=self.fixture(tmp,competition='cellular')
   with self.assertRaisesRegex(ValueError,'competition mismatch'):
    validate_discovery_semantics(Path(tmp),rows,audit)

 def test_candidate_must_match_audit_row(self):
  with tempfile.TemporaryDirectory() as tmp:
   rows,audit=self.fixture(tmp)
   rows[0]['mapping_status']='not_requested'
   with self.assertRaisesRegex(ValueError,'Candidate/audit row mismatch'):
    validate_discovery_semantics(Path(tmp),rows,audit)

 def test_invalid_domain_threshold_fails(self):
  with tempfile.TemporaryDirectory() as tmp:
   rows,audit=self.fixture(tmp)
   domain={'query_id':'q1','orf_id':'orf1','model':'RdRP','accession':'PF00000','i_evalue':'1','hmm_coverage':'0.1','aa_start':'1','aa_end':'100'}
   self.write_tsv(Path(tmp)/'domain_evidence.tsv',list(domain),[domain])
   with self.assertRaisesRegex(ValueError,'domain_evidence.*threshold'):
    validate_discovery_semantics(Path(tmp),rows,audit)

 def test_invalid_competitive_evalue_fails(self):
  with tempfile.TemporaryDirectory() as tmp:
   rows,audit=self.fixture(tmp,competitive_evalue='1')
   with self.assertRaisesRegex(ValueError,'competitive_annotated.*threshold'):
    validate_discovery_semantics(Path(tmp),rows,audit)

 def test_nonretained_contig_with_qualifying_hit_fails(self):
  with tempfile.TemporaryDirectory() as tmp:
   rows,audit=self.fixture(tmp)
   extra=dict(audit[0]);extra['query_id']='q2';extra['retained']='0';extra['evidence']=''
   audit.append(extra)
   nt=['q2','AB999999','1000','1000','90','400','40','40','40','0','0','1','400','1','400','1e-20','300']
   with (Path(tmp)/'blastn.tsv').open('a',encoding='utf-8') as handle:
    handle.write('\t'.join(nt)+'\n')
   with self.assertRaisesRegex(ValueError,'Audit evidence mismatch'):
    validate_discovery_semantics(Path(tmp),rows,audit)


class GradeRequiresNucleotideTests(unittest.TestCase):
 def test_protein_plus_rdrp_domain_is_not_A(self):
  # BLASTX and the RdRp HMM usually hit the same polymerase region.
  self.assertEqual(grade('protein;RdRP_domain','viral',True)[0],'B_supported')

 def test_nt_plus_protein_is_still_A(self):
  self.assertEqual(grade('nt;protein','viral',True)[0],'A_supported')

 def test_nt_plus_domain_without_mapping_is_B(self):
  self.assertEqual(grade('nt;RdRP_domain','viral',False)[0],'B_supported')

 def test_conflict_still_overrides_new_rule(self):
  self.assertEqual(grade('nt;protein','cellular',True)[0],'D_conflicting')


class VmrNameCleaningTests(unittest.TestCase):
 def test_complete_genome_suffix(self):
  text='NC_009087.2 | Chrysanthemum virus B, complete genome'
  self.assertEqual(legacy.clean_reference_name(text),'Chrysanthemum virus B')

 def test_isolate_and_complete_cds(self):
  text='NC_076803.1 | Plasmopara viticola lesion associated ourmia-like virus 74 isolate DMG-C_28482 RNA dependent RNA polymerase (QKR79_gp1) gene, complete cds'
  self.assertEqual(legacy.clean_reference_name(text),'Plasmopara viticola lesion associated ourmia-like virus 74')

 def test_segment_isolate_combination(self):
  text='NC_017989.1 | Persimmon cryptic virus segment 1, isolate SSPI, complete genome'
  self.assertEqual(legacy.clean_reference_name(text),'Persimmon cryptic virus')

 def test_strain_and_segment(self):
  text='NC_028385.1 | Wuhan aphid virus 1 strain WHYC-1 segment 4, complete sequence'
  self.assertEqual(legacy.clean_reference_name(text),'Wuhan aphid virus 1')

 def test_parenthesized_strain(self):
  self.assertEqual(legacy.clean_reference_name('Tobacco mosaic virus (strain O)'),'Tobacco mosaic virus')

 def test_word_isolate_in_real_name_is_not_truncated(self):
  self.assertEqual(legacy.clean_reference_name('Grapevine isolate virus'),'Grapevine isolate virus')

 def test_uniprot_fragment_stays_unresolved(self):
  text='tr|A0A8S5TBU1|A0A8S5TBU1_9CAUD 30S ribosomal protein S4'
  self.assertEqual(legacy.extract_accessions(text),[])
  maps=({},{},{},{},{},{},{},{})
  with mock.patch.object(legacy,'load_vmr_genome_maps',return_value=maps):
   rows=[{'best_nt_hit':text,'best_protein_hit':'','putative_virus':'Unknown'}]
   legacy.annotate_with_vmr(rows,'unused.xlsx')
  self.assertEqual(rows[0]['genome_composition'],'')

 def test_cleaned_nt_name_resolves_vmr(self):
  key=legacy.normalize_name('Chrysanthemum virus B')
  maps=({},{key:{'ssRNA(+)'}},{},{},{},{key:{'plants'}},{},{})
  with mock.patch.object(legacy,'load_vmr_genome_maps',return_value=maps):
   rows=[{'best_nt_hit':'NC_009087.2 | Chrysanthemum virus B, complete genome','best_protein_hit':'','putative_virus':'Unknown'}]
   legacy.annotate_with_vmr(rows,'unused.xlsx')
  self.assertEqual(rows[0]['genome_composition'],'ssRNA(+)')
  self.assertEqual(rows[0]['baltimore_group'],'IV')
  self.assertEqual(rows[0]['plant_virus_status'],'plant virus')

 def test_ambiguous_cleaned_name_is_not_used(self):
  key=legacy.normalize_name('Wuhan aphid virus 1')
  maps=({},{key:{'ssRNA(+)','dsRNA'}},{},{},{},{},{},{})
  with mock.patch.object(legacy,'load_vmr_genome_maps',return_value=maps):
   rows=[{'best_nt_hit':'NC_028385.1 | Wuhan aphid virus 1 strain WHYC-1 segment 4, complete sequence','best_protein_hit':'','putative_virus':'Unknown'}]
   legacy.annotate_with_vmr(rows,'unused.xlsx')
  self.assertEqual(rows[0]['genome_composition'],'')


class GenomeTypeMappingTests(unittest.TestCase):
 def test_ambisense_ssdna_type_and_baltimore(self):
  for value in ('ssDNA(+/-)','ssDNA(-/+)'):
   self.assertEqual(legacy.genome_type_cn(value),'单链DNA病毒')
   self.assertEqual(legacy.baltimore_group(value),'II')

 def test_existing_mappings_unchanged(self):
  self.assertEqual(legacy.baltimore_group('dsRNA'),'III')
  self.assertEqual(legacy.baltimore_group('dsDNA-RT'),'VII')
  self.assertEqual(legacy.baltimore_group('unknown genome'),'')


class CandidateLabelTests(unittest.TestCase):
 def test_nt_only_is_nucleotide_supported(self):
  row={'evidence':'nt','best_nt_hit':'NC_007346.1 | Emiliania huxleyi virus 86, complete genome','best_protein_hit':'','blastx_anno':''}
  self.assertEqual(legacy.infer_segment_or_protein(row),'Nucleotide-supported candidate')

 def test_rdrp_domain_only_is_labelled_correctly(self):
  row={'evidence':'RdRP_domain','best_nt_hit':'','best_protein_hit':'','blastx_anno':''}
  self.assertEqual(legacy.infer_segment_or_protein(row),'RdRP-domain-only candidate')

 def test_protein_keyword_takes_priority(self):
  row={'evidence':'protein;RdRP_domain','best_protein_hit':'tr|X|Y','blastx_anno':'RNA-directed RNA polymerase OS=Example virus OX=1'}
  self.assertEqual(legacy.infer_segment_or_protein(row),'RdRP/replicase')

 def test_protein_hit_without_keyword(self):
  row={'evidence':'protein','best_protein_hit':'tr|X|Y','blastx_anno':'Uncharacterized protein OS=Example virus OX=1'}
  self.assertEqual(legacy.infer_segment_or_protein(row),'Best protein-supported candidate')

 def test_legacy_row_without_evidence_keeps_old_fallback(self):
  self.assertEqual(legacy.infer_segment_or_protein({'blastx_anno':'','blastn_hit':''}),'Best protein-supported candidate')

 def test_legacy_nt_hit_wins_over_non_keyword_protein_description(self):
  row={'blastx_anno':'hypothetical protein','blastn_hit':'NC_000000','best_protein_hit':''}
  self.assertEqual(legacy.infer_segment_or_protein(row),'Best nucleotide-supported candidate')

 def test_empty_evidence_is_unresolved(self):
  self.assertEqual(legacy.infer_segment_or_protein({'evidence':'','best_nt_hit':'','blastx_anno':''}),'Unresolved viral feature')


class ReferenceConcordanceTests(unittest.TestCase):
 def test_concordant_strain_variant(self):
  self.assertEqual(legacy.reference_name_concordance('Tobacco mosaic virus','Tobacco mosaic virus (strain O)'),'concordant')

 def test_discordant_references(self):
  self.assertEqual(legacy.reference_name_concordance('Chrysanthemum virus B','Chrysanthemum indicum carlavirus'),'discordant')
  self.assertEqual(legacy.reference_name_concordance('Petunia vein clearing virus','Aristotelia chilensis virus 1'),'discordant')

 def test_single_and_missing_references(self):
  self.assertEqual(legacy.reference_name_concordance('','Tobacco mosaic virus'),'single_reference')
  self.assertEqual(legacy.reference_name_concordance('',''),'unresolved')


class ZeroDepthCoverageTests(unittest.TestCase):
 def test_explicit_zero_creates_full_zero_series(self):
  values=coverage_batch.explicit_zero_depth_series({'length':'5000','mapq20_mean_depth':'0.0','mapq20_breadth_1x':'0.0'})
  self.assertEqual(len(values),5000)
  self.assertTrue(all(value==0 for value in values))

 def test_nonzero_or_blank_never_synthesizes(self):
  self.assertIsNone(coverage_batch.explicit_zero_depth_series({'length':'5000','mapq20_mean_depth':'12.0','mapq20_breadth_1x':'0.9'}))
  self.assertIsNone(coverage_batch.explicit_zero_depth_series({'length':'5000','mapq20_mean_depth':'','mapq20_breadth_1x':''}))
  self.assertIsNone(coverage_batch.explicit_zero_depth_series({'length':'0','mapq20_mean_depth':'0.0','mapq20_breadth_1x':'0.0'}))
  self.assertIsNone(coverage_batch.explicit_zero_depth_series({'length':'nan','mapq20_mean_depth':'0.0','mapq20_breadth_1x':'0.0'}))

 def test_zero_depth_contig_gets_figure_and_manifest_status(self):
  root=Path(__file__).resolve().parent
  with tempfile.TemporaryDirectory() as tmp:
   tmp=Path(tmp)
   fields=['query_id','contig_id','length','mapq20_mean_depth','mapq20_breadth_1x','evidence_grade','putative_virus','blastn_anno','blastx_anno']
   rows=[
    {'query_id':'zero_contig','contig_id':'zero_contig','length':'5000','mapq20_mean_depth':'0.0','mapq20_breadth_1x':'0.0','evidence_grade':'B_supported','putative_virus':'Tobacco mosaic virus','blastn_anno':'','blastx_anno':''},
    {'query_id':'explicit_rows_contig','contig_id':'explicit_rows_contig','length':'200','mapq20_mean_depth':'0.0','mapq20_breadth_1x':'0.0','evidence_grade':'B_supported','putative_virus':'Example virus','blastn_anno':'','blastx_anno':''},
    {'query_id':'missing_contig','contig_id':'missing_contig','length':'3000','mapq20_mean_depth':'8.0','mapq20_breadth_1x':'0.9','evidence_grade':'A_supported','putative_virus':'Example virus','blastn_anno':'','blastx_anno':''},
   ]
   with (tmp/'plot_candidates.tsv').open('w',encoding='utf-8',newline='') as handle:
    writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t')
    writer.writeheader();writer.writerows(rows)
   (tmp/'depth_mapq20.tsv').write_text(''.join('explicit_rows_contig\t%d\t0\n'%i for i in range(1,201)),encoding='utf-8')
   argv=['batch_visualize_coverage.py','--final-tsv',str(tmp/'plot_candidates.tsv'),'--depth',str(tmp/'depth_mapq20.tsv'),'--module',str(root/'legacy'/'visualize_coverage.py'),'--outdir',str(tmp/'coverage_plots'),'--no-pdf','--bins','10']
   with mock.patch('sys.argv',argv):
    coverage_batch.main()
   with (tmp/'coverage_plots'/'coverage_graphic_summary_manifest.tsv').open(encoding='utf-8') as handle:
    manifest=list(csv.DictReader(handle,delimiter='\t'))
   status={r['contig_id']:r['status'] for r in manifest}
   self.assertEqual(status,{'zero_contig':'ok_zero_depth','explicit_rows_contig':'ok_zero_depth','missing_contig':'skipped_no_depth'})
   html=(tmp/'coverage_plots'/'zero_contig.coverage_summary.html').read_text(encoding='utf-8')
   self.assertIn('zero depth',html)


if __name__=='__main__':unittest.main()
