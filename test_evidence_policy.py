import unittest
from enhanced_v1.evidence_policy import allowed_model, regional_competition, discovery_priority
from validate_candidates import coverage_metrics

def hit(cls,score,start=1,end=300,sid='ref'):
    return dict(query_id='q',subject_id=sid,classification=cls,bitscore=str(score),qstart=str(start),qend=str(end))

class EvidencePolicyTests(unittest.TestCase):
    def test_cellular_models_rejected(self):
        for acc,name in [('PF24572.2','RBD_RDR6'),('PF24577.2','RDR6_2nd'),('PF25358.1','PH_fung_RdRP')]:
            self.assertFalse(allowed_model(acc,name))
        self.assertTrue(allowed_model('PF00680.25','RdRP_1'))
        self.assertFalse(allowed_model('PF00680.25','RDR6'))
    def test_ties_and_near_ties(self):
        for score in [100,105,110]:
            rows=[hit('viral',100),hit('cellular',score)]
            self.assertEqual(regional_competition(rows)[0]['q'],'ambiguous')
            self.assertEqual(regional_competition(rows[::-1])[0]['q'],'ambiguous')
    def test_same_region_clear_winner(self):
        self.assertEqual(regional_competition([hit('viral',200),hit('cellular',100)])[0]['q'],'viral')
    def test_separate_regions_are_mixed(self):
        self.assertEqual(regional_competition([hit('viral',100),hit('cellular',900,500,800)])[0]['q'],'mixed')
    def test_opposite_strands_do_not_compete(self):
        self.assertEqual(regional_competition([hit('viral',100),hit('cellular',900,300,1)])[0]['q'],'mixed')
    def test_novel_domain_candidate_retained_for_review(self):
        self.assertEqual(discovery_priority('RdRP_domain','no_hit',True,True)[0],'novelty_review_high')
        self.assertEqual(discovery_priority('RdRP_domain','no_hit',True,False)[0],'homology_review')
        self.assertEqual(discovery_priority('RdRP_domain','mixed',True,True)[0],'review_conflict')
    def test_gaps_and_spikes(self):
        result=coverage_metrics(100,{1:1000,2:1,100:1})
        self.assertGreater(result['mean_depth'],5)
        self.assertEqual(result['longest_zero_run'],97)
        self.assertEqual(result['breadth_5x'],0.01)
        self.assertEqual(coverage_metrics(100,{})['longest_zero_run'],100)

if __name__=='__main__':unittest.main()
