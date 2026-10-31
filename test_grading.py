import unittest
from report import grade
class GradingTests(unittest.TestCase):
 def test_conflict_overrides_support(self):self.assertEqual(grade('nt;protein','cellular',True)[0],'D_conflicting')
 def test_mapping_required_for_A(self):self.assertEqual(grade('nt;protein','viral',False)[0],'B_supported')
 def test_A(self):self.assertEqual(grade('nt;protein','viral',True)[0],'A_supported')
 def test_domain_only(self):self.assertEqual(grade('RdRP_domain','no_hit',False)[0],'C_provisional')
 def test_unclassified(self):self.assertEqual(grade('nt;protein','unclassified',True)[0],'B_supported')
if __name__=='__main__':unittest.main()
