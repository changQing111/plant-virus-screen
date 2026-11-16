import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'legacy'))
import visualize_summary as summary
import render_common as rc

class SummaryLayoutTests(unittest.TestCase):
    def test_long_names_remain_outside_bar_column(self):
        fig=rc.Figure(1400,1000)
        name='Long reference-related group '+('unbroken_identifier_'*14)
        groups=[dict(reference_related_group=name,contig_count='10',grade_counts='{"A_supported":10}')]
        bottom=summary._panel_group_stack(fig,groups,355,15,24,1330)
        labels=[op for op in fig.ops if op['op']=='text' and op.get('title')==name]
        self.assertGreater(len(labels),1)
        self.assertTrue(all(op['x']+rc.text_width(op['s'],op['size'])<=416 for op in labels))
        self.assertGreater(bottom,max(op['y'] for op in labels))

    def test_dense_top_points_have_tooltips_and_separate_details(self):
        rows=[dict(contig_id='candidate_'+str(i),evidence_grade='A_supported',
                   mapq20_mean_depth='100',mapq20_breadth_1x='1') for i in range(30)]
        fig=rc.Figure(1400,1000)
        bottom=summary._panel_coverage_scatter(fig,rows,900,1260,24,1330,0.8,5)
        points=[op for op in fig.ops if op['op']=='circle' and op.get('title')]
        labels=[op for op in fig.ops if op['op']=='text' and op['s'].startswith('candidate_')]
        self.assertEqual(len(points),30)
        self.assertEqual(len(labels),30)
        self.assertTrue(all(op['y']>1260 for op in labels))
        self.assertGreater(bottom,max(op['y'] for op in labels))
        self.assertTrue(all(970<op['cy']<1260 for op in points))

    def test_empty_scatter_returns_layout_boundary(self):
        fig=rc.Figure(1400,1000)
        self.assertGreater(summary._panel_coverage_scatter(fig,[],900,1260,24,1330,0.8,5),900)

if __name__=='__main__':unittest.main()
