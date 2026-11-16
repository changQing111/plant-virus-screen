import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "legacy"))
import render_common as rc
import visualize_blastn_reference as vbn
import visualize_blastx_protein as vbx


class AlignmentPlotsLayoutTests(unittest.TestCase):
    def test_pident_color_ramp_is_monotonic_in_luminance(self):
        luminances = []
        for pident in (60.0, 70.0, 80.0, 90.0, 100.0):
            hex_col = rc.pident_to_color(pident, 60.0, 100.0)
            r, g, b = rc.hex_to_rgb01(hex_col)
            lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
            luminances.append(lum)
        for idx in range(len(luminances) - 1):
            self.assertGreater(
                luminances[idx],
                luminances[idx + 1],
                "Sequential identity ramp should darken monotonically for grayscale printability",
            )

    def test_blastn_long_ids_and_dense_hsps_do_not_cross_columns_or_min_font(self):
        long_id = "spades_NODE_99999_length_8989_cov_6230.698183_g39_i2_extra_long_identifier_segment"
        rows = [
            [long_id, "NC_009087.2", "8989", "8855", "88.5", "4800", "90.0", "54.0", "54.0", "100", "2", "37", "4840", "4040", "8843", "0.0", "4990.0"],
            [long_id, "NC_009087.2", "8989", "8855", "79.2", "1040", "90.0", "12.0", "12.0", "50", "1", "7873", "8915", "1", "1043", "0.0", "684.0"],
            [long_id, "NC_009087.2", "8989", "8855", "82.1", "900", "90.0", "10.0", "10.0", "40", "1", "7900", "8800", "50", "950", "0.0", "610.0"],
        ]
        hits = [vbn.BlastnHit(r) for r in rows]
        grouped = vbn.group_hits_by_contig(hits)
        doc = vbn.build_html("NC_009087.2", "NC_009087.2 | Chrysanthemum virus B, complete genome", grouped)
        fig = doc.figure
        text_ops = [op for op in fig.ops if op["op"] == "text"]
        # All PDF text >= 5 pt
        self.assertTrue(all(op["size"] * rc.PT_PER_PX >= 5.0 for op in text_ops))
        # Left column row labels stay strictly left of track x0 = 420
        section_titles = {
            "a",
            "b",
            "Reference genome & aligned HSP intervals",
            "Assembled contig alignments (ordered by reference position & bitscore)",
        }
        left_ops = [op for op in text_ops if op["x"] < 400 and op["y"] > 120 and op["s"] not in section_titles]
        self.assertGreater(len(left_ops), 0)
        for op in left_ops:
            self.assertLessEqual(op["x"] + rc.text_width(op["s"], op["size"], op["bold"]), 400.0)

    def test_blastx_two_panel_layout_keeps_columns_separated(self):
        long_id = "spades_NODE_433_length_6366_cov_1033.720960_g44_i3_unbroken_suffix_identifier"
        rows = [
            [long_id, "A9CMG0", "6366", "2094", "91.5", "2094", "98.0", "40", "0", "75", "6356", "1", "2094", "0.0", "3706.0"],
            [long_id, "A9CMG0", "6366", "2094", "86.2", "250", "12.0", "20", "0", "481", "1230", "1", "250", "0.0", "487.0"],
        ]
        hits = [vbx.BlastxHit(r) for r in rows]
        grouped = vbx.group_hits_by_contig(hits)
        doc = vbx.build_html("A9CMG0", "tr|A9CMG0|A9CMG0_CVB | Replicase OS=Chrysanthemum virus B", grouped)
        fig = doc.figure
        text_ops = [op for op in fig.ops if op["op"] == "text"]
        self.assertTrue(all(op["size"] * rc.PT_PER_PX >= 5.0 for op in text_ops))
        left_ops = [op for op in text_ops if op["x"] < 380 and op["y"] > 120]
        self.assertGreater(len(left_ops), 0)
        for op in left_ops:
            self.assertLessEqual(op["x"] + rc.text_width(op["s"], op["size"], op["bold"]), 380.0)


if __name__ == "__main__":
    unittest.main()
