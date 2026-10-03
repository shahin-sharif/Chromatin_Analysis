"""Independent interval/bin oracles and CLI regression cases."""

import contextlib
import io
import json
from pathlib import Path
import random
import sys
import tempfile
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import ChromHMMTools as tool


class ChromTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.wt = self.write("wt.bed", "chr1\t0\t400\tE1\nchr1\t400\t800\tE2\n")
        self.mt = self.write("mt.bed", "chr1\t0\t200\tE1\nchr1\t200\t800\tE2\n")
        self.out = self.root / "out"

    def write(self, name, text):
        p = self.root / name
        p.write_text(text)
        return p

    def run_tool(self, *options, command="run", expected=0):
        args = [command, "--outdir", str(self.out)]
        if command != "report":
            args += [
                "--wt",
                str(self.wt),
                "--mt",
                str(self.mt),
                "--state-mode",
                "shared",
            ]
        else:
            args += ["--indir", str(self.out)]
        if command != "compare":
            args += ["--no-plots"]
        with contextlib.redirect_stderr(io.StringIO()):
            result = tool.main(args + list(map(str, options)))
        self.assertEqual(result, expected)

    def test_known_counts_and_report(self):
        self.run_tool("--write-change-bed")
        c, _, _ = tool.read_matrix(self.out / "overlap.bp.tsv")
        np.testing.assert_array_equal(c, [[200, 200], [0, 400]])
        c, _, _ = tool.read_matrix(self.out / "WT.trans.counts.tsv")
        np.testing.assert_array_equal(c, [[1, 1], [0, 1]])
        self.assertIn("75.00%", (self.out / "report.html").read_text())
        self.assertEqual(
            (self.out / "state_changes.bed").read_text(),
            "chr1\t200\t400\tE1->E2\t0\t.\n",
        )

    def test_gap_does_not_create_transition(self):
        c, bins = tool.transitions(
            {"x": [(0, 400, "a"), (600, 1000, "b")]}, ["a", "b"], 200
        )
        np.testing.assert_array_equal(c, [[1, 0], [0, 1]])
        self.assertEqual(bins, 4)

    def test_partial_bins_excluded(self):
        c, bins = tool.transitions(
            {"x": [(1, 399, "a"), (399, 801, "b")]}, ["a", "b"], 200
        )
        self.assertEqual(bins, 2)
        np.testing.assert_array_equal(c, [[0, 0], [0, 1]])

    def test_random_bin_oracle(self):
        rng = random.Random(91)
        for _ in range(100):
            labels = [rng.choice(["a", "b", None]) for _ in range(100)]
            rows = []
            for i, label in enumerate(labels):
                if label is None:
                    continue
                if rows and rows[-1][1] == i and rows[-1][2] == label:
                    rows[-1] = (rows[-1][0], i + 1, label)
                else:
                    rows.append((i, i + 1, label))
            bins = [
                (
                    labels[i]
                    if labels[i] is not None and len(set(labels[i : i + 5])) == 1
                    else None
                )
                for i in range(0, 100, 5)
            ]
            expected = np.zeros((2, 2), dtype=int)
            for a, b in zip(bins, bins[1:]):
                if a and b:
                    expected["ab".index(a), "ab".index(b)] += 1
            actual, used = tool.transitions({"x": rows}, ["a", "b"], 5)
            np.testing.assert_array_equal(actual, expected)
            self.assertEqual(used, sum(x is not None for x in bins))

    def test_normalization_zero_support(self):
        p = tool.normalize(np.array([[0, 0], [1, 3]]))
        self.assertTrue(np.isnan(p[0]).all())
        np.testing.assert_allclose(p[1], [0.25, 0.75])

    def test_annotation_partition_conserves_bp(self):
        anno = self.write(
            "anno.bed", "chr1\t100\t600\tgenic\nchr1\t300\t500\tpromoter\n"
        )
        self.run_tool("--anno-bed", anno)
        info = json.loads((self.out / "RUNINFO.json").read_text())
        sums = {
            key: int(tool.read_matrix(self.out / value)[0].sum())
            for key, value in info["annotation_files"].items()
        }
        self.assertEqual(sums, {"genic": 300, "intergenic": 300, "promoter": 200})

    def test_gtf_strands_and_zero_based(self):
        gtf = self.write(
            "g.gtf",
            'chr1\tx\tgene\t101\t200\t.\t+\t.\tgene_id "g";\nchr2\tx\tgene\t101\t200\t.\t-\t.\tgene_id "h";\n',
        )
        args = tool.parser().parse_args(
            [
                "compare",
                "--wt",
                str(self.wt),
                "--mt",
                str(self.mt),
                "--state-mode",
                "shared",
                "--outdir",
                str(self.out),
                "--gtf",
                str(gtf),
                "--promoter-up",
                "10",
                "--promoter-down",
                "2",
            ]
        )
        a = tool.annotations(args)
        self.assertIn((90, 103, "promoter"), a["chr1"])
        self.assertIn((197, 210, "promoter"), a["chr2"])
        self.assertIn((100, 200, "genic"), a["chr1"])

    def test_empty_include_rejects(self):
        inc = self.write("empty.bed", "")
        self.run_tool("--include-bed", inc, expected=2)
        self.assertFalse((self.out / "RUNINFO.json").exists())

    def test_masks_exact_intervals(self):
        data = {"x": [(0, 100, "a")]}
        self.assertEqual(
            tool.clip_regions(
                data, {"x": [(10, 90)]}, {"x": [(20, 30), (25, 40), (50, 60)]}
            ),
            {"x": [(10, 20, "a"), (40, 50, "a"), (60, 90, "a")]},
        )

    def test_distinct_rectangular_no_conservation(self):
        self.mt.write_text("chr1\t0\t800\tX\n")
        self.run_tool("--state-mode", "distinct")
        c, r, col = tool.read_matrix(self.out / "overlap.bp.tsv")
        self.assertEqual(c.shape, (2, 1))
        self.assertEqual(col, ["X"])
        self.assertNotIn(
            "Matching-label retention:", (self.out / "report.html").read_text()
        )

    def test_explicit_mapping(self):
        self.mt.write_text("chr1\t0\t400\tX\nchr1\t400\t800\tY\n")
        mapping = self.write("map.tsv", "X\tE1\nY\tE2\n")
        self.run_tool("--state-mode", "mapped", "--state-map", mapping)
        self.assertIn("100.00%", (self.out / "report.html").read_text())

    def test_nonbijective_mapping_rejected(self):
        mapping = self.write("map.tsv", "E1\tE1\nE2\tE1\n")
        self.run_tool("--state-mode", "mapped", "--state-map", mapping, expected=2)

    def test_overlapping_input_rejected(self):
        self.wt.write_text("chr1\t0\t500\tE1\nchr1\t400\t800\tE2\n")
        self.run_tool(expected=2)

    def test_disjoint_coverage_rejected(self):
        self.mt.write_text("chr2\t0\t800\tE1\n")
        self.run_tool(expected=2)

    def test_report_same_dir_force(self):
        self.run_tool()
        self.run_tool("--force", command="report")
        self.assertTrue((self.out / "report.html").is_file())

    def test_no_overwrite_and_inputs_unchanged(self):
        self.run_tool()
        before = (self.out / "RUNINFO.json").read_bytes()
        self.run_tool(expected=2)
        self.assertEqual((self.out / "RUNINFO.json").read_bytes(), before)
        self.assertEqual(self.wt.read_text(), "chr1\t0\t400\tE1\nchr1\t400\t800\tE2\n")

    def test_zero_support_state_and_no_zero_ranked_pairs(self):
        self.wt.write_text(self.wt.read_text() + "chr2\t0\t200\tE3\n")
        self.run_tool()
        self.assertIn("E3\t0\tNA\tNA\tNA", (self.out / "state_summary.tsv").read_text())
        self.assertEqual(len((self.out / "top_pairs.tsv").read_text().splitlines()), 2)

    def test_label_aligned_conservation(self):
        self.run_tool(command="compare")
        tool.write_matrix(
            self.out / "overlap.bp.tsv",
            np.array([[200, 200], [400, 0]]),
            ["E1", "E2"],
            ["E2", "E1"],
        )
        self.run_tool(command="report")
        self.assertIn("75.00%", (self.out / "report.html").read_text())

    def test_min_change_len(self):
        self.run_tool("--write-change-bed", "--min-change-len", "201")
        self.assertEqual((self.out / "state_changes.bed").read_text(), "")

    def test_html_escapes_state_names(self):
        self.mt.write_text("chr1\t0\t800\t<script>\n")
        self.run_tool("--state-mode", "distinct")
        h = (self.out / "report.html").read_text()
        self.assertIn("&lt;script&gt;", h)
        self.assertNotIn("<script>", h)


if __name__ == "__main__":
    unittest.main()
