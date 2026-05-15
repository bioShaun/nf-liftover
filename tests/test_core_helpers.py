from pathlib import Path
import importlib.util
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class DeriveChromPairsTests(unittest.TestCase):
    def test_derives_pairs_by_shared_numeric_suffix_when_names_differ(self):
        module = load_module("derive_chrom_pairs", REPO_ROOT / "bin" / "derive_chrom_pairs.py")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            ref_fai = tmp / "ref.fa.fai"
            query_fai = tmp / "query.fa.fai"
            out = tmp / "chrom_pairs.tsv"

            ref_fai.write_text(
                "SL4.0ch01\t1000\t0\t80\t81\n"
                "SL4.0ch02\t2000\t0\t80\t81\n",
                encoding="utf-8",
            )
            query_fai.write_text(
                "la2093.chr02\t2100\t0\t80\t81\n"
                "la2093.chr01\t1100\t0\t80\t81\n",
                encoding="utf-8",
            )

            module.derive_chrom_pairs(ref_fai, query_fai, out, strategy="suffix")

            self.assertEqual(
                out.read_text(encoding="utf-8"),
                "SL4.0ch01\tla2093.chr01\nSL4.0ch02\tla2093.chr02\n",
            )

    def test_user_mapping_is_validated_and_copied(self):
        module = load_module("derive_chrom_pairs", REPO_ROOT / "bin" / "derive_chrom_pairs.py")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            ref_fai = tmp / "ref.fa.fai"
            query_fai = tmp / "query.fa.fai"
            mapping = tmp / "mapping.tsv"
            out = tmp / "chrom_pairs.tsv"

            ref_fai.write_text("A\t100\t0\t80\t81\n", encoding="utf-8")
            query_fai.write_text("B\t100\t0\t80\t81\n", encoding="utf-8")
            mapping.write_text("A\tB\n", encoding="utf-8")

            module.derive_chrom_pairs(ref_fai, query_fai, out, mapping=mapping)

            self.assertEqual(out.read_text(encoding="utf-8"), "A\tB\n")


class RestoreSplitPafTests(unittest.TestCase):
    def test_restores_query_coordinates_and_original_query_length(self):
        module = load_module("restore_split_paf", REPO_ROOT / "bin" / "restore_split_paf.py")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            split_paf = tmp / "split.paf"
            fai = tmp / "query.fa.fai"
            out = tmp / "restored.paf"

            fai.write_text("SL4.0ch01\t100000\t0\t80\t81\n", encoding="utf-8")
            split_paf.write_text(
                "SL4.0ch01_sliding:10001-20000\t10000\t5\t20\t+\tla2093.chr01\t120000\t30\t45\t15\t15\t60\tcs:Z::15\n",
                encoding="utf-8",
            )

            module.restore_split_paf(split_paf, fai, out)

            fields = out.read_text(encoding="utf-8").strip().split("\t")
            self.assertEqual(fields[:4], ["SL4.0ch01", "100000", "10005", "10020"])
            self.assertEqual(fields[4:], ["+", "la2093.chr01", "120000", "30", "45", "15", "15", "60", "cs:Z::15"])


class LiftoverByIdTests(unittest.TestCase):
    def test_parses_chrom_pos_id_from_last_underscore(self):
        module = load_module("liftover_by_id", REPO_ROOT / "bin" / "liftover_by_id.py")

        self.assertEqual(module.parse_id_to_chrom_pos("SL4.0ch01_123"), ("SL4.0ch01", 123))
        self.assertEqual(module.parse_id_to_chrom_pos("chr_name_with_underscores_456"), ("chr_name_with_underscores", 456))

    def test_sorts_bed_records_by_fai_order_then_start(self):
        module = load_module("liftover_by_id", REPO_ROOT / "bin" / "liftover_by_id.py")

        bed_df = module.pd.DataFrame(
            [
                {"chrom": "chr2", "start": 8, "pos": 9},
                {"chrom": "chr1", "start": 20, "pos": 21},
                {"chrom": "chr1", "start": 10, "pos": 11},
            ]
        )
        chrom_df = module.pd.DataFrame(
            [
                {"chrom": "chr1", "chrom_size": 100},
                {"chrom": "chr2", "chrom_size": 100},
            ]
        )

        sorted_df = module.sort_bed_by_fai(bed_df, chrom_df)

        self.assertEqual(sorted_df["chrom"].astype(str).tolist(), ["chr1", "chr1", "chr2"])
        self.assertEqual(sorted_df["start"].tolist(), [10, 20, 8])

    def test_slop_and_merge_clamps_to_chromosome_bounds(self):
        module = load_module("liftover_by_id", REPO_ROOT / "bin" / "liftover_by_id.py")

        bed_df = module.pd.DataFrame(
            [
                {"chrom": "chr1", "start": 2, "pos": 3},
                {"chrom": "chr1", "start": 8, "pos": 9},
                {"chrom": "chr1", "start": 97, "pos": 98},
                {"chrom": "chr2", "start": 10, "pos": 11},
            ]
        )

        merged_df = module.slop_and_merge(bed_df, {"chr1": 100, "chr2": 15}, flank=5)

        self.assertEqual(
            merged_df.to_dict(orient="records"),
            [
                {"chrom": "chr1", "start": 0, "end": 14},
                {"chrom": "chr1", "start": 92, "end": 100},
                {"chrom": "chr2", "start": 5, "end": 15},
            ],
        )


if __name__ == "__main__":
    unittest.main()
