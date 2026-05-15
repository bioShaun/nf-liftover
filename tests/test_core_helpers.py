from pathlib import Path
import importlib.util
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
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


if __name__ == "__main__":
    unittest.main()
