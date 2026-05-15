from pathlib import Path
import importlib.util
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
TOMATO_SMOKE = REPO_ROOT / "tests" / "data" / "tomato-smoke"


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


class TomatoSmokeFixtureTests(unittest.TestCase):
    def test_tomato_smoke_fixture_has_pairable_fastas_and_ids(self):
        derive_module = load_module("derive_chrom_pairs_fixture", REPO_ROOT / "bin" / "derive_chrom_pairs.py")
        liftover_module = load_module("liftover_by_id_fixture", REPO_ROOT / "bin" / "liftover_by_id.py")

        ref_fai = TOMATO_SMOKE / "SL4.0ch01.10kb.fa.fai"
        query_fai = TOMATO_SMOKE / "LA2093.chr01.10kb.fa.fai"
        ids = TOMATO_SMOKE / "tomato-smoke.id"
        mapping = TOMATO_SMOKE / "chrom_pairs.tsv"

        for fixture in (ref_fai, query_fai, ids, mapping):
            self.assertTrue(fixture.is_file(), f"Missing fixture: {fixture}")

        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "derived.tsv"
            derive_module.derive_chrom_pairs(ref_fai, query_fai, out, mapping=mapping)
            self.assertEqual(out.read_text(encoding="utf-8"), "SL4.0ch01\tchr01\n")

        parsed_ids = [
            liftover_module.parse_id_to_chrom_pos(line.strip())
            for line in ids.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        self.assertEqual(parsed_ids, [("SL4.0ch01", 1000), ("SL4.0ch01", 5000), ("SL4.0ch01", 9000)])


class NextflowTypedMigrationTests(unittest.TestCase):
    def test_prepare_genomes_uses_records_for_high_risk_channel_shapes(self):
        main_nf = (REPO_ROOT / "main.nf").read_text(encoding="utf-8")
        prepare_nf = (REPO_ROOT / "subworkflows" / "prepare_genomes.nf").read_text(encoding="utf-8")

        self.assertIn("record(role: 'ref'", main_nf)
        self.assertIn("record(role: 'query'", main_nf)
        self.assertIn("record(has_mapping: true", main_nf)
        self.assertIn("record(has_mapping: false", main_nf)
        self.assertIn("record(", prepare_nf)
        self.assertNotIn("tuple val(role), val(source_fasta), path(fasta), val(fai_hint)", prepare_nf)
        self.assertNotIn("tuple val(has_mapping), path(mapping_file)", prepare_nf)

    def test_liftover_process_uses_typed_inputs_and_outputs(self):
        liftover_nf = (REPO_ROOT / "subworkflows" / "liftover.nf").read_text(encoding="utf-8")

        self.assertIn("nextflow.enable.types = true", liftover_nf)
        self.assertIn("id_file: Path", liftover_nf)
        self.assertIn("chain: Path", liftover_nf)
        self.assertIn("files = files('out/*')", liftover_nf)
        self.assertNotIn("path id_file", liftover_nf)
        self.assertNotIn("path 'out/*', emit: files", liftover_nf)

    def test_align_and_chain_process_uses_typed_inputs_and_outputs(self):
        align_nf = (REPO_ROOT / "subworkflows" / "align_and_chain.nf").read_text(encoding="utf-8")

        self.assertIn("nextflow.enable.types = true", align_nf)
        self.assertIn("ref_fa: Path", align_nf)
        self.assertIn("query_fa: Path", align_nf)
        self.assertIn("paf: Path = file('all.paf')", align_nf)
        self.assertIn("chain: Path = file('all.chain')", align_nf)
        self.assertNotIn("path ref_fa", align_nf)
        self.assertNotIn("path 'all.paf', emit: paf", align_nf)

    def test_main_workflow_uses_explicit_call_outputs(self):
        main_nf = (REPO_ROOT / "main.nf").read_text(encoding="utf-8")

        self.assertIn("prepared = PREPARE_GENOMES(", main_nf)
        self.assertIn("aligned = ALIGN_AND_CHAIN(", main_nf)
        self.assertNotIn("PREPARE_GENOMES.out.", main_nf)
        self.assertNotIn("ALIGN_AND_CHAIN.out.", main_nf)


if __name__ == "__main__":
    unittest.main()
