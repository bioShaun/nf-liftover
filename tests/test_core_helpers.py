from pathlib import Path
import gzip
import importlib.util
import json
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
            ref_fai = tmp / "ref.fa.fai"
            out = tmp / "restored.paf"

            ref_fai.write_text("SL4.0ch01\t100000\t0\t80\t81\n", encoding="utf-8")
            split_paf.write_text(
                "SL4.0ch01_sliding:10001-20000\t10000\t5\t20\t+\tla2093.chr01\t120000\t30\t45\t15\t15\t60\tcs:Z::15\n",
                encoding="utf-8",
            )

            module.restore_split_paf(split_paf, ref_fai, out)

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

    def test_split_bed_helpers_convert_bed_outputs_only(self):
        module = load_module("liftover_by_id_split", REPO_ROOT / "bin" / "liftover_by_id.py")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            split_bed = tmp / "split.bed"
            split_fai = tmp / "genome.fa.fai"
            outdir = tmp / "out"
            outdir.mkdir()

            split_bed.write_text(
                "chr01\t0\t3000\tchr01_part1\n"
                "chr01\t3000\t10000\tchr01_part2\n",
                encoding="utf-8",
            )
            split_fai.write_text(
                "chr01_part1\t3000\t0\t80\t81\n"
                "chr01_part2\t7000\t0\t80\t81\n",
                encoding="utf-8",
            )

            sorted_bed = module.pd.DataFrame(
                [
                    {
                        "chrom": "chr01",
                        "start": 999,
                        "pos": 1000,
                        "id": "SL4.0ch01_1000",
                        "pos_id": "chr01_1000",
                    },
                    {
                        "chrom": "chr01",
                        "start": 4999,
                        "pos": 5000,
                        "id": "SL4.0ch01_5000",
                        "pos_id": "chr01_5000",
                    },
                ]
            )
            snpcalling_sorted = module.pd.DataFrame(
                [
                    {"chrom": "chr01", "start": 899, "end": 1100},
                    {"chrom": "chr01", "start": 4899, "end": 5100},
                ]
            )

            module.write_liftover_outputs(
                sorted_bed=sorted_bed,
                snpcalling_sorted=snpcalling_sorted,
                outdir=outdir,
                probe_name="probe",
                split_bed=split_bed,
                split_genome_fai=split_fai,
            )

            self.assertEqual((outdir / "probe.id").read_text(encoding="utf-8"), "chr01_1000\nchr01_5000\n")
            self.assertEqual(
                (outdir / "probe.pos.tsv").read_text(encoding="utf-8"),
                "chr01\t1000\tSL4.0ch01_1000\nchr01\t5000\tSL4.0ch01_5000\n",
            )
            self.assertEqual(
                (outdir / "probe.bed").read_text(encoding="utf-8"),
                "chr01_part1\t999\t1000\nchr01_part2\t1999\t2000\n",
            )
            self.assertEqual(
                (outdir / "probe.snpcalling.bed").read_text(encoding="utf-8"),
                "chr01_part1\t899\t1100\nchr01_part2\t1899\t2100\n",
            )

    def test_split_bed_helpers_infer_genome_fai_from_split_bed_directory(self):
        module = load_module("liftover_by_id_split_infer", REPO_ROOT / "bin" / "liftover_by_id.py")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            split_bed = tmp / "split.bed"
            inferred_fai = tmp / "genome.fa.fai"
            split_bed.write_text("chr01\t0\t10\tchr01_part1\n", encoding="utf-8")
            inferred_fai.write_text("chr01_part1\t10\t0\t80\t81\n", encoding="utf-8")

            self.assertEqual(module.resolve_split_genome_fai(split_bed, None), inferred_fai)

    def test_split_bed_dataframe_handles_boundary_overlaps_with_clamping(self):
        module = load_module("liftover_by_id_split_overlap", REPO_ROOT / "bin" / "liftover_by_id.py")

        split_bed_df = module.pd.DataFrame(
            [
                {"chrom": "chr01", "split_start": 0, "split_end": 100, "new_chrom": "chr01_part1"},
                {"chrom": "chr01", "split_start": 100, "split_end": 200, "new_chrom": "chr01_part2"},
            ]
        )

        bed_df = module.pd.DataFrame(
            [
                {"chrom": "chr01", "start": 90, "pos": 110},
            ]
        )

        split_df = module.split_bed_dataframe(bed_df, split_bed_df, start_col="start", end_col="pos")
        records = split_df.to_dict(orient="records")

        self.assertEqual(len(records), 2)
        self.assertEqual(records[0], {"new_chrom": "chr01_part1", "new_start": 90, "new_end": 100})
        self.assertEqual(records[1], {"new_chrom": "chr01_part2", "new_start": 0, "new_end": 10})


class MergeBlastRescueTests(unittest.TestCase):
    def test_merges_only_strict_blast_rescue_and_reports_status(self):
        module = load_module("merge_blast_rescue", REPO_ROOT / "bin" / "merge_blast_rescue.py")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            id_file = tmp / "panel.id"
            transanno_pos = tmp / "panel.pos.tsv"
            rejected_vcf = tmp / "rejected.panel.id.vcf.gz"
            selection = tmp / "blast.selection.tsv"
            mapping = tmp / "mapping.tsv"
            query_fai = tmp / "query.fa.fai"
            rescue_out = tmp / "blast_rescue"

            id_file.write_text("src1_10\nsrc1_20\nsrc1_30\nsrc1_40\nsrc1_50\n", encoding="utf-8")
            transanno_pos.write_text("chr1\t101\tsrc1_10\n", encoding="utf-8")
            with gzip.open(rejected_vcf, "wt", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.2\n")
                for marker_id, reason in [
                    ("src1_20", "NO_CHAIN"),
                    ("src1_30", "NO_CHAIN"),
                    ("src1_40", "MULTIMAP"),
                    ("src1_50", "UNEXPECTED_REF"),
                ]:
                    handle.write(f"src1\t1\t{marker_id}\tA\t.\t.\t.\tFAILED_REASON={reason}\n")
            selection.write_text(
                "id\tnew_id\tchrom\tpos\tstrand\talleles\trank\tbitscore\tmatch_ratio\talign_len\t"
                "informative_len\tn_count\tmismatches\tgap_opens\tchr_map_status\tselection_reason\n"
                "src1_20\tchr1_200\tchr1\t200\t+\t-/-\t1\t100\t0.97\t100\t100\t0\t0\t0\tmatched\tok\n"
                "src1_30\tchr1_300\tchr1\t300\t+\t-/-\t1\t90\t0.90\t90\t100\t0\t0\t0\tmatched\tlow\n"
                "src1_40\tchr1_400\tchr1\t400\t+\t-/-\t1\t80\t0.99\t99\t100\t0\t0\t0\tmatched\ttie1\n"
                "src1_40\tchr1_401\tchr1\t401\t+\t-/-\t2\t80\t0.99\t99\t100\t0\t0\t0\tmatched\ttie2\n"
                "src1_50\tchr2_500\tchr2\t500\t+\t-/-\t1\t100\t0.99\t99\t100\t0\t0\t0\tfallback\tfallback\n",
                encoding="utf-8",
            )
            mapping.write_text("src1\tchr1\n", encoding="utf-8")
            query_fai.write_text("chr1\t1000\t0\t80\t81\nchr2\t1000\t0\t80\t81\n", encoding="utf-8")

            module.merge_blast_rescue(
                id_file=id_file,
                transanno_pos=transanno_pos,
                rejected_vcf=rejected_vcf,
                blast_selection=selection,
                mapping_tsv=mapping,
                query_fai=query_fai,
                outdir=rescue_out,
                panel="panel",
                flank=10,
            )

            self.assertEqual((rescue_out / "panel.id").read_text(encoding="utf-8"), "chr1_101\nchr1_200\n")
            self.assertEqual(
                (rescue_out / "panel.bed").read_text(encoding="utf-8"),
                "chr1\t100\t101\nchr1\t199\t200\n",
            )
            self.assertEqual(
                (rescue_out / "panel.snpcalling.bed").read_text(encoding="utf-8"),
                "chr1\t90\t111\nchr1\t189\t210\n",
            )

            summary = module.pd.read_table(rescue_out / "panel.summary.tsv")
            self.assertEqual(summary.set_index("metric")["count"].to_dict(), {
                "transanno_success": 1,
                "blast_rescue": 1,
                "failed": 3,
            })
            failed_reasons = (rescue_out / "panel.failed-reasons.tsv").read_text(encoding="utf-8")
            self.assertIn("NO_CHAIN;BLAST_LOW_MATCH_RATIO", failed_reasons)
            self.assertIn("MULTIMAP;BLAST_TIED_TOP_HIT", failed_reasons)
            self.assertIn("UNEXPECTED_REF;BLAST_CHR_MAP_NOT_MATCHED+BLAST_TARGET_CHROM_NOT_IN_MAPPING", failed_reasons)

    def test_split_output_keeps_panel_filenames_in_rescue_outdir(self):
        module = load_module("merge_blast_rescue_split", REPO_ROOT / "bin" / "merge_blast_rescue.py")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            outdir = tmp / "rescue"
            query_fai = tmp / "query.fa.fai"
            split_bed = tmp / "split.bed"
            split_fai = tmp / "split.fa.fai"

            query_fai.write_text("chr1\t1000\t0\t80\t81\n", encoding="utf-8")
            split_bed.write_text("chr1\t0\t150\tchr1a\nchr1\t150\t1000\tchr1b\n", encoding="utf-8")
            split_fai.write_text("chr1a\t150\t0\t80\t81\nchr1b\t850\t0\t80\t81\n", encoding="utf-8")
            combined = module.pd.DataFrame(
                [
                    {"chrom": "chr1", "pos": 101, "id": "src1_10", "method": "transanno"},
                    {"chrom": "chr1", "pos": 200, "id": "src1_20", "method": "blast_rescue"},
                ]
            )

            module.write_combined_outputs(
                combined,
                query_fai,
                outdir,
                "panel",
                flank=10,
                split_bed=split_bed,
                split_genome_fai=split_fai,
            )

            self.assertEqual((outdir / "panel.id").read_text(encoding="utf-8"), "chr1_101\nchr1_200\n")
            self.assertEqual((outdir / "panel.bed").read_text(encoding="utf-8"), "chr1a\t100\t101\nchr1b\t49\t50\n")
            self.assertEqual(
                (outdir / "panel.snpcalling.bed").read_text(encoding="utf-8"),
                "chr1a\t90\t111\nchr1b\t39\t60\n",
            )
            self.assertTrue((outdir / "panel.method.tsv").is_file())



class ChainMetaTests(unittest.TestCase):
    def test_validate_accepts_matching_chain_headers(self):
        module = load_module("chain_meta", REPO_ROOT / "bin" / "chain_meta.py")
        blocks = [
            {
                "t_name": "SL4.0ch01",
                "t_size": 10000,
                "t_strand": "+",
                "q_name": "chr01",
                "q_size": 10000,
                "q_strand": "+",
            }
        ]
        ref_len = {"SL4.0ch01": 10000}
        query_len = {"chr01": 10000}
        module.validate_chain_against_genomes(blocks, ref_len, query_len)

    def test_validate_rejects_reversed_chain(self):
        module = load_module("chain_meta", REPO_ROOT / "bin" / "chain_meta.py")
        blocks = [
            {
                "t_name": "chr01",
                "t_size": 10000,
                "t_strand": "+",
                "q_name": "SL4.0ch01",
                "q_size": 10000,
                "q_strand": "+",
            }
        ]
        ref_len = {"SL4.0ch01": 10000}
        query_len = {"chr01": 10000}
        with self.assertRaises(ValueError) as ctx:
            module.validate_chain_against_genomes(blocks, ref_len, query_len)
        self.assertIn("reversed", str(ctx.exception).lower())

    def test_validate_rejects_size_mismatch(self):
        module = load_module("chain_meta", REPO_ROOT / "bin" / "chain_meta.py")
        blocks = [
            {
                "t_name": "SL4.0ch01",
                "t_size": 9999,
                "t_strand": "+",
                "q_name": "chr01",
                "q_size": 10000,
                "q_strand": "+",
            }
        ]
        ref_len = {"SL4.0ch01": 10000}
        query_len = {"chr01": 10000}
        with self.assertRaises(ValueError) as ctx:
            module.validate_chain_against_genomes(blocks, ref_len, query_len)
        self.assertIn("size", str(ctx.exception).lower())

    def test_sidecar_sha_mismatch_fails(self):
        module = load_module("chain_meta_side", REPO_ROOT / "bin" / "chain_meta.py")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            side = tmp / "chain_meta.yml"
            side.write_text(
                "source_fasta:\n"
                "  sha256: \"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\"\n"
                "target_fasta:\n"
                "  sha256: \"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\"\n"
                "chain:\n"
                "  sha256: \"cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc\"\n",
                encoding="utf-8",
            )
            with self.assertRaises(ValueError):
                module.validate_sidecar_required(
                    str(side),
                    "dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd",
                    "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
                    "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc",
                )

    def test_yaml_escaping_for_special_paths(self):
        module = load_module("chain_meta_escape", REPO_ROOT / "bin" / "chain_meta.py")
        lines = module.build_chain_meta_doc(
            mode="reuse",
            ref_path=r"/data/ref$weird/name's.fa",
            query_path=r"/data/query`x`.fa",
            chain_path=r"/data/all.chain",
            ref_sha="abc",
            ref_size=1,
            query_sha="def",
            query_size=2,
            chain_sha="ghi",
            chain_size=3,
            blocks=[{"t_name": "a", "q_name": "b"}],
            pair_lines=["a\tb"],
            pair_strategy="suffix",
            align_mode="auto",
            split_threshold=1,
            split_size=1,
            minimap2_args="-cx asm5 --cs",
            sidecar_info=None,
        )
        text = "\n".join(lines)
        self.assertIn(json.dumps(r"/data/ref$weird/name's.fa"), text)
        self.assertIn(json.dumps(r"/data/query`x`.fa"), text)


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
        prepare_nf = (
            (REPO_ROOT / "subworkflows" / "prepare_genomes.nf").read_text(encoding="utf-8") + "\n" +
            (REPO_ROOT / "modules" / "local" / "maybe_faidx.nf").read_text(encoding="utf-8") + "\n" +
            (REPO_ROOT / "modules" / "local" / "derive_chrom_pairs.nf").read_text(encoding="utf-8")
        )

        self.assertIn("record(role: 'ref'", main_nf)
        self.assertIn("record(role: 'query'", main_nf)
        self.assertIn("record(mapping_file: file(params.mapping))", main_nf)
        self.assertIn("record(mapping_file: null)", main_nf)
        self.assertIn("record(", prepare_nf)
        self.assertNotIn("record GenomeInput", prepare_nf)
        self.assertNotIn("record ChromMapping", prepare_nf)
        self.assertNotIn("tuple val(role), val(source_fasta), path(fasta), val(fai_hint)", prepare_nf)
        self.assertNotIn("tuple val(has_mapping), path(mapping_file)", prepare_nf)

    def test_liftover_process_uses_typed_inputs_and_outputs(self):
        liftover_nf = (
            (REPO_ROOT / "subworkflows" / "liftover.nf").read_text(encoding="utf-8") + "\n" +
            (REPO_ROOT / "modules" / "local" / "liftover_by_id.nf").read_text(encoding="utf-8")
        )

        self.assertIn("nextflow.enable.types = true", liftover_nf)
        self.assertIn("id_file: Path", liftover_nf)
        self.assertIn("chain: Path", liftover_nf)
        self.assertIn("ref_fai: Path", liftover_nf)
        self.assertIn("files = files('out/*')", liftover_nf)
        self.assertNotIn("path id_file", liftover_nf)
        self.assertNotIn("path 'out/*', emit: files", liftover_nf)

    def test_align_and_chain_process_uses_typed_inputs_and_outputs(self):
        modules_content = ""
        for name in ["build_align_pairs", "extract_chrom_fastas", "build_query_mmi", "align_whole_chromosome", "split_ref_chromosome", "align_split_window", "combine_split_pafs", "combine_all_pafs", "paf_to_chain"]:
            modules_content += (REPO_ROOT / "modules" / "local" / f"{name}.nf").read_text(encoding="utf-8") + "\n"
        align_nf = (REPO_ROOT / "subworkflows" / "align_and_chain.nf").read_text(encoding="utf-8") + "\n" + modules_content

        self.assertIn("nextflow.enable.types = true", align_nf)
        self.assertIn("ref_fa: Path", align_nf)
        self.assertIn("query_fa: Path", align_nf)
        self.assertIn("paf: Path = file('all.paf')", align_nf)
        self.assertIn("chain: Path = file('all.chain')", align_nf)
        self.assertIn("process ALIGN_WHOLE_CHROMOSOME", align_nf)
        self.assertIn("process SPLIT_REF_CHROMOSOME", align_nf)
        self.assertIn("process ALIGN_SPLIT_WINDOW", align_nf)
        self.assertIn("process COMBINE_SPLIT_PAFS", align_nf)
        self.assertIn("process COMBINE_ALL_PAFS", align_nf)
        self.assertNotIn("while IFS=", align_nf)
        self.assertNotIn("done <", align_nf)
        self.assertNotIn("path ref_fa", align_nf)
        self.assertNotIn("path 'all.paf', emit: paf", align_nf)

    def test_liftover_process_passes_split_bed_options(self):
        liftover_nf = (
            (REPO_ROOT / "subworkflows" / "liftover.nf").read_text(encoding="utf-8") + "\n" +
            (REPO_ROOT / "modules" / "local" / "liftover_by_id.nf").read_text(encoding="utf-8")
        )

        self.assertIn("record SplitLiftoverOptions", liftover_nf)
        self.assertIn("split_bed: Path?", liftover_nf)
        self.assertIn("split_genome_fai: Path?", liftover_nf)
        self.assertIn("stageAs ref_fa, 'liftover_ref.fa'", liftover_nf)
        self.assertIn("stageAs query_fa, 'liftover_query.fa'", liftover_nf)
        self.assertIn("stageAs ref_fai, 'liftover_ref.fa.fai'", liftover_nf)
        self.assertIn("stageAs query_fai, 'liftover_query.fa.fai'", liftover_nf)
        self.assertIn("stageAs split_bed, 'split_liftover.bed'", liftover_nf)
        self.assertIn("stageAs split_genome_fai, 'split_liftover.genome.fai'", liftover_nf)
        self.assertIn("--split-bed", liftover_nf)
        self.assertIn("--split-genome-fai", liftover_nf)

    def test_vcf_liftover_process_is_available(self):
        main_nf = (REPO_ROOT / "main.nf").read_text(encoding="utf-8")
        vcf_nf = (REPO_ROOT / "modules" / "local" / "liftover_vcf.nf").read_text(encoding="utf-8")

        self.assertIn("params.vcf", main_nf)
        self.assertIn("include { LIFTOVER_VCF", main_nf)
        self.assertIn("process LIFTOVER_VCF", vcf_nf)
        self.assertIn("vcf_file: Path", vcf_nf)
        self.assertIn("transanno liftvcf", vcf_nf)
        self.assertIn("tabix -f -p vcf", vcf_nf)
        self.assertIn("stageAs ref_fai, 'liftover_ref.fa.fai'", vcf_nf)
        self.assertIn("stageAs query_fai, 'liftover_query.fa.fai'", vcf_nf)

    def test_processes_use_stage_aliases_for_same_named_inputs(self):
        prepare_nf = (
            (REPO_ROOT / "subworkflows" / "prepare_genomes.nf").read_text(encoding="utf-8") + "\n" +
            (REPO_ROOT / "modules" / "local" / "maybe_faidx.nf").read_text(encoding="utf-8") + "\n" +
            (REPO_ROOT / "modules" / "local" / "derive_chrom_pairs.nf").read_text(encoding="utf-8")
        )
        modules_content = ""
        for name in ["build_align_pairs", "extract_chrom_fastas", "build_query_mmi", "align_whole_chromosome", "split_ref_chromosome", "align_split_window", "combine_split_pafs", "combine_all_pafs", "paf_to_chain"]:
            modules_content += (REPO_ROOT / "modules" / "local" / f"{name}.nf").read_text(encoding="utf-8") + "\n"
        align_nf = (REPO_ROOT / "subworkflows" / "align_and_chain.nf").read_text(encoding="utf-8") + "\n" + modules_content

        self.assertIn("stageAs ref_fai, 'ref.fa.fai'", prepare_nf)
        self.assertIn("stageAs query_fai, 'query.fa.fai'", prepare_nf)
        self.assertIn("stageAs chrom_pairs, 'chrom_pairs.tsv'", align_nf)
        self.assertIn("stageAs ref_fai, 'ref.fa.fai'", align_nf)
        self.assertIn("stageAs ref_fa, 'ref_genome.fa'", align_nf)
        self.assertIn("stageAs query_fa, 'query_genome.fa'", align_nf)

    def test_schema_defines_alignment_and_split_bed_parameters(self):
        schema = json.loads((REPO_ROOT / "nextflow_schema.json").read_text(encoding="utf-8"))
        properties = schema["properties"]

        self.assertEqual(properties["align_mode"]["enum"], ["auto", "whole", "split"])
        self.assertEqual(properties["align_mode"]["default"], "auto")
        self.assertEqual(properties["split_bed"]["type"], ["string", "null"])
        self.assertEqual(properties["split_genome_fai"]["type"], ["string", "null"])
        self.assertEqual(properties["id"]["type"], ["string", "null"])
        self.assertEqual(properties["vcf"]["type"], ["string", "null"])

    def test_main_workflow_uses_explicit_call_outputs(self):
        main_nf = (REPO_ROOT / "main.nf").read_text(encoding="utf-8")

        self.assertIn("prepared = PREPARE_GENOMES(", main_nf)
        self.assertIn("aligned = ALIGN_AND_CHAIN(", main_nf)
        self.assertNotIn("PREPARE_GENOMES.out.", main_nf)
        self.assertNotIn("ALIGN_AND_CHAIN.out.", main_nf)

    def test_nextflow_scripts_do_not_use_legacy_out_property(self):
        scripts = [REPO_ROOT / "main.nf", *sorted((REPO_ROOT / "subworkflows").glob("*.nf"))]
        offenders = [
            str(script.relative_to(REPO_ROOT))
            for script in scripts
            if ".out." in script.read_text(encoding="utf-8")
        ]

        self.assertEqual(offenders, [])


if __name__ == "__main__":
    unittest.main()
