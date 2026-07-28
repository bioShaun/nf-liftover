#!/usr/bin/env nextflow

nextflow.enable.types = true

include { PREPARE_GENOMES } from './subworkflows/prepare_genomes'
include { ALIGN_AND_CHAIN  } from './subworkflows/align_and_chain'
include { LIFTOVER         } from './subworkflows/liftover'
include { LIFTOVER_VCF     } from './modules/local/liftover_vcf'
include { COLLATE_VERSIONS } from './modules/local/collate_versions'
include { WRITE_CHAIN_META } from './modules/local/write_chain_meta'
include { WRITE_RUN_META   } from './modules/local/write_run_meta'

def helpMessage() {
    return """
    nf-liftover  —  end-to-end genome chain generation and coordinate liftover

    USAGE
      nextflow run . --ref_fa REF.fa --query_fa QUERY.fa --id SITES.id [options]
      nextflow run . --ref_fa REF.fa --query_fa QUERY.fa --vcf INPUT.vcf.gz [options]
      nextflow run . --help

    REQUIRED
      --ref_fa FILE          Source/reference FASTA (coordinates of input sites)
      --query_fa FILE        Target/query FASTA (liftover destination)
      --id FILE | --vcf FILE Exactly one input mode

    OPTIONAL INPUTS
      --ref_fai FILE         FASTA index for --ref_fa (else generated)
      --query_fai FILE       FASTA index for --query_fa (else generated)
      --mapping FILE         Two-column TSV: source_chr  target_chr
      --chain FILE           Existing chain; skip minimap2 / PAF-to-chain
      --chain_meta FILE      Optional chain_meta.yml sidecar to verify FASTA checksums
      --split_bed FILE       ID mode: emit BED in split-genome coordinates
      --split_genome_fai FILE  FAI for sorting split BED (default: sibling genome.fa.fai)

    ALIGNMENT
      --pair_strategy STR    Auto chrom pairing: suffix | order  [${params.pair_strategy}]
      --aligner STR          minimap2 | mm2plus  [${params.aligner}]
      --align_mode STR       auto | whole | split  [${params.align_mode}]
      --split_threshold N    auto: split chroms with length >= N  [${params.split_threshold}]
      --split_size N         Sliding window size in split mode  [${params.split_size}]

    LIFTOVER / OUTPUT
      --flank N              snpcalling.bed flank bases  [${params.flank}]
      --outdir DIR           Results directory  [${params.outdir}]
      --publish_paf BOOL     Publish all.paf to outdir/chain  [${params.publish_paf}]

    RESOURCES (per-task caps; also clamp retry growth)
      --max_cpus N           [${params.max_cpus}]
      --max_memory MEM       [${params.max_memory}]
      --max_time DURATION    [${params.max_time}]

    ENVIRONMENT
      --conda_dir DIR        Conda/mamba root  [${params.conda_dir}]

    Full schema: nextflow_schema.json
    Docs: README.md, docs/usage.md
    """.stripIndent()
}

def isTruthyHelp(value) {
    if (value instanceof Boolean) {
        return value
    }
    if (value == null) {
        return false
    }
    def text = value.toString().trim().toLowerCase()
    return text in ['true', 't', '1', 'yes', 'y', 'help']
}

def requireExistingFile(pathValue, name) {
    if (!pathValue) {
        return
    }
    def f = file(pathValue)
    if (!f.exists() || f.isDirectory()) {
        throw new IllegalArgumentException("Parameter --${name} does not exist or is not a file: ${pathValue}")
    }
}

def requirePositiveInt(value, name) {
    if (value == null) {
        throw new IllegalArgumentException("Parameter --${name} must be a positive integer")
    }
    def n = value as Integer
    if (n < 1) {
        throw new IllegalArgumentException("Parameter --${name} must be >= 1 (got ${value})")
    }
}

def requireNonNegativeInt(value, name) {
    if (value == null) {
        throw new IllegalArgumentException("Parameter --${name} must be a non-negative integer")
    }
    def n = value as Integer
    if (n < 0) {
        throw new IllegalArgumentException("Parameter --${name} must be >= 0 (got ${value})")
    }
}


def requireMemory(value, name) {
    if (value == null) {
        throw new IllegalArgumentException("Parameter --${name} must be a memory value like 60.GB")
    }
    def text = value.toString().trim()
    def matcher = (text =~ /^(?i)(\d+(?:\.\d+)?)\s*\.?(B|KB|MB|GB|TB)$/)
    if (!matcher.matches()) {
        throw new IllegalArgumentException("Invalid --${name} '${value}'. Expected form like 8.GB or 512.MB")
    }
    def amount = matcher[0][1] as BigDecimal
    if (amount <= 0) {
        throw new IllegalArgumentException("Parameter --${name} must be > 0 (got ${value})")
    }
}

def requireDuration(value, name) {
    if (value == null) {
        throw new IllegalArgumentException("Parameter --${name} must be a duration like 48.h")
    }
    def text = value.toString().trim()
    def matcher = (text =~ /^(?i)(\d+(?:\.\d+)?)\s*\.?(ms|s|m|h|d)$/)
    if (!matcher.matches()) {
        throw new IllegalArgumentException("Invalid --${name} '${value}'. Expected form like 48.h or 30.m")
    }
    def amount = matcher[0][1] as BigDecimal
    if (amount <= 0) {
        throw new IllegalArgumentException("Parameter --${name} must be > 0 (got ${value})")
    }
}

def gitProvenance(projectDir) {
    def commit = 'unknown'
    def dirty = 'unknown'
    try {
        def pb = new ProcessBuilder('git', '-C', projectDir.toString(), 'rev-parse', 'HEAD')
        pb.redirectErrorStream(true)
        def proc = pb.start()
        def out = proc.inputStream.text.trim()
        proc.waitFor()
        if (proc.exitValue() == 0 && out) {
            commit = out
        }
    } catch (Exception ignored) {
        // leave unknown
    }
    try {
        def pb = new ProcessBuilder('git', '-C', projectDir.toString(), 'status', '--porcelain')
        pb.redirectErrorStream(true)
        def proc = pb.start()
        def out = proc.inputStream.text
        proc.waitFor()
        if (proc.exitValue() == 0) {
            dirty = out.trim() ? 'true' : 'false'
        }
    } catch (Exception ignored) {
        // leave unknown
    }
    return [commit: commit, dirty: dirty]
}

def validateParameters() {
    if (isTruthyHelp(params.help)) {
        return
    }

    ['ref_fa', 'query_fa'].each { name ->
        if (!params[name]) {
            throw new IllegalArgumentException("Missing required parameter --${name}")
        }
    }

    if (!params.id && !params.vcf) {
        throw new IllegalArgumentException("Missing required input: provide --id or --vcf")
    }
    if (params.id && params.vcf) {
        throw new IllegalArgumentException("Provide only one input type: --id or --vcf")
    }

    def alignerEnvs = params.aligner_envs ?: [:]
    def aligner = params.aligner ?: 'minimap2'
    def supportedAligners = ['minimap2', 'mm2plus']
    if (!(aligner in supportedAligners)) {
        throw new IllegalArgumentException(
            "Unsupported --aligner '${params.aligner}'; supported: minimap2, mm2plus.\n" +
            "If you upgraded an existing nextflow.config, copy the `aligner_envs` block from nextflow.config.example.")
    }
    // Bare minimap2 is allowed without a complete aligner_envs map (legacy configs).
    // Any other supported aligner needs a non-empty env directory name.
    if (aligner != 'minimap2') {
        def envName = alignerEnvs instanceof Map ? alignerEnvs.get(aligner) : null
        if (!(envName instanceof CharSequence) || !envName.toString().trim()) {
            throw new IllegalArgumentException(
                "Missing aligner_envs['${aligner}'] for --aligner '${aligner}'; supported: minimap2, mm2plus.\n" +
                "If you upgraded an existing nextflow.config, copy the `aligner_envs` block from nextflow.config.example.")
        }
    }

    def allowedPair = ['order', 'suffix'] as Set
    if (!(params.pair_strategy in allowedPair)) {
        throw new IllegalArgumentException("Invalid --pair_strategy '${params.pair_strategy}'. Allowed: ${allowedPair.join(', ')}")
    }

    def allowedAlign = ['auto', 'whole', 'split'] as Set
    if (!(params.align_mode in allowedAlign)) {
        throw new IllegalArgumentException("Invalid --align_mode '${params.align_mode}'. Allowed: ${allowedAlign.join(', ')}")
    }

    requirePositiveInt(params.split_threshold, 'split_threshold')
    requirePositiveInt(params.split_size, 'split_size')
    requireNonNegativeInt(params.flank, 'flank')
    requirePositiveInt(params.max_cpus, 'max_cpus')
    requireMemory(params.max_memory, 'max_memory')
    requireDuration(params.max_time, 'max_time')

    if (params.vcf && (params.split_bed || params.split_genome_fai)) {
        throw new IllegalArgumentException("--split_bed / --split_genome_fai apply only to ID mode, not --vcf")
    }
    if (params.split_genome_fai && !params.split_bed) {
        throw new IllegalArgumentException("--split_genome_fai requires --split_bed")
    }
    if (params.chain_meta && !params.chain) {
        throw new IllegalArgumentException("--chain_meta requires --chain")
    }

    if (params.chain && params.align_mode != 'auto') {
        log.warn "Ignoring --align_mode=${params.align_mode} because --chain was provided (alignment is skipped)"
    }

    // File existence (fail-fast before process submission)
    requireExistingFile(params.ref_fa, 'ref_fa')
    requireExistingFile(params.query_fa, 'query_fa')
    requireExistingFile(params.ref_fai, 'ref_fai')
    requireExistingFile(params.query_fai, 'query_fai')
    requireExistingFile(params.mapping, 'mapping')
    requireExistingFile(params.chain, 'chain')
    requireExistingFile(params.chain_meta, 'chain_meta')
    requireExistingFile(params.id, 'id')
    requireExistingFile(params.vcf, 'vcf')
    requireExistingFile(params.split_bed, 'split_bed')
    requireExistingFile(params.split_genome_fai, 'split_genome_fai')
}

def chainMetaPayload(mode) {
    return [
        mode            : mode,
        ref_fa          : params.ref_fa?.toString(),
        query_fa        : params.query_fa?.toString(),
        chain_path      : params.chain ? params.chain.toString() : 'all.chain',
        chain_meta      : params.chain_meta?.toString(),
        pair_strategy   : params.pair_strategy?.toString(),
        aligner         : (params.aligner ?: 'minimap2').toString(),
        align_mode      : params.align_mode?.toString(),
        split_threshold : params.split_threshold?.toString(),
        split_size      : params.split_size?.toString(),
        minimap2_args   : (params.minimap2_args ?: '-cx asm5 --cs').toString(),
    ]
}

workflow {
    if (isTruthyHelp(params.help)) {
        log.info helpMessage()
        return
    }

    validateParameters()

    mapping_ch = params.mapping
        ? channel.value(record(mapping_file: file(params.mapping)))
        : channel.value(record(mapping_file: null))

    split_liftover_ch = channel.value(record(
        split_bed: params.split_bed ? file(params.split_bed) : null,
        split_genome_fai: params.split_genome_fai ? file(params.split_genome_fai) : null
    ))

    prepared = PREPARE_GENOMES(
        channel.of(record(role: 'ref', fasta: file(params.ref_fa), fai: params.ref_fai ? file(params.ref_fai) : null)),
        channel.of(record(role: 'query', fasta: file(params.query_fa), fai: params.query_fai ? file(params.query_fai) : null)),
        mapping_ch,
        channel.value(params.pair_strategy)
    )

    chain_meta_opt_ch = channel.value(record(
        chain_meta_file: params.chain_meta ? file(params.chain_meta) : null
    ))

    if (params.chain) {
        // Reuse an existing chain; skip expensive genome alignment.
        // WRITE_CHAIN_META validates chain headers against ref/query FAI (names, lengths, orientation).
        chain_ch = channel.value(file(params.chain))
        def reuseMetaJson = groovy.json.JsonOutput.toJson(chainMetaPayload('reuse'))
        chain_meta = WRITE_CHAIN_META(
            prepared.ref_fa,
            prepared.query_fa,
            prepared.ref_fai,
            prepared.query_fai,
            prepared.chrom_pairs,
            chain_ch,
            channel.value('reuse'),
            channel.value(reuseMetaJson),
            chain_meta_opt_ch
        )
        aligned_versions = chain_meta.versions
        // Prefer the published/validated copy for liftover so paths are consistent.
        chain_for_lift = chain_meta.chain_out
    } else {
        aligned = ALIGN_AND_CHAIN(
            prepared.ref_fa,
            prepared.query_fa,
            prepared.ref_fai,
            prepared.query_fai,
            prepared.chrom_pairs
        )
        def builtMetaJson = groovy.json.JsonOutput.toJson(chainMetaPayload('built'))
        chain_meta = WRITE_CHAIN_META(
            prepared.ref_fa,
            prepared.query_fa,
            prepared.ref_fai,
            prepared.query_fai,
            prepared.chrom_pairs,
            aligned.chain,
            channel.value('built'),
            channel.value(builtMetaJson),
            channel.value(record(chain_meta_file: null))
        )
        aligned_versions = aligned.versions.mix(chain_meta.versions)
        chain_for_lift = aligned.chain
    }

    if (params.vcf) {
        lifted = LIFTOVER_VCF(
            channel.value(file(params.vcf)),
            chain_for_lift,
            prepared.ref_fa,
            prepared.query_fa,
            prepared.ref_fai,
            prepared.query_fai
        )
    } else {
        lifted = LIFTOVER(
            channel.value(file(params.id)),
            chain_for_lift,
            prepared.ref_fa,
            prepared.query_fa,
            prepared.ref_fai,
            prepared.query_fai,
            split_liftover_ch
        )
    }

    ch_versions = prepared.versions
        .mix(aligned_versions)
        .mix(lifted.versions)
        .collect()

    // Single-output process: call returns the versions path channel directly.
    versions_yml_ch = COLLATE_VERSIONS(workflow.nextflow.version.toString(), ch_versions)

    def gitMeta = gitProvenance(workflow.projectDir)
    def runMeta = [
        id               : params.id?.toString(),
        vcf              : params.vcf?.toString(),
        ref_fa           : params.ref_fa?.toString(),
        query_fa         : params.query_fa?.toString(),
        mapping          : params.mapping?.toString(),
        chain            : params.chain?.toString(),
        chain_meta       : params.chain_meta?.toString(),
        pair_strategy    : params.pair_strategy?.toString(),
        aligner          : (params.aligner ?: 'minimap2').toString(),
        align_mode       : params.align_mode?.toString(),
        split_threshold  : params.split_threshold?.toString(),
        split_size       : params.split_size?.toString(),
        flank            : params.flank?.toString(),
        publish_paf      : params.publish_paf?.toString(),
        outdir           : params.outdir?.toString(),
        profile          : workflow.profile?.toString(),
        manifest_version : workflow.manifest.version?.toString() ?: 'unknown',
        git_commit       : workflow.commitId?.toString() ?: gitMeta.commit,
        git_dirty        : gitMeta.dirty,
        revision         : workflow.revision?.toString() ?: workflow.manifest.version?.toString() ?: 'unknown',
        run_name         : workflow.runName?.toString(),
        command_line     : workflow.commandLine?.toString(),
    ]
    def runMetaJson = groovy.json.JsonOutput.toJson(runMeta)

    WRITE_RUN_META(
        versions_yml_ch,
        channel.value(runMetaJson)
    )
}
