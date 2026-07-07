nextflow.enable.types = true

include { LIFTOVER_BY_ID } from '../modules/local/liftover_by_id'

record SplitLiftoverOptions {
    split_bed: Path?
    split_genome_fai: Path?
}

workflow LIFTOVER {
    take:
    id_file
    chain
    ref_fa
    query_fa
    ref_fai
    query_fai
    split_options

    main:
    lifted_out = LIFTOVER_BY_ID(id_file, chain, ref_fa, query_fa, ref_fai, query_fai, split_options)

    emit:
    files    = lifted_out.files
    versions = lifted_out.versions
}
