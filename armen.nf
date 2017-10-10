#!/usr/bin/env nextflow
/************************************************
 * Armen - Antibiotic resistance in metagenomes
 * Copyright (c) Authors 2017 
 * Authors:
 *  Fredrik Boulund <fredrik.boulund@ki.se>
************************************************/

params.input_reads = '' // Must be specified on the command line
Channel
    .fromFilePairs(params.input_reads)
    .ifEmpty{ exit 1, "Found no input reads, did you specify --input_reads? I got: '${params.input_reads}'"}
    .into {input_reads_average_genome_size,
           input_reads_qa}


process average_genome_size {
    tag {pair_id}
    publishDir "${params.outdir}/microbecensus", mode: 'copy'

    input:
    set pair_id, file(reads) from input_reads_average_genome_size

    output:
    file "${pair_id}.microbecensus.txt"
    file "${pair_id}.microbecensus.stdout"

    """
    run_microbe_census.py \
        -t ${task.cpus} \
        ${reads[0]},${reads[1]} \
        ${pair_id}.microbecensus.txt \
        > ${pair_id}.microbecensus.stdout
    """
}


process qa_reads {
    tag {pair_id}
    publishDir "${params.outdir}/qa_reads", mode: 'copy'

    input:
    set pair_id, file(reads) from input_reads_qa

    output:
    set pair_id, file("${pair_id}_{1,2}.trimmed.fq.gz") into input_remove_human
    file "${pair_id}.stats.txt"
    file "${pair_id}.bhist.txt"
    file "${pair_id}.qhist.txt"
    file "${pair_id}.qchist.txt"
    file "${pair_id}.aqhist.txt"
    file "${pair_id}.bqhist.txt"
    file "${pair_id}.lhist.txt"
    file "${pair_id}.gchist.txt"

    """
    bbduk.sh \
        in1=${reads[0]} \
        in2=${reads[1]} \
        ref=${params.bbduk_ref} \
        out1=${pair_id}_1.trimmed.fq \
        out2=${pair_id}_2.trimmed.fq \
        stats=${pair_id}.stats.txt \
        bhist=${pair_id}.bhist.txt \
        qhist=${pair_id}.qhist.txt \
        qchist=${pair_id}.qchist.txt \
        aqhist=${pair_id}.aqhist.txt \
        bqhist=${pair_id}.bqhist.txt \
        lhist=${pair_id}.lhist.txt \
        gchist=${pair_id}.gchist.txt \
        minlen=${params.bbduk_minlen} \
        qtrim=${params.bbduk_qtrim} \
        trimq=${params.bbduk_trimq} \
        ktrim=${params.bbduk_ktrim} \
        k=${params.bbduk_k} \
        mink=${params.bbduk_mink} \
        hdist=${params.bbduk_hdist} \
    """
} 


process remove_human {
    tag {pair_id}
    publishDir "${params.outdir}/human_filtered_reads", mode: 'copy'

    input:
    set pair_id, file(reads) from input_remove_human

    output:
    set pair_id, file("${pair_id}_{1,2}.human.fq.gz") into input_reads_megares
    file "${pair_id}.human.fq.gz"
    file "${pair_id}.statsfile.txt"

    """
    bbmap.sh \
        threads=${task.cpus} \
        in1=${reads[0]} \
        in2=${reads[1]} \
        path=${params.human_db} \
        outu1=${pair_id}_1.fq.gz \
        outu2=${pair_id}_2.fq.gz \
        outm=${pair_id}.human.fq.gz \
        statsfile=${pair_id}.statsfile.txt \
        minid=${params.removehuman_minid} \
        maxindel=${params.removehuman_maxindel} \
        minhits=${params.removehuman_minhits} \
        bandwidthratio=${params.removehuman_bandwidthratio} \
        bandwidth=${params.removehuman_bandwidth} \
        qtrim=${params.removehuman_qtrim} \
        trimq=${params.removehuman_trimq} \
        ${params.removehuman_quickmatch} \
        ${params.removehuman_fast} \
        ${params.removehuman_untrim} \
    """
}


process megares {
    tag {sample_id}
    publishDir "${params.outdir}/megares", mode: 'copy'

    input:
    set sample_id, file(reads) from input_reads_megares

    output:
    set sample_id, file("${sample_id}.sam.gz") into output_megares_mapped
    file "${sample_id}.scafstats.txt"
    file "${sample_id}.statsfile.txt"
    file "${sample_id}.covstats.txt"
    file "${sample_id}.rpkm.txt"

    """
    bbmap.sh \
        minid=${params.bbmap_minid} \
        threads=${task.cpus} \
        path=${params.megares_path} \
        in1=${reads[0]} \
        in2=${reads[1]} \
        out=${sample_id}.sam.gz \
        scafstats=${sample_id}.scafstats.txt \
        statsfile=${sample_id}.statsfile.txt \
        covstats=${sample_id}.covstats.txt \
        rpkm=${sample_id}.rpkm.txt \
    """
} 
