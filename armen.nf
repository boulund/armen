#!/usr/bin/env nextflow
/****************************************
 * Armen -- 
 * Antibiotic Resistance of Metagenomes
 * Copyright (c) Authors 2017 
 * Authors:
 *  Fredrik Boulund <fredrik.boulund@ki.se>
 ****************************************/

params.input_reads = '' // Specify on command line

Channel
    .fromFilePairs(params.input_reads)
    .ifEmpty{ exit 1, "Found no input reads, did you specify --input_reads? I got: '${params.input_reads}'"}
    .into {input_reads_megares}


process megares {
    tag {sample_id}
    publishDir "${params.outdir}/megares", mode: 'copy'

    input:
    set sample_id, file(reads) from input_reads_megares

    output:
    set sample_id, file("${sample_id}.sam.gz") into output_megares_mapped
    file "${sample_id}.scafstats.txt.gz"
    file "${sample_id}.statsfile.txt.gz"
    file "${sample_id}.covstats.txt.gz"
    file "${sample_id}.rpkm.txt.gz"

    """
    bbmap.sh \
		minid=${params.bbmap_minid} \
		threads=${task.cpus} \
		path=${params.megares_path} \
        in1=${reads[0]} \
        in2=${reads[1]} \
		out=${sample_id}.sam.gz \
		scafstats=${sample_id}.scafstats.txt.gz \
		statsfile=${sample_id}.statsfile.txt.gz \
		covstats=${sample_id}.covstats.txt.gz \
		rpkm=${sample_id}.rpkm.txt.gz \
    """
} 

