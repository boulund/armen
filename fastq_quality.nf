#!/usr/bin/env nextflow
/****************************************
 * Metagenomics read QA workflow
 * Copyright (c) Authors 2017 
 * Authors:
 *  Fredrik Boulund <fredrik.boulund@ki.se>
 ****************************************/

params.input_reads = '' // Specify on command line

Channel
    .fromFilePairs(params.input_reads)
    .ifEmpty{ exit 1, "Found no input reads, did you specify --input_reads? I got: '${params.input_reads}'"}
    .into {input_reads_bbduk}


process qa_reads {
    tag {pair_id}
    publishDir "${params.outdir}/qa_reads", mode: 'copy'

    input:
    set pair_id, file(reads) from input_reads_bbduk

    output:
	set pair_id, file("{$outfiles}.trimmed.fq.gz") into input_remove_human
    file "${reads[0].baseName}.fq.gz" 
    file "${reads[1].baseName}.fq.gz" 
    file "${pair_id}.stats.txt.gz"
    file "${pair_id}.bhist.txt.gz"
    file "${pair_id}.qhist.txt.gz"
    file "${pair_id}.qchist.txt.gz"
    file "${pair_id}.aqhist.txt.gz"
    file "${pair_id}.bqhist.txt.gz"
    file "${pair_id}.lhist.txt.gz"
    file "${pair_id}.gchist.txt.gz"

	script:
	outfiles = reads.collect{ it.baseName }.join(',')
    """
    bbduk.sh \
        in1=${reads[0]} \
        in2=${reads[1]} \
        ref=${params.bbduk_ref} \
        out1=${reads[0].baseName}.trimmed.fq.gz \
        out2=${reads[1].baseName}.trimmed.fq.gz \
        stats=${pair_id}.stats.txt.gz \
        bhist=${pair_id}.bhist.txt.gz \
        qhist=${pair_id}.qhist.txt.gz \
        qchist=${pair_id}.qchist.txt.gz \
        aqhist=${pair_id}.aqhist.txt.gz \
        bqhist=${pair_id}.bqhist.txt.gz \
        lhist=${pair_id}.lhist.txt.gz \
        gchist=${pair_id}.gchist.txt.gz \
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
	file "${reads[0].baseName}.fq.gz"
	file "${reads[1].baseName}.fq.gz"
	file "${reads[0].baseName}.human.fq.gz"

	"""
	bbmap.sh \
		in1=${reads[0]} \
		in2=${reads[1]} \
		outu1=${reads[0].baseName}.fq.gz \
		outu2=${reads[1].baseName}.fq.gz \
		outm=${reads[0].baseName}.human.fq.gz \
		path=${params.human_seq} \
		minid=0.95 \
		maxindel=3 \
		minhits=2 \
		bandwithratio=0.16 \
		bandwith=12 \
		quickmatch \
		fast \
		qtrim=rl \
		trimq=10 \
		untrim \
	"""
}
