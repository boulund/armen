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
    .into {input_reads_qa}


process qa_reads {
    tag {pair_id}
    publishDir "${params.outdir}/qa_reads", mode: 'copy'

    input:
    set pair_id, file(reads) from input_reads_qa

    output:
	set pair_id, file("${pair_id}_{1,2}.trimmed.fq.gz") into input_remove_human
    file "${pair_id}.stats.txt.gz"
    file "${pair_id}.bhist.txt.gz"
    file "${pair_id}.qhist.txt.gz"
    file "${pair_id}.qchist.txt.gz"
    file "${pair_id}.aqhist.txt.gz"
    file "${pair_id}.bqhist.txt.gz"
    file "${pair_id}.lhist.txt.gz"
    file "${pair_id}.gchist.txt.gz"

    """
    bbduk.sh \
        in1=${reads[0]} \
        in2=${reads[1]} \
        ref=${params.bbduk_ref} \
        out1=${pair_id}_1.trimmed.fq.gz \
        out2=${pair_id}_2.trimmed.fq.gz \
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
	file "${pair_id}_1.fq.gz"
	file "${pair_id}_2.fq.gz"
	file "${pair_id}.human.fq.gz"
	file "${pair_id}.statsfile.txt.gz"

	"""
	bbmap.sh \
		in1=${reads[0]} \
		in2=${reads[1]} \
		path=${params.human_db} \
		outu1=${pair_id}_1.fq.gz \
		outu2=${pair_id}_2.fq.gz \
		outm=${pair_id}.human.fq.gz \
		statsfile=${pair_id}.statsfile.txt.gz \
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
