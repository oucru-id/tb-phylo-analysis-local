nextflow.enable.dsl = 2

process FHIR_ANALYSIS {
    publishDir "${params.results_dir}/phylo", mode: 'copy'
    
    input:
    path fhir_files
    path reference
    path anchor_files

    output:
    path "distance_matrix.tsv", emit: matrix
    path "phylo_tree.nwk",      emit: tree
    path "metadata.tsv",        emit: metadata
    path "consensus.fasta",     emit: fasta 

    script:
    def anchor_arg = anchor_files ? "--anchors ${anchor_files}" : ""
    """
    python3 $baseDir/scripts/fhir_phylo.py \\
        --inputs ${fhir_files} \\
        --reference ${reference} \\
        $anchor_arg
    """
}

workflow PHYLO_ANALYSIS {
    take:
    fhir_files
    reference
    anchors

    main:
    FHIR_ANALYSIS(fhir_files.collect(), reference, anchors)

    emit:
    matrix   = FHIR_ANALYSIS.out.matrix
    tree     = FHIR_ANALYSIS.out.tree
    metadata = FHIR_ANALYSIS.out.metadata
    fasta    = FHIR_ANALYSIS.out.fasta 
}
