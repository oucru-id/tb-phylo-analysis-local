nextflow.enable.dsl = 2

process FHIR_ANALYSIS {
    publishDir "${params.results_dir}/phylo", mode: 'copy'
    
    input:
    path fhir_files
    path reference

    output:
    path "distance_matrix.tsv", emit: matrix
    path "phylo_tree.nwk",      emit: tree
    path "metadata.tsv",        emit: metadata
    path "consensus.fasta",     emit: fasta 

    script:
    """
    python3 $baseDir/scripts/fhir_phylo.py \\
        --inputs ${fhir_files} \\
        --reference ${reference}
    """
}

workflow PHYLO_ANALYSIS {
    take:
    fhir_files
    reference

    main:
    FHIR_ANALYSIS(fhir_files.collect(), reference)

    emit:
    matrix   = FHIR_ANALYSIS.out.matrix
    tree     = FHIR_ANALYSIS.out.tree
    metadata = FHIR_ANALYSIS.out.metadata
    fasta    = FHIR_ANALYSIS.out.fasta 
}
