nextflow.enable.dsl = 2

process UPLOAD_MATRIX_TO_FHIR {
    publishDir "${params.results_dir}/phylo", mode: 'copy'

    input:
    path matrix
    val  url
    path token_file

    output:
    path "upload_result.json", emit: upload_result

    script:
    """
    python3 $baseDir/scripts/upload_distance_matrix.py \\
        --matrix ${matrix} \\
        --url "${url}" \\
        --token-file "${token_file}"
    """
}

workflow UPLOAD_FHIR {
    take:
    matrix
    url
    token_file

    main:
    UPLOAD_MATRIX_TO_FHIR(matrix, url, token_file)

    emit:
    upload_result = UPLOAD_MATRIX_TO_FHIR.out.upload_result
}
