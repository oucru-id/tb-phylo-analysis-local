nextflow.enable.dsl = 2

process VERSIONS {
    publishDir "${params.results_dir}", mode: 'copy'

    output:
    path "software_versions.yml"

    script:
    """
    echo "pipeline:" > software_versions.yml
    echo "  name: phylo_tb_analysis" >> software_versions.yml
    echo "  version: ${params.version}" >> software_versions.yml
    echo "  nextflow: $nextflow.version" >> software_versions.yml
    
    python3 $baseDir/scripts/get_versions.py >> software_versions.yml
    """
}