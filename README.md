# TB FHIR Phylogeny Analysis Pipeline (Local)

This pipeline processes FHIR bundle JSON files containing Mycobacterium tuberculosis genomics data to generate SNP distance matrices, phylogenetic trees, transmission network visualizations. Please refer to our full documentation: https://tb-pipeline-docs.readthedocs.io/en/latest/index.html

## Features

- **FHIR Genomics Data Gateway:** Directly fetch patient FHIR genomic bundles from FHIR servers.
- **Alternative local FHIR Genomics input:** Support for local processing of FHIR Bundles.
- **Phylogenetic Analysis:** Generates SNP distance matrices and phylogenetic trees (Rectangular, Circular, Unrooted).
- **Transmission Network:** Interactive Graph visualization and statistical plots (histogram, heatmap, violin plot).
- **Federated Analysis:** Generates pooled coefficient input for phylogenetic federated analytics (global model).
- **Clinical Metadata Integration:** Extracts metadata (location, lineage, patient ID) from FHIR resources.

## Usage

### Requirements

- [Nextflow](https://www.nextflow.io/)
- Python 3.8+
- Python packages: `biopython`, `pandas`, `networkx`, `pyvis`, `matplotlib`, `seaborn`, `numpy`, `requests`

Install Python dependencies:
```bash
pip install biopython pandas networkx pyvis matplotlib seaborn numpy requests
```

### Run the Pipeline

```bash
nextflow run main.nf
```

### Input

- FHIR Bundle Genomics JSON files in `data/JSON/` (local) 
- FHIR Bundle Genomics from FHIR server (params.fhir_server_url and params.fhir_server_auth)
- Reference genome FASTA in `data/H37Rv.fasta` and genomic anchor in `data/anchor`
