# CodeHandover

A comprehensive collection of code and analytical pipelines developed during my tenure at Bayer, focusing on fungicide mode-of-action characterization through transcriptomics analysis.

## Repository Structure

This repository is organized into four key projects:

### 1. BRIDGE_pipeline

A complete Nextflow pipeline for processing BRIDGE RNA-seq data that:
- Aligns reads to specified reference genomes
- Performs normalization using DESeq
- Generates visualization plots
- Automatically compiles results into PowerPoint presentations

#### Key Files

| File/Directory | Description |
|----------------|-------------|
| `BarcodeSetups/` | Barcode information files |
| `Configs/` | Pipeline configuration files |
| `BRIDGE.sh` | Wrapper bash script for the Nextflow pipeline |
| `bridge_local.nf` | Nextflow pipeline containing all analysis steps |
| `bridge_library.py` | Core Python library with classes and functions |
| `differential_expression_analysis.R` | R script for DE analysis |
| `make_powerpoint.py` | Creates presentations from analysis results |
| `setup_R.R` | Sets up required R packages |
| `setup_bash.sh` | Configures the Ubuntu environment |

#### Configuration

Before running any analysis, you must create a config file with these parameters:

| Parameter | Description |
|-----------|-------------|
| `accession` | Name of your data folder under Experiments |
| `plate_format` | Number of wells (typically 96 or 384) |
| `meta_data` | Path to metadata spreadsheet |
| `read_pooling` | Substrings to group fastq files (e.g., "plate1,plate2") |
| `treatment_sheets` | Treatment information for each read pool |
| `barcode_setup` | Location of the barcode setup file |
| `gtf` | Genome annotation file path |
| `genome` | Reference genome fasta file path |
| `genome_dir` | Directory containing STAR genome index |
| `trim_mode` | Trimming strictness (1-4, with 4 being strictest) |
| `email` | Notification email address |

#### Usage

```bash
./BRIDGE.sh Configs/<your_config_file>
```

#### Pipeline Workflow

1. **Data Download**: Retrieves read data and experiment metadata from the datapond
2. **Quality Control & Trimming**: Processes reads based on selected trim_mode:
    - Mode 1: Basic read length trimming with cutadapt
    - Mode 2: TrimGalore
    - Mode 3: TrimGalore followed by FastP (default)
    - Mode 4: TrimGalore + Trimmomatic + FastP
3. **Genome Alignment & Demultiplexing**: Uses STARsolo to align reads and demultiplex libraries
4. **Summary & Results**: Cleans up intermediate files and generates a PowerPoint summary

#### Outputs

- Main results: `<experiment_name>_BRIDGE_summary.pptx`
- Additional expression plots: `DeSeq/Plots/`
- Quality assessment: `QC/`

### 2. accessory_BRIDGE_scripts

Specialized analytical notebooks extending the primary pipeline:

| Notebook | Description |
|----------|-------------|
| `clustering_analysis_for_publication_normalised_deseq.ipynb` | Compares Alithea vs bulk RNA-seq; MoA clustering; PCA/UMAP; Mahalanobis distances |
| `shap_analysis.ipynb` | SHAP values for gene expression marker identification |
| `go_analysis.ipynb` | Gene Ontology enrichment analysis |
| `fill_genestack_table.ipynb` | Converts matrices to Genestack format |
| `plotting4tabea.ipynb` | Alithea sequencing parameter visualization |
| `differential_expression_analysis.ipynb` | DEG analysis with multiple testing correction |
| `pathway_enrichment.ipynb` | GSEA and pathway analysis |

### 3. multiomics

Analysis of paired cell painting and transcriptomic datasets:

- `picassoXbridge.ipynb`: Multi-omics analyses including MoFA
- `multi_omics_explore.ipynb`: Exploratory multi-omics analysis of basic parameters

### 4. cell_painting

Cell painting data processing with a focus on quality control:

- `prep_files.ipynb`: Preparation for QC analysis
- `qc_cell_painting.ipynb`: Analysis of batch effects (plate position, timepoints, lenses)
- `ATPase.ipynb`: Statistical analysis of CellProfiler columns related to ATPase MoA
- `MaMel.ipynb`: Statistical analysis of CellProfiler columns related to MaMel MoA

## Transcriptomics Data Structure

1. **Count matrices**: Genes (columns) × Samples (rows)
2. **Metadata columns**: Prefixed with "Metadata_", including:
    - Compound identifiers (`Metadata_treatments`)
    - Concentration (`Metadata_concentration`)
    - Time points (`Metadata_Timepoints`)
    - Mode of action (`Metadata_Mode_of_Action`)
    - FRAC subgroups (`Metadata_FRAC_sub_group`)
    - Chemical groups (`Metadata_Name_Chemical_group`)
    - Biological replicates (`Metadata_replicate`)
    - Batch information (`Metadata_batch`)
    - Sequencing depth (`Metadata_sequencing_depth`)
    - Species information (`Metadata_species`)

## Related Work

For transcriptomics integration with MOTHER platform, see the [MOTHER repository](https://github.com/bayer-int/smol-cls-mother/tree/master), which includes:
- Multiple normalization techniques
- Data simulation capabilities
- Feature selection methods
- End-to-end ML workflows

Example notebook available [here](https://github.com/bayer-int/smol-cls-mother/blob/master/examples/notebooks/example_rna_preprocessing.ipynb).

## Dependencies

Key requirements (see `bridge_env.yml` for full list):
- pandas (≥1.3.0), numpy (≥1.20.0)
- scanpy (≥1.8.0), anndata (≥0.7.6)
- matplotlib (≥3.4.0), seaborn (≥0.11.0), plotly (≥5.0.0)
- scikit-learn (≥1.0.0), statsmodels (≥0.12.0)
- DESeq2 (via rpy2), umap-learn (≥0.5.1)
- shap (≥0.39.0), pptx (≥0.6.18)
- biopython (≥1.79), goatools (≥1.0.15)

## Contact

For questions or support, please contact [n.arning@gmx.de].
