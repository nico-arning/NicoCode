# CodeHandover

A comprehensive collection of code and analytical pipelines developed during my tenure at Bayer.

## Repository Structure

This repository is organized into several key projects:

### 1. BRIDGE_pipeline
A complete pipeline for processing BRIDGE data that:
- Aligns reads to specified reference genomes
- Performs normalization using DESeq
- Generates visualization plots
- Automatically compiles results into PowerPoint presentations

### 2. accessory_BRIDGE_scripts
Specialized analytical notebooks extending beyond the primary pipeline:
- MoA disentanglement analysis
- Bulk vs. Alithea data comparison
- SHAP value analysis for feature importance
- Utility scripts for converting internal metadata to Genestack format

### 3. multiomics
Analysis of paired cell painting and transcriptomic datasets:
- Basic exploratory data analysis
- Machine learning implementations
- Multi-Omics Factor Analysis (MoFA) workflows

### 4. cell_painting
Cell painting data processing and analysis:
- Quality control procedures
- Data preparation and formatting workflows

## Related Work

For transcriptomics integration into the MOTHER platform, please refer to the [MOTHER repository](https://github.com/bayer-int/smol-cls-mother/tree/master). This codebase includes:
- Multiple normalization techniques
- Normalization object export functionality
- RNA data simulation
- Feature selection methods
- End-to-end machine learning workflows for transcriptomics data

## Usage

Each project folder contains dedicated documentation and example notebooks demonstrating typical workflows.
