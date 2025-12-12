# GNU Make Pipelines Collection

Welcome! This repository contains a personal collection of small, modular pipelines built using **GNU Make**.

Each subfolder in this repository corresponds to a self-contained pipeline, typically structured around a specific task or analysis. These pipelines are designed to be simple, reproducible, and easy to extend or adapt.

## Structure

Each pipeline folder includes:

* A `Makefile` & `Makefile_toplevel` to orchestrate tasks and manage dependencies
* Required input/output structure and reference files (on request)
* Singularity containers (on request) for full reproducibility
* Documentation or usage notes specific to the pipeline

## Example

The `facets_pipeline/` folder is an example of a complete setup for somatic copy number analysis using the FACETS algorithm. See its README for detailed usage instructions.

## Requirements

* GNU Make (v4+ recommended)
* Singularity (in order to use containerized tools)
* Bash and standard UNIX tools

## Usage

To run a specific pipeline:

```bash
cp your_pipeline_directory/Makefile_toplevel Makefile
make snp-pileup NORMAL TUMOR
```

Refer to the individual folder's README for detailed setup and run instructions.
