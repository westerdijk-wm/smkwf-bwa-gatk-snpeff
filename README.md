# Snakemake workflow: smkwf-bwa-gatk-snpeff

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/westerdijk-wm/smkwf-bwa-gatk-snpeff/actions/workflows/main.yaml/badge.svg?branch=main)](https://github.com/westerdijk-wm/smkwf-bwa-gatk-snpeff/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/westerdijk-wm/smkwf-bwa-gatk-snpeff)

A Snakemake workflow for for mapping, variant-calling and annotating variants using SNPeff

- [Snakemake workflow: smkwf-bwa-gatk-snpeff](#snakemake-workflow-smkwf-bwa-gatk-snpeff)
  - [Usage](#usage)
  - [Deployment options](#deployment-options)
  - [Workflow profiles](#workflow-profiles)
  - [Authors](#authors)
  - [References](#references)
  - [TODO](#todo)

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/westerdijk-wm/smkwf-bwa-gatk-snpeff).
This includes a visualization of the workflow diagram and a table with all workflow parameters.

Detailed information about input data and workflow configuration can also be found in the [`config/README.md`](config/README.md).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository or its DOI.

## Deployment options

To run the workflow from command line, change the working directory.

```bash
cd path/to/smkwf-bwa-gatk-snpeff
```

Adjust options in the default config file `config/config.yaml`.
Before running the complete workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

To run the workflow with test files using **conda**:

```bash
snakemake --cores 2 --sdm conda --directory .test
```

To run the workflow with **apptainer** / **singularity**, add a link to a container registry in the `Snakefile`, for example `container: "oras://ghcr.io/<user>/<repository>:<version>"` for Github's container registry.
Run the workflow with:

```bash
snakemake --cores 2 --sdm conda apptainer --directory .test
```

## Workflow profiles

The `profiles/` directory can contain any number of [workflow-specific profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) that users can choose from.
The [profiles `README.md`](profiles/README.md) provides more details.

## Authors

- Balazs Brankovics
  - [Westerdijk Fungal Biodiversity Institute - KNAW](https://wi.knaw.nl/Balazs_Brankovics)
  - [https://orcid.org/0000-0003-0536-7787](https://orcid.org/0000-0003-0536-7787)
- Manuel Leeuwerik
  - [Westerdijk Fungal Biodiversity Institute - KNAW](https://wi.knaw.nl/Manuel_Leeuwerik)
  - [https://orcid.org/0009-0004-9658-1117](https://orcid.org/0009-0004-9658-1117)

## References


## TODO

- [x] Replace `<owner>` and `<repo>` everywhere in the template with the correct user name/organization, and the repository name. The workflow will be automatically added to the [snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/index.html) once it is publicly available on Github.
- [x] Replace `<name>` with the workflow name (can be the same as `<repo>`).
- [x] Replace `<description>` with a description of what the workflow does.
- [x] Update the [deployment](#deployment-options), [authors](#authors) and [references](#references) sections.
- [x] Update the `README.md` badges. Add or remove badges for `conda`/`singularity`/`apptainer` usage depending on the workflow's [deployment](#deployment-options) options.
- Do not forget to also adjust the configuration-specific `config/README.md` file.
