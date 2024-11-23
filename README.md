

# catg/nanotax

[![GitHub Actions CI Status](https://github.com/catg/nanotax/actions/workflows/ci.yml/badge.svg)](https://github.com/catg/nanotax/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/catg/nanotax/actions/workflows/linting.yml/badge.svg)](https://github.com/catg/nanotax/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**catg/nanotax** is a bioinformatics pipeline for the analysis of 16S rRNA gene sequencing data obtained by Nanopore sequencing. It takes a samplesheet with POD5 or FastQ files and barcodes (optional) and groups (optional) as input and performs basecalling and demultiplexing, quality control (QC), taxonomic assignment with databases, functional prediction and alpha diversity metrics and produces tables and plots with all the results. All the process are optional except taxonomic assignment.

The pipeline then:
1. Basecalling and demultiplexing with ([`Dorado`](https://github.com/nanoporetech/dorado)).
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) and ([`NanoQ`](https://github.com/esteinig/nanoq)).
3. Quality and length filter with ([`NanoQ`](https://github.com/esteinig/nanoq)).
4. Sampling by quality with ([`Filtlong`](https://github.com/rrwick/Filtlong))
5. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))
6. Assigns taxonomy to reads using ([`MMSeqs2`](https://github.com/soedinglab/MMseqs2)) with ([`Genbank`](https://www.ncbi.nlm.nih.gov/refseq/targetedloci/16S_process/)) or ([`SILVA`](https://www.arb-silva.de/)) database.
7. (optionally) alpha diversity metrics with ([`Vegan`](https://cran.r-project.org/web/packages/vegan/vegan.pdf)) (R).
8. (optionally) functional prediction with ([`PICRUSt2`](https://github.com/picrust/picrust2)) .
9. (optionally) differential expression for functional prediction with([`LEfSe`](https://huttenhower.sph.harvard.edu/lefse/)). 

All plots and tables are generated using Python, through either ([`pandas`](https://github.com/pandas-dev/pandas)) or ([`polars`](https://github.com/pola-rs/polars)). 

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

The columns in the sample file will vary depending on the analyses to be performed. Each row represents a sample, specifying either its associated `FASTQ` file or the `barcode` used during sequencing.

By default, the pipeline performs quality control and taxonomic assignment. For these tasks, the sample file must include the columns "samples" and "fastq".

```csv
sample,fastq
sample_1,sample_1.fastq.gz
sample_2,sample_2.fastq.gz
```

If you wish to start from the basecalling step, the `fastq` column should be replaced with `barcode` , and the directory containing the POD5 files must be specified using the `--basecalling.pod5_dir`  parameter.
```csv
sample,barcode
sample_1,barcode01
sample_2,barcode02
```

To perform diversity analyses or differential expression of metabolic pathways, the sample file must include a groups column.
```csv
sample,fastq,group
sample_1,sample_1.fastq.gz,G1
sample_2,sample_2.fastq.gz,G2
```


Now, you can run the pipeline using:

```bash
nextflow run catg/nanotax \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

catg/nanotax was originally written by JacquelineAldridge.

<!-- We thank the following people for their extensive assistance in the development of this pipeline: -->

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use catg/nanotax for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
