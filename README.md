# demultiplextagseq

## Introduction

**demultiplextagseq** is a bioinformatics pipeline to demultiplex 3'-tag-seq data. It performs the following steps:

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) multiplexed sample (paired-end)
2. Extract barcodes [`(UMI-tools extract)`](https://umi-tools.readthedocs.io/en/latest/reference/extract.html)
3. Demultiplex samples ([`fqtk`](https://github.com/fulcrumgenomics/fqtk))
4. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) demultiplex samples (single-end)
5. Plot selected quality metrics across 96-well plate
6. Combine and visualize ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
SAMPLEID1_PAIRED_END,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLEID2_PAIRED_END,HEG587A1_S1_L002_R1_001.fastq.gz,HEG587A1_S1_L002_R2_001.fastq.gz
```

Each row represents a pair of fastq files (paired end). This should come from the 3'tag-seq protocol. In many cases there will be only one row, but if several plates have been prepared, they can be demultiplexed in parallel.

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run berger demultiplextagseq \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Parameters

The pipeline takes the following optional pipeline specific paramters:

| Parameter     | Description                                                                                                                  |
| ------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| --split_fastq | Should the multiplexed input fastq file(s) be split<br />into smaller chuncks to speed up pipleline? <br />(Default: true) |
| --split_size  | Size of the fastq chunks.<br />(Default 50.000.000)                                                                          |

In addition, nextflow specific parameters can be used, such as: -**resume**, **-w** etc. See [Nextflow Docs](https://www.nextflow.io/docs/latest/reference/cli.html#run) for more details.

## Output

The pipeline outputs a folder named after the -outdir parameter. This folder contain the following subfolders:

* fastqc: fastqc report for multiplexed and demultiplexed samples
* fqtk: **demultiplexed fastq files to be used for analysis** named by well
* multiqc: **multiqc report** summarizing QC and methods
* pipeline_info: info about the pipeline run
* umitools: umitool logs

## Credits

demultiplextagseq was written by Elin Axelsson-Ekker, strongly inspired by Yoav Voicheck's original python script

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
