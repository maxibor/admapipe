<img src="img/logo.png" height="150">

* * *

## Introduction

Admapipe (Ancient DNA Metagenomics Analysis PIPEline) is a [Nextflow](https://nextflow.io) pipeline designed to reliably identify species in an ancient DNA metagenomics Illumina dataset. It makes use of an ensemblist approach and uses four different metagenomics softwares.

## Dependancies

### Softwares

-   [Nextflow > 0.30](https://www.nextflow.io/)
-   [Conda](https://conda.io/miniconda.html)
-   [MALT](http://ab.inf.uni-tuebingen.de/software/malt/)
-   [BASTA](https://github.com/timkahlke/BASTA)

### Databases

-   [Kraken database](https://ccb.jhu.edu/software/kraken/)
-   [Blast database](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
-   MALT database: to [build](http://ab.inf.uni-tuebingen.de/data/software/malt/download/manual.pdf) from reference sequence file.
-   [BASTA database](https://github.com/timkahlke/BASTA/wiki/2.-Initial-Setup)

## The Pipeline

**Pipeline overview:**

-   **1:**     FASTQC
-   **2:**     AdapterRemoval: Adapter trimming, quality filtering, and read merging
-   **3:**     Megahit: Metagenome de novo assembly
-   **4:**     Building bowtie index with contigs
-   **5:**     Aligning reads back on the contigs
-   **6:**     Indexing BAM alignment
-   **7:**     Converting BAM to SAM
-   **8:**     Position specific coverage for contigs
-   **9.1:**   Contig filtering on 10th percentile coverage
-   **9.2:**   Fasta file filtering with filtered contigs on coverage
-   **10:**    Fasta file filtering on contig length
-   **11.1:**  Kraken metagenome taxonomic classification
-   **11.2:**  Kraken report generation
-   **12:**    Metaphlan metagenome analyis
-   **13.1:**  MegaBlast contig analysis
-   **13.2:**  LCA with BASTA from MegaBlast result
-   **14.1:**  MALT metagenome analysis
-   **14.2:**  Conversion of MALT results to standard blast output file
-   **14.3:**  LCA with BASTA from MALT result
-   **15:**    Result summary
-   **\*:**     MultiQC

## Documentation

    $ nextflow run maxibor/admapipe --help

    N E X T F L O W  ~  version 0.30.1
    Launching `main.nf` [modest_kalam] - revision: 513b4a15b9

    =========================================
    Admapipe: Ancient DNA Metagenomics Analysis PIPEline
    Homepage / Documentation: <https://github.com/maxibor/admapipe>
    Author: Maxime Borry <mailto:maxime.borry@gmail.com>
    Version 0.2

    # Last updated on June 13th, 2018

    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run maxibor/admapipe --reads '\*\_R{1,2}.fastq.gz'
    Mandatory arguments:
    --reads                       Path to input data (must be surrounded with quotes)

    Options:
    --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to 64
    --minimum_coverage            Specifies the minimum 10th percentile coverage to consider a contig. Defaults to 2
    --min_length                  Specifies the minimum length of a contig to be considered. Defaults to 300
    --trimmingCPU                 Specifies the number of CPU used to trimming/cleaning by AdapterRemoval. Defaults to 16
    --megahitCPU                  Specifies the number of CPU used for Assembly by megahit. Defaults to 16
    --bowtieCPU                   Specifies the number of CPU used by bowtie2 aligner. Defaults to 20
    --krakenCPU                   Specifies the number of CPU used by Kraken taxonomic classifier. Defaults to 20
    --metaphlanCPU                Specifies the number of CPU used by Metaphlan. Defaults to 20
    --megablastCPU                Specifies the number of CPU used by MegaBlast. Defaults to 20
    --maltCPU                     Specifies the number of CPU used by MALT taxonomic classifier. Default to 16

    References: (files and directories must exist if used)
    --krakendb                    Path to Kraken index database. Defaults to /home/dist/maxime.borry/db/kraken/minikraken_20171101_8GB_dustmasked
    --blastdb                     Path to Blast index database. Defaults to /home/dist/maxime.borry/db/nt/blast/blast/nt
    --maltdb                      Path to MALT index database. Defaults to /home/dist/maxime.borry/db/malt/refseq_bacteria_step2
    --bastadb                     Path to BASTA (LCA) database. Defaults to /home/dist/maxime.borry/BASTA/taxonomy

    Other options:
    --results                     Name of result directory. Defaults to ./admapipe_results
    --help  --h                   Shows this help page
