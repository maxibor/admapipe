<img src="./img/logo.png" height="300">

* * *

IN DEVELOPMENT

## Documentation
```
betsy:admapipe Maxime$ nextflow run main.nf --help
N E X T F L O W  ~  version 0.29.1
Launching `main.nf` [stupefied_hopper] - revision: 3c445188f9

=========================================
 Admapipe: Ancient DNA Metagenomics Analysis PIPEline
 Homepage / Documentation: https://github.com/maxibor/admapipe
 Author: Maxime Borry <maxime.borry@gmail.com>
 Version 0.1
 Last updated on May 31th, 2018
=========================================
Usage:
The typical command for running the pipeline is as follows:
nextflow run maxibor/admapipe --reads '*_R{1,2}.fastq.gz'
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
  --help                        Shows this help page
  --h                           Shows this help page
```
