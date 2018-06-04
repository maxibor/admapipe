#!/usr/bin/env nextflow

/*
========================================================================================
                                      Admapipe
========================================================================================
 Admapipe: Ancient DNA Metagenomics Analysis PIPEline
#### Homepage / Documentation
https://github.com/maxibor/admapipe
#### Authors
 Maxime Borry <maxime.borry@gmail.com>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Pipeline overview:
 - 1:     FASTQC
 - 2:     AdapterRemoval: Adapter trimming, quality filtering, and read merging
 - 3:     Megahit: Metagenome de novo assembly
 - 4:     Building bowtie index with contigs
 - 5:     Aligning reads back on the contigs
 - 6:     Indexing BAM alignment
 - 7:     Converting BAM to SAM
 - 8:     Position specific coverage for contigs
 - 9.1:   Contig filtering on 10th percentile coverage
 - 9.2:   SAM file filtering with filtered contigs on coverage
 - 9.3:   Fasta file filtering with filtered contigs on coverage
 - 10:    Fasta file filtering on contig length
 - 11.1:  Kraken metagenome taxonomic classification
 - 11.2:  Kraken report generation
 - 12:    Metaphlan metagenome analyis
 - 13.1:  MegaBlast contig analysis
 - 13.2:  LCA with BASTA from MegaBlast result
 - 14.1:  MALT metagenome analysis
 - 14.2:  Conversion of MALT results to standard blast output file
 - 14.3:  LCA with BASTA from MALT result
 - 15:    Result summary
 - 16:    Plot UpSet
 - *:     MultiQC

 ----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     Admapipe: Ancient DNA Metagenomics Analysis PIPEline
     Homepage / Documentation: https://github.com/maxibor/admapipe
     Author: Maxime Borry <maxime.borry@gmail.com>
     Version ${version}
     Last updated on ${version_date}
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run maxibor/admapipe --reads '*_R{1,2}.fastq.gz'
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)

    Options:
      --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to ${params.phred}
      --minimum_coverage            Specifies the minimum 10th percentile coverage to consider a contig. Defaults to ${params.min_coverage}
      --min_length                  Specifies the minimum length of a contig to be considered. Defaults to ${params.min_length}
      --trimmingCPU                 Specifies the number of CPU used to trimming/cleaning by AdapterRemoval. Defaults to ${params.trimmingCPU}
      --megahitCPU                  Specifies the number of CPU used for Assembly by megahit. Defaults to ${params.megahitCPU}
      --bowtieCPU                   Specifies the number of CPU used by bowtie2 aligner. Defaults to ${params.bowtieCPU}
      --krakenCPU                   Specifies the number of CPU used by Kraken taxonomic classifier. Defaults to ${params.krakenCPU}
      --metaphlanCPU                Specifies the number of CPU used by Metaphlan. Defaults to ${params.metaphlanCPU}
      --megablastCPU                Specifies the number of CPU used by MegaBlast. Defaults to ${params.megablastCPU}
      --maltCPU                     Specifies the number of CPU used by MALT taxonomic classifier. Default to ${params.maltCPU}

    References: (files and directories must exist if used)
      --krakendb                    Path to Kraken index database. Defaults to ${params.krakendb}
      --blastdb                     Path to Blast index database. Defaults to ${params.blastdb}
      --maltdb                      Path to MALT index database. Defaults to ${params.maltdb}
      --bastadb                     Path to BASTA (LCA) database. Defaults to ${params.bastadb}

    Other options:
      --results                     Name of result directory. Defaults to ${params.results}
      --help  --h                   Shows this help page

    """.stripIndent()
}

version = "0.1"
version_date = "June 1st, 2018"

params.results = "./admapipe_results"
params.reads = "/home/maxime/Documents/data/mock_metagenome_deaminated/*.{1,2}.fastq"



params.phred = 64
params.min_coverage = 2
params.min_length = 300

params.krakendb = "/home/dist/maxime.borry/db/kraken/minikraken_20171101_8GB_dustmasked"

params.blastdb = "/home/dist/maxime.borry/db/nt/blast/blast/nt"

params.bastadb = "/home/dist/maxime.borry/BASTA/taxonomy"

params.maltdb = "/home/dist/maxime.borry/db/malt/refseq_bacteria_step2"


params.trimmingCPU = 16
params.megahitCPU = 16
params.bowtieCPU = 20
params.krakenCPU = 20
params.metaphlanCPU = 20
params.megablastCPU = 20
params.maltCPU = 16



basta = "/home/dist/maxime.borry/BASTA/bin/basta"
malt = "/home/dist/maxime.borry/malt/malt-run"

multiqc_conf = "$baseDir/conf/.multiqc_config.yaml"

// Show help message
params.help = false
params.h = false
if (params.help || params.h){
    helpMessage()
    exit 0
}

//Logging parameters
log.info "================================================================"
log.info " Admapipe: Ancient DNA Metagenomics Analysis PIPEline"
log.info " Homepage / Documentation: https://github.com/maxibor/admapipe"
log.info " Author: Maxime Borry <maxime.borry@gmail.com>"
log.info " Version ${version}"
log.info " Last updated on ${version_date}"
log.info "================================================================"
def summary = [:]
summary['Reads'] = params.reads
summary["Result directory"] = params.results
summary['phred quality'] = params.phred
summary["Minimum 10th percentile contig coverage"] = params.min_coverage
summary["Minimum contig length"] = params.min_length
summary["Kraken database"] = params.krakendb
summary["Blast database"] = params.blastdb
summary["MALT database"] = params.maltdb
summary["BASTA (LCA) database"] = params.bastadb
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

// Creating reads channel
Channel
    .fromFilePairs( params.reads, size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n" }
	  .into { reads_fastqc; reads_to_trim }

// Step 1 - FASTQC
process fastqc {
    tag "$name"

    label 'normal'

    cpus 1

    publishDir "${params.results}/fastqc", mode: 'copy'

    errorStrategy 'ignore'

    beforeScript "set +u; source activate fastqc"
    afterScript "source deactivate"

    input:
        set val(name), file(reads) from reads_fastqc

    output:
        file '*_fastqc.{zip,html}' into fastqc_results
    script:
        """
        fastqc -q $reads
        """
}

// Step 2 - AdapterRemoval: Adapter trimming, quality filtering, and read merging
process adapter_removal_ancient_dna_PE {
    tag "$name"

    label 'normal'

    cpus params.trimmingCPU

    publishDir "${params.results}/trimmed", mode: 'copy'

    beforeScript "set +u; source activate adapterremoval"
    afterScript "source deactivate"


    input:
        set val(name), file(reads) from reads_to_trim

    output:
        set val(name), file('*.collapsed.fastq') into trimmed_reads_assembly, trimmed_reads_mapping, trimmed_reads_kraken, trimmed_reads_metaphlan, trimmed_reads_malt
        set val(name), file("*.settings") into adapter_removal_results

    script:
        out1 = name+".pair1.discarded.fastq"
        out2 = name+".pair2.discarded.fastq"
        col_out = name+".collapsed.fastq"
        """
        AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --collapse --output1 $out1 --output2 $out2 --outputcollapsed $col_out --threads ${task.cpus} --qualitybase ${params.phred}
        """

}

// Step 3 - Megahit: Metagenome de novo assembly
process megahit_assembly{
    tag "$name"

    // label 'big_mem'
    label 'normal'

    cpus params.megahitCPU

    publishDir "${params.results}/assembly", mode: 'copy'

    beforeScript "set +u; source activate megahit"
    afterScript "source deactivate"

    input:
        set val(name), file(merged) from trimmed_reads_assembly
    output:
        set val(name), file ("megahit_out/*.contigs.fa") into contigs_mapping, contigs_filter
    script:
        """
        megahit -r $merged -t ${task.cpus} --out-prefix $name
        """
}

// Step 4 - Building bowtie index with contigs
process bowtie_index_contigs{
    tag "$name"

    label 'normal'

    cpus params.bowtieCPU

    publishDir "${params.results}/bowtie_index", mode: 'copy'

    beforeScript "set +u; source activate bowtie"
    afterScript "source deactivate"

    input:
        set val(name), file(contig) from contigs_mapping
    output:
        set val(name), file("*.bt2") into bt_index
    script:
        """
        bowtie2-build --threads ${task.cpus} $contig $name
        """
}

// Step 5 - Aligning reads back on the contigs
process align_reads_to_contigs{
    tag "$name"

    label 'normal'

    cpus params.bowtieCPU

    publishDir "${params.results}/alignment", mode: 'copy'

    beforeScript "set +u; source activate bowtie"
    afterScript "source deactivate"

    input:
        set val(name), file(reads) from trimmed_reads_mapping
        set val(name), file(contig) from bt_index
    output:
        set val(name), file("*.sorted.bam") into alignment_to_index, alignment_to_coverage, alignment_to_sam
    script:
        outfile = name+".sorted.bam"
        """
        bowtie2 -x $name -U $reads --very-fast --threads ${task.cpus} | samtools view -S -b -F 4 - | samtools sort - > $outfile
        """

}

// Step 6 - Indexing BAM alignment
process bam_index {
    tag "$name"

    label 'normal'

    cpus 1

    publishDir "${params.results}/alignment", mode: 'copy'

    beforeScript "set +u; source activate samtools"
    afterScript "source deactivate"

    input:
        set val(name), file(bam) from alignment_to_index
    output:
        set val(name), file("*.bam.bai")
    script:
        """
        samtools index $bam
        """
}

// Step 7 - Converting BAM to SAM
process bam2sam {
    tag "$name"

    label 'normal'

    cpus 1

    beforeScript "set +u; source activate samtools"
    afterScript "source deactivate"

    input:
        set val(name), file(bam) from alignment_to_sam
    output:
        set val(name), file("*.sam") into alignment_sam
    script:
        outfile = name+".sam"
        """
        samtools view -sB -F 4 $bam > $outfile
        """

}

// Step 8 - Position specific coverage for contigs
process bedtools_genomecov {
    tag "$name"

    label 'normal'

    cpus 1

    publishDir "${params.results}/coverage", mode: 'copy'

    beforeScript "set +u; source activate bedtools"
    afterScript "source deactivate"

    input:
        set val(name), file(bam) from alignment_to_coverage
    output:
        set val(name), file("*.bed") into bedfile
    script:
        outfile = name+".bed"
        """
        bedtools genomecov -ibam $bam -d > $outfile
        """
}

// Step 9.1 - Contig filtering on 10th percentile coverage
process filter_contigs {
    tag "$name"

    label 'normal'

    cpus 1

    publishDir "${params.results}/coverage", mode: 'copy', pattern: "*.contig_coverage.txt"

    input:
        set val(name), file(bed) from bedfile
    output:
        set val(name), file("*.contig_coverage.txt") into contig_coverage
        set val(name), file("*.filtered_contigs_list.txt") into filtered_contigs_sam, filtered_contigs_fasta
    script:
        outbed = name+".contig_coverage.txt"
        contig_list = name+".filtered_contigs_list.txt"
        """
        poscov2featurecov $bed -m ${params.min_coverage} > $outbed
        awk '{print \$1}' $outbed > $contig_list
        """
}

// Step 9.2 - SAM file filtering with filtered contigs on coverage
/*
process filter_sam_coverage {
    tag "$name"

    label 'normal'

    cpus 1

    publishDir "${params.results}/alignment", mode: 'copy'

    errorStrategy 'ignore'

    input:
        set val(name), file(sam) from alignment_sam
        set val(name), file(contigs) from filtered_contigs_sam
    output:
        set val(name), file("*.filtered.sam") into filtered_alignment
    script:
        outfile = name+".filtered.sam"
        """
        grep -f $contigs -w $sam > $outfile
        """
}
*/

// Step 9.3 - Fasta file filtering with filtered contigs on coverage
process filter_fasta_coverage{
    tag "$name"

    label 'normal'

    cpus 1

    input:
        set val(name), file(contig) from contigs_filter
        set val(name), file(filt_contig) from filtered_contigs_fasta
    output:
        set val(name), file("*.filtered_contigs.fa") into filtered_fasta_coverage
    script:
        outfile = name+".filtered_contigs.fa"
        """
        filterFastaByName $contig $filt_contig -o $outfile
        """
}

// Step 10 - Fasta file filtering on contig length
process filter_fasta_length {
    tag "$name"

    label 'normal'

    cpus 1

    publishDir "${params.results}/assembly", mode: 'copy'

    input:
        set val(name), file(fasta) from filtered_fasta_coverage
    output:
        set val(name), file("*.filtered.fa") into fasta2megablast
    script:
        """
        filterFastaByLength -min ${params.min_length} $fasta
        """
}

// Step 11.1 - Kraken metagenome taxonomic classification
process kraken {
    tag "$name"

    label 'normal'

    cpus params.krakenCPU

    publishDir "${params.results}/kraken", mode: 'copy'

    beforeScript "set +u; source activate kraken"
    afterScript "source deactivate kraken"

    input:
        set val(name), file(fasta) from trimmed_reads_kraken
    output:
        set val(name), file("*.kraken") into kraken_output
    script:
        outfile = name+".kraken"
        """
        kraken --db ${params.krakendb} --threads ${task.cpus} --fastq-input $fasta > $outfile
        """
}


// Step 11.2 - Kraken report generation
process kraken_report {
    tag "$name"

    label 'normal'

    cpus 1

    publishDir "${params.results}/kraken", mode: 'copy'

    beforeScript "set +u; source activate kraken"
    afterScript "source deactivate kraken"

    input:
        set val(name), file(kraken_out) from kraken_output
    output:
        set val(name), file("*.kraken.report") into kraken_res
    script:
        outfile = name+".kraken.report"
        """
        kraken-report --db ${params.krakendb} $kraken_out > $outfile
        """
}

// Step 12 - Metaphlan metagenome analyis
process metaphlan {
    tag "$name"

    label 'normal'

    cpus params.metaphlanCPU

    publishDir "${params.results}/kraken", mode: 'copy'

    beforeScript "set +u; source activate metaphlan"
    afterScript "source deactivate metaphlan"

    input:
        set val(name), file(fasta) from trimmed_reads_metaphlan
    output:
        set val(name), file("*.metaphlan.out") into metaphlan_res
    script:
        outfile = name+".metaphlan.out"
        """
        metaphlan2.py --nproc ${params.metaphlanCPU} --input_type fastq --output_file $outfile $fasta
        """
}

// Step 13.1 - MegaBlast contig analysis
process megablast {
    tag "$name"

    label 'normal'

    cpus params.megablastCPU

    publishDir "${params.results}/megablast", mode: 'copy'

    beforeScript "set +u; source activate blast"
    afterScript "source deactivate"

    input:
        set val(name), file(fasta) from fasta2megablast
    output:
        set val(name), file("*.blast") into blast_output
    script:
        outfile = name+".blast"
        """
        blastn -task megablast -num_threads ${task.cpus} -outfmt 6 -db ${params.blastdb} -query $fasta -out $outfile
        """
}

// Step 13.2 - LCA with BASTA from MegaBlast result
process basta_from_blast {
    tag "$name"

    errorStrategy 'ignore'

    label 'normal'

    cpus 1

    publishDir "${params.results}/megablast", mode: 'copy'

    beforeScript "set +u; source activate basta"
    afterScript "source deactivate"

    input:
        set val(name), file(blast_out) from blast_output
    output:
        set val(name), file("*.blast.basta") into blast_res
    script:
        outfile = name+".blast.basta"
        basta_max = 30
        basta_id = 97
        basta_min = 1
        basta_len = 300
        """
        $basta sequence -t all -n $basta_max -m $basta_min -i $basta_id -l $basta_len -x True -d ${params.bastadb} $blast_out $outfile  gb
        """
}

// Step 14.1 - MALT metagenome analysis
process malt {
    tag "$name"

    label 'very_big_mem'

    cpus params.maltCPU

    publishDir "${params.results}/malt", mode: 'copy'

    input:
        set val(name), file(fasta) from trimmed_reads_malt
    output:
        set val(name), file("*.aligned.malt") into malt_output
        set val(name), file("*.rma") into malt_rma
    script:
        outfile = name+".aligned.malt"
        outrma = name+".rma"
        """
        $malt --mode BlastN --alignments $outfile --gzipAlignments false --format Tab --index ${params.maltdb} --numThreads ${task.cpus} --inFile $fasta --output $outrma
        """
}

// Step 14.2 - Conversion of MALT results to standard blast output file
process malt_convert {
    tag "$name"

    errorStrategy 'ignore'

    label 'normal'

    cpus 1

    publishDir "${params.results}/malt", mode: 'copy'

    input:
        set val(name), file(malt_out) from malt_output
    output:
        set val(name), file("*.blast_converted.malt") into malt_converted
    script:
        outfile = name+".blast_converted.malt"
        """
        sed -e 's/|tax|\\([0-9]\\+\\)|//g' $malt_out > $outfile
        """
}

// Step 14.3 - LCA with BASTA from MALT result
process basta_from_malt {
    tag "$name"

    errorStrategy 'ignore'

    label 'normal'

    cpus 1

    publishDir "${params.results}/malt", mode: 'copy'

    beforeScript "set +u; source activate basta"
    afterScript "source deactivate"

    input:
        set val(name), file(malt_out) from malt_converted
    output:
        set val(name), file("*.malt.basta") into malt_res
    script:
        outfile = name+".malt.basta"
        basta_max = 30
        basta_id = 97
        basta_min = 1
        basta_len = 100
        """
        $basta sequence -t all -n $basta_max -m $basta_min -i $basta_id -l $basta_len -x True -d ${params.bastadb} $malt_out $outfile gb
        """
}

// Step 15 - Result summary
process summarize_results {
    tag "$name"

    errorStrategy 'ignore'

    label 'normal'

    cpus 1

    publishDir "${params.results}/summary", mode: 'copy'

    input:
        set val(name), file(metaphlan_reads) from metaphlan_res
        set val(name), file(kraken_reads) from kraken_res
        set val(name), file(malt_reads) from malt_res
        set val(name), file(blast_contigs) from blast_res
    output:
        set val(name), file("*.csv") into summary_result

    script:
        outfile = name+".summary.csv"
        """
        summarize_all_methods -mr $metaphlan_reads -krr $kraken_reads -mar $malt_reads -bc $blast_contigs -o $outfile
        """

}

// Step 16 - Plot UpSet

process upset {
    tag "$name"

    errorStrategy 'ignore'

    label 'normal'

    cpus 1

    publishDir "${params.results}/summary", mode: 'copy'

    beforeScript "set +u; source activate upset"
    afterScript "source deactivate"

    input:
        set val(name), file(table) from summary_result
    output:
        set val(name), file("*.upsetplot.png") into upset_plot
    script:
        outfile = name+".upsetplot.png"
        """
        plot_upset $table -o $outfile
        """
}

// MultiQC
process multiqc {
    tag "$prefix"

    errorStrategy 'ignore'

    label 'normal'

    cpus 1

    publishDir "${params.results}/MultiQC", mode: 'copy'

    beforeScript "set +u; source activate multiqc"
    afterScript "source deactivate"


    input:
        file (fastqc:'fastqc/*') from fastqc_results.collect()
        file ('adapter_removal/*') from adapter_removal_results.collect()

    output:
        file '*multiqc_report.html' into multiqc_report
        file '*_data' into multiqc_data

    script:
        prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
        """
        multiqc -f -d fastqc adapter_removal -c $multiqc_conf
        """
}
