# ChipSeqPipeline
Bash pipeline for Chip-Seq data processing and analysis.

Authors: 
  * Antonio Álvarez Gómez
  * Gabriel Ferreras Garrucho
  * Helena Victoria Cotán

## Introduction

The pipeline hereby presented was designed in order to ChIP-seq (Chromatin ImmunopreciPitation sequencing) data processing and analysis. This material is a task of the course "Omics technologies and Bioinformatics" given by Francisco J. Romero-Campero and Ignacio Pérez Hurtado de Mendoza in the Biochemistry degree at the University of Seville.

## Pipeline summary

1. **ChipSeqPipeline.sh**
   * Read parameters
   * Prepare workspace
   * Copy the data
   * Create a genome index ([`bowtie2-build`](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml))
   * Sample processing
 2. **SampleProcessing.sh**
    * Quality control ([`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
    * Map to genome ([`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml))
    * Generate sorted bam file ([`samtools`](http://www.htslib.org))
    * Peak determination
 3. **PeakDetermination.sh**
    * Peak determination ([`masc2 callpeak`](https://github.com/macs3-project/MACS))
    * Peak annotation
    * Homer-motifs-finding ([`findMotifsGenome.pl`](http://homer.ucsd.edu/homer/ngs/peakMotifs.html))
    * Experiment global analysis
 4. **PeakAnnotation.R**
    * Install packages
    * Read arguments
    * Read peak file
    * Definition of promoter region
    * Peak annotation
    * Convert annotation to target genes
 5. **ExperimentGlobalAnalysis.R**
    * Install packages
    * Read arguments
    * Group the replicates
    * Extraction of overlapping replicates
    * Gene ontology enrichment
    * KEGG pathway enrichment

## What does it do?

  * Quality control reads using fastQC 
  * Map reads to genome using Bowtie2
  * Call peaks for multiple combinations of samples using MACS2
  * Motifs-finding using findMotifsGenome.pl
  * A differential analysis using R

## What does it output?

  * QC reports
  * Genome index
  * Bam files
  * Peak annotation files
  * Differential analysis report

## Dependencies

For running this pipeline, you need to assure that the following softwares are installed in your computer:

  * [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  * [`SAMTOOLS`](https://sourceforge.net/projects/samtools/files/samtools/)
  * [`MACS2`](https://github.com/macs3-project/MACS)
  * [`R`](https://www.r-project.org)
    - [`BiocManager`](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)
    - [`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
    - [`TxDb.Athaliana.BioMart.plantsmart28`](https://bioconductor.org/packages/release/data/annotation/html/TxDb.Athaliana.BioMart.plantsmart28.html)
    - [`DO.db`](http://bioconductor.org/packages/release/data/annotation/html/DO.db.html)
    - [`org.At.tair.db`](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)
    - [`VennDiagram`](https://cran.r-project.org/web/packages/VennDiagram/VennDiagram.pdf)
    - [`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
    - [`pathview`](https://bioconductor.org/packages/release/bioc/html/pathview.html)
  * [`Homer`](http://homer.ucsd.edu/homer/download.html)

## Installation

This pipeline does not require any installation apart from those mentionated above in Dependencies section. Usage just does require download the repository from GitHub. If you do not have Git install, consult https://github.com/git-guides/install-git.

The steps you need to follow are these described:

  1. Open a terminal
  2. Go to the working directory you desire to save the repository
  3. Introduce the following command

      ```bash
      git clone https://github.com/gabriel-ferreras/ChipSeqPipeline.git
      ```

## Running pipeline
  
To run the pipeline on your experimental data, first enter the necessary parameters in the spreadsheet file (see following section), and then from the terminal type:

```sh
bash ChipSeqPipeline test_params.txt
```

### Input parameters

The pipeline requires a file as input to specify the samples and the design of the analysis. Within the repository you can find an example of parameters file (test_params.txt). The parameters you need to include are:

 * **General parameters**:
    - **analysis_name**: is the name you desire to give to the analysis. The folder created during the process will carry this name.
    - **number_samples**: the number of samples that are going to be processed. Count the pair Chip/Control just as one sample. Does not matter wether if the samples are paired or not.
    - **number_experiments**: the number of experiments that are going to be processed
    - **experimental_design**: specify the samples that correspond to each experiment. For example, if there are two experiments and each experiment has two samples, indicate it as "(1, 1, 2, 2 )"
 
 * **Directories**:
    - **working_directory**: the directory where all the data is going to be placed during the analysis
    - **installation_directory**: the directory where the repository is placed
    - **path_annotation**: the directory where the annotation is placed
    - **path_genome**: the directory where the genome is placed
    - **path_chip_1**: the directory where the chip samples are placed
    - **path_control_1**: the directory where the control samples are placed  
    
 * **Eligible parameters**:
    - **paired**: indicate whether the samples are from a paired-end sequencing (True/1) or a single-end sequencing (False/0).
       > [`Paired-end sequencing`](https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html) allows users to sequence both ends of a fragment and generate high-quality, alignable sequence data. It produces twice the bumber of reads in teh same time and detect insertion-deletion (indel) variants, which is not possible with single-read data.
       > [`Single-end sequencing`](https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html) involves sequencing DNA from only one end, delivering large volumes of high-quality data, rapidly and economically.
    - **broad**: indicate whether, in the peak analysis, you prefer to indentify broad peaks (True/1) or narrow peaks (False/0). 
       > [`MACS2`](https://github.com/macs3-project/MACS) has two different specific methods in searching for enrichment regions. [`Broad peaks`](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html ) or broad domains (i.e. histone modifications that cover entire gene bodies) or [`narrow peaks`](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html ) (i.e. a transcription factor binding). Narrow peaks are easier to detect as we are looking for regions that have higher amplitude and are easier to distinguish from the background, compared to broad or dispersed marks.
  
 * **Specific parameters**:
    - **upstream_limit**: indicate the window length you desire to define as a promoter region or TSS (transcription start site) upstream the peak identified, in the analysis of the possible target genes. For more information, visit[`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html).
    - **downstream_limit**: indicate the window length you desire to define as a promoter region or TSS (transcription start site) downstream the peak identified, in the analysis of the possible target genes. For more information, visit [`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html).
    - **motif_length**: indicate the length of motifs to be found. The length of time it takes to find motifs increases greatly with increasing size.  In general, it is recommended to start searching for enrichment with shorter lengths (i.e. less than 15) before trying longer lengths. For more information, visit [`Homer`](http://homer.ucsd.edu/homer/download.html)
    - **motif_size**: indicate the size of the region analyze when searching for motifs. For example, if it you are aimed at register ChIP-Seq peaks from a transcription factor, it is recommended a value between 50-200, whereas if you are interested in histone marked regions, 500-1000 would be suitable. For more information, visit [`Homer`](http://homer.ucsd.edu/homer/download.html)
    - **chromosome**: indicate the chromosome you are referring to

### Output and how to interpret the results

This section describes the output produced by the pipeline

#### Pipeline overview

The pipeline is built using [`GNU nano`](https://www.nano-editor.org), a command line text editor for Unix and Linux operating systems. See Pipeline summary for a condensed overview of the steps in the pipeline, and the bioinformatics tools used at each step.
For more information regarding the ChIP-seq protocol, and for an extensive list of publications, visit [`abcam website`](https://www.abcam.com/protocols/cross-linking-chromatin-immunoprecipitation-x-chip-protocol)

The directories listed below will be created in the output directory after the pipeline has finished. All paths are relative to the top-level results directory.

#### QC Reports

[`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) analyses the quality metrics of your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%A/C/G/T). You get information about adapter contamination and other overrepresented sequences.

<details markdown="1">
    <summary>Output files</summary>

```bash 
<working_directory>/<analysis_name>/samples/sample_<#>
```
  
`*_fastqc.html`: FastQC report containing quality metrics for read 1 (and read2 if paired-end).

`*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

#### Genome index

[`Bowtie2-build`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) builds a Bowtie index from a set of DNA sequences. The output files, listed below, together constitute the index: they are all that is needed to align reads to that reference. The original sequence FASTA files are no longer used by Bowtie 2 once the index is built.

<details markdown="1">
    <summary>Output files</summary>

```bash 
<working_directory>/<analysis_name>/genome
```

`index.1.bt2`,`index.2.bt2`,`index.3.bt2`,`index.4.bt2`,`index.rev.1.bt2`

</details>

#### Sample processing

[`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is used for mapping the samples to the genome, which has been previously indexed.

[`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) outputs alignments in [`SAM`](http://samtools.github.io/hts-specs/SAMv1.pdf) (Sequence Alignment/Map) format, enabling interoperation with a large number of other tools, including ([`samtools`](http://www.htslib.org). 

   > Note that *.sam files will be deleted once they are processed, so they will not appear in the final results.

([`Samtools`](http://www.htslib.org) is a set of utilities that manipulate alignments in the [`BAM`](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm) (Binary Alignment/ Map) format. It imports from and exports to the [`SAM`](http://samtools.github.io/hts-specs/SAMv1.pdf) (Sequence Alignment/Map) format, does sorting, merging and indexing, and allows to retrieve reads in any regions swiftly.

The pipeline uses in the first place [`samtools sort`](http://www.htslib.org/doc/samtools-sort.html), sorting and converting the *.sam file generated previously into a *.bam file. A BAM file is the compressed binary version of a SAM file. As the *.bam file is still too large to work with, it is later index by using [`samtools index`](http://www.htslib.org/doc/samtools-index.html).

<details markdown="1">
    <summary>Output files</summary>
    
```bash
<working_directory>/<analysis_name>/samples/sample_<#>
```
    
`*.bam`: file generated by samtools sort, sorting and converting *.sam file generated by bowtie2 into a *.bam file
`*.bam.bai`: file generated by samtools index where *.bam file is indexed.

</details>

#### Call peaks

[`MACS`](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html) (Model-based Analysis of ChIP-Seq) is one of the most common algorithm for identifying transcript factor binding sites, so called peak-calling. [`MACS`](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html) evaluate the significance of enriched ChIP regions and find out the binding sites through combining the information of both sequencing tag position and orientation.

As indicated previously in the description of Input parameters, calling can be done in a broad manner (True/1) or as in narrow manner (False/0). In this sense, if broad-peaks are called, the algorithm will put nearby highly enriched regions into a broad region with loose cutoff

<details markdown="1">
    <summary>Output files</summary>
    
```bash
<working_directory>/<analysis_name>/results
```
    
`*peaks.xls`: is a tabular file which contains information about called peaks. You can open it in excel and sort/filter using excel functions. It does includes:

   * chromosome name
   * start position of peak
   * end position of peak
   * length of peak region
   * absolute peak summit position
   * pileup height at peak summit
   * -log10(pvalue) for the peak summit (e.g. pvalue =1e-10, then this value should be 10)
   * fold enrichment for this peak summit against random Poisson distribution with local lambda
   * -log10(qvalue) at peak summit
   
`*.peaks.narrowPeak`: contains the peak locations together with peak summit, p-value, and q-value. You can load it to the UCSC genome browser. Definition of some specific columns are:

   * 5th: integer score for display. It's calculated as int(-10*log10pvalue) or int(-10*log10qvalue) depending on whether -p (pvalue) or -q (qvalue) is used as score cutoff. Please note that currently this value might be out of the [0-1000] range defined in UCSC ENCODE narrowPeak format. You can let the value saturated at 1000 (i.e. p/q-value = 10^-100) by using the following 1-liner awk: awk -v OFS="\t" '{$5=$5>1000?1000:$5} {print}' NAME_peaks.narrowPeak
   * 7th: fold-change at peak summit
   * 8th: -log10pvalue at peak summit
   * 9th: -log10qvalue at peak summit
   * 10th: relative summit position to peak start

If the pipeline is run in a broad peak mode, this file will be called *.peaks.broadPeak, offering same piece of information, except for 10th column, as this mode will not register peak summits.

`*peaks.summits.bed`: file in BED format containing the peak summits locations for every peak. The file can be loaded directly to the UCSC genome browser. This file is useful if you are aimed to find the motifs at the binding sites.

`*model.r`: R script which you can use to produce a PDF image of the model based on your data. Once the script is run, the PDF will automatically appear in the current directory.

</details>

#### Motif finding

[`HOMER`](http://homer.ucsd.edu/homer/index.html) (Hypergeometric Optimization of Motif EnRichment) is a suite of tools for Motif Discovery and next-generation sequencing analysis. It contains a novel motif discovery algorithm that was designed for regulatory element analysis in genomics applications. It is a differential motif discovery algorithm, which means that it takes two sets of sequences and tries to identify the regulatory elements that are specifically enriched in on set relative to the other.

There are several workflows for running motif analysis with [`HOMER`](http://homer.ucsd.edu/homer/index.html). The one used in this pipeline is [`findMotifsGenome.pl`](http://homer.ucsd.edu/homer/ngs/peakMotifs.html), which manage all the steps for discovering motifs in genomic regions. By default, this will perform de novo motif discovery as well as check the enrichment of known motifs.

<details markdown="1">
    <summary>Output files</summary>
    
```bash
<working_directory>/<analysis_name>/results/motifs_sample<#>
```
     
`homerMotifs.motifs<#>`: output files from the de novo motif finding, separated by motif length, and represent separate runs of the algorithm.

`homerMotifs.all.motifs`: file composed of all the homerMotifs.motifs<#> files.

`knownResults.html`: formatted output of known motif finding.

`knownResults/ directory`: contains files for the knownResults.html webpage, including known<#>.motif files for use in finding specific instance of each motif.

`knownResults.txt`: text file containing statistics about known motif enrichment. By default is opened in Text.editor. For optimized view, open manually in Excel.

`homerResults.html`: formatted output of de novo motif finding.

`homerResults/ directory`: contains files for the homerResults.html webpage, including motif<#>.motif files for use in finding specific instance of each motif.

`seq.autonorm.tsv`: autonormalization statistics for lower-order oligo normalization.

`motifFindingParameters.txt`: command used to execute findMotifsGenome.pl

</details>

## Contribution and Support

If you enjoy what we have done and want to support us, please, give us a 10 in the subject. Thank you so much for using ChipSeqPipeline!

