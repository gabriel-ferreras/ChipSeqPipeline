# ChipSeqPipeline
Bash pipeline for *Arabidopsis thaliana* Chip-Seq data processing and analysis.

Authors: 
  * Antonio Álvarez Gómez
  * Gabriel Ferreras Garrucho
  * Helena Victoria Cotán

## Introduction

The pipeline hereby presented was designed in order to perform a complete processing and analysis of ChIP-seq (Chromatin ImmunopreciPitation sequencing) data from the plant model *Arabidopsis thaliana*. This material is part of the course "Omics technologies and Bioinformatics" given by Francisco J. Romero Campero and Ignacio Pérez Hurtado de Mendoza in the Biochemistry degree at the University of Seville.

## Pipeline summary

1. **ChipSeqPipeline.sh**
   * Read parameters
   * Prepare workspace
   * Copy the data
   * Create a genome index ([`bowtie2-build`](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml))
   * Submit chip_sample_processing.sh and control_sample_processing for each sample
 2. **chip_sample_processing.sh** and **control_sample_processing.sh** 
    * Quality control ([`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
    * Map to reference genome ([`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml))
    * Generate sorted bam file ([`samtools`](http://www.htslib.org))
    * Submit peak_determination.sh for each sample
 3. **peak_determination.sh**
    * Peak determination ([`masc2 callpeak`](https://github.com/macs3-project/MACS))
    * Peak annotation by submitting *target_genes.R* for each sample
    * Homer-motifs-finding ([`findMotifsGenome.pl`](http://homer.ucsd.edu/homer/ngs/peakMotifs.html))
    * Overlapping target genes for replicates, GO and KEGG enrichment by submitting *exp_analysis.R* for each experiment
 4. **target_genes.R**
    * Install packages if neccesary ([`BiocManager`](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html))
    * Read arguments 
    * Read peak file ([`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html))
    * Definition of promoter region ([`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html))
    * Peak annotation ([`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html),[`DO.db`](http://bioconductor.org/packages/release/data/annotation/html/DO.db.html))
    * Extract target genes from annotation
 5. **exp_analysis.R**
    * Install packages if neccesary ([`BiocManager`](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html))
    * Read arguments
    * Read target gene files for each replicate
    * Extraction of overlapping genes in all replicates ([`VennDiagram`](https://cran.r-project.org/web/packages/VennDiagram/VennDiagram.pdf))
    * Gene ontology enrichment ([`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html),[`TxDb.Athaliana.BioMart.plantsmart28`](https://bioconductor.org/packages/release/data/annotation/html/TxDb.Athaliana.BioMart.plantsmart28.html),[`pathview`](https://bioconductor.org/packages/release/bioc/html/pathview.html),[`org.At.tair.db`](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html))
    * KEGG pathway enrichment ([`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html),[`TxDb.Athaliana.BioMart.plantsmart28`](https://bioconductor.org/packages/release/data/annotation/html/TxDb.Athaliana.BioMart.plantsmart28.html),[`pathview`](https://bioconductor.org/packages/release/bioc/html/pathview.html),[`org.At.tair.db`](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html))

## What does it do?

  * Quality control reads using fastQC 
  * Map reads to reference genome using Bowtie2
  * Call peaks for each sample (chip-control pair) using MACS2
  * Motifs-finding using findMotifsGenome.pl for each sample (chip-control pair)
  * Asociation of the peaks to target genes using ChIPseeker
  * Overlap of target genes from samples using VennDiagram
  * Analysis of enrichment in terms of gene ontology (GO) and Kyoto Encycopedia of Genes and Genomes (KEGG) using clusterProfiler

## What does it output?

  * QC reports
  * Genome index
  * Bam and bam.bai files
  * Peak annotation files
  * Motifs-finding analysis reports
  * Analysis of enrichment in terms of gene ontology (GO) and Kyoto Encycopedia of Genes and Genomes (KEGG) analysis plots

## Dependencies

For running this pipeline, you need to assure that the following softwares are installed in your computer:

  * [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  * [`SAMTOOLS`](https://sourceforge.net/projects/samtools/files/samtools/)
  * [`MACS2`](https://github.com/macs3-project/MACS)
  * [`R`](https://www.r-project.org), with the following packages:
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

This pipeline does not require any installation apart from those mentionated above in Dependencies section. Usage just does require download the repository from GitHub. If you do not have Git installed, consult https://github.com/git-guides/install-git.

The steps you need to follow are these described:

  1. Open a terminal
  2. Go to the working directory you desire to save the repository
  3. Introduce the following command

      ```bash
      git clone https://github.com/gabriel-ferreras/ChipSeqPipeline.git
      ```

## Running pipeline
  
To run the pipeline on your experimental data, first enter the necessary parameters in the spreadsheet file (see Input parameters section), and then from the terminal type:

```sh
bash ChipSeqPipeline test_params.txt
```

### Input parameters

The pipeline requires a file as an input to specify the samples and the design of the analysis. Within the repository you can find an example of parameters file (test_params.txt). The parameters you need to include are:

 * **General parameters**:
    - **analysis_name**: is the name you desire to give to the analysis. The folder created during the process will be named after it.
    - **number_samples**: the number of samples that are going to be processed, that is, the total number of Chip/Control pairs. Does not matter whether the samples are paired or not.
    - **number_experiments**: the number of experiments (ChipSeq analyses for different transcription factors or histone modifications) that are going to be processed.
    - **experimental_design**: specify the samples that correspond to each experiment, that is, the replicates for a same ChipSeq analysis. For example, if there are two experiments and each experiment has two samples (replicates), indicate it as "(1,1,2,2)".
 
 * **Directories**:
    - **working_directory**: the directory where all the data and results are going to be placed during the analysis.
    - **installation_directory**: the directory where the cloned github repository is placed.
    
 * **Paths**:
    - **path_annotation**: the path to where the annotation (in gtf format) is placed, including the file.
    - **path_genome**: the path to where the reference genome (in fasta format) is placed, including the file.
    - **path_chip_1**: the path to where the chip samples (in fastq.gz format) are placed, including the file. If there are more than one samples, you must also indicate the path to those by adding a parameter called "path_chip_i", where "i" stands for the number of the sample.
    - **path_control_1**: the path to where the control samples (in fastq.gz format) are placed, including the file. If there are more than one sampls, you must also indicate the path to those by adding a parameter called "path_control_i", where "i" stands for the number of the sample.
    
 * **Eligible parameters**:
    - **paired**: indicate whether the samples are from a paired-end sequencing (True/1) or a single-end sequencing (False/0).
       > [`Paired-end sequencing`](https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html) allows users to sequence both ends of a fragment and generate high-quality, alignable sequence data. It produces twice the bumber of reads in the same time and detect insertion-deletion (indel) variants, which is not possible with single-read data.
       > [`Single-end sequencing`](https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html) involves sequencing DNA from only one end, delivering large volumes of high-quality data, rapidly and economically.
    - **broad**: indicate whether, in the peak analysis, you prefer to identify broad peaks (True/1) or narrow peaks (False/0). 
       > [`MACS2`](https://github.com/macs3-project/MACS) has two different specific methods in searching for enrichment regions. [`Broad peaks`](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html ) or broad domains (i.e. histone modifications that cover entire gene bodies) or [`narrow peaks`](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html ) (i.e. a transcription factor binding). Narrow peaks are easier to detect as we are looking for regions that have higher amplitude and are easier to distinguish from the background, compared to broad or dispersed marks.
  
 * **Specific parameters**:
    - **upstream_limit**: indicate the window length you desire to define as a promoter region upstream of the TSS (transcription start site) in the analysis of the possible target genes. For more information, visit [`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html).
    - **downstream_limit**: indicate the window length you desire to define as a promoter region downstream of the TSS (transcription start site) in the analysis of the possible target genes. For more information, visit [`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html).
    - **motif_length**: indicate the length of motifs to be found. The length of time it takes to find motifs increases greatly with increasing size.  In general, it is recommended to start searching for enrichment with shorter lengths (i.e. less than 15) before trying longer lengths. For more information, visit [`Homer`](http://homer.ucsd.edu/homer/download.html).
    - **motif_size**: indicate the size of the region analyze when searching for motifs. For example, if it you are aimed at register ChIP-Seq peaks from a transcription factor, it is recommended a value between 50-200, whereas if you are interested in histone marked regions, 500-1000 would be suitable. For more information, visit [`Homer`](http://homer.ucsd.edu/homer/download.html).
    - **chromosome**: indicate the universe you are refering to when performing [`Gene Ontology (GO)`](http://geneontology.org/docs/go-enrichment-analysis/) enrichment. This parameter can take several values: `ALL`, if you are referring to the whole genome; `C`, if you are referring to chloroplast DNA; `M`, if you are referring to mitochondrial DNA; or an `integer` if you are referring to a certain chromosome.

### Output and how to interpret the results

This section describes the output produced by the pipeline.

#### Pipeline overview

The pipeline is built using [`GNU nano`](https://www.nano-editor.org), a command line text editor for Unix and Linux operating systems. See the Pipeline summary section for a condensed overview of the steps in the pipeline, and the bioinformatics tools used at each step.

For more information regarding the ChIP-seq protocol, and for an extensive list of publications, visit [`abcam website`](https://www.abcam.com/protocols/cross-linking-chromatin-immunoprecipitation-x-chip-protocol).

The directories listed below will be created in the output directory after the pipeline has finished. All paths are relative to the top-level results directory.

#### QC Reports

[`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) analyses the quality metrics of your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%A/C/G/T), adapter contamination and other overrepresented sequences, among other issues.

<details markdown="1">
    <summary>Output files</summary>

```bash 
<working_directory>/<analysis_name>/samples/sample_<i>
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

   > Note that *.sam* files, as well as the initial *.fastq.gz* files, will be deleted once they are processed, so they will not appear in the final results.

([`Samtools`](http://www.htslib.org) is a set of utilities that manipulate alignments in the [`BAM`](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm) (Binary Alignment/ Map) format. It imports from and exports to the [`SAM`](http://samtools.github.io/hts-specs/SAMv1.pdf) (Sequence Alignment/Map) format, does sorting, merging and indexing, and allows to retrieve reads in any regions swiftly.

The pipeline uses in the first place [`samtools sort`](http://www.htslib.org/doc/samtools-sort.html), sorting and converting the *.sam* file generated previously into a *.bam file*. A BAM file is the compressed binary version of a SAM file. As the *.bam file* is still too large to work with, it is later indexed by using [`samtools index`](http://www.htslib.org/doc/samtools-index.html), generating a *bam.bai* file.

<details markdown="1">
    <summary>Output files</summary>
    
```bash
<working_directory>/<analysis_name>/samples/sample_<i>
```
    
`*.bam`: file generated by samtools sort, sorting and converting *.sam* file generated by bowtie2 into a *.bam* file.

`*.bam.bai`: file generated by samtools index by indexing *.bam* file.

</details>

#### Call peaks

[`MACS`](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html) (Model-based Analysis of ChIP-Seq) is one of the most common algorithm for identifying transcript factor binding sites, so called peak-calling, evaluate the significance of enriched ChIP regions and find out the binding sites through combining the information of both sequencing tag position and orientation. Peaks are therefore define as regiones in the genome where there is an accumulation of lectures, indicating the binding of a transcript factor.

As indicated previously in Input parameters section, calling can be done in a broad manner (True/1) or as in narrow manner (False/0). In this sense, if broad-peaks are called, the algorithm will put nearby highly enriched regions into a broad region with loose cutoff.

<details markdown="1">
    <summary>Output files</summary>
    
```bash
<working_directory>/<analysis_name>/results/sample_<i>_result
```
    
`<analysis_name>_sample_<i>_peaks.xls`: is a tabular file which contains information about called peaks. You can open it in excel and sort/filter using excel functions. It does includes:

   * chromosome name
   * start position of peak
   * end position of peak
   * length of peak region
   * absolute peak summit position
   * pileup height at peak summit
   * -log10(pvalue) for the peak summit (e.g. pvalue =1e-10, then this value should be 10)
   * fold enrichment for this peak summit against random Poisson distribution with local lambda
   * -log10(qvalue) at peak summit
   
`<analysis_name>_sample_<i>_peaks.narrowPeak`: contains the peak locations together with peak summit, p-value, and q-value. You can load it to the UCSC genome browser. Definition of some specific columns are:

   * 5th: integer score for display. It's calculated as int(-10*log10pvalue) or int(-10*log10qvalue) depending on whether -p (pvalue) or -q (qvalue) is used as score cutoff. Please note that currently this value might be out of the [0-1000] range defined in UCSC ENCODE narrowPeak format.
   * 7th: fold-change at peak summit
   * 8th: -log10pvalue at peak summit
   * 9th: -log10qvalue at peak summit
   * 10th: relative summit position to peak start

If the pipeline is run in a broad peak mode, this file will be called *<analysis_name>_sample_<i>_peaks.broadPeak*, offering same piece of information, except for 10th column, as this mode will not register peak summits.

`<analysis_name>_sample_<i>_summits.bed`: file in BED format containing the peak summits locations for every peak. The file can be loaded directly to the UCSC genome browser. This file is useful if you are aimed to find the motifs at the binding sites.

`<analysis_name>_sample_<i>_model.r`: R script which you can use to produce a PDF image of the model based on your data. Once the script is run, the PDF will automatically appear in the current directory.

</details>

#### Motifs finding

[`HOMER`](http://homer.ucsd.edu/homer/index.html) (Hypergeometric Optimization of Motif EnRichment) is a suite of tools for Motif Discovery and next-generation sequencing analysis. It contains a novel motif discovery algorithm that was designed for regulatory element analysis in genomics applications. It is a differential motif discovery algorithm, which means that it takes two sets of sequences and tries to identify the regulatory elements that are specifically enriched in on set relative to the other.

There are several workflows for running motif analysis with [`HOMER`](http://homer.ucsd.edu/homer/index.html). The one used in this pipeline is [`findMotifsGenome.pl`](http://homer.ucsd.edu/homer/ngs/peakMotifs.html), which manage all the steps for discovering motifs in genomic regions. By default, this will perform de novo motif discovery as well as check the enrichment of known motifs.

<details markdown="1">
    <summary>Output files</summary>
    
```bash
<working_directory>/<analysis_name>/results/sample_<i>_result/motifs_sample<i>
```
     
`homerMotifs.motifs<motif_length>`: output files from the de novo motif finding, separated by motif length, and represent separate runs of the algorithm.

`homerMotifs.all.motifs`: file composed of all the homerMotifs.motifs<motif_length> files.

`knownResults.html`: formatted output of known motif finding.

`knownResults/ directory`: contains files for the knownResults.html webpage, including known<#>.motif files for use in finding specific instance of each motif.

`knownResults.txt`: text file containing statistics about known motif enrichment. By default is opened in Text.editor. For optimized view, open manually in Excel.

`homerResults.html`: formatted output of de novo motif finding.

`homerResults/ directory`: contains files for the homerResults.html webpage, including motif<#>.motif files for use in finding specific instance of each motif.

`seq.autonorm.tsv`: autonormalization statistics for lower-order oligo normalization.

`motifFindingParameters.txt`: command used to execute findMotifsGenome.pl.

</details>

#### Association of the peaks to target genes

[`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html) is a package that implements the functions necessary for retrieving the nearest genes around the peak, annotate genomic region of the peak and estimate the significance of overlap among ChIP peak data sets, among other functions. 

The criteria selected for the determination of target genes (that is, the potential regulome of the transcription factor) is Nearest Downstream Gene (NDG), associating to each peak its nearest downstream gene. This is based on the assumption that transcription factors usually bind to promoters, that is, 5' of the regulated gene. Among the two genes on each strand, the closest one is chosen as the target gene. 

This peak annotation is applied both to the peaks (either *.narrowPeak* or *.broadPeak* files) and the summits (*.bed* files). Apart from the determination of target genes, several visualization functions are implemented to summarize the finding of the analysis.

<details markdown="1">
    <summary>Output files</summary>
    
```bash
<working_directory>/<analysis_name>/results/sample_<i>_result
```

`<analysis_name>_sample_<i>_peak_targetgenes.txt`: text file including the genes of *Arabidopsis thaliana* (TAIR ID) identified as target genes from peaks (*.narrowPeak* or *.broadPeak* files).

`<analysis_name>_sample_<i>_summit_targetgenes.txt`: text file including the genes of *Arabidopsis thaliana* (TAIR ID) identified as target genes from summits (*.bed* files).

`R.plots`: a pdf file including all the plots generated during the sample processing. These plots are, always in this order:

  1. Pie plot of the distribution of peaks in different genetic elements, generated by *plotAnnoPie()*.
  2. Bar plot of the distribution of peaks in different genetic elements, generated by *plotAnnoBar()*.
  3. Distance plot of the distribution of peaks regarding their distance to the TSS, generated by *plotDistToTSS()*.
  4. Pie plot of the distribution of summits in different genetic elements, generated by *plotAnnoPie()*.
  5. Bar plot of the distribution of summits in different genetic elements, generated by *plotAnnoBar()*.
  6. Distance plot of the distribution of summits regarding their distance to the TSS, generated by *plotDistToTSS()*.

</details>

#### Determination of the potentital regulome of the transcription factor

[`VennDiagram`](https://www.rdocumentation.org/packages/VennDiagram/versions/1.6.20) is a package that implement a series of functions in order to generate high-resolution Venn and Euler plots. It will be used to overlap the target genes (from peaks, not summits) identified in all samples from each experiment, that is, the different replicates of a same ChipSeq analysis, therefore obtaining the a more reliable potential regulome of the transcription factor studied than the ones obtained from each individual sample. 

The regulome is defined as the global set of genes regulated by a transcription factor under certain circumstances, and it can be inferred from the cistrome, the set of genetic positions where the transcription factor binds under those certain conditions, by the association of ChipSeq peaks to target genes via the NDG criteria. However, this method of determination of the regulome is highly is problematic as it results in a high percentage of false positives (genes to whose promoters the transcription factor is binding, but is not regulating). In this pipeline this lack of reliability of the method is diminished by the overlap of the individual regulomes determined for each replicate of the same ChipSeq experiment. It should be mentioned that an even more efective and reliable methodology would be the overlap of these ChipSeq-determined potential regulome with differentially-expressed genes (DEGs) obtained via transcriptomic analysis such as RNA-Seq. 

<details markdown="1">
    <summary>Output files</summary>
    
```bash
<working_directory>/<analysis_name>/results/exp_<j>_result
```

`<analysis_name>_experiment_<j>_overlappinggenes.txt`: text file containing the overlapping genes that constitute the determined potential regulome of the transcription factor.

</details>

#### Gene ontology enrichment (GO) and Kyoto Encycopedia of Genes and Genomes (KEGG)

[`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) is a package that implements methods to analyze and visualize functional profiles of gene and gene clusters, these are, gene ontology enrichment (GO) and Kyoto Encycopedia of Genes and Genomes (KEGG). 

The [`Gene Ontology (GO)`](http://geneontology.org/docs/go-enrichment-analysis/) refers to a structured and specific terminology to describe biological processes, cellular components or molecular functions in a a systematic (and therefore faster) and unequivocal way. As the intereset is to find those items that, as a consequence of the treatment or mutation, an enrichment is performed. An enrichment is a mathematical-computational method that determines if there are GO terms in a set of genes that appear with statistical significance with respect to the genome of the organism.

The [`Kyoto Encyclopedia of Genes and Genomes (KEGG)`](https://www.genome.jp/kegg/) is a database intended for the study of High-level biological functions and utilities from large-scale molecular information from genomic or transcriptomic analysis. Therefore, performing a KEGG analysis is very similar to that perform in GO terms, but at a higher functional level, of complete metabolic or regulatory pathways, rather than at the level of specific biological functions or processes.

These enrichment analyses, both of GO and KEGG, are perfomed using *clusterProfiler* for the lists of overlapping target genes from all the replicates of each experiment, and will enable the association of possible functions to the transcription factor studied.

<details markdown="1">
    <summary>Output files</summary>
    
```bash
<working_directory>/<analysis_name>/results/exp_<j>_result
```
`<analysis_name>_experiment_<j>_BP_GOs.txt`: text file containing the enriched biological process GO terms found for the overlapping target gene list. Could be empty if no enriched terms are found.
`<analysis_name>_experiment_<j>_MF_GOs.txt`: text file containing the enriched molecular function GO terms found for the overlapping target gene list. Could be empty if no enriched terms are found.
`<analysis_name>_experiment_<j>_CC_GOs.txt`: text file containing the enriched celullar component GO terms found for the overlapping target gene list. Could be empty if no enriched terms are found.
`<analysis_name>_experiment_<j>_KEGGs.txt`: text file containing the enriched KEGG terms found for the overlapping target gene list. Could be empty if no enriched terms are found.
`R.plots`: a pdf file including all the plots generated during the enrichment analysis. These plots are, always in this order:
  1. Bar plot of enriched biological process GO terms, generated by *batplot()*, showing 20 categories.
  2. Dot plot of enriched biological process GO terms, generated by *dotplot()*, showing 20 categories.
  3. Cnet plot of enriched biological process GO terms, generated by *cnetplot()*.
  4. Ema plot of enriched biological process GO terms, generated by *emapplot()*.
  5. Go plot of enriched biological process GO terms, generated by *goplot()*.
  6. Bar plot of enriched molecular function GO terms, generated by *batplot()*, showing 20 categories.
  7. Dot plot of enriched molecular function GO terms, generated by *dotplot()*, showing 20 categories.
  8. Cnet plot of enriched molecular function GO terms, generated by *cnetplot()*.
  9. Ema plot of enriched molecular function GO terms, generated by *emapplot()*.
  10. Go plot of enriched molecular function GO terms, generated by *goplot()*.
  11. Bar plot of enriched celullar component GO terms, generated by *batplot()*, showing 20 categories.
  12. Dot plot of enriched celullar component GO terms, generated by *dotplot()*, showing 20 categories.
  13. Cnet plot of enriched celullar component GO terms, generated by *cnetplot()*.
  14. Ema plot of enriched celullar component GO terms, generated by *emapplot()*.
  15. Go plot of enriched celullar component GO terms, generated by *goplot()*.
  11. Bar plot of enriched KEGG terms, generated by *batplot()*, showing 10 categories.
  12. Dot plot of enriched KEGG GO terms, generated by *dotplot()*, showing 10 categories.
 
</details>

## Contribution and Support

If you enjoy what we have done and want to support us, please, give us a 10 in the subject. Thank you so much for using ChipSeqPipeline!

