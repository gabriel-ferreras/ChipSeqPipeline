############################################################################
#             PEAK AND SUMMIT ASSOCIATION TO TARGET GENES                  # 
#                          CHIP-SEQ ANALYSIS                               #
############################################################################
## Authors:
    # Antonio Álvarez Gómez (alvarezgomezantonio@gmail.com)
    # Gabriel Ferreras Garrucho (gabrifg10@gmail.com)
    # Helena Victoria Cotán (hevico99@gmail.com)

## Install the following packages if not already installed:

  #if (!requireNamespace("BiocManager", quietly = TRUE))
  #  install.packages("BiocManager")
  #BiocManager::install()
  #BiocManager::install("ChIPseeker")
  #BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")
  #BiocManager::install("DO.db")
  #BiocManager::install("org.At.tair.db")

  library(DO.db)
  library(ChIPseeker)
  library(TxDb.Athaliana.BioMart.plantsmart28)
  txdb <- TxDb.Athaliana.BioMart.plantsmart28
  library(org.At.tair.db)
  
## Reading arguments:
  args = commandArgs(trailingOnly=TRUE)
  
  peak_file <- args[1]
  summit_file<-args[2]
  up_limit <- as.numeric(args[3])
  down_limit <- as.numeric(args[4])
  peak_output <- args[5]
  summit_output <- args[6]

## Reading peak file:
  
  peaks <- readPeakFile(peakfile = peak_file,header=FALSE)
  summits <- readPeakFile(peakfile = summit_file,header=FALSE)
  
## Definition of promoter region:
  
  promoter <- getPromoters(TxDb=txdb, upstream=up_limit, downstream=down_limit)
  tagMatrix <- getTagMatrix(peaks, windows=promoter)

## Peak annotation:

  peakAnno <- annotatePeak(peak = peaks, tssRegion=c(-up_limit, down_limit), 
                         TxDb=txdb, annoDb = "org.At.tair.db")
  summitAnno <- annotatePeak(peak = summits, tssRegion=c(-up_limit, down_limit), 
                           TxDb=txdb, annoDb = "org.At.tair.db")
  plotAnnoPie(peakAnno)
  plotAnnoBar(peakAnno)
  plotDistToTSS(peakAnno,
                title="Distribution of genomic loci relative to TSS",
                ylab = "Genomic Loci (%) (5' -> 3')")
  plotAnnoPie(summitAnno)
  plotAnnoBar(summitAnno)
  plotDistToTSS(summitAnno,
                title="Distribution of genomic loci relative to TSS",
                ylab = "Genomic Loci (%) (5' -> 3')")

## Converting annotation to data frame and exporting:
  peakannotation <- as.data.frame(peakAnno)
  peak.target.genes <- peakannotation$geneId[grepl("Promoter", peakannotation$annotation)]
  write(x = peak.target.genes, file = peak_output)

  summitannotation <- as.data.frame(summitAnno)
  summit.target.genes <- summitannotation$geneId[grepl("Promoter", summitannotation$annotation)]
  write(x = summit.target.genes, file = summit_output)

