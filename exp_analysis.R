############################################################################
#                         CHIP-SEQ EXPERIMENT                              # 
#                     TARGET GENES GLOBAL ANALYSIS                         #
############################################################################

## Authors:

# Antonio Álvarez Gómez (alvarezgomezantonio@gmail.com)
# Gabriel Ferreras Garrucho (gabrifg10@gmail.com)
# Helena Victoria Cotán (hevico99@gmail.com)


## Install the following packages if not already installed:

 # install.packages("VennDiagram")
 # if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
 # BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")

 library(VennDiagram)
 library(TxDb.Athaliana.BioMart.plantsmart28)
 txdb <- TxDb.Athaliana.BioMart.plantsmart28


## Reading arguments:

 args <- commandArgs(trailingOnly = TRUE)

 experiment <- args[1]
 exp_design <- args[2]
 num_samples <- args[3]
 analysis_name <- args[4]
 
 exp_design_vector <- as.numeric(unlist(strsplit(substr(exp_design, num_samples, nchar(hola)-1), ",")))

 
## Grouping the replicates of the experiment:
 
 samples_in_exp <- which(exp_design_vector == experiment)
 exp_all_peaks <- list()
 for (i in samples_in_exp) 
 {
  exp_all_peaks[i] <- as.vector(read.table(file = paste(analysis_name, "sample", i, "peaks_targetgenes.txt", sep = "_"), header = FALSE))
 }

 
## Extracting the overlapping peaks, that is, the regulome of the transcription factor:
 
 experiment_overlap <- calculate.overlap(exp_all_peaks)
 if (length(samples_in_exp) == 1) 
 {
   regulome <- experiment_overlap[[1]]
 } else if (length(samples_in_exp) == 2) 
 {
   regulome <- experiment_overlap$a3
 } else if (length(samples_in_exp) == 3)
 {
   regulome <- experiment_overlap$a5
 } else if (length(samples_in_exp) == 3)
 {
   regulome <- experiment_overlap$a6
 } else if (length(samples_in_exp) == 5)
 {
   regulome <- experiment_overlap$a31
 }

 
## Gene ontology enrichment:
 
 genes_atha <- as.data.frame(genes(txdb))

 my_universe <- genes_atha$gene_id

 BP <- enrichGO(gene = regulome,
               OrgDb = "org.At.tair.db",
               keyType = "TAIR", 
               ont = "BP", 
               universe = my_universe)

 barplot(BP, showCategory = 20)
 dotplot(BP, showCategory = 20)
 cnetplot(BP)
 emapplot(BP)

 MF <- enrichGO(gene = regulome,
                OrgDb = "org.At.tair.db",
                keyType = "TAIR", 
                ont = "MF", 
                universe = my_universe)

 barplot(MF, showCategory = 20)
 dotplot(MF, showCategory = 20)
 cnetplot(MF)
 emapplot(MF)

 CC <- enrichGO(gene = regulome,
                OrgDb = "org.At.tair.db",
                keyType = "TAIR", 
                ont = "CC", 
                universe = my_universe)

 barplot(CC, showCategory = 20)
 dotplot(CC, showCategory = 20)
 cnetplot(CC)
 emapplot(CC)


## KEGG pathway enrichment:
 
 KEGG <- enrichKEGG(gene = regulome,
                    organism = "ath",
                    pvalueCutoff = 0.05)

 browseKEGG(KEGG, "ath00010")

 pathview(gene.data  = regulome,
          pathway.id = "ath00010",
          species    = "ath")
