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
 library(clusterProfiler)
 library(org.At.tair.db)
 library(pathview)
 
## Reading arguments:

 args <- commandArgs(trailingOnly = TRUE)

 experiment <- as.numeric(args[1])
 exp_design <- args[2]
 num_samples <- as.numeric(args[3])
 analysis_name <- args[4]
 chromosome <- args[5]
 
 exp_design_vector <- as.numeric(unlist(strsplit(substr(exp_design, 2, nchar(exp_design)-1), ",")))

 
## Grouping the replicates of the experiment:
 
 samples_in_exp <- which(exp_design_vector == experiment)
 exp_all_peaks <- list()
 for (i in samples_in_exp) 
 {
  exp_all_peaks[i] <- as.vector(read.table(file = paste(analysis_name, "sample", i, "peaks_targetgenes.txt", sep = "_"), header = FALSE))
 }

 
## Extracting the overlapping peaks, that is, the potential regulome of the transcription factor:
 
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
 } else if (length(samples_in_exp) == 4)
 {
   regulome <- experiment_overlap$a6
 } else if (length(samples_in_exp) == 5)
 {
   regulome <- experiment_overlap$a31
 }
 name<-paste(analysis_name, "experiment", experiment, "overlappinggenes.txt", sep = "_")
 write(x = regulome, file = name)
 
## Gene ontology enrichment:
 
 genes_atha <- as.data.frame(genes(txdb))
 if(chromosome == "ALL")
 {
   genes_selection<-genes_atha
 }
 
 if(chromosome != "ALL")
 {
   chromosome<-as.numeric(chromosome)
   genes_selection<-subset(genes_atha, seqnames == chromosome)
 }
 
 my_universe <- genes_selection$gene_id
 
 BP <- enrichGO(gene = regulome,
               OrgDb = "org.At.tair.db",
               keyType = "TAIR", 
               ont = "BP", 
               universe = my_universe)
 
 BP_name<-paste(analysis_name, "experiment", experiment, "BP_GOs.txt", sep = "_")
 write.table(BP, file = BP_name, sep = ",")
 
 if (nrow(as.data.frame(BP)) > 0)
 {
    print("")
    print("BIOLOGICAL PROCESS ENRICHMENT FOUND")
    print("")
    barplot(BP, showCategory = 20)
 }
 if (nrow(as.data.frame(BP)) > 0)
 {
    dotplot(BP, showCategory = 20)
 }
 if (nrow(as.data.frame(BP)) > 0)
 {
    cnetplot(BP)
 }
 if (nrow(as.data.frame(BP)) > 1)
 {
    emapplot(BP)
 }
 if (nrow(as.data.frame(BP)) > 1)
 {
    goplot(BP)
 }
 if (nrow(as.data.frame(BP)) == 0)
 {
    print("")
    print("NO BIOLOGICAL PROCESS ENRICHMENT FOUND")
    print("")
 }
 
 MF <- enrichGO(gene = regulome,
                OrgDb = "org.At.tair.db",
                keyType = "TAIR", 
                ont = "MF", 
                universe = my_universe)
 
 MF_name<-paste(analysis_name, "experiment", experiment, "MF_GOs.txt", sep = "_")
 write.table(MF, file = MF_name, sep = ",")
 
 if (nrow(as.data.frame(MF)) > 0)
 {
    print("")
    print("MOLECULAR FUNCTION ENRICHMENT FOUND")
    print("")
    barplot(MF, showCategory = 20)
 }
 if (nrow(as.data.frame(MF)) > 0)
 {
    dotplot(MF, showCategory = 20)
 }
 if (nrow(as.data.frame(MF)) > 0)
 {
    cnetplot(MF)
 }
 if (nrow(as.data.frame(MF)) > 1)
 {
    emapplot(MF)
 }
 if (nrow(as.data.frame(MF)) > 1)
 {
    goplot(MF)
 }
 if (nrow(as.data.frame(MF)) == 0)
 {
    print("")
    print("NO MOLECULAR FUNCTION ENRICHMENT FOUND")
    print("")
 }
 
 CC <- enrichGO(gene = regulome,
                OrgDb = "org.At.tair.db",
                keyType = "TAIR", 
                ont = "CC", 
                universe = my_universe)
 
 CC_name<-paste(analysis_name, "experiment", experiment, "CC_GOs.txt", sep = "_")
 write.table(CC, file = CC_name, sep = ",")
 
 if (nrow(as.data.frame(CC)) > 0)
 {
    print("")
    print("CELLULAR COMPONENT ENRICHMENT FOUND")
    print("")
    barplot(CC, showCategory = 20)
 }
 if (nrow(as.data.frame(CC)) > 0)
 {
    dotplot(CC, showCategory = 20)
 }
 if (nrow(as.data.frame(CC)) > 0)
 {
    cnetplot(CC)
 }
 if (nrow(as.data.frame(CC)) > 1)
 {
    emapplot(CC)
 }
 if (nrow(as.data.frame(CC)) > 1)
 {
    goplot(CC)
 }
 if (nrow(as.data.frame(CC)) == 0)
 {
    print("")
    print("NO CELLULAR COMPONENT ENRICHMENT FOUND")
    print("")
 }

## KEGG pathway enrichment:
 
 kk <- enrichKEGG(gene = regulome,
                  universe = my_universe, 
                  organism = "ath",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05)

 kk_name<-paste(analysis_name, "experiment", experiment, "KEGGs.txt", sep = "_")
 write.table(kk, file = kk_name, sep = ",")
 
 if (nrow(as.data.frame(kk)) > 0)
 {
    print("")
    print("KEGG ENRICHMENT FOUND")
    print("")
    for (i in 1:nrow(as.data.frame(kk)))
    {
       pathview(gene.data  = regulome, 
                gene.idtype = "KEGG", 
                pathway.id = kk[i]$ID,
                species    = "ath")
    }
 }
 if (nrow(as.data.frame(kk)) > 0)
 {
    barplot(kk, showCategory = 10)
 }
 if (nrow(as.data.frame(kk)) > 0)
 {
    dotplot(kk, showCategory = 10)
 }

 if (nrow(as.data.frame(kk)) == 0)
 {
    print("")
    print("NO KEGG ENRICHMENT FOUND")
    print("")
 }
 
