############################################################################
#                         CHIP-SEQ EXPERIMENT                              # 
#                     TARGET GENES GLOBAL ANALYSIS                         #
############################################################################

## Authors:

# Antonio Álvarez Gómez (alvarezgomezantonio@gmail.com)
# Gabriel Ferreras Garrucho (gabrifg10@gmail.com)
# Helena Victoria Cotán (hevico99@gmail.com)


## Installing the following packages if not already installed:
install.packages("VennDiagram")

library(VennDiagram)


## Reading arguments:

args <- commandArgs(trailingOnly = TRUE)

experiment <- args[1]
exp_design <- args[2]
num_samples <- args[3]
analysis_name <- args[4]


exp_design_vector <- as.numeric(unlist(strsplit(substr(exp_design, num_samples, nchar(hola)-1), ",")))

samples_in_exp <- which(exp_design_vector == experiment)

exp_all_peaks <- list()

for (i in samples_in_exp) 
{
  exp_all_peaks[i] <- as.vector(read.table(file = paste(analysis_name, "sample", i, "peaks_targetgenes.txt", sep = "_"), header = FALSE))
}

experiment_overlap <- calculate.overlap(exp_all_peaks)

## Al parecer no podemos extraerlo así:
enriched_target_genes <- experiment_overlap[[1]]

