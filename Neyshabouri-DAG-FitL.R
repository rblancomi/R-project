#################################################################

library(OncoSimulR)

# Neyshabouri et al. (2020) propose a probabilistic model of mutually exclusive 
# linearly ordered driver pathways and analyze two large datasets of colorectal
# adenocarcinoma and GBM from IntOGen-mutations database

### Colorectal adenocarcinoma (COADREAD) dataset

## DAG construction under fitness specifications (satisfied and no satisfied) 
## based on the MCMC method inferred model from Neyshabouri et al.

# The MCMC method generates a linear model, thus generating an oncogenic tree

COADREAD_fitness_df <- data.frame(parent = c("Root", "A", "B", "C", "D",
                                             "E", "F"),
                                  child = c("A", "B", "C", "D", "E", "F",
                                            "G"),
                                  s = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
                                  sh = -1,
                                  typeDep = "MN")

COADREAD_fitness_df

COADREAD_DAG <- allFitnessEffects(COADREAD_fitness_df, 
                                  geneToModule = c( "Root" = "Root",
                                                    "A" = "APC",
                                                    "B" = "TP53",
                                                    "C" = "KRAS",
                                                    "D" = "PIK3CA, NRAS, LRP1B",
                                                    "E" = "TCF7L2, FBXW7, FAT4, ARID1A",
                                                    "F" = "ATM, SMAD2, ERBB3, MTOR, CTNNB1",
                                                    "G" = "SOX9, SMAD4"))

# Posible bug? Si pones los genes de un mismo módulo en líneas distintas,
# salta un error
plot(COADREAD_DAG)
dev.off()

COADREAD_DAG_FL <- evalAllGenotypes(COADREAD_DAG, max = 131072)

# Due to memory problems, parallelize

# install.packages("parallel")
# library(parallel)

# mclapply(COADREAD_DAG_FL, plotFitnessLandscape(), mc.cores = detectCores())

plotFitnessLandscape(COADREAD_DAG_FL) # Not working. Error: vector memory exhausted (limit reached?)
dev.off()


### GMB dataset

## DAG construction under fitness specifications (satisfied and no satisfied) 
## based on the MCMC method inferred model from Neyshabouri et al.

GMB_fitness_df <- data.frame(parent = c("Root", "A", "B", "C", "D",
                                             "E", "F"),
                                  child = c("A", "B", "C", "D", "E", "F",
                                            "G"),
                                  s = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
                                  sh = -1,
                                  typeDep = "MN")

GMB_fitness_df

GMB_DAG <- allFitnessEffects(GMB_fitness_df, 
                                  geneToModule = c( "Root" = "Root",
                                                    "A" = "NF1, EGFR",
                                                    "B" = "PTEN",
                                                    "C" = "TP53",
                                                    "D" = "PIK3R1, PIK3CA",
                                                    "E" = "RB1, IDH1, STAG2",
                                                    "F" = "ATRX, LZTR1",
                                                    "G" = "BCOR, DCAF12L2"))

plot(GMB_DAG)
dev.off()

GMB_DAG_FL <- evalAllGenotypes(GMB_DAG, max = 8192)

plotFitnessLandscape(GMB_DAG_FL)
dev.off()

