#################################################################

library(OncoSimulR)

# Neyshabouri et al. (2020) propose a probabilistic model of mutually exclusive 
# linearly ordered driver pathways and analyze two large datasets of colorectal
# adenocarcinoma and GBM from IntOGen-mutations database. Their model assumes:
# 1. Driver genes are over-represented among those mutated across a large tumor 
# collection and, thus, they can be identified in terms of frequency
# 2. Driver genes participating of the same pathway are mutated in a mutually 
# exclusive manner because more than one mutation in a pathway does not give any 
# selective advantage to the clone

# The analysis and re-evaluation of Neyshabouri et al. model will be perform in
# the following steps:
# 1. Map their model to an evolutionary model (fitness landscape) using OncoSimulR accounting for
# their inferred order of restrictions, mutual exclusivity relationships and the
# rest of assumptions.
# 2. Map their model to an evolutionary model using OncoSimulR accounting only for
# their inferred order of restrictions and mutual exclusivity relationships. This
# is to illustrate how their inferred order of restrictions fits many more 
# evolutionary landscapes than the one they assume.
# 3. Simulate tumor progression and evolutionary trajectories under the different
# fitness landscapes to illustrate how differently tumors can progress even though
# the restrictions in the order of mutations are always maintained

#################################################################################
############# Colorectal adenocarcinoma (COADREAD) dataset ######################
#################################################################################

# Due to memory problems, the following genes from the dataset have been removed: 
# FAT4 (module E), CTNNB1 (module F), TCF7L2 (module E)

## MODEL 1: model from Neyshabouri et al + equal fitness effect for all relationships (s)
## + equal fitness effect for all non-satisfied relationships (sh) forcing absolute
## compliance with the DAG

# Dataframe with specification of restrictions
COADREAD_rT_M1 <- data.frame(parent = c("Root", "A", "B", "C", "D",
                                             "E", "F"),
                                  child = c("A", "B", "C", "D", "E", "F",
                                            "G"),
                                  s = 0.5,
                                  sh = -1,
                                  typeDep = "MN")

COADREAD_rT_M1

# Create fitness specifications from DAG and modules 
COADREAD_fitness_M1 <- allFitnessEffects(COADREAD_rT_M1, 
                                  geneToModule = c( "Root" = "Root",
                                                    "A" = "APC",
                                                    "B" = "TP53",
                                                    "C" = "KRAS",
                                                    "D" = "PIK3CA, NRAS, LRP1B",
                                                    "E" = "FBXW7, ARID1A",
                                                    "F" = "ATM, SMAD2, ERBB3, MTOR",
                                                    "G" = "SOX9, SMAD4"))

COADREAD_fitness_M1

# Plot DAG
plot(COADREAD_fitness_M1, expandModules = TRUE, autofit = TRUE)
dev.off()

# Evaluation of all possible genotypes fitness under the previous fitness specifications
COADREAD_FL_M1 <- evalAllGenotypes(COADREAD_fitness_M1, max = 131072)
COADREAD_FL_M1

# Plot fitness landscape
plotFitnessLandscape(COADREAD_FL_M1) 
dev.off()

# Simulate tumor progression under this fitness landscape (one simulation)

set.seed(2)

COADREAD_Simul_M1 <- oncoSimulIndiv(COADREAD_fitness_M1,
                                    model = "McFL", 
                                    mu = 1e-5,
                                    sampleEvery = 0.025,
                                    keepEvery = 1,
                                    initSize = 2000, 
                                    finalTime = 20000,
                                    keepPhylog = TRUE, 
                                    onlyCancer = FALSE)

# Plot clones diversity (genotypes)

plot(COADREAD_Simul_M1, show = "genotypes", type = "line", plotDiversity = TRUE, 
     legend.ncols = 5, thinData = TRUE, thinData.keep = 0.1)
plot(COADREAD_Simul_M1, show = "genotypes", type = "stacked", plotDiversity = TRUE,
     legend.ncols = 5, thinData = TRUE, thinData.keep = 0.1)
plot(COADREAD_Simul_M1, show = "genotypes", type = "stream", plotDiversity = TRUE,
     legend.ncols = 3, thinData = TRUE, thinData.keep = 0.1)

## MODEL 2: model from Neyshabouri et al. + equal fitness effect for all relationships (s)
## + equal fitness effect for all non-satisfied relationships (sh) without forcing absolute
## compliance with the DAG

# Less penalize deviations from the DAG of restrictions are based on the confidence
# in the pathway score implemented by Neyshabouri et al. in their model

# Dataframe with specification of restrictions
COADREAD_rT_M2 <- data.frame(parent = c("Root", "A", "B", "C", "D",
                                        "E", "F"),
                             child = c("A", "B", "C", "D", "E", "F",
                                       "G"),
                             s = 0.5,
                             sh = c(rep(-1, 4), rep(-.5, 2), -.2),
                             typeDep = "MN")

COADREAD_rT_M2

# Create fitness specifications from DAG and modules 
COADREAD_fitness_M2 <- allFitnessEffects(COADREAD_rT_M2, 
                                         geneToModule = c( "Root" = "Root",
                                                           "A" = "APC",
                                                           "B" = "TP53",
                                                           "C" = "KRAS",
                                                           "D" = "PIK3CA, NRAS, LRP1B",
                                                           "E" = "FBXW7, ARID1A",
                                                           "F" = "ATM, SMAD2, ERBB3, MTOR",
                                                           "G" = "SOX9, SMAD4"))

COADREAD_fitness_M2

# Plot DAG
plot(COADREAD_fitness_M2, expandModules = TRUE, autofit = TRUE)
dev.off()

# Evaluation of all possible genotypes fitness under the previous fitness specifications
COADREAD_FL_M1 <- evalAllGenotypes(COADREAD_fitness_M1, max = 131072)
COADREAD_FL_M1

# Plot fitness landscape
plotFitnessLandscape(COADREAD_FL_M1) 
dev.off()

# Simulate tumor progression under this fitness landscape (one simulation)

set.seed(2)

COADREAD_Simul_M2 <- oncoSimulIndiv(COADREAD_fitness_M2,
                                    model = "McFL", 
                                    mu = 1e-5,
                                    sampleEvery = 0.025,
                                    keepEvery = 1,
                                    initSize = 2000, 
                                    finalTime = 20000,
                                    keepPhylog = TRUE, 
                                    onlyCancer = FALSE)

# Plot clones diversity (genotypes)

plot(COADREAD_Simul_M2, show = "genotypes", type = "line", plotDiversity = TRUE, 
     legend.ncols = 5, thinData = TRUE, thinData.keep = 0.1)
plot(COADREAD_Simul_M2, show = "genotypes", type = "stacked", plotDiversity = TRUE,
     legend.ncols = 5, thinData = TRUE, thinData.keep = 0.1)
plot(COADREAD_Simul_M2, show = "genotypes", type = "stream", plotDiversity = TRUE,
     legend.ncols = 3, thinData = TRUE, thinData.keep = 0.1)

plotClonePhylog(COADREAD_Simul_M2, N = 1, keepEvents = TRUE)

#################################################################################
###############################  GMB dataset ####################################
#################################################################################

## MODEL 1: model from Neyshabouri et al + equal fitness effect for all relationships (s)
## + equal fitness effect for all non-satisfied relationships (sh) forcing absolute
## compliance with the DAG

# Dataframe with specification of restrictions
GMB_rT_M1 <- data.frame(parent = c("Root", "A", "B", "C", "D",
                                        "E", "F"),
                             child = c("A", "B", "C", "D", "E", "F",
                                       "G"),
                             s = 0.5,
                             sh = -1,
                             typeDep = "MN")

GMB_rT_M1

# Create fitness specifications from DAG and modules 
GMB_fitness_M1 <- allFitnessEffects(GMB_rT_M1, 
                             geneToModule = c( "Root" = "Root",
                                               "A" = "NF1, EGFR",
                                               "B" = "PTEN",
                                               "C" = "TP53",
                                               "D" = "PIK3R1, PIK3CA",
                                               "E" = "RB1, IDH1, STAG2",
                                               "F" = "ATRX, LZTR1",
                                               "G" = "BCOR, DCAF12L2"))

GMB_fitness_M1

plot(GMB_fitness_M1,  expandModules = TRUE, autofit = TRUE)
dev.off()

GMB_FL_M1 <- evalAllGenotypes(GMB_fitness_M1, max = 8192)
GMB_FL_M1

plotFitnessLandscape(GMB_FL_M1)
dev.off()

## MODEL 2: model from Neyshabouri et al. + equal fitness effect for all relationships (s)
## + equal fitness effect for all non-satisfied relationships (sh) without forcing absolute
## compliance with the DAG

# Less penalize deviations from the DAG of restrictions are based on the confidence
# in the pathway score implemented by Neyshabouri et al. in their model

# Dataframe with specification of restrictions
GMB_rT_M2 <- data.frame(parent = c("Root", "A", "B", "C", "D",
                                   "E", "F"),
                        child = c("A", "B", "C", "D", "E", "F",
                                  "G"),
                        s = 0.5,
                        sh = c(rep(-1, 6), -.2),
                        typeDep = "MN")

GMB_rT_M2

# Create fitness specifications from DAG and modules 
GMB_fitness_M2 <- allFitnessEffects(GMB_rT_M2, 
                                    geneToModule = c( "Root" = "Root",
                                                      "A" = "NF1, EGFR",
                                                      "B" = "PTEN",
                                                      "C" = "TP53",
                                                      "D" = "PIK3R1, PIK3CA",
                                                      "E" = "RB1, IDH1, STAG2",
                                                      "F" = "ATRX, LZTR1",
                                                      "G" = "BCOR, DCAF12L2"))

GMB_fitness_M2

plot(GMB_fitness_M2,  expandModules = TRUE, autofit = TRUE)
dev.off()

GMB_FL_M2 <- evalAllGenotypes(GMB_fitness_M2, max = 8192)
GMB_FL_M2

plotFitnessLandscape(GMB_FL_M2)
dev.off()
