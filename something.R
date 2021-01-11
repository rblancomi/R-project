# Load some libraries
library(OncoSimulR)
library(parallel)
library(graph)
library(igraph)

# Set the working directory
setwd("/home/henry/Downloads/R_programming/Programming_exercise/Ex-2")

# DAGS derived from Raphael and Vandin, 2015: 10.1089/cmb.2014.0161
# In the paper, the authors alerady used modules so no relationships
# of mutual exclusivity will be defined for individual genes!

# Define the dataframe that will contain the modules relationships, fitness
# effects (s: respect the restrictions imposed in the DAG and sh: deviations
# from the restrictions in the graph. Positive values means an advantage in 
# in fitness, whereas negative values means disadvantage in fitness) and type 
# of restrictions in the order of mutations.

# We will use the function allFitnessEffects to set the restrictions in the
# poset/DAG and the function evalAllGenotypes to see the fitness effects of 
# genotypes.

# Additionally we will plot the DAG as a graph using the default plot function
# for DAG and the from the igraph(library).

## Colorectal cancer study from Wood et al., 2007 for eight mutated genes.

# Pilot: as in the vignette of the package
CRC_W <- allFitnessEffects(data.frame(parent = c("Root", "A", "B", "KRAS"),
                                      child = c("A", "B", "KRAS", "FBXW7"),
                                      s = 0.1,
                                      sh = -0.01,
                                      typeDep = "MN"), 
                           geneToModule = c("Root" = "Root",
                                            "A" = "APC, EPHA3, TCF7L2",
                                            "B" = "EVC2, PIK3CA, TP53", 
                                            "KRAS" = "KRAS",
                                            "FBXW7" = "FBXW7"))


plot(CRC_W, expandModules = TRUE, autofit = TRUE)

#plot(CRC_W, "igraph", layout = layout.reingold.tilford, expandModules = TRUE)

# Map to fitness
# Without order of restrictions
CRC_F <- evalAllGenotypes(CRC_W, order = FALSE, addwt = TRUE)

# Plot fitness landscapes

# Using directly the output from allFitnessEffects
plotFitnessLandscape(CRC_W, use_ggrepel = TRUE) #The same
# Using the output from evalAllGenotypes
plot(CRC_F, use_ggrepel = TRUE) #The same, but we will use this since you
# can specify the order of restrictions

# With order of restrictions from evalAllGenotypes
CRC_FO <- evalAllGenotypes(CRC_W, order = TRUE, addwt = TRUE, max = 109600)
# No difference observed since the order of restrictions is not specified 
# in allFitnessEffects as a vector

# Some simulations:
# For simulations we will only use the McFarland model (continuous-time,
# logistic-like, and death rate depends on population size), since is the more
# realistic. Also, initial population sizes were set at 2000 (this value ensures
# that most simulations reach cancer, see Diaz-Uriarte 2017). Ony one mutation 
# rate will be used: 1e-5, this value is in range to others previously used.
# IDEAL VALUES BUT I do not have the computational power to test this!!

# For reproducibility
set.seed(14)

CRC_W_S <- oncoSimulIndiv(CRC_W,
                          model = "McFL",
                          mu = 1e-5,
                          sampleEvery = 0.025, #small since we are using McFL
                          keepEvery = 1, #since we want a fast execution time
                          initSize = 2000, 
                          finalTime = 20000, #time to reach cancer
                          keepPhylog = TRUE, #keep genealogy of clones
                          onlyCancer = FALSE)

# ERROR if onlyCancer = TRUE: Hitted maxtries. Exiting.
# Hitting max tries is regarded as an error.

# Summary of simulation
summary(CRC_W_S)

# We will plot the information of the matrix that contains the abundances of all
# the clones (or genotypes) at each of the sampling periods.

par(mfrow = c(7, 1))
# Drivers
plot(CRC_W_S, type = "line", addtot = TRUE, lwdClone = 0.9, log = "", 
     thinData = TRUE, thinData.keep = 0.1, plotDiversity = TRUE,
     legend.ncols = 2)
plot(CRC_W_S, type = "stacked", thinData = TRUE, thinData.keep = 0.1,
     plotDiversity = TRUE, legend.ncols = 5)
plot(b1, type = "stream", legend.ncols = 5)

# Genotypes
plot(CRC_W_S, show = "genotypes", type = "line", plotDiversity = TRUE, 
     legend.ncols = 5, thinData = TRUE, thinData.keep = 0.1)
plot(CRC_W_S, show = "genotypes", type = "stacked", plotDiversity = TRUE,
     legend.ncols = 5, thinData = TRUE, thinData.keep = 0.1)
plot(CRC_W_S, show = "genotypes", type = "stream", plotDiversity = TRUE,
     legend.ncols = 3, thinData = TRUE, thinData.keep = 0.1)

# Plot phylogenetic relationships and of clones that exist, this is N = 1
plotClonePhylog(CRC_W_S, N = 1, keepEvents = TRUE)

##FIXATION FOR EPISTASIS? (ALTHOUGH BETTER FOR NOT CANCER DATASETS)
#############################Order of effects###################################

# Adding order of effects to modules wit arbitrary values

CRC_WOE <- allFitnessEffects(data.frame(parent = c("Root", "A", "B", "KRAS"),
                                      child = c("A", "B", "KRAS", "FBXW7"),
                                      s = 0.1,
                                      sh = -0.01,
                                      typeDep = "MN"), 
                             orderEffects = c(
                               "FBXW7 > KRAS > B > A" = -0.5,
                               "KRAS > FBXW7 > B > A" = -0.3,
                               "KRAS > B > FBXW7 > A" = -0.3,                               
                               "KRAS > B > A > FBXW7" = -0.4,
                               "KRAS > B > A" = -0.2,
                               "B > KRAS > A" = -0.2,
                               "B > A > KRAS" = -0.2,
                               "B > A" = -0.1,
                               "A > B" = 0.05),
                           geneToModule = c("Root" = "Root",
                                            "A" = "APC, EPHA3, TCF7L2",
                                            "B" = "EVC2, PIK3CA, TP53", 
                                            "KRAS" = "KRAS",
                                            "FBXW7" = "FBXW7"))

# THe plot of the DAG does not change

#Fitness landscape with order of effects

CRC_WOE_F <- evalAllGenotypes(CRC_WOE, order = TRUE, addwt = TRUE, max = 109600)
# Now there are different values for fitness

plot(CRC_WOE_F, use_ggrepel = TRUE) #OncoSimulR cannot plot order effects yet
## Limitation!! Note that in the fitness landscape without order of effects
# some genotypes are something like: "FBXW7, B, KRAS, A" This genotypes
# does not make sens since they do not respect the restriction of the DAG
# so there is a penalty in fitness, but then why are they local maxima??

# Simulation

set.seed(3)
CRC_WOE_S <- oncoSimulIndiv(CRC_WOE,
                          model = "McFL",
                          mu = 1e-5,
                          sampleEvery = 0.02, #small since we are using McFL
                          keepEvery = 1, #since we want a fast execution time
                          initSize = 2000, 
                          finalTime = 20000, #time to reach cancer
                          keepPhylog = TRUE, #keep genealogy of clones
                          onlyCancer = FALSE)

# Plots

# Drivers
plot(CRC_WOE_S, type = "line", addtot = TRUE, lwdClone = 0.9, log = "", 
     thinData = TRUE, thinData.keep = 0.1, plotDiversity = TRUE,
     legend.ncols = 2)

#Genotypes
plot(CRC_WOE_S, show = "genotypes", type = "line", plotDiversity = TRUE, 
     legend.ncols = 10, thinData = TRUE, thinData.keep = 0.1)

CRC_WOE_S <- oncoSimulIndiv(CRC_WOE,
                            model = "McFL",
                            mu = 1e-5,
                            sampleEvery = 0.02, #small since we are using McFL
                            keepEvery = 1, #since we want a fast execution time
                            initSize = 10, 
                            finalTime = 100, #time to reach cancer
                            keepPhylog = TRUE, #keep genealogy of clones
                            onlyCancer = TRUE) #ERROR

#############################Using other fitness values#########################
# We will use other fitness values to represent sign epistasis
# (see Diad-Uriarte, 2017: https://doi.org/10.1093/bioinformatics/btx663)
# Moreover, the fitness values assigned will be of a small magnitude, since 
# cancer can be explained as a large number of mutations with a small fitness
# increases that drive tumor progression. 
# I will increase the fitness for each modules as we go deeper in the DAG
# this idea assumes that the clonal population is always under some selective
# pressure, therefore fitness increase is necessary. Although this is not 
# entirely true since a new selective pressure may benefit previous mutations
# more than recently arisen mutations (This is just an idea!). 
# Furthermore, the genes in the modules are known as superdrivers (e.g. 
# TP53, KRAS, APC, see: https://doi.org/10.1371/journal.pone.0027136). 
# One possible limitation could be that not all the genes of a modules
# are superdrivers, therefore, their fitness is not the same! For example, 
# according to Sjöblom,. et al 2006: genes that are mutated in more than 10%
# of the study are: TP53, APC, KRAS, FBXW7; whereas these genes have a low
# mutation prevalence: EPHA3, TCF7L2. (However, this problem can be overcome
# if we do not use modules and specify epistasis for those genes that are
# mutually exclusive).
# Also, the penalty for genotypes not possible under the DAG is set to a higher 
# value.

CRC_W1 <- allFitnessEffects(data.frame(parent = c("Root", "A", "B", "KRAS"),
                                        child = c("A", "B", "KRAS", "FBXW7"),
                                        s = c(rep(0.5, 2), 0.55, 0.6),
                                        sh = -0.9,
                                        typeDep = "MN"),
                             geneToModule = c("Root" = "Root",
                                              "A" = "APC, EPHA3, TCF7L2",
                                              "B" = "EVC2, PIK3CA, TP53", 
                                              "KRAS" = "KRAS",
                                              "FBXW7" = "FBXW7"))
                         
# DAG plot
plot(CRC_W1, expandModules = TRUE, autofit = TRUE)

# Fitness landscape
CRC_W1_F <- evalAllGenotypes(CRC_W1, order = FALSE, addwt = TRUE, max = 1000)

# Plot the fitness landscape and save the plot as pdf (Does not work for me :c)
#pdf(file = "FL_1.pdf", height = 8, width = 6)
plot(CRC_W1_F, use_ggrepel = TRUE)
#dev.off

# Note how the fitness laandscape changes if we use a lower penalty value
# for the deviations from DAG

CRC_W2 <- allFitnessEffects(data.frame(parent = c("Root", "A", "B", "KRAS"),
                                       child = c("A", "B", "KRAS", "FBXW7"),
                                       s = c(rep(0.5, 2), 0.55, 0.6),
                                       sh = -0.1,
                                       typeDep = "MN"),
                            geneToModule = c("Root" = "Root",
                                             "A" = "APC, EPHA3, TCF7L2",
                                             "B" = "EVC2, PIK3CA, TP53", 
                                             "KRAS" = "KRAS",
                                             "FBXW7" = "FBXW7")) 

CRC_W2_F <- evalAllGenotypes(CRC_W2, order = FALSE, addwt = TRUE, max = 1000)
plot(CRC_W2_F, use_ggrepel = TRUE) #Not a visible change

CRC_W3 <- allFitnessEffects(data.frame(parent = c("Root", "A", "B", "KRAS"),
                                    child = c("A", "B", "KRAS", "FBXW7"),
                                    s = c(rep(0.5, 2), 0.55, 0.6),
                                    sh = -0.05,
                                    typeDep = "MN"),
                         geneToModule = c("Root" = "Root",
                                          "A" = "APC, EPHA3, TCF7L2",
                                          "B" = "EVC2, PIK3CA, TP53", 
                                          "KRAS" = "KRAS",
                                          "FBXW7" = "FBXW7"),
                         drvNames = c("TP53", "APC", "KRAS", "FBXW7"))

CRC_W3_F <- evalAllGenotypes(CRC_W3, order = FALSE, addwt = TRUE, max = 1000)
plot(CRC_W3_F, use_ggrepel = TRUE)

# Although the changes, visually, are subtle. Note that if we use a 
# higher value as penalty for deviations lead to a higher number of genotypes
# near the sinks of the landscape. So we can imagine a landscape with 
# multiple peaks and deeper valleys? Does this make sense biologically?
# On the other hand, if we use a smaller number for deviations from monotonicity
# then, we obtain some genotypes that are between local and minima. We can say
# that the landscape is more "smooth" although with a high degree of 
# ruggedness??

# Some simulations

set.seed(32)
CRC_W3_S <- oncoSimulIndiv(CRC_W3,
                            model = "McFL",
                            mu = 1e-5,
                            sampleEvery = 0.02, #small since we are using McFL
                            keepEvery = 1, #since we want a fast execution time
                            initSize = 2000, 
                            finalTime = 20000, #time to reach cancer
                            keepPhylog = TRUE, #keep genealogy of clones
                            onlyCancer = FALSE, 
                           detectionDrivers = 4)

# Plots

# Drivers
plot(CRC_W3_S, type = "line", addtot = TRUE, lwdClone = 0.9, log = "", 
     thinData = TRUE, thinData.keep = 0.1, plotDiversity = TRUE,
     legend.ncols = 2)

#Genotypes
plot(CRC_W3_S, show = "genotypes", type = "line", plotDiversity = TRUE, 
     legend.ncols = 10, thinData = TRUE, thinData.keep = 0.1)

set.seed(33)

CRC_W3_S1 <- oncoSimulIndiv(CRC_W3,
                           model = "McFL",
                           mu = 1e-5,
                           sampleEvery = 0.02, #small since we are using McFL
                           keepEvery = 1, #since we want a fast execution time
                           initSize = 2000, 
                           finalTime = 20000, #time to reach cancer
                           keepPhylog = TRUE, #keep genealogy of clones
                           onlyCancer = TRUE, 
                           detectionDrivers = 4)

# Plots

# Drivers
plot(CRC_W3_S1, type = "line", addtot = TRUE, lwdClone = 0.9, log = "", 
     thinData = TRUE, thinData.keep = 0.1, plotDiversity = TRUE,
     legend.ncols = 2)

#Genotypes
plot(CRC_W3_S1, show = "genotypes", type = "line", plotDiversity = TRUE, 
     legend.ncols = 10, thinData = TRUE, thinData.keep = 0.1)

############################Some epistasis######################################

# As explained in the vignette of OncoSimulR, the advantage of using modules
# is that genes in the same modules are mutually exclusive. This is, mutating
# mutating a two genes of the same module does not increase the fitness
# significantly. 

# As mentioned in Diaz-Uriarte, 2017 (10.1093/bioinformatics/btx663), reciprocal
# sign epistasis or synthetic lethality is important in cancer, although it cannot
# be properly represented in cancer progression models and DAGs.
# Also, Reciprocal sign epistasis also increases the ruggedness of the landscape
# (lead to multiple peaks). (MAYBE ASK RAMON OR MALEY IF THIS MAKES SENSE:
# adding synthetic lethality between modules who are produced only by a single
# parent (vertical relationship)?). It could that some of the no tsuperdriver
# genes of a modules when produced with another gene lead to synthetic lethality
# Nevertheless, let's do a silly example with synthetic lethality
# at the end of the DAG. We will use a lower penalty value
# to see if there is a change in ruggedness. 

CRC_W4 <- allFitnessEffects(data.frame(parent = c("Root", "A", "B", "KRAS"),
                                       child = c("A", "B", "KRAS", "FBXW7"),
                                       s = c(rep(0.5, 2), 0.55, 0.6),
                                       sh = -0.05,
                                       typeDep = "MN"),
                            epistasis = c("-KRAS : FBXW7" = 0.2,
                                          "KRAS : -FBXW7" = 0.3,
                                          "KRAS : FBXW7" = -0.1),
                            geneToModule = c("Root" = "Root",
                                             "A" = "APC, EPHA3, TCF7L2",
                                             "B" = "EVC2, PIK3CA, TP53", 
                                             "KRAS" = "KRAS",
                                             "FBXW7" = "FBXW7"),
                            drvNames = c("TP53", "APC", "KRAS", "FBXW7")) 

CRC_W4_F <- evalAllGenotypes(CRC_W4, order = FALSE, addwt = TRUE, max = 1000)

# Plot the fitness landscape
plot(CRC_W4_F, use_ggrepel = TRUE)

# What if we increase the value for synthetic lethality?

CRC_W5 <- allFitnessEffects(data.frame(parent = c("Root", "A", "B", "KRAS"),
                                       child = c("A", "B", "KRAS", "FBXW7"),
                                       s = c(rep(0.5, 2), 0.55, 0.6),
                                       sh = -0.05,
                                       typeDep = "MN"),
                            epistasis = c("-KRAS : FBXW7" = 0.2,
                                          "KRAS : -FBXW7" = 0.3,
                                          "KRAS : FBXW7" = -0.9),
                            geneToModule = c("Root" = "Root",
                                             "A" = "APC, EPHA3, TCF7L2",
                                             "B" = "EVC2, PIK3CA, TP53", 
                                             "KRAS" = "KRAS",
                                             "FBXW7" = "FBXW7")) 


CRC_W5_F <- evalAllGenotypes(CRC_W5, order = FALSE, addwt = TRUE, max = 1000)

# Plot the fitness landscape
plot(CRC_W5_F, use_ggrepel = TRUE)

# Synthetic lethality indeed increases the ruggedness of the landscape. 
# Also, if we use a lower value for the there are a low number of local max
# and some pathways are accessible from local min. Moreover, there are several
# neutral changes with intermediate fitness between the local max and min.
# On the other, hand, if we use a higher value for synthetic lethality,
# there are a more local maxima than any of the previous models. Also, note that
# some genotypes have a fitness lower than the local min, but they are not local
# min. DOES THIS MAKE SENSE? I do not think so!

# Some simulations

set.seed(34)

CRC_W4_S <- oncoSimulIndiv(CRC_W4,
                            model = "McFL",
                            mu = 1e-5,
                            sampleEvery = 0.02, #small since we are using McFL
                            keepEvery = 1, #since we want a fast execution time
                            initSize = 2000, 
                            finalTime = 20000, #time to reach cancer
                            keepPhylog = TRUE, #keep genealogy of clones
                            onlyCancer = FALSE, 
                            detectionDrivers = 4)

# Plots

# Drivers
plot(CRC_W4_S, type = "line", addtot = TRUE, lwdClone = 0.9, log = "", 
     thinData = TRUE, thinData.keep = 0.1, plotDiversity = TRUE,
     legend.ncols = 2)

#Genotypes
plot(CRC_W4_S, show = "genotypes", type = "line", plotDiversity = TRUE, 
     legend.ncols = 10, thinData = TRUE, thinData.keep = 0.1)



###################Fitness effects with an evolutionary model###################

# For this, I will use data without epistasis and order of effects: CRC_W3

# THe bozic model is a type of exponential growth model were birth rate is fixed
# and only death rate is allowed to change
CRC_W3_FEM <- evalAllGenotypes(CRC_W3, order = FALSE, addwt = TRUE,
                               model = "Bozic")

# Plot fitness landscape

plot(CRC_W3_FEM, use_ggrepel = TRUE) #ERROR!
# Error in data.frame: arguments imply differing number of rows: 256, 0

# There are obvious changes in the fitness vaues of genotyoes between CRC_W3_F
# and CRC_W3_FEM. Nevertheless, it seems that fitness landscapes cannot be plot!
# THe model used in evalAllGenotypes is McFarland (kind of logistic growth).

################################################################################
############ Colorectal data from TCGA with 14 genes ###########################
################################################################################

CRC_TCGA <- allFitnessEffects(data.frame(parent = c("Root", "1", "2", "3", "4"),
                                         child = c("1", "2", "3", "4", "ELF3"),
                                         s = 0.1,
                                         sh = -0.01,
                                         typeDep = "MN"),
                              geneToModule = c("Root" = "Root",
                                               "1" = "APC, FBXW7",
                                               "2" = "ACVR2A, FAM123B, PIK3CA, TCF7L2, TP53",
                                               "3" = "BRAF, KRAS, NRAS",
                                               "4" = "SMAD2, SMAD4, SOX9",
                                               "ELF3" = "ELF3"))

plot(CRC_TCGA, expandModules = TRUE, autofit = TRUE)

plot(CRC_TCGA, "igraph", layout = layout.reingold.tilford, expandModules = TRUE)

# Map to fitness

CRC_TCGA_F <- evalAllGenotypes(CRC_TCGA, order = FALSE, addwt = TRUE, max = 16384)

# Plot fitness landscapes (Run this only if you have enough computational power)
plot(CRC_TCGA_F, use_ggrepel = TRUE)


################################################################################
############ Glioblastoma from TCGA with 24 genes ##############################
################################################################################

GBC <- allFitnessEffects(data.frame(parent = c("Root", "1", "2",
                                                  "3", "4", "5"), 
                                       child = c("1", "2", "3", "4", "5", "FGFR1"),
                                       s = 0.1,
                                       sh = -0.01,
                                       typeDep = "MN"),
                            geneToModule = c("Root" = "Root",
                                             "1" = "BRAF, NF1, EGFR, PDGFRA
                                             , IDH1, CDKN2C, PIK3CDG",
                                             "2" = "CDKN2A, CDK4, FGFR3, RB1",
                                             "3" = "CDK6, PIK3CA, PIK3CB, PIK3CG,
                                             PIK3R1, PIK3R2, PTEN",
                                             "4" = "MDM2, MDM4, PIK3C2A, TP53",
                                             "5" = "ATRX, FGFR2, MET, PIK3C2B",
                                             "FGFR1" = "FGFR1"))

plot(GBC, expandModules = TRUE, autofit = TRUE)

plot(GBC, "igraph", layout = layout.reingold.tilford, expandModules = TRUE)

# Map to fitness

evalAllGenotypes(GBC_14, order = FALSE, addwt = TRUE, max = 268435456)

# Plot fitness landscapes

#Almost impossible, if not impossible, to plot with my computational resources
plotFitnessLandscape(GBC_14, use_ggrepel = TRUE)

################################################################################
# As explained in Diaz-Uriarte, 2017: "seven genes is probably close to the 
# upper limits of fitness landscapes that can be easily visualized and related
# to their true DAG. Also, with modules, the number of genotypes
# grows extremely fast: 2^m, where m is the number of genes. EVen more genotypes
# exist if order of effects are taken into aacount.
# Thus, only the colorectal data set from Raphael and Vandin, 2015 was used for
# simulations
################################################################################

