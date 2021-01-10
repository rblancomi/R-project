#################################################################

library(OncoSimulR)
library(graph)
library(igraph)

# Colorectal cancer dataset from  https://doi.org/10.1089/cmb.2016.0171

### First, we construct the DAG graph derived from pathTiMEx method
### using dataset from Wood et al. 2007

# Dataframe containing all genes and modules (collection of genes)
# in addition to fitness effect that applies if the relationship is satisfied (s)
# and fitness effect that applies if the relationship is not satisifed (sh)

# I suppose all relationships between genes and modules are monotonic.

colcancer_monoto <- data.frame(parent = c(rep("Root",3), "A", "B", "C"),
                 child = c("A", "B", "D", "C", "E", "E"),
                 s = c(rep(0.5, 3), rep(0.05, 3)),
                 sh = -0.5,
                 typeDep = "MN"
)


## Function that returns fitness specifications for each genotype
# geneToModule is used to specify genes and modules

colcancer_monoto_efec <- allFitnessEffects(colcancer_monoto,
                             geneToModule = c("Root" = "Root",
                                              "A" = "APC",
                                              "B" = "TP53, EVC2",
                                              "C" = "KRAS",
                                              "D" = "PI3KCA, EPHA",
                                              "E" = "FBXW7, TCF7L2"),
                           drvNames = c("APC", "TP53", "EVC2", "KRAS",
                                        "PI3KCA", "EPHA", "FBXW7", "TCF7L2"))

# DAG reprsentation
plot(colcancer_monoto_efec, expandModules = TRUE, autofit = TRUE)

## Evaluation of different genotypes based on their genotype
colcancer_monoto_efec_FL <- evalAllGenotypes(DAG_1, max = 110000)

## Object with all possible genotypes
colcancer_monoto_efec_FL

# Plot of fitness landscape (busy landscape !!!!!)
plotFitnessLandscape(colcancer_monoto_efec_FL)


#### -------- #######

## Simplified version of Wood et al. colorectal cancer model derived from https://doi.org/10.1089/cmb.2016.0171
## Mutually exclusivity means that the mutation of a gene in the module
## is enough to affect fitness. Mutation of another gene from the same pathway
## does not improve or reduce fitness of the genotype.

## Thus, modules are not specified in function "evalAllGenotypes"
## In this case, modules are not expanded

## First, monotonic model.
Wood_colo <- allFitnessEffects(colcancer_monoto)

## DAG visualization
plot(Wood_colo, autofit = TRUE)

## Obtain all genotypes from the fitnessEffect
Wood_colo_geno <- evalAllGenotypes(Wood_colo)

## Plot the fitness landscape.
plotFitnessLandscape(Wood_colo_geno,
                     use_ggrepel = TRUE
                     )

## In this case, we can visualize better the fitness landscape, allowing to 
## know how does a genotype really affect in the cancer progression.

set.seed(257) ## Fix the seed, so we can repeat it

## A simulation is performed using function "oncoSimulIndiv"
ep <- oncoSimulIndiv(Wood_colo, ## A fitnessEffects object
                      model = "McFL", ## Model used
                      mu = 1e-4, ## Mutation rate
                      sampleEvery = 0.02, ## How often the whole population is sampled
                      keepEvery = 1,
                      initSize = 400, ## Initial population time
                      finalTime = 220,
                      keepPhylog = TRUE, ## Allow to see parent-child relationships
                      onlyCancer = FALSE
)

### All three posibble representation of the model are used

## Plot of simulation
plot(ep, ## OncoSimulIndv model
     show = "genotypes",
     type = "stacked"
     #plotDiversity = TRUE ## Show a small plot on top with Shannon's diversity index
)

## Plot of simulation
plot(ep, ## OncoSimulIndv model
     show = "genotypes",
     type = "line"
     #plotDiversity = TRUE ## Show a small plot on top with Shannon's diversity index
)

## Plot of simulation
plot(ep, ## OncoSimulIndv model
     show = "genotypes",
     type = "stream"
     #plotDiversity = TRUE ## Show a small plot on top with Shannon's diversity index
)

## Parent-child relationship derived from simulation
## This could be used as the DAG derived from the implementation of the model
plotClonePhylog(ep,
                #fixOverlap = TRUE,
                N = 0, ## Specify clones that exist
                keepEvents = TRUE ## Arrows showing how many times each clones appeared
)



###############################################################################
##### Slightly different model ############

## From the restrictions (DAG) from https://doi.org/10.1089/cmb.2016.0171 (previous DAG)
# I suppose all relationships between genes and modules are monotonic
# except for the last relationship between C and B with E.
colcancer_semi <- data.frame(parent = c(rep("Root",3), "A", "B", "C"),
                               child = c("A", "B", "D", "C", "E", "E"),
                               s = c(rep(0.5, 3), rep(0.05, 3)),
                               sh = -0.5,
                               typeDep = c(rep("MN", 4), rep("SM", 2))
)

## Simplified version of Wood et al. colorectal cancer.
Wood_colo <- allFitnessEffects(colcancer_semi)

## DAG visualization
plot(Wood_colo)

## Obtain all genotypes from the fitnessEffect
Wood_colo_geno <- evalAllGenotypes(Wood_colo)

## Plot the fitness landscape.
plotFitnessLandscape(Wood_colo_geno,
                     use_ggrepel = TRUE
)

set.seed(257) ## Fix the seed, so we can repeat it

## A simulation is performed using function "oncoSimulIndiv"
ep <- oncoSimulIndiv(Wood_colo, ## A fitnessEffects object
                     model = "McFL", ## Model used
                     mu = 1e-4, ## Mutation rate
                     sampleEvery = 0.02, ## How often the whole population is sampled
                     keepEvery = 1,
                     initSize = 400, ## Initial population time
                     finalTime = 220,
                     keepPhylog = TRUE, ## Allow to see parent-child relationships
                     onlyCancer = FALSE
)

## Plot of simulation
plot(ep, ## OncoSimulIndv model
     show = "genotypes",
     type = "stacked"
     #plotDiversity = TRUE ## Show a small plot on top with Shannon's diversity index
)

## Plot of simulation
plot(ep, ## OncoSimulIndv model
     show = "genotypes",
     type = "line"
     #plotDiversity = TRUE ## Show a small plot on top with Shannon's diversity index
)

## Plot of simulation
plot(ep, ## OncoSimulIndv model
     show = "genotypes",
     type = "stream"
     #plotDiversity = TRUE ## Show a small plot on top with Shannon's diversity index
)

## Parent-child relationship derived from simulation
plotClonePhylog(ep,
                #fixOverlap = TRUE,
                N = 0, ## Specify clones that exist
                keepEvents = TRUE ## Arrows showing how many times each clones appeared
)


####################################################################################################################
