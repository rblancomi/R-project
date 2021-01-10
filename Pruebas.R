#################################################################

## Packages requiered

library(OncoSimulR)
library(graph)
library(igraph)

###################################################################

##### IMPORTANT ##### READ THIS BEFORE #####

###############################################################################
#### Model derived from Gerstung et al. (2011)
#### https://doi.org/10.1371/journal.pone.0027136

# Specification of restrictions in the order of accumulation of mutations 
# using a DAG with colorectal cancer.

# The main advantage of this work is that the give an aproximate value
# to the fitness of each genotype, based on the frequency of apparence
# in cancer progression. Thus, it could give us an idea about the values to
# include when defining fitness, rather than ranomd values as I was using.

### In the paper is written the following...

# Under the assumption of
# an identical mutation rate per gene, an accumulation rate of 1 per
# year corresponds approximately to a fitness surplus of 2.6%
# (Methods Section). Therefore the high accumulation rates of
# TP53, KRAS, and APC can be explained by a fitness effect on the
# order of a few percent, compared to other mutations with lower
# fitness gains of order 10^-3 or 10^-4
# Thus, these critical genes,
# which also form the mountains in the mutation landscape
# , may act as ‘superdrivers’ that provide a higher fitness
# gain than other genes.

## Based on this values and the model propossed in figure 2A from the paper
## I construct the folloing DAG, fitness landscape and simulation from the
# fitness effects

## 1. Fitness effects
colrect <- allFitnessEffects(
  data.frame(parent = c(rep("Root", 3), "APC", rep("KRAS", 2),
                         "PIK3CA", "EVC2", "FBXW7", "TCF7L2"),
             child = c("APC", "TP53", "PIK3CA", "KRAS", "FBXW7",
                       rep("EVC2", 2), "TCF7L2", rep("EPHA3", 2)),
             s = c(rep(0.26, 2), 0.001, 0.26, rep(0.0001, 6)),
             sh = c(rep(-0.2, 2), -0.05, -0.2, rep(-0.02, 6)),
             typeDep = "MN"),
            drvNames = c("APC", "TP53", "KRAS", "PIK3CA", 
                                    "EVC2", "FBXW7", "TCF7L2", "EPHA3")
  )


## Plot the DAG of the fitnessEffects object
plot(colrect)

## Obtain the fitness effects of all genotypes
colrect_FL <- evalAllGenotypes(colrect)

## Plot the fitness landscape of all genotypes
plotFitnessLandscape(colrect_FL, use_ggrepel = TRUE, show_labels = TRUE
                      )

## 2. Simulate from it. We change several possible options. 

set.seed(154) ## Fix the seed, so we can repeat it

ep2 <- oncoSimulIndiv(colrect, model = "McFL",
                      mu = 9e-4,
                      sampleEvery = 0.02,
                      keepEvery = 1,
                      initSize = 500,
                      finalTime = 50000,
                      detectionDrivers = 3,
                      keepPhylog = TRUE,
                      onlyCancer = TRUE)

## Plot the model derived
plot(ep2, show = "genotypes", type = "stacked",
     plotDiversity = TRUE)

### Another simulation
set.seed(149) ## Fix the seed, so we can repeat it

ep3 <- oncoSimulIndiv(colrect, ## A fitnessEffects object
                      model = "McFL", ## Model used
                      mu = 7e-5, ## Mutation rate
                      sampleEvery = 0.01, ## How often the whole population is sampled
                      keepEvery = 1,
                      initSize = 100, ## Initial population time
                      finalTime = 1000,
                      detectionDrivers = 3,
                      keepPhylog = TRUE, ## Allow to see parent-child relationships
                      onlyCancer = FALSE
                      )

## Parent-child relationship derived from simulation
plotClonePhylog(ep3,
                #fixOverlap = TRUE,
                N = 0, ## Specify clones that exist
                keepEvents = TRUE ## Arrows showing how many times each clones appeared
                )
## Model derived from the simulation
plot(ep3, ## OncoSimulIndv model
     show = "genotypes",
     type = "stacked"
     #plotDiversity = TRUE ## Show a small plot on top with Shannon's diversity index
     )

########################################################################################

############################## pathTiMEx dataset #######################################

########################################################################################


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

################ Alternative simulation of the monotonic model #############
############################################################################

## vector with the mutation rate associated to each gene
## mutation rates are based on mutation frequency in cross-sectional data
## Module A, B and C contains APC, TP53 and KRAS respectively
## It's frequency in data is higher than the others genes. Thus,
## It's mutation frequency is higher. In contrast, the frequency
## of appearance in tumours of modules D and E is lower and therefore,
## its mutation rate

muvect <- c("A" = 1e-4, "B" = 1e-4, "C" = 1e-4, "D" = 1e-5, "E" = 1e-5) 


set.seed(113) ## Fix the seed, so we can repeat it
ep <- oncoSimulIndiv(Wood_colo, ## A fitnessEffects object
                     model = "McFL", ## Model used
                     mu = muvect, ## Mutation rate
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


####### ########Alternative simulation of the semimonotonic model #############
############################################################################

## We assume the same mutation rate vector as in the previous case
muvect <- c("A" = 1e-4, "B" = 1e-4, "C" = 1e-4, "D" = 1e-5, "E" = 1e-5) 


set.seed(113) ## Fix the seed, so we can repeat it
ep <- oncoSimulIndiv(Wood_colo, ## A fitnessEffects object
                     model = "McFL", ## Model used
                     mu = muvect, ## Mutation rate
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



###############################################################################
################################# Order Effects ###############################
###############################################################################

## To assess how does order affect in fitness, we will specify different fitness
## to a specific order in gene mutation. 
## For this purpose, the simplify model of pathTiMEx of colorectal cancer is used

## 1. Construct the fitness object
## pathTiMEx model includes a parameter called (lambda) which refers to the estimates
## of its waiting time rate. It gives information about when an event occurs. 
## A higher lambda value means an early event in cancer progression.

## From pathTiMEx model, it is obtained that: "APC (module A) is
## an early event (k = 9:5), followed by KRAS (module B) (k = 2:89) 
## and TP53 (module C) (k = 1:92).

## Thus, a higher fitness values is done to this order in cancer progression.
## While the other combination will get a lower value.

## As a first approach, only the three previous named genes are used

cc <- data.frame(parent = c(rep("Root", 2), "A"),
                 child = c("A", "C", "B"),
                 typeDep = "MN")

## Fitness object
cc_order <- allFitnessEffects(
  orderEffects = c("A > B > C" = 0.5, "B > A > C" = 0.2,
                   "B > C > A" = 0.1,
                   "B > C" = 0.2,
                   "C > B" = 0.1,
                   "B > A" = 0.1,
                   "A > B" = 0.3),

  geneToModule =
    c("A" = "A",
      "B" = "B",
      "C" = "C") )

## All genotypes derived
cc_order_geno <- evalAllGenotypes(cc_order, order = TRUE)

## To visualize the order
cc_order_geno

## Fitness landscape derived from genotypes ordered
plotFitnessLandscape(cc_order_geno) ## Error in to_Fitness_Matrix(x, max_num_genotypes = max_num_genotypes) : 
                                    ## We cannot deal with order effects yet.


##### #################### An alternative is to use epistasis ##################

## As before, the number of genes is reduced to A, B and C.
## Fitness values are given according to the restrictions derived from
## the generative model.

## 1 Fitness object defined using epistasias
cc_epi <- allFitnessEffects(epistasis = 
                              c("A: -B: -C" = 0.4,
                                "-A: B: -C" = -0.4, ## Presence of B in absence of A is penalized
                                "-A: -B: C" = 0.3,
                                "A: B: -C" = 0.8,
                                "A : B: C" = 1.4,
                                "-A: B: C" = 0.1, ## Presence of B in absence of A is penalized
                                "A : -B: C" = 0.5
                                ))

## 2 Genotypes derived from fitness defined with epistasia relationships
cc_epi_geno <- evalAllGenotypes(cc_epi)

## Visualization of fitness
cc_epi_geno

## Fitness landscape from this relationships
plot(cc_epi_geno)

## Now, this fitness landscape is used to model cancer progression

set.seed(2) ## Fix the seed, so we can repeat it
ep <- oncoSimulIndiv(cc_epi, ## A fitnessEffects object
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
## The only way to construct the DAG from epistasias
plotClonePhylog(ep,
                #fixOverlap = TRUE,
                N = 0, ## Specify clones that exist
                keepEvents = TRUE ## Arrows showing how many times each clones appeared
)




####################################################################################################################
