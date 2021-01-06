#################################################################

library(OncoSimulR)

# Colorectal cancer dataset from  https://doi.org/10.1089/cmb.2016.0171

### First, we construct the DAG graph derived from pathTiMEx method

# Dataframe containing all genes and modules (collection of genes)
# in addition to fitness effect that applies if the relationship is satisfied (s)
# and fitness effect that applies if the relationship is not satisifed (sh)

# I suppose all relationships between genes and modules are monotonic.
# In this first model I am not covering mutual exclusivity in modules.

df <- data.frame(parent = c(rep("Root",3), "A", "B", "C"),
                 child = c("A", "B", "D", "C", "E", "E"),
                 s = c(0.15, 0.20, 0.50, 0.30, 0.17, 0.17),
                 sh = c(rep(0,3), -.2, rep(-.4, 2)),
                 typeDep = c(rep("--", 3), rep("MN",3))
)


## Function that returns fitness specifications for each genotype
# geneToModule is used to specify genes and modules

DAG_1 <- allFitnessEffects(df,
                             geneToModule = c("Root" = "Root",
                                              "A" = "APC",
                                              "B" = "TP53, EVC2",
                                              "C" = "KRAS",
                                              "D" = "PI3KCA, EPHA",
                                              "E" = "FBXW7, TCF7L2"))

# DAG reprsentation
plot(DAG_1, expandModules = T)

## Evaluation of different genotypes based on their genotype

DAG_1_FL <- evalAllGenotypes(DAG_1, max = 110000)

# Plot of fitness landscape (busy landscape !!!!! uggly hugh hugh)

plotFitnessLandscape(DAG_1_FL)


####### Second attempt, dividing modules in specific genes to show mutually 
# exclusivity between them.
# First, only with module "B" 


df2 <- data.frame(parent = c(rep("Root",4), "A", "B", "C", "E", "B"),
                  child = c("A", "B", "C", "D", "E", "F", "F", "F", "C"),
                  s = c(0.15, 0.20, 0.50, 0.30, 0.21, rep(0.17, 3), 0.50),
                  sh = -1,
                  typeDep = c("MN", "XMPN", "XMPN","MN", "MN", "MN", "MN", "MN", "XMPN")
)

# c("MN", "XMPN", "XMPN","MN", "MN", "MN", "MN", "MN", "XMPN", "XMPN")

DAG_2 <- allFitnessEffects(df2,
                           geneToModule = c("Root" = "Root",
                                            "A" = "APC",
                                            "B" = "TP53",
                                            "C" = "EVC2",
                                            "D" = "PI3KCA, EPHA",
                                            "E" = "KRAS",
                                            "F" = "FBXW7, TCF7L2"))
plot(DAG_2, expandModules = T)
#With 8 nodes, is not possible to execute 
### plot(DAG_2, expandModules = T)
### Error in `*tmp*`[[i]] : subíndice fuera de  los límites

DAG_2_FL <- evalAllGenotypes(DAG_2, max = 2500 )

plotFitnessLandscape(DAG_2_FL)
