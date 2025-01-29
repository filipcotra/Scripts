# BIOL702 Final Project
# Author: Filip Cotra
# ----------------------------- Phylogeny Workflow --------------------------- #
# Main Research Question: Is there a link between the evolutionary distance of
# a gene and prionic behaviour?
# Rationale: Evolutionary distance can indicate the genetic stability of a gene,
# with those that mutate more often having greater distance in a phylogeny. 
# While prions tend to have high rates of specific amino acids, their
# fundamentally unstructured and nonfunctional tendencies might make them more
# flexible. On the other hand, its possible that due to the dangerous nature of
# priogenicity, there is a distinct evolutionary pressure to limit changes 
# which would further 
# Main Steps:
# 1. Environment Set-Up
# 2. Handling Input
# 3. Multiple Sequence Alignment
# 4. Model Fitting and Selection
# 5. Phylogenetic Tree Construction
# 6. Producing Plots
# 7. Evolutionary Distance Analysis
# ------------------------- 1. Environment Set-Up ---------------------------- #
# This step  involves installing any necessary packages (if necessary), loading
# necessary libraries, and setting up the working directory.

# a) Checking for and installing required packages.
# BiocManager - This package is necessary to install some subsequent packages.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Biostring - This package will allow for DNA sequences to be read in a better
# format than strings, useful for multiple sequence alignments.
if (!require("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")
# msa - This package will be used to perform mulitple sequence alignments.
if (!require("msa", quietly = TRUE))
  BiocManager::install("msa")
# ape - This package will be used along with phangorn for tree building.
if (!require("ape", quietly = TRUE))
  install.packages("ape")
# phangorn - This package will be used for tree building. 
if (!require("phangorn", quietly = TRUE))
  install.packages("phangorn")
# rstudioapi - This package will be used to access the RStudio API for the
# purpose of identifying the source directory.
if (!require("rstudioapi", quietly = TRUE))
  install.packages("rstudioapi")
# stringr - This package will be used for the integration of regular 
# expressions.
if (!require("stringr", quietly = TRUE))
  install.packages("stringr")
# tidyr - This package will be used for data formatting ahead of statistical
# testing.
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
# car - This package will be used for some statistical testing.
if (!require("car", quietly = TRUE))
  install.packages("car")

# b) Loading libraries.
library(Biostrings)
library(msa)
library(ape)
library(phangorn)
library(rstudioapi)
library(stringr)
library(tidyr)
library(car)

# c) Setting up working directory to wherever the R script is installed. This
# is where any output will be saved after the analysis is done and where input
# will be assumed to be stored.
setwd(dirname(getSourceEditorContext()$path))

# --------------------------- 2. Handling Input ------------------------------ #
# The input will be several orthologues of priogenic and non-priogenic proteins.
# The prion-like group includes PRNP, MAVS, and RIP3. Thee control group is a 
# set of housekeeping genes, including GAPDH, ACTB, and RPLP0. Each orthology 
# consists of the same 14 mammals, including Homo sapiens, and Caretta caretta 
# (the loggerhead sea turtle), which will be as the outgroup in every tree.

# a) Reading orthologues from files.
# Defining a vector containing all the file paths. This will also serve as the
# names for the information in all susbequent lists, so it will be named as
# such.
ortho_names <- c("OrthoSeqs/PRNP_refseq_transcript.fasta",
                 "OrthoSeqs/MAVS_refseq_transcript.fasta",
                 "OrthoSeqs/RIPK3_refseq_transcript.fasta",
                 "OrthoSeqs/GAPDH_refseq_transcript.fasta",
                 "OrthoSeqs/ACTB_refseq_transcript.fasta",
                 "OrthoSeqs/RPLP0_refseq_transcript.fasta")
# Reading the files.
ortho_stringSet <- lapply(ortho_names, readDNAStringSet)

# b) Defining a function to edit the names of the string sets.
editLabels <- function(stringSet){
  currNames <- names(stringSet)
  # This lapply command will extract elements of each current name depending on
  # some conditionals. This code is specifically designed around the example
  # orthologies used, and would require editing to specifically tailor naming
  # conventions to a different set of species orthologues. These conditionals
  # just guide which indices to extract in order to retrieve the full scientific
  # name for each organism.
  newNames <- lapply(strsplit(currNames, " "), function(x) {
    if ("PREDICTED:" %in% x){
      if ("Gorilla" %in% x){
        paste(x[3:5], collapse = " ")
      }
      else{
        paste(x[3:4], collapse = " ")
      }
    }
    else{
      paste(x[2:3], collapse = " ")
    }
  })
  names(stringSet) <- newNames
  return(stringSet)
}

# c) Editing the names of the string sets. This will be beneficial for the 
# legibility of the eventual phylogenetic tree.
ortho_stringSet_renamed <- lapply(ortho_stringSet, editLabels)

# ---------------------- 3. Multiple Sequence Alignment ---------------------- #
# These alignments are being performed in preparation for building the 
# phylogenetic tree - it is simply a necessary step in tree construction. An
# alignment will be produced for every set of orthologues.

# a) Making a method to produce alignments, as this block of code be called 
# repeatedly to produce an MSA for each gene orthology.
produceMSA <- function(ortho){
  # Non-default parameters:
  # .method - MUSCLE, due to its speed and accuracy, relative to other tools
  # .cluster - Neighbour joining, as upgma clustering assumes rates of evolution
  # are equal, which defeats the purpose of this analysis.
  # .type - DNA, as the sequences to be analyzed are DNA.
  ortho_msa <- msa(inputSeqs = ortho, method = "Muscle", 
                     cluster = "neighborjoining", type = "DNA")
  return(ortho_msa)
}

# b) Producing MSAs for each orthology. This will take a while.
ortho_MSA <- lapply(ortho_stringSet_renamed, produceMSA)

# ---------------------- 4. Model Fitting and Selection ---------------------- #
# There are various models that can be used to describe the aligned sequences
# change over time, and picking the one that best fits with the data will ensure
# that their representation on the phylogenetic tree are as accurate as 
# possible. This process involves testing the models and selecting the best one
# according to various metrics.

# a) Converting MSA objects to phyDat. This is necessary for use of the 
# modelTest function in phangorn.
ortho_phyDat <- lapply(ortho_MSA, as.phyDat)

# b) Testing each phangorn model on the phyDat objects. This will take a while.
ortho_modelEval <- lapply(ortho_phyDat, modelTest, model = "all")

# c) Defining a set of functions to identify the best model for each orthology. 
# First, defining a vector containing the metrics used to rank the models.
evalMetrics <- c("AIC", "BIC", "AICc", "AICw", "AICcw")
# The first function will rank the metrics for each model. Lower AIC, BIC, and 
# AICc are better, while higher AICw and AICcw are better. No metric is
# definitively best, so looking at all of them should give the clearest picture.
rankModels <- function(modelInfo){
  # Looping over each metric to rank them.
  for (metric in evalMetrics){
    # If w is in par, it is a weight parameter and so largest is best. Thus, we
    # want to sort in decreasing order.
    decr_flag = grepl(pattern = "w", x = metric) 
    modelRanks <- modelInfo[order(modelInfo[,metric], decreasing = decr_flag),]
    modelRanks[metric] <- 1:nrow(modelRanks)
  }
  return(modelRanks)
}
# The second function will extract the model with the highest average rank.
extractBest <- function(modelEval){
  # Creating a data frame containing the information for each model.
  modelInfo <- data.frame(AIC = modelEval$AIC, BIC = modelEval$BIC,
                          AICc = modelEval$AICc, AICw = modelEval$AICw,
                          AICcw = modelEval$AICcw, row.names = modelEval$Model)
  # Ranking each model according to the different metrics.
  modelRanks <- rankModels(modelInfo)
  # Now, finding the average rank for each model across the five parameters
  # and selecting the model with the best.
  modelRanks$AVE <- rowMeans(modelRanks)
  bestModel <- rownames(modelRanks[which("AVE" == max("AVE")),])
  return(bestModel)
}

# d) Finding the best model for each orthology.
ortho_bestModels <- lapply(ortho_modelEval, extractBest)

# --------------------- 5. Phylogenetic Tree Construction -------------------- #
# This step will construct bootstrapped phylogenetic trees with branch lengths
# representing evolutionary distance. One tree will be produced for each 
# orthology using the respective best models determined in the previous step.

# a) Defining a set of functions to fit the models to the alignments, construct
# the phylogenetic trees, and convert the phylogenetic trees to midpoint 
# rooted trees.
fitModel <- function(orthoName){
  nameInd <- which(ortho_names == orthoName)
  # Fitting the model to the MSA.
  modelFit <- pml_bb(ortho_phyDat[[nameInd]], 
                     model = ortho_bestModels[[nameInd]])
  return(modelFit)
}
produceTree <- function(modelFit){
  # Creating the tree. Bootstrapping is set to 100, meaning that the tree will
  # be resampled 100 times to ensure only the most likely branches appear.
  phyloTree <- bootstrap.pml(modelFit, bs = 100)
  return(phyloTree)
}
# The last function just converts trees to midpoint rooted trees, which is 
# important for interpreting the tree. Midpoint rooting assumes that the two
# most distant taxa represent the evolutionary extremes of the tree, and so 
# their midpoint serves as a viable rooting point.
convertMidpoint <- function(fitModel){
  midpointTree <- midpoint(fitModel$tree)
  return(midpointTree)
}

# b) Fitting the best models to the multiple sequence alignments. This will
# take a while. The orthology names are the same as the file paths, so we can
# use that initial vector to vectorize this operation.
ortho_fitModel <- lapply(ortho_names, fitModel)

# c) Retrieving midpoint rooted trees.
ortho_midpointTree <- lapply(ortho_fitModel, convertMidpoint)

# d) Producing phylogenetic trees. This will take a while.
ortho_phyloTree <- lapply(ortho_fitModel, produceTree)

# --------------------------- 6. Producing Plots ----------------------------- #
# This step will just involve making finishing touches and exporting plots as 
# PNGs.

# a) Defining a set of functions: one to create phylogenetic tree plots and 
# export them to PDFs and one to change the colours of the Homo sapiens
# tree branch
exportPlot <- function(orthoName){
  # This extracts the name of the orthologous gene using a regular expression.
  geneName <- str_extract(orthoName, "[A-Z]*[0-9]*(?=_)")
  # This creates the file path with the corresponding gene name, placing it in
  # a subdirectory of the working directory.
  filePath <- sub(pattern = "%f", replacement = geneName, 
                  x = "phyloPlots/%f.png")
  nameInd <- which(ortho_names == orthoName)
  phyloTree <- ortho_phyloTree[[nameInd]]
  fitTree <- ortho_midpointTree[[nameInd]]
  branchColours <- ortho_branchColours[[nameInd]]
  # Exporting the plot to a file.
  png(file = filePath)
  # Creating the plot.
  plotBS(fitTree, phyloTree, p = 50, type = "p", digits = 3,
         edge.color = branchColours, main = geneName)
  add.scale.bar(length = 0.1)
  dev.off()
}
# This function will identify the Homo sapiens branch and colour it red for 
# legibility.
colourBranches <- function(tree){
  # Setting default colour for branches.
  branchColours <- rep("black", nrow(tree$edge))
  # Identifying the Homo sapiens branch. This code will fail if Homo sapiens is
  # not included in a set of orthologues.
  humanIndex <- which(tree$edge[, 2] == 
                        which(tree$tip.label == "Homo sapiens"))
  branchColours[humanIndex] <- "red"
  return(branchColours)
}

# c) Retrieving the branch colours for each midpoint tree.
ortho_branchColours <- lapply(ortho_midpointTree, colourBranches)

# d) Plotting the trees and saving the plots as PDFs. Representing the trees as 
# phylograms as their genetic distance is important for this analysis.
# Prion-like group.
lapply(ortho_names, exportPlot)

# ------------------ 7. Evolutionary Distance Analysis ----------------------- #
# This step will simply involve calculating the ratio of human-specific
# relative to the entire tree to evaluate the genetic stability of the human
# gene groups. Statistical analyses will be performed to evaluate significance 
# and draw conclusions about the results. This step is the least transferable
# to other steps; depending on the question one is hoping to answer, they can
# use the constructed phylogenetic data to perform other analyses. 

# a) Defining a functions to extract the relative evolutionary distance of the 
# human branch. First, it will identify the human branch indices by referencing
# the branch colours, where the index of "red" is the human branch, and then
# find the relative length of the human branch to the cumulative length of the
# tree.
getEvolutionaryDistance <- function(orthoName){
  nameInd <- which(ortho_names == orthoName)
  currColours <- ortho_branchColours[[nameInd]]
  humanInd <- which(currColours == "red")
  currFit <- ortho_fitModel[[nameInd]]
  humanBranch <- currFit$tree$edge.length[humanInd]
  cumLength <- sum(currFit$tree$edge.length)
  return(humanBranch/cumLength)
}

# b) Getting the evolutionary distance of each human branch.
ortho_evolutionaryDistances <- lapply(ortho_names, getEvolutionaryDistance)

# c) Grouping the evolutionary distances by gene group. This step is very 
# specific to the amount of groups used and their order in the initial ortho 
# file list, so any re-use of this code would require editing.
prion_results <- unlist(ortho_evolutionaryDistances[1:3])
hk_results <- unlist(ortho_evolutionaryDistances[4:6])
evolutionaryDistance_results <- data.frame(Prion.Like = prion_results, 
                                           Housekeeping.Control = hk_results)
evolutionaryDistance_df <- pivot_longer(data = evolutionaryDistance_results,
                                        cols = c(Prion.Like, Housekeeping.Control),
                                        names_to = "Experimental.Group",
                                        values_to = "Evolutionary.Distance")


# d) Assessing significance between the groups using a t-test.
# Doing levene's test to check for equality of variances.
levene_result <- leveneTest(Evolutionary.Distance ~ Experimental.Group, 
                            data = evolutionaryDistance_df)
# Since the p-value is greater than 0.05, variance can be assumed to be equal
# between the two groups and the students' t-test can be used.
# Now, performing the T-test.
t.test_result <- t.test(prion_results, hk_results)
# It appears that there are no significant differences in the data set. This 
# implies that prion-like genes are not significantly more stable or unstable
# than housekeeping genes. 

# d) Producing and exporting a barplot for the t-test results.
result_means <- c(mean(hk_results), mean(prion_results))
result_sds <- c(sd(hk_results), sd(prion_results))
png(file = "evolutionaryDistance.png")
ave_barplot <- barplot(result_means, beside = TRUE, 
                        ylim = c(0, max(result_means + result_sds) * 1.2), 
                        col = c("skyblue", "orange"), 
                        names.arg = c("Housekeeping", "Prion-Like"), 
                        ylab = "Mean Evolutionary Distance", xlab = "Experimental Group")
arrows(x0 = ave_barplot, y0 = result_means - result_sds, x1 = ave_barplot, 
        y1 = result_means + result_sds, code = 3, angle = 90, length = 0.1, 
        col = "black"
) 
dev.off()
