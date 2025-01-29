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
# 5. Phylogenetc Tree Construction
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
  BiocManager::install("msa)")
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

# b) Loading libraries.
library(Biostrings)
library(msa)
library(ape)
library(phangorn)
library(rstudioapi)

# c) Setting up working directory to wherever the R script is installed. This
# is where any output will be saved after the analysis is done and where input
# will be assumed to be stored.
setwd(dirname(getSourceEditorContext()$path))

# --------------------------- 2. Handling Input ------------------------------ #
# The input will be several orthologues of priogenic and non-priogenic proteins.
# The prion-like group includes PRNP, MAVS, and RIP3. There are two control 
# groups: a set of cancer genes, including APC, EGFR, and BRCA1, and a set of
# housekeeping genes, including GAPDH, ACTB, and RPLP0. Each orthology consists
# of the same 14 mammals, including Homo sapiens, and Caretta caretta 
# (the loggerhead sea turtle), which will be as the outgroup in every tree.

# a) Reading orthologues from files.
# Prion-like group. This is essentially the experimental group.
PRNP_ortho <- readDNAStringSet("OrthoSeqs\\PRNP_refseq_transcript.fasta")
MAVS_ortho <- readDNAStringSet("OrthoSeqs\\MAVS_refseq_transcript.fasta")
RIP3_ortho <- readDNAStringSet("OrthoSeqs\\RIPK3_refseq_transcript.fasta")
# Cancer control group. This group will act as a positive control, as 
# cancer genes are generally highly mutagenic and therefore genetically
# unstable.
APC_ortho <- readDNAStringSet("OrthoSeqs\\APC_refseq_transcript.fasta")
EGFR_ortho <- readDNAStringSet("OrthoSeqs\\EGFR_refseq_transcript.fasta")
BRCA1_ortho <- readDNAStringSet("OrthoSeqs\\BRCA1_refseq_transcript.fasta")
# Housekeeping control group. This group will act as a negative control, as
# housekeeping genes are generally very genetically stable due to their 
# essential nature.
GAPDH_ortho <- readDNAStringSet("OrthoSeqs\\GAPDH_refseq_transcript.fasta")
ACTB_ortho <- readDNAStringSet("OrthoSeqs\\ACTB_refseq_transcript.fasta")
RPLP0_ortho <- readDNAStringSet("OrthoSeqs\\RPLP0_refseq_transcript.fasta")

# ---------------------- 3. Multiple Sequence Alignment ---------------------- #
# These alignments are being performed in preparation for building the 
# phylogenetic tree - it is simply a necessary step in tree construction. An
# alignment will be produced for every set of orthologues.

# a) Making a method to produce alignments, as this blcok of code be called 
# repeatedly to produce an MSA for each gene orthology.
produceMSA <- function(ortho){
  # Non-default parameters:
  # .method - MUSCLE, due to its speed and accuracy, relative to other tools
  # .cluster - Neighbour joining, as upgma clustering assumeds rates of evolution
  # are equal, which defeats the purpose of this analysis.
  # .type - DNA, as the sequences to be analyzed are DNA.
  ortho_msa <- msa(inputSeqs = ortho, method = "Muscle", 
                     cluster = "neighborjoining", type = "DNA")
  return(ortho_msa)
}
# b) Producing MSAs for each orthology. This may take a while.
# Prion-like group.
PRNP_MSA <- produceMSA(PRNP_ortho)
MAVS_MSA <- produceMSA(MAVS_ortho)
RIP3_MSA <- produceMSA(RIP3_ortho)
# Cancer control group.
APC_MSA <- produceMSA(APC_ortho)
EGFR_MSA <- produceMSA(EGFR_ortho)
BRCA1_MSA <- produceMSA(BRCA1_ortho)
# Housekeeping control group.
GADPH_MSA <- produceMSA(GADPH_ortho)
ACTB_MSA <- produceMSA(ACTB_ortho)
RPLP0_MSA <- produceMSA(RPLP0_ortho)

# ---------------------- 4. Model Fitting and Selection ---------------------- #

# --------------------- 5. Phylogenetic Tree Construction -------------------- #


# ----------------- 4. Phylogenetic Analysis --------------------------------- #
# a) Fit Model
# Evaluating each model to see which will be best for fitting.
msa_alignment <- as.phyDat(protein_msa)
modelEval <- modelTest(object = msa_alignment, model = "all")
# Now that we have results from testing, select the best with which to proceed.
# Doing so by looking at the 5 parameters: AIC, BIC, AICc, AICw, and AICcw. 
# In general, lowest AIC, BIC, and AICc are better, while higher AICw and AICcw
# are better. Since there is really no way of knowing which is best, the
# workflow will pick the model with the best average ranking.
model_ranks <- data.frame(AIC = modelEval$AIC, BIC = modelEval$BIC,
                           AICc = modelEval$AICc, AICw = modelEval$AICw,
                           AICcw = modelEval$AICcw, row.names = modelEval$Model)
# Defining a method to convert data frame values to ranks. It will directly
# change model_ranks to allow easier vectorization.
rankParams <- function(par){
  # If w is in par, it is a weight parameter and so largest is best. Thus, we
  # want to sort in decreasing order.
  decr_flag = grepl(pattern = "w", x = par)
  model_ranks <<- model_ranks[order(model_ranks[,par], decreasing = decr_flag),]
  model_ranks[par] <<- 1:nrow(model_ranks)
}
params <- c("AIC", "BIC", "AICc", "AICw", "AICcw")
lapply(params, rankParams)
# Now, finding the average rank for each model across the five parameters
# and selecting the model with the best.
model_ranks$AVE <- rowMeans(model_ranks)
best_model <- rownames(model_ranks[which("AVE" == max("AVE")),])
# b) Construct Tree
# Now fitting the model.
best_fit <- pml_bb(msa_alignment, model = best_model)
# Using bootstrapping of 100.
tree_bs <- bootstrap.pml(best_tree, bs = 100, optNni = TRUE, 
                            control = pml.control(trace = 0))
# Choosing to represent the tree as a phylogram so that branch lengths have
# meaning. Only branches supported by 50% of bootstrapping runs are being
# plotted.
plotBS(midpoint(best_fit$tree), tree_bs, p = 50, type="p",
       main="Standard bootstrap")

# ------------------------- 5. Analyses -------------------------------------- # 
# Getting gene identifiers from protein sequences. These will  be used for
# protein ontology analysis using the database provided with the R script.
getUniprotID <- function(protein_info){
  protein_info[2]
}
protein_split <- strsplit(names(protein_sequences), split = "|", fixed = TRUE) 
protein_id <- sapply(X = protein_split, FUN = getUniprotID)
