# CSCI 5352 Clauset
# Independent Project
# Henry Li

# code file 1
# 1_get_parameters.R
# import empirical food webs and obtain parameters needed to generate synthetic webs

# set working directories and load packages
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
source("./0_niche.R")
source("./0_sbm.R")
edge_files <- Sys.glob("../empirical_webs/*edges.csv")

#-----------------------------------------------------------------------------#
# build adjacency matrix for each empirical web
for (i in 1:length(edge_files)) { # loop over all webs
  edge_list <- read.csv(edge_files[i]) # import edge list
  edge_list <- as.matrix(edge_list[edge_list$Type == "Feeding", c(3,2)]) # keep just trophic interactions, keep just prey and consumer columns
  adj_mat <- matrix(0, max(edge_list), max(edge_list)) # initialize adjacency matrix of all 0s
  adj_mat[edge_list] <- 1 # assign 1 for all interactions in edge list
  isolated_species_to_remove <- which(colSums(adj_mat[]) + rowSums(adj_mat[]) == 0) # get vector of isolated species
  if (length(isolated_species_to_remove) > 0) {
    adj_mat <- adj_mat[-isolated_species_to_remove, -isolated_species_to_remove] # remove isolated species from adjacency matrix
  }
  identical_species_to_remove <- identicals(adj_mat) # get vector of trophically identical species
  while(length(identical_species_to_remove) > 0) { # continue removing while there are identical species
    adj_mat <- adj_mat[-identical_species_to_remove[1], -identical_species_to_remove[1]] # remove first trophically identical species from adjacency matrix
    identical_species_to_remove <- identicals(adj_mat)# update vector of trophically identical species
  }
  matrix_name <- paste("web_", i, "_adj_mat", sep = "") # rename the adjacency matrix for storage
  assign(matrix_name, adj_mat)
  write.csv(adj_mat, paste("../adjacency_matrices/", matrix_name, ".csv", sep = "")) # export adjacency matrix as .csv
}

# store all adjacency matrices in a list
adj_mat_list <- list(web_1_adj_mat,
                     web_2_adj_mat,
                     web_3_adj_mat,
                     web_4_adj_mat,
                     web_5_adj_mat,
                     web_6_adj_mat)

#-----------------------------------------------------------------------------#
# get parameters for NM
nm_pars <- data.frame(matrix(ncol = 3, nrow = length(adj_mat_list)))
colnames(nm_pars) <- c("Web ID", "Size", "Connectance")

for (i in 1:nrow(nm_pars)) {
  nm_pars$'Web ID'[i] <- i
  nm_pars$Size[i] <- nrow(adj_mat_list[[i]])
  nm_pars$Connectance[i] <- sum(adj_mat_list[[i]]) / (nrow(adj_mat_list[[i]]))^2
}

#-----------------------------------------------------------------------------#
# get parameters for SBM

# helper function: return a vector of trophic levels for each species given an adjacency matrix
make_trophic_levels <- function(adjacency_matrix) {
  A <- adjacency_matrix
  A <- A / rowSums(A)
  A[is.nan(A)] <- 0
  S <- nrow(A)
  tl <- solve(diag(S) - A, rep(1, S))
  return(tl)
}

# helper function: return a mixing matrix given an adjacency matrix and vector of trophic levels
make_mixing_matrix <- function(adjacency_matrix, tl_vector) {
  M <- as.data.frame(matrix(0, ncol = length(unique(tl_vector)), nrow = length(unique(tl_vector)))) # inialize mixing matrix
  for (i in 1:nrow(M)) { # loop over all rows of mixing matrix
    for (j in 1:nrow(M)) { # loop over all columns of mixing matrix 
      consumer_species <- which(tl_vector == i) # get vector of all species in trophic level i
      prey_species <- which(tl_vector == j) # get vector of all species in trophic level j
      ij_interactions <- sum(adjacency_matrix[consumer_species, prey_species]) # get number of feeding interactions of a species in trophic level i eating a species in trophic level j
      i_interactions <- sum(adjacency_matrix[consumer_species,]) # get number of all feeding interactions of a species in trophic level i
      ifelse(i_interactions == 0, M[i,j] <- 0, M[i,j] <- ij_interactions/i_interactions) # set M_ij as the fraction of these two numbers
    }
  }
  return(M)
}

sbm_pars <- data.frame(matrix(ncol = 4, nrow = length(adj_mat_list)))
colnames(sbm_pars) <- c("Web ID", "Communities", "Community Label Vector", "Mixing Matrix")

for (i in 1:nrow(sbm_pars)) {
  trophic_levels <- round(make_trophic_levels(adj_mat_list[[i]]), digits = 0)
  sbm_pars$'Web ID'[i] <- i
  sbm_pars$Communities[i] <- length(unique(trophic_levels))
  sbm_pars$'Community Label Vector'[i] <- list(trophic_levels)
  sbm_pars$'Mixing Matrix'[i] <- list(make_mixing_matrix(adj_mat_list[[i]], trophic_levels))
}

#-----------------------------------------------------------------------------#
# save parameter tables for NM and SBM for synthetic network generation
saveRDS(nm_pars, "../model_pars/nm_pars.rds")
saveRDS(sbm_pars, "../model_pars/sbm_pars.rds")