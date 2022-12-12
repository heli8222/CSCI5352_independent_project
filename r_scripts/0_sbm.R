# CSCI 5352 Clauset
# Independent Project
# Henry Li

# code file 0b
# 0_sbm.R
# generate a synthetic food web using the stochastic block model

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
source("./0_niche.R")

# make_sbm_web(): master function that generates SBM network
make_sbm_web <- function(z, M) {
  adjacency_matrix <- matrix(0, nrow = length(z), ncol = length(z)) # initialize adjacency matrix with all 0s
  for (i in 1:nrow(adjacency_matrix)) { # loop over all species in adjacency matrix
    for (j in 1:nrow(adjacency_matrix)) { # loop over all possible prey of species i in adjacency matrix
      i_group <- z[i] # get the trophic level of species i
      j_group <- z[j] # get the trophic level of species j
      if (runif(1, min = 0, max = 1) < M[i_group,j_group]) {adjacency_matrix[i,j] <- 1} # i eats j with the probability defined in the mixing matrix
    }
  }
  adjacency_matrix <- DFS(adjacency_matrix) # remove all cycles from the adjacency matrix
  adjacency_matrix <- lower_triangular(adjacency_matrix) # make the adjacency matrix lower triangular
  
  return(adjacency_matrix)
}

# correct_M(): function that scales M to generate a SBM with targeted connectance C
correct_M <- function(M, z, C) {
  n <- rep(0, max(z)) # initialize vector to store the number of species in each community
  for (i in 1:length(n)) { # loop over all communities in network
    n[i] <- length(which(z == i)) # get number of species in community
  }
  S <- sum(n) # get the total number of species in network
  M <- as.matrix(M) # convert M from a dataframe to a matrix
  alpha <- (C * S^2) / as.numeric(n %*% M %*% n) # calculate the correcting coefficient for M
  corrected_M <- alpha * M
  
  return(corrected_M)
}