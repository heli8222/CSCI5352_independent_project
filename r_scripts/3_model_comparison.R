# CSCI 5352 Clauset
# Independent Project
# Henry Li

# code file 3
# 3_model_comparison.R
# compare synthetic networks according to network properties
# see Williams and Martinez (2000) for network property definitions
# Williams and Martinez (2000): 10.1038/35004572

# set working directories and load packages
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
source("./0_niche.R")
source("./0_sbm.R")
web_1_emp <- as.matrix(read.csv("../adjacency_matrices/web_1_adj_mat.csv", row.names = 1))
web_2_emp <- as.matrix(read.csv("../adjacency_matrices/web_2_adj_mat.csv", row.names = 1))
web_3_emp <- as.matrix(read.csv("../adjacency_matrices/web_3_adj_mat.csv", row.names = 1))
web_4_emp <- as.matrix(read.csv("../adjacency_matrices/web_4_adj_mat.csv", row.names = 1))
web_5_emp <- as.matrix(read.csv("../adjacency_matrices/web_5_adj_mat.csv", row.names = 1))
web_6_emp <- as.matrix(read.csv("../adjacency_matrices/web_6_adj_mat.csv", row.names = 1))
web_1_nm <- readRDS("../synthetic_webs/web_1_nm.rds")
web_2_nm <- readRDS("../synthetic_webs/web_2_nm.rds")
web_3_nm <- readRDS("../synthetic_webs/web_3_nm.rds")
web_4_nm <- readRDS("../synthetic_webs/web_4_nm.rds")
web_5_nm <- readRDS("../synthetic_webs/web_5_nm.rds")
web_6_nm <- readRDS("../synthetic_webs/web_6_nm.rds")
web_1_sbm <- readRDS("../synthetic_webs/web_1_sbm.rds")
web_2_sbm <- readRDS("../synthetic_webs/web_2_sbm.rds")
web_3_sbm <- readRDS("../synthetic_webs/web_3_sbm.rds")
web_4_sbm <- readRDS("../synthetic_webs/web_4_sbm.rds")
web_5_sbm <- readRDS("../synthetic_webs/web_5_sbm.rds")
web_6_sbm <- readRDS("../synthetic_webs/web_6_sbm.rds")
web_1_sbm_c <- readRDS("../synthetic_webs/web_1_sbm_c.rds")
web_2_sbm_c <- readRDS("../synthetic_webs/web_2_sbm_c.rds")
web_3_sbm_c <- readRDS("../synthetic_webs/web_3_sbm_c.rds")
web_4_sbm_c <- readRDS("../synthetic_webs/web_4_sbm_c.rds")
web_5_sbm_c <- readRDS("../synthetic_webs/web_5_sbm_c.rds")
web_6_sbm_c <- readRDS("../synthetic_webs/web_6_sbm_c.rds")

#-----------------------------------------------------------------------------#
# define helper functions to calculate web properties

# T(): function to calculate the fraction of top species within a web
# input: adjacency matrix
# output: non-negative real number
T <- function(adjacency_matrix) {
  S <- nrow(adjacency_matrix) # store the number of species
  counter <- 0 # initialize counter for number of top species
  for (i in 1:S) { # loop through all species in web
    if (sum(adjacency_matrix[,i]) == 0) { # if no one eats species i, add it to counter
      counter <- counter + 1
    }
  }
  return(counter / S)
}

# I(): function to calculate the fraction of intermediate species within a web
# input: adjacency matrix
# output: non-negative real number
I <- function(adjacency_matrix) {
  S <- nrow(adjacency_matrix) # store the number of species
  counter <- 0 # initialize counter for number of intermediate species
  for (i in 1:S) { # loop through all species in web
    if (sum(adjacency_matrix[i,]) > 0 & sum(adjacency_matrix[,i]) > 0) { # if species i eats someone and someone eats species i, add it to counter
      counter <- counter + 1
    }
  }
  return(counter / S)
}

# B(): function to calculate the fraction of basal species within a web
# input: adjacency matrix
# output: non-negative real number
B <- function(adjacency_matrix) {
  S <- nrow(adjacency_matrix) # store the number of species
  counter <- 0 # initialize counter for number of basal species
  for (i in 1:S) { # loop through all species in web
    if (sum(adjacency_matrix[i,]) == 0) { # if species i eats no one, add it to counter
      counter <- counter + 1
    }
  }
  return(counter / S)
}

# GenSD(): function to calculate the SD of generality
# input: adjacency matrix
# output: non-negative real number
GenSD <- function(adjacency_matrix) {
  S <- nrow(adjacency_matrix) # store the number of species
  L <- sum(adjacency_matrix) # store the number of links in web
  genSD_vector <- rep(NA, S) # initialize vector to store generalities of species
  for (i in 1:S) { # loop through all species in web
    genSD_vector[i] <- S / L * sum(adjacency_matrix[i,]) # calculate normalized generality of species
  }
  return(sd(genSD_vector))
}

# VulSD(): function to calculate the SD of vulnerability
# input: adjacency matrix
# output: non-negative real number
VulSD <- function(adjacency_matrix) {
  S <- nrow(adjacency_matrix) # store the number of species
  L <- sum(adjacency_matrix) # store the number of links in web
  vulSD_vector <- rep(NA, S) # initialize vector to store generalities of species
  for (i in 1:S) { # loop through all species in web
    vulSD_vector[i] <- S / L * sum(adjacency_matrix[,i]) # calculate normalized vulnerability of species
  }
  return(sd(vulSD_vector))
}

# MxSim(): function to calculate the mean maximum similarity
# input: adjacency matrix
# output: non-negative real number
MxSim <- function(adjacency_matrix) {
  S <- nrow(adjacency_matrix) # store the number of species
  sim_matrix <- matrix(0, S, S) # initialize matrix to store pairwise similarities
  for (i in 1:S) { # loop through all species in web
    i_eats <- which(adjacency_matrix[i,] == 1) # initialize a vector to store all species that species i eats
    i_eaten <- which(adjacency_matrix[,i] == 1) # initialize a vector to store all species that eat species i
    i_neighbors <- unique(c(i_eats, i_eaten)) # initialize a vector to store all species that neighbor species i
    for (j in 1:S) {
      j_eats <- which(adjacency_matrix[j,] == 1) # initialize a vector to store all species that species j eats
      j_eaten <- which(adjacency_matrix[,j] == 1) # initialize a vector to store all species that eat species j
      j_neighbors <- unique(c(j_eats, j_eaten)) # initialize a vector to store all species that neighbor species j
      if (j == i) {
        sim_matrix[i,j] <- 0
      }
      else {
        sim_matrix[i,j] <- length(intersect(i_neighbors, j_neighbors)) / length(union(i_neighbors, j_neighbors)) # divide the number of common neighbors by the number of total neighbors
      }
    }
  }
  sim_vector <- rep(0, S) # initialize vector to store max pairwise similarities
  for (i in 1:S) { # loop through all species in web
    sim_vector[i] <- max(sim_matrix[i,]) # store the max pairwise similarity for each species
  }
  return(mean(sim_vector))
}

#-----------------------------------------------------------------------------#
# initialize tables for storing network properties
properties <- c("T", "I", "B", "GenSD", "VulSD", "MxSim")
web_column <- rep(0, 6)
col_names <- c("Property", "Web 1", "Web 2", "Web 3", "Web 4", "Web 5", "Web 6")
col_names_WM <- c("Property", "Web 2", "Web 5", "Web 6")

web_2_WM <- c(0.14, 0.59, 0.27, 1.13, 0.72, 0.49)
web_5_WM <- c(0.043, 0.85, 0.11, 1.08, 0.58, 0.69)
web_6_WM <- c(0.078, 0.75, 0.17, 1.09, 0.62, 0.62)

stat_table_nm <- data.frame(properties, web_column, web_column, web_column, web_column, web_column, web_column)
colnames(stat_table_nm) <- col_names

stat_table_sbm <- data.frame(properties, web_column, web_column, web_column, web_column, web_column, web_column)
colnames(stat_table_sbm) <- col_names

stat_table_sbm_c <- data.frame(properties, web_column, web_column, web_column, web_column, web_column, web_column)
colnames(stat_table_sbm_c) <- col_names

stat_table_emp <- data.frame(properties, web_column, web_column, web_column, web_column, web_column, web_column)
colnames(stat_table_emp) <- col_names

stat_table_WM <- data.frame(properties, web_2_WM, web_5_WM, web_6_WM)
colnames(stat_table_WM) <- col_names_WM

#------------------------------------------------------------------------------#
# calculate empirical network properties

# store webs in a list
emp_webs <- list(web_1_emp,
                 web_2_emp,
                 web_3_emp,
                 web_4_emp,
                 web_5_emp,
                 web_6_emp)

for (i in 1:length(emp_webs)) { # calculate properties for all webs
  stat_table_emp[1, i+1] <- T(emp_webs[[i]])
  stat_table_emp[2, i+1] <- I(emp_webs[[i]])
  stat_table_emp[3, i+1] <- B(emp_webs[[i]])
  stat_table_emp[4, i+1] <- GenSD(emp_webs[[i]])
  stat_table_emp[5, i+1] <- VulSD(emp_webs[[i]])
  stat_table_emp[6, i+1] <- MxSim(emp_webs[[i]])
}

#------------------------------------------------------------------------------#
# calculate synthetic network properties

# store webs in a list
nm_webs <- list(web_1_nm,
                web_2_nm,
                web_3_nm,
                web_4_nm,
                web_5_nm,
                web_6_nm)
sbm_webs <- list(web_1_sbm,
                 web_2_sbm,
                 web_3_sbm,
                 web_4_sbm,
                 web_5_sbm,
                 web_6_sbm)
sbm_c_webs <- list(web_1_sbm_c,
                   web_2_sbm_c,
                   web_3_sbm_c,
                   web_4_sbm_c,
                   web_5_sbm_c,
                   web_6_sbm_c)

for (i in 1:length(nm_webs)) { # calculate properties for all webs
  stat_table_nm[[1, i+1]] <- list(sapply(nm_webs[[i]], T)) # calculate T for all NM webs
  stat_table_sbm[[1, i+1]] <- list(sapply(sbm_webs[[i]], T)) # calculate T for all SBM webs
  stat_table_sbm_c[[1, i+1]] <- list(sapply(sbm_c_webs[[i]], T)) # calculate T for all SBM_c webs
  
  stat_table_nm[[2, i+1]] <- list(sapply(nm_webs[[i]], I)) # calculate I for all NM webs
  stat_table_sbm[[2, i+1]] <- list(sapply(sbm_webs[[i]], I)) # calculate I for all SBM webs
  stat_table_sbm_c[[2, i+1]] <- list(sapply(sbm_c_webs[[i]], I)) # calculate I for all SBM_c webs
  
  stat_table_nm[[3, i+1]] <- list(sapply(nm_webs[[i]], B)) # calculate B for all NM webs
  stat_table_sbm[[3, i+1]] <- list(sapply(sbm_webs[[i]], B)) # calculate B for all SBM webs
  stat_table_sbm_c[[3, i+1]] <- list(sapply(sbm_c_webs[[i]], B)) # calculate B for all SBM_c webs
  
  stat_table_nm[[4, i+1]] <- list(sapply(nm_webs[[i]], GenSD)) # calculate GenSD for all NM webs
  stat_table_sbm[[4, i+1]] <- list(sapply(sbm_webs[[i]], GenSD)) # calculate GenSD for all SBM webs
  stat_table_sbm_c[[4, i+1]] <- list(sapply(sbm_c_webs[[i]], GenSD)) # calculate GenSD for all SBM_c webs
  
  stat_table_nm[[5, i+1]] <- list(sapply(nm_webs[[i]], VulSD)) # calculate VulSD for all NM webs
  stat_table_sbm[[5, i+1]] <- list(sapply(sbm_webs[[i]], VulSD)) # calculate VulSD for all SBM webs
  stat_table_sbm_c[[5, i+1]] <- list(sapply(sbm_c_webs[[i]], VulSD)) # calculate VulSD for all SBM_c webs
  
  stat_table_nm[[6, i+1]] <- list(sapply(nm_webs[[i]], MxSim)) # calculate MxSim for all NM webs
  stat_table_sbm[[6, i+1]] <- list(sapply(sbm_webs[[i]], MxSim)) # calculate MxSim for all SBM webs
  stat_table_sbm_c[[6, i+1]] <- list(sapply(sbm_c_webs[[i]], MxSim)) # calculate MxSim for all SBM_c webs
}

#-----------------------------------------------------------------------------#
# save network property tables and statistics
saveRDS(stat_table_emp, "../web_properties/empirical_values.rds")
saveRDS(stat_table_WM, "../web_properties/WM_values.rds")
saveRDS(stat_table_nm, "../web_properties/nm_properties.rds")
saveRDS(stat_table_sbm, "../web_properties/sbm_properties.rds")
saveRDS(stat_table_sbm_c, "../web_properties/sbm_c_properties.rds")