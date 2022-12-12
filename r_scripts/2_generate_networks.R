# CSCI 5352 Clauset
# Independent Project
# Henry Li

# code file 2
# 2_generate_networks.R
# generate synthetic webs for analysis using empirical parameters

# set working directories and load packages
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
source("./0_niche.R")
source("./0_sbm.R")
nm_pars <- readRDS("../model_pars/nm_pars.rds")
sbm_pars <- readRDS("../model_pars/sbm_pars.rds")

set.seed(0)

#-----------------------------------------------------------------------------#
# generate web 1 synthetic NM webs
web_1_nm <- list() # initialize list to store synthetic webs
S <- nm_pars$Size[1] # get size of web 1
C <- nm_pars$Connectance[1] # get connectance of web 1
while(length(web_1_nm) < 1000) { # generate 1000 synthetic webs of web 1
  new_web <- make_nm_web(S, 1, C)[[2]] # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0 ) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_nm_web(S, 1, C)[[2]]
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_1_nm[[length(web_1_nm) + 1]] <- new_web # append candidate web to list of synthetic webs
}

# generate web 2 synthetic NM webs
web_2_nm <- list()
S <- nm_pars$Size[2]
C <- nm_pars$Connectance[2]
while(length(web_2_nm) < 1000) { # generate 1000 synthetic webs of web 2
  if (length(web_2_nm) %% 100 == 0) {print(length(web_2_nm))}
  new_web <- make_nm_web(S, 1, C)[[2]] # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0 ) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_nm_web(S, 1, C)[[2]]
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_2_nm[[length(web_2_nm) + 1]] <- new_web # append candidate web to list of synthetic webs
}

# generate web 3 synthetic NM webs
web_3_nm <- list()
S <- nm_pars$Size[3]
C <- nm_pars$Connectance[3]
while(length(web_3_nm) < 1000) { # generate 1000 synthetic webs of web 3
  if (length(web_3_nm) %% 100 == 0) {print(length(web_3_nm))}
  new_web <- make_nm_web(S, 1, C)[[2]] # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0 ) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_nm_web(S, 1, C)[[2]]
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_3_nm[[length(web_3_nm) + 1]] <- new_web # append candidate web to list of synthetic webs
}

# generate web 4 synthetic NM webs
web_4_nm <- list()
S <- nm_pars$Size[4]
C <- nm_pars$Connectance[4]
while(length(web_4_nm) < 1000) { # generate 1000 synthetic webs of web 4
  if (length(web_4_nm) %% 100 == 0) {print(length(web_4_nm))}
  new_web <- make_nm_web(S, 1, C)[[2]] # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0 ) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_nm_web(S, 1, C)[[2]]
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_4_nm[[length(web_4_nm) + 1]] <- new_web # append candidate web to list of synthetic webs
}

# generate web 5 synthetic NM webs
web_5_nm <- list()
S <- nm_pars$Size[5]
C <- nm_pars$Connectance[5]
while(length(web_5_nm) < 1000) { # generate 1000 synthetic webs of web 5
  if (length(web_5_nm) %% 100 == 0) {print(length(web_5_nm))}
  new_web <- make_nm_web(S, 1, C)[[2]] # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0 ) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_nm_web(S, 1, C)[[2]]
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_5_nm[[length(web_5_nm) + 1]] <- new_web # append candidate web to list of synthetic webs
}

# generate web 6 synthetic NM webs
web_6_nm <- list()
S <- nm_pars$Size[6]
C <- nm_pars$Connectance[6]
while(length(web_6_nm) < 1000) { # generate 1000 synthetic webs of web 6
  if (length(web_6_nm) %% 100 == 0) {print(length(web_6_nm))}
  new_web <- make_nm_web(S, 1, C)[[2]] # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0 ) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_nm_web(S, 1, C)[[2]]
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_6_nm[[length(web_6_nm) + 1]] <- new_web # append candidate web to list of synthetic webs
}

#-----------------------------------------------------------------------------#
# generate web 1 synthetic SBM webs
web_1_sbm <- list()
z <- sbm_pars$`Community Label Vector`[[1]]
M <- sbm_pars$`Mixing Matrix`[[1]]
while(length(web_1_sbm) < 1000) { # generate 1000 synthetic webs of web 1
  if (length(web_1_sbm) %% 100 == 0) {print(length(web_1_sbm))}
  new_web <- make_sbm_web(z, M) # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_sbm_web(z, M)
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_1_sbm[[length(web_1_sbm) + 1]] <- new_web # append candidate web to list of synthetic webs
}

# generate web 2 synthetic SBM webs
web_2_sbm <- list()
z <- sbm_pars$`Community Label Vector`[[2]]
M <- sbm_pars$`Mixing Matrix`[[2]]
while(length(web_2_sbm) < 1000) { # generate 1000 synthetic webs of web 2
  if (length(web_2_sbm) %% 100 == 0) {print(length(web_2_sbm))}
  new_web <- make_sbm_web(z, M) # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_sbm_web(z, M)
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_2_sbm[[length(web_2_sbm) + 1]] <- new_web # append candidate web to list of synthetic webs
}

# generate web 3 synthetic SBM webs
web_3_sbm <- list()
z <- sbm_pars$`Community Label Vector`[[3]]
M <- sbm_pars$`Mixing Matrix`[[3]]
while(length(web_3_sbm) < 1000) { # generate 1000 synthetic webs of web 3
  if (length(web_3_sbm) %% 100 == 0) {print(length(web_3_sbm))}
  new_web <- make_sbm_web(z, M) # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_sbm_web(z, M)
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_3_sbm[[length(web_3_sbm) + 1]] <- new_web # append candidate web to list of synthetic webs
}

# generate web 4 synthetic SBM webs
web_4_sbm <- list()
z <- sbm_pars$`Community Label Vector`[[4]]
M <- sbm_pars$`Mixing Matrix`[[4]]
while(length(web_4_sbm) < 1000) { # generate 1000 synthetic webs of web 4
  if (length(web_4_sbm) %% 100 == 0) {print(length(web_4_sbm))}
  new_web <- make_sbm_web(z, M) # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_sbm_web(z, M)
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_4_sbm[[length(web_4_sbm) + 1]] <- new_web # append candidate web to list of synthetic webs
}

# generate web 5 synthetic SBM webs
web_5_sbm <- list()
z <- sbm_pars$`Community Label Vector`[[5]]
M <- sbm_pars$`Mixing Matrix`[[5]]
while(length(web_5_sbm) < 1000) { # generate 1000 synthetic webs of web 5
  if (length(web_5_sbm) %% 100 == 0) {print(length(web_5_sbm))}
  new_web <- make_sbm_web(z, M) # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_sbm_web(z, M)
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_5_sbm[[length(web_5_sbm) + 1]] <- new_web # append candidate web to list of synthetic webs
}

# generate web 6 synthetic SBM webs
web_6_sbm <- list()
z <- sbm_pars$`Community Label Vector`[[6]]
M <- sbm_pars$`Mixing Matrix`[[6]]
while(length(web_6_sbm) < 1000) { # generate 1000 synthetic webs of web 6
  if (length(web_6_sbm) %% 100 == 0) {print(length(web_6_sbm))}
  new_web <- make_sbm_web(z, M) # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_sbm_web(z, M)
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_6_sbm[[length(web_6_sbm) + 1]] <- new_web # append candidate web to list of synthetic webs
}

#-----------------------------------------------------------------------------#
# generate web 1 connectance-corrected synthetic SBM webs
web_1_sbm_c <- list()
target_C <- nm_pars$Connectance[1]
z <- sbm_pars$`Community Label Vector`[[1]]
M <- correct_M(sbm_pars$`Mixing Matrix`[[1]], z, target_C)
while(length(web_1_sbm_c) < 1000) { # generate 1000 synthetic webs of web 1
  if (length(web_1_sbm_c) %% 100 == 0) {print(length(web_1_sbm_c))}
  new_web <- make_sbm_web(z, M) # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_sbm_web(z, M)
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_1_sbm_c[[length(web_1_sbm_c) + 1]] <- new_web # append candidate web to list of synthetic webs
}

# generate web 2 connectance-corrected synthetic SBM webs
web_2_sbm_c <- list()
target_C <- nm_pars$Connectance[2]
z <- sbm_pars$`Community Label Vector`[[2]]
M <- correct_M(sbm_pars$`Mixing Matrix`[[2]], z, target_C)
while(length(web_2_sbm_c) < 1000) { # generate 1000 synthetic webs of web 2
  if (length(web_2_sbm_c) %% 100 == 0) {print(length(web_2_sbm_c))}
  new_web <- make_sbm_web(z, M) # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_sbm_web(z, M)
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_2_sbm_c[[length(web_2_sbm_c) + 1]] <- new_web # append candidate web to list of synthetic webs
}
# generate web 3 connectance-corrected synthetic SBM webs
web_3_sbm_c <- list()
target_C <- nm_pars$Connectance[3]
z <- sbm_pars$`Community Label Vector`[[3]]
M <- correct_M(sbm_pars$`Mixing Matrix`[[3]], z, target_C)
while(length(web_3_sbm_c) < 1000) { # generate 1000 synthetic webs of web 3
  if (length(web_3_sbm_c) %% 100 == 0) {print(length(web_3_sbm_c))}
  new_web <- make_sbm_web(z, M) # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_sbm_web(z, M)
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_3_sbm_c[[length(web_3_sbm_c) + 1]] <- new_web # append candidate web to list of synthetic webs
}
# generate web 4 connectance-corrected synthetic SBM webs
web_4_sbm_c <- list()
target_C <- nm_pars$Connectance[4]
z <- sbm_pars$`Community Label Vector`[[4]]
M <- correct_M(sbm_pars$`Mixing Matrix`[[4]], z, target_C)
while(length(web_4_sbm_c) < 1000) { # generate 1000 synthetic webs of web 4
  if (length(web_4_sbm_c) %% 100 == 0) {print(length(web_4_sbm_c))}
  new_web <- make_sbm_web(z, M) # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_sbm_web(z, M)
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_4_sbm_c[[length(web_4_sbm_c) + 1]] <- new_web # append candidate web to list of synthetic webs
}
# generate web 5 connectance-corrected synthetic SBM webs
web_5_sbm_c <- list()
target_C <- nm_pars$Connectance[5]
z <- sbm_pars$`Community Label Vector`[[5]]
M <- correct_M(sbm_pars$`Mixing Matrix`[[5]], z, target_C)
while(length(web_5_sbm_c) < 1000) { # generate 1000 synthetic webs of web 5
  if (length(web_5_sbm_c) %% 100 == 0) {print(length(web_5_sbm_c))}
  new_web <- make_sbm_web(z, M) # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_sbm_web(z, M)
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_5_sbm_c[[length(web_5_sbm_c) + 1]] <- new_web # append candidate web to list of synthetic webs
}
# generate web 6 connectance-corrected synthetic SBM webs
web_6_sbm_c <- list()
target_C <- nm_pars$Connectance[6]
z <- sbm_pars$`Community Label Vector`[[6]]
M <- correct_M(sbm_pars$`Mixing Matrix`[[6]], z, target_C)
while(length(web_6_sbm_c) < 1000) { # generate 1000 synthetic webs of web 6
  if (length(web_6_sbm_c) %% 100 == 0) {print(length(web_6_sbm_c))}
  new_web <- make_sbm_web(z, M) # generate a candidate web
  problem_species <- c(isolates(new_web), identicals(new_web)) # get vector of problem species in candidate web
  while (length(problem_species) > 0) { # keep generating new candidate webs as long as problem species are present
    new_web <- make_sbm_web(z, M)
    problem_species <- c(isolates(new_web), identicals(new_web))
  }
  web_6_sbm_c[[length(web_6_sbm_c) + 1]] <- new_web # append candidate web to list of synthetic webs
}

#-----------------------------------------------------------------------------#
# save synthetic webs for model comparison
saveRDS(web_1_nm, "../synthetic_webs/web_1_nm.rds")
saveRDS(web_2_nm, "../synthetic_webs/web_2_nm.rds")
saveRDS(web_3_nm, "../synthetic_webs/web_3_nm.rds")
saveRDS(web_4_nm, "../synthetic_webs/web_4_nm.rds")
saveRDS(web_5_nm, "../synthetic_webs/web_5_nm.rds")
saveRDS(web_6_nm, "../synthetic_webs/web_6_nm.rds")
saveRDS(web_1_sbm, "../synthetic_webs/web_1_sbm.rds")
saveRDS(web_2_sbm, "../synthetic_webs/web_2_sbm.rds")
saveRDS(web_3_sbm, "../synthetic_webs/web_3_sbm.rds")
saveRDS(web_4_sbm, "../synthetic_webs/web_4_sbm.rds")
saveRDS(web_5_sbm, "../synthetic_webs/web_5_sbm.rds")
saveRDS(web_6_sbm, "../synthetic_webs/web_6_sbm.rds")
saveRDS(web_1_sbm_c, "../synthetic_webs/web_1_sbm_c.rds")
saveRDS(web_2_sbm_c, "../synthetic_webs/web_2_sbm_c.rds")
saveRDS(web_3_sbm_c, "../synthetic_webs/web_3_sbm_c.rds")
saveRDS(web_4_sbm_c, "../synthetic_webs/web_4_sbm_c.rds")
saveRDS(web_5_sbm_c, "../synthetic_webs/web_5_sbm_c.rds")
saveRDS(web_6_sbm_c, "../synthetic_webs/web_6_sbm_c.rds")