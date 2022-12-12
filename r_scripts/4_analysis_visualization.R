# CSCI 5352 Clauset
# Independent Project
# Henry Li

# code file 4
# 4_analysis_visualization.R
# analyze and visualize data

# set working directories and load packages
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(igraph)
library(reshape2)
library(grid)
library(gridExtra)
source("../r_scripts/0_sbm.R")
emp_dat <- readRDS("../web_properties/empirical_values.rds")
WM_dat <- readRDS("../web_properties/WM_values.rds")
nm_dat <- readRDS("../web_properties/nm_properties.rds")
sbm_dat <- readRDS("../web_properties/sbm_properties.rds")
sbm_c_dat <- readRDS("../web_properties/sbm_c_properties.rds")
nm_pars <- readRDS("../model_pars/nm_pars.rds")
sbm_c_pars <- readRDS("../model_pars/sbm_pars.rds")
emp_web <- t(as.matrix(read.csv("../adjacency_matrices/web_2_adj_mat.csv", row.names = 1)))
colnames(emp_web) <- c(1:31)
nm_web <- readRDS("../synthetic_webs/web_2_nm.rds")[[5]]
sbm_c_web <- readRDS("../synthetic_webs/web_2_sbm_c.rds")[[1]]

#-----------------------------------------------------------------------------#
# prepare data for plotting and analysis

# convert Property to a factor in all datasets
emp_dat$Property <- factor(emp_dat$Property, levels = c("T", "I", "B", "GenSD", "VulSD", "MxSim"))
WM_dat$Property <- factor(WM_dat$Property, levels = c("T", "I", "B", "GenSD", "VulSD", "MxSim"))
nm_dat$Property <- factor(nm_dat$Property, levels = c("T", "I", "B", "GenSD", "VulSD", "MxSim"))
sbm_dat$Property <- factor(sbm_dat$Property, levels = c("T", "I", "B", "GenSD", "VulSD", "MxSim"))
sbm_c_dat$Property <- factor(sbm_c_dat$Property, levels = c("T", "I", "B", "GenSD", "VulSD", "MxSim"))

# unlist lists of web properties
nm_dat <- unnest(nm_dat, cols = c(`Web 1`, `Web 2`, `Web 3`, `Web 4`, `Web 5`, `Web 6`))
sbm_dat <- unnest(sbm_dat, cols = c(`Web 1`, `Web 2`, `Web 3`, `Web 4`, `Web 5`, `Web 6`))
sbm_c_dat <- unnest(sbm_c_dat, cols = c(`Web 1`, `Web 2`, `Web 3`, `Web 4`, `Web 5`, `Web 6`))

# unnest vectors of web properties
nm_dat <- unnest(nm_dat, cols = c(`Web 1`, `Web 2`, `Web 3`, `Web 4`, `Web 5`, `Web 6`))
sbm_dat <- unnest(sbm_dat, cols = c(`Web 1`, `Web 2`, `Web 3`, `Web 4`, `Web 5`, `Web 6`))
sbm_c_dat <- unnest(sbm_c_dat, cols = c(`Web 1`, `Web 2`, `Web 3`, `Web 4`, `Web 5`, `Web 6`))

# pivot longer
emp_dat <- pivot_longer(emp_dat,
                        cols = c(`Web 1`, `Web 2`, `Web 3`, `Web 4`, `Web 5`, `Web 6`),
                        names_to = 'Web',
                        values_to = 'Value')
WM_dat <- pivot_longer(WM_dat,
                       cols = c(`Web 2`, `Web 5`, `Web 6`),
                       names_to = 'Web',
                       values_to = 'Value')
nm_dat <- pivot_longer(nm_dat,
                       cols = c(`Web 1`, `Web 2`, `Web 3`, `Web 4`, `Web 5`, `Web 6`),
                       names_to = 'Web',
                       values_to = 'Value')
sbm_dat <- pivot_longer(sbm_dat,
                        cols = c(`Web 1`, `Web 2`, `Web 3`, `Web 4`, `Web 5`, `Web 6`),
                        names_to = 'Web',
                        values_to = 'Value')
sbm_c_dat <- pivot_longer(sbm_c_dat,
                          cols = c(`Web 1`, `Web 2`, `Web 3`, `Web 4`, `Web 5`, `Web 6`),
                          names_to = 'Web',
                          values_to = 'Value')

#-----------------------------------------------------------------------------#
# Question 1: Does my niche model replicate the results of the niche model described in Williams and Martinez 2000?

# Web 2
q1_w2 <- ggplot() +
  geom_boxplot(data = nm_dat[nm_dat$Web == "Web 2", ], aes(x = Property, y = Value)) +
  geom_point(data = WM_dat[WM_dat$Web == "Web 2", ], aes(x = Property, y = Value), color = 'red') +
  labs(title = "Validation of Niche Model Generated Webs", subtitle = "Chesapeake Bay", x = "Web Property", y = "Value") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

# Web 5
q1_w5 <- ggplot() +
  geom_boxplot(data = nm_dat[nm_dat$Web == "Web 5", ], aes(x = Property, y = Value)) +
  geom_point(data = WM_dat[WM_dat$Web == "Web 5", ], aes(x = Property, y = Value), color = 'red') +
  labs(title = "Validation of Niche Model Generated Webs", subtitle = "Little Rock Lake", x = "Web Property", y = "Value") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

# Web 6
q1_w6 <- ggplot() +
  geom_boxplot(data = nm_dat[nm_dat$Web == "Web 6", ], aes(x = Property, y = Value)) +
  geom_point(data = WM_dat[WM_dat$Web == "Web 6", ], aes(x = Property, y = Value), color = 'red') +
  labs(title = "Validation of Niche Model Generated Webs", subtitle = "St Martin Island", x = "Web Property", y = "Value") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

#-----------------------------------------------------------------------------#
# Question 2: Do the NM and connectance-corrected SBM generate network properties that differ in distribution?

# Kolmogorov-Smirnov 2 sample test
Web <- c("Web 1", "Web 2", "Web 3", "Web 4", "Web 5", "Web 6")
Property <- c("T", "I", "B", "GenSD", "VulSD", "MxSim")
q2_df <- as.data.frame(matrix(0, nrow = 6, ncol = 6))
colnames(q2_df) <- Web
rownames(q2_df) <- Property
for (i in 1:nrow(q2_df)) {
  for (j in 1:ncol(q2_df)) {
    nm_vector <- nm_dat[nm_dat$Web == Web[j] & nm_dat$Property == Property[i],]$Value
    sbm_c_vector <- sbm_c_dat[sbm_c_dat$Web == Web[j] & sbm_c_dat$Property == Property[i],]$Value
    q2_df[i,j] <- ks.test(nm_vector, sbm_c_vector)$statistic
  }
}
grid.table(q2_df)

# Web 1
q2_w1 <- ggplot() +
  geom_density(data = nm_dat[nm_dat$Web == "Web 1",], aes(x = Value, fill = "a"), alpha = 0.5) +
  geom_density(data = sbm_c_dat[sbm_c_dat$Web == "Web 1",], aes(x = Value, fill = "b"), alpha = 0.5) +
  geom_vline(data = emp_dat[emp_dat$Web == "Web 1",], aes(xintercept = Value), linetype = "dashed", color = "black", size = 1) +
  scale_fill_manual(name = "Legend", values = c("a" = "red", "b" = "blue", "c" = "black"), labels = c("a" = "NM", "b" = "SBM", "c" = "Empirical")) +
  facet_wrap(~Property, ncol = 3, scales = "free") +
  labs(title = "Comparison of Network Property Distributions", subtitle = "Bahia de Falsa San Quintin", x = "Value", y = "Density") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.text.y = element_blank())

# Web 2
q2_w2 <- ggplot() +
  geom_density(data = nm_dat[nm_dat$Web == "Web 2",], aes(x = Value, fill = "a"), alpha = 0.5) +
  geom_density(data = sbm_c_dat[sbm_c_dat$Web == "Web 2",], aes(x = Value, fill = "b"), alpha = 0.5) +
  geom_vline(data = emp_dat[emp_dat$Web == "Web 2",], aes(xintercept = Value), linetype = "dashed", color = "black", size = 1) +
  scale_fill_manual(name = "Legend", values = c("a" = "red", "b" = "blue", "c" = "black"), labels = c("a" = "NM", "b" = "SBM", "c" = "Empirical")) +
  facet_wrap(~Property, ncol = 3, scales = "free") +
  labs(title = "Comparison of Network Property Distributions", subtitle = "Chesapeake Bay", x = "Value", y = "Density") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.text.y = element_blank())

# Web 3
q2_w3 <- ggplot() +
  geom_density(data = nm_dat[nm_dat$Web == "Web 3",], aes(x = Value, fill = "a"), alpha = 0.5) +
  geom_density(data = sbm_c_dat[sbm_c_dat$Web == "Web 3",], aes(x = Value, fill = "b"), alpha = 0.5) +
  geom_vline(data = emp_dat[emp_dat$Web == "Web 3",], aes(xintercept = Value), linetype = "dashed", color = "black", size = 1) +
  scale_fill_manual(name = "Legend", values = c("a" = "red", "b" = "blue", "c" = "black"), labels = c("a" = "NM", "b" = "SBM", "c" = "Empirical")) +
  facet_wrap(~Property, ncol = 3, scales = "free") +
  labs(title = "Comparison of Network Property Distributions", subtitle = "Carpinteria Salt Marsh", x = "Value", y = "Density") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.text.y = element_blank())

# Web 4
q2_w4 <- ggplot() +
  geom_density(data = nm_dat[nm_dat$Web == "Web 4",], aes(x = Value, fill = "a"), alpha = 0.5) +
  geom_density(data = sbm_c_dat[sbm_c_dat$Web == "Web 4",], aes(x = Value, fill = "b"), alpha = 0.5) +
  geom_vline(data = emp_dat[emp_dat$Web == "Web 4",], aes(xintercept = Value), linetype = "dashed", color = "black", size = 1) +
  scale_fill_manual(name = "Legend", values = c("a" = "red", "b" = "blue", "c" = "black"), labels = c("a" = "NM", "b" = "SBM", "c" = "Empirical")) +
  facet_wrap(~Property, ncol = 3, scales = "free") +
  labs(title = "Comparison of Network Property Distributions", subtitle = "Estero de Punta Banda", x = "Value", y = "Density") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.text.y = element_blank())

# Web 5
q2_w5 <- ggplot() +
  geom_density(data = nm_dat[nm_dat$Web == "Web 5",], aes(x = Value, fill = "a"), alpha = 0.5) +
  geom_density(data = sbm_c_dat[sbm_c_dat$Web == "Web 5",], aes(x = Value, fill = "b"), alpha = 0.5) +
  geom_vline(data = emp_dat[emp_dat$Web == "Web 5",], aes(xintercept = Value), linetype = "dashed", color = "black", size = 1) +
  scale_fill_manual(name = "Legend", values = c("a" = "red", "b" = "blue", "c" = "black"), labels = c("a" = "NM", "b" = "SBM", "c" = "Empirical")) +
  facet_wrap(~Property, ncol = 3, scales = "free") +
  labs(title = "Comparison of Network Property Distributions", subtitle = "Little Rock Lake", x = "Value", y = "Density") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.text.y = element_blank())

# Web 6
q2_w6 <- ggplot() +
  geom_density(data = nm_dat[nm_dat$Web == "Web 6",], aes(x = Value, fill = "a"), alpha = 0.5) +
  geom_density(data = sbm_c_dat[sbm_c_dat$Web == "Web 6",], aes(x = Value, fill = "b"), alpha = 0.5) +
  geom_vline(data = emp_dat[emp_dat$Web == "Web 6",], aes(xintercept = Value), linetype = "dashed", color = "black", size = 1) +
  scale_fill_manual(name = "Legend", values = c("a" = "red", "b" = "blue", "c" = "black"), labels = c("a" = "NM", "b" = "SBM", "c" = "Empirical")) +
  facet_wrap(~Property, ncol = 3, scales = "free") +
  labs(title = "Comparison of Network Property Distributions", subtitle = "St Martin Island", x = "Value", y = "Density") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.text.y = element_blank())

#-----------------------------------------------------------------------------#
# Question 3: Which of the two models generates network properties that are closer to the empirically measured values?

Web <- c("Web 1", "Web 2", "Web 3", "Web 4", "Web 5", "Web 6")
Property <- c("T", "I", "B", "GenSD", "VulSD", "MxSim")
Model <- c("Niche Model", "Stochastic Block Model")
q3_df <- expand_grid(Web, Property, Model)
q3_df$Model <- factor(q3_df$Model, levels = c("Niche Model", "Stochastic Block Model"))
q3_df$z <- NA

for (i in 1:nrow(q3_df)) {
  ifelse(q3_df$Model[i] == "Niche Model", df <- nm_dat, df <- sbm_c_dat)
  value <- mean(df[df$Property == q3_df$Property[i] & df$Web == q3_df$Web[i],]$Value)
  mean <- emp_dat[emp_dat$Property == q3_df$Property[i] & emp_dat$Web == q3_df$Web[i],]$Value
  sd <- sd(df[df$Property == q3_df$Property[i] & df$Web == q3_df$Web[i],]$Value)
  ifelse (sd == 0, q3_df$z[i] <- 0, q3_df$z[i] <- (value - mean) / sd)
}

q3_nm_df <- data.frame(matrix(0, nrow = 6, ncol = 6))
rownames(q3_nm_df) <- Property
colnames(q3_nm_df) <- Web
for (i in 1:nrow(q3_nm_df)) {
  for (j in 1:ncol(q3_nm_df)) {
    q3_nm_df[i,j] <- q3_df[q3_df$Web == Web[j] & q3_df$Property == Property[i] & q3_df$Model == "Niche Model",]$z
  }
}
grid.table(q3_nm_df)

q3_sbm_c_df <- data.frame(matrix(0, nrow = 6, ncol = 6))
rownames(q3_sbm_c_df) <- Property
colnames(q3_sbm_c_df) <- Web
for (i in 1:nrow(q3_sbm_c_df)) {
  for (j in 1:ncol(q3_sbm_c_df)) {
    q3_sbm_c_df[i,j] <- q3_df[q3_df$Web == Web[j] & q3_df$Property == Property[i] & q3_df$Model == "Stochastic Block Model",]$z
  }
}
grid.table(q3_sbm_c_df)

# plot together
q3_hist <- ggplot(data = q3_df, aes(x = z)) +
  geom_histogram(color = "black", fill = "white", bins = 50) +
  facet_wrap(~Model, ncol = 1) +
  scale_y_continuous(breaks = c(2,4,6,8,10)) +
  labs(title = "Normalized Errors for Network Properties", x = "Normalized Error", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))

#-----------------------------------------------------------------------------#
# food web visualizations

# empirical web
emp_graph <- graph_from_adjacency_matrix(emp_web, mode = "directed")

# NM web
nm_graph <- graph_from_adjacency_matrix(t(nm_web), mode = "directed")

# connectance-corrected SBM web
sbm_c_graph <- graph_from_adjacency_matrix(t(sbm_c_web), mode = "directed")

layout <- layout_with_gem(sbm_c_graph)
par(mfrow = c(1,3))
par(mar = c(1,1,1,1))
plot(emp_graph, vertex.size = 12, edge.arrow.size = 0.15, layout = layout)
title("Empirical", cex.main = 1, col.main = "black", line = -12)
plot(nm_graph, vertex.size = 12, edge.arrow.size = 0.15, layout = layout)
title("Niche Model", cex.main = 1, col.main = "black", line = -12)
plot(sbm_c_graph, vertex.size = 12, edge.arrow.size = 0.15, layout = layout)
title("Stochastic Block Model", cex.main = 1, col.main = "black", line = -12)
title("Chesapeake Bay", cex.main = 2, col.main = "black", line = -10, outer = TRUE)

#-----------------------------------------------------------------------------#
# connectance-corrected SBM mixing matrices

# calculate corrected mixing matrices
M_df <- as.data.frame(matrix(NA, nrow = 6, ncol = 2))
colnames(M_df) <- c("Web ID", "Corrected Mixing Matrix")
for (i in 1:nrow(M_df)) {
  M_df$`Web ID`[[i]] <- i
  M_df$`Corrected Mixing Matrix`[i] <- list(as.data.frame(correct_M(sbm_c_pars$`Mixing Matrix`[[i]], sbm_c_pars$`Community Label Vector`[[i]], nm_pars$Connectance[i])))
}

# Web 1
w1_M <- as.matrix(M_df$`Corrected Mixing Matrix`[[1]])
colnames(w1_M) <- c("TL 1", "TL2", "TL3", "TL4")
rownames(w1_M) <- c("TL 1", "TL2", "TL3", "TL4")
w1_M <- melt(w1_M)
colnames(w1_M) <- c("Consumer", "Prey", "Value")
MM_w1 <- ggplot(data = w1_M, aes(x = Prey, y = Consumer, fill = Value)) +
  geom_tile() +
  labs(title = "BSQ") +
  theme(plot.title = element_text(hjust = 0.5, size = 10), axis.title.x = element_blank())

# Web 2
w2_M <- as.matrix(M_df$`Corrected Mixing Matrix`[[2]])
colnames(w2_M) <- c("TL 1", "TL2", "TL3", "TL4")
rownames(w2_M) <- c("TL 1", "TL2", "TL3", "TL4")
w2_M <- melt(w2_M)
colnames(w2_M) <- c("Consumer", "Prey", "Value")

MM_w2 <- ggplot(data = w2_M, aes(x = Prey, y = Consumer, fill = Value)) +
  geom_tile() +
  labs(title = "CHB") +
  theme(plot.title = element_text(hjust = 0.5, size = 10), axis.title.x = element_blank(), axis.title.y = element_blank())

# Web 3
w3_M <- as.matrix(M_df$`Corrected Mixing Matrix`[[3]])
colnames(w3_M) <- c("TL 1", "TL2", "TL3", "TL4")
rownames(w3_M) <- c("TL 1", "TL2", "TL3", "TL4")
w3_M <- melt(w3_M)
colnames(w3_M) <- c("Consumer", "Prey", "Value")

MM_w3 <- ggplot(data = w3_M, aes(x = Prey, y = Consumer, fill = Value)) +
  geom_tile() +
  labs(title = "CSM") +
  theme(plot.title = element_text(hjust = 0.5, size = 10), axis.title.x = element_blank(), axis.title.y = element_blank())

# web 4
w4_M <- as.matrix(M_df$`Corrected Mixing Matrix`[[4]])
colnames(w4_M) <- c("TL 1", "TL2", "TL3", "TL4", "TL5")
rownames(w4_M) <- c("TL 1", "TL2", "TL3", "TL4", "TL5")
w4_M <- melt(w4_M)
colnames(w4_M) <- c("Consumer", "Prey", "Value")

MM_w4 <- ggplot(data = w4_M, aes(x = Prey, y = Consumer, fill = Value)) +
  geom_tile() +
  labs(title = "EPB") +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

# Web 5
w5_M <- as.matrix(M_df$`Corrected Mixing Matrix`[[5]])
colnames(w5_M) <- c("TL 1", "TL2", "TL3", "TL4", "TL5")
rownames(w5_M) <- c("TL 1", "TL2", "TL3", "TL4", "TL5")
w5_M <- melt(w5_M)
colnames(w5_M) <- c("Consumer", "Prey", "Value")

MM_w5 <- ggplot(data = w5_M, aes(x = Prey, y = Consumer, fill = Value)) +
  geom_tile() +
  labs(title = "LRL") +
  theme(plot.title = element_text(hjust = 0.5, size = 10), axis.title.y = element_blank())

# Web 6
w6_M <- as.matrix(M_df$`Corrected Mixing Matrix`[[6]])
colnames(w6_M) <- c("TL 1", "TL2", "TL3", "TL4", "TL5")
rownames(w6_M) <- c("TL 1", "TL2", "TL3", "TL4", "TL5")
w6_M <- melt(w6_M)
colnames(w6_M) <- c("Consumer", "Prey", "Value")
MM_w6 <- ggplot(data = w6_M, aes(x = Prey, y = Consumer, fill = Value)) +
  geom_tile() +
  labs(title = "STM") +
  theme(plot.title = element_text(hjust = 0.5, size = 10), axis.title.y = element_blank())

SBM_c_pars <- grid.arrange(MM_w1, MM_w2, MM_w3, MM_w4, MM_w5, MM_w6,
                   nrow = 2,
                   top = textGrob("Estimated Mixing Matrices for Stochastic Block Model",
                   gp= gpar(fontsize = 15)))

#-----------------------------------------------------------------------------#
# NM parameter table
nm_pars$`Web Name` <- c("Bahia Falsa de San Quintin",
                        "Chesapeake Bay",
                        "Carpinteria Salt Marsh",
                        "Estero de Punta Banda",
                        "Little Rock Lake",
                        "St Martin Island")
nm_pars <- nm_pars[,c(1,4,2,3)]

grid.table(nm_pars)