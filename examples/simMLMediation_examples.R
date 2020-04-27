###########################################################################
# This script contains examples of how to use the simMLMediation function #
# to generate data according to the: 2-1-1, 2-2-1, and 1-1-1 models.      #
###########################################################################

# Load packages that simMLMediation depends on 
library(magrittr)
library(dplyr)

# Load lme4 package to estimate linear mixed models using the simulated data
library(lme4)

# Source the simMLMMediation.R script 
source("scripts/simMLMMediation.R")

# Set seed to make results reproducible 
set.seed(1)

# Simulate data where the DGP is a 2-1-1 model 
data211 <- 
  simMLMMediation(
    mlm.mediation.model = "2-1-1", 
    total.sample.size = 6000, 
    number.clusters = 300,
    model.one.fixed.effect = c(M.X = .40), 
    model.two.fixed.effect = c(Y.M_CLUSTER_MEAN = .14, Y.M_CWC = -.39),
    conditional.icc = c(ICC_M = .15, ICC_Y = .10)
  )

m1_211 <- lmer(M ~ X + (1|CLUSTER), data211)
m2_211 <- lmer(Y ~ M_CLUSTER_MEAN + M_CWC + (1|CLUSTER), data211)

# Simulate data where the DGP is a 2-2-1 model
data221 <- 
  simMLMMediation(
    mlm.mediation.model = "2-2-1", 
    total.sample.size = 6000, 
    number.clusters = 300,
    model.one.fixed.effect = c(M.X = .40), 
    model.two.fixed.effect = c(Y.M = .30),
    conditional.icc = c(ICC_Y = .10)
  )

# m1_221 is a lm model, so we don't need to repeat observations within
# a cluster.
m1_221 <- lm(M ~ X, data221[!duplicated(data221$CLUSTER), ]) 
m2_221 <- lmer(Y ~ M + (1|CLUSTER), data221)

# Simulate data where the DGP is a 1-1-1 model
data111 <- 
  simMLMMediation(
    mlm.mediation.model = "1-1-1", 
    total.sample.size = 6000, 
    number.clusters = 300,
    model.one.fixed.effect = c(M.X_CLUSTER_MEAN = .40, M.X_CWC = .10), 
    model.two.fixed.effect = c(Y.M_CLUSTER_MEAN = .30, Y.M_CWC = .20),
    conditional.icc = c(ICC_X = .15, ICC_M = .10, ICC_Y = .10)
  )

m1_111 <- lmer(M ~ X_CLUSTER_MEAN + X_CWC + (1|CLUSTER), data111)
m2_111 <- lmer(Y ~ M_CLUSTER_MEAN + M_CWC + (1|CLUSTER), data111)




