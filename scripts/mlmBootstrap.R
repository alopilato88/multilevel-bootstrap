#######################################################
# This script contains the code needed to create the  #
# mlmBootstrap function.                              #  
#######################################################
test <- mlmBootstrap(
  data = data,
  nsims = 100,
  type = "parametric",
  func = indirectEffect,
  m1Formula = as.formula(M ~ X_CLUSTER_MEAN + X_CWC + (1|CLUSTER)),
  m2Formula = as.formula(Y ~ M_CLUSTER_MEAN + M_CWC + (1|CLUSTER)),
  iv = "X_CWC",
  med = "M_CWC",
  cluster = "CLUSTER"
)

mlmBootstrap <- function(
  data, 
  nsims = 1000,
  type,
  func,
  m1Formula,
  m2Formula,
  iv,
  med,
  cluster
) {
  
  if(!("magrittr" %in% .packages())) { library(magrittr) }
  if(!("dplyr" %in% .packages())) { library(dplyr) }
  if(!("lme4" %in% .packages())) { library(lme4) }
  
  if(type == "parametric") {
    
    m1 <- lme4::lmer(m1Formula, data = data)
    m2 <- lme4::lmer(m2Formula, data = data)
    
    # Predict mean response from m1 and m2 
    m1PredMean <- predict(m1)
    m2PredMean <- predict(m2)
    
    # Identify predictor names 
    m1PredictorNames <-
      row.names(summary(m1)$coefficients)[-1]
    
    m2PredictorNames <- 
      row.names(summary(m2)$coefficients)[-1]
    
    predictorNames <- c(
      m1PredictorNames,
      m2PredictorNames
    ) %>%
      unique()
    
    # Determine number of clusters
    clusterId <- data[, cluster] %>%
      unique()
    
    # Create a data frame to merge bootstrapped data with 
    bootData <- 
      data[, c(predictorNames, cluster)]
    
    bootData$M1_PRED_MEAN <- m1PredMean
    bootData$M2_PRED_MEAN <- m2PredMean
    
    # Save the random-effects covariance matrix from both models 
    m1Level2Cov <- summary(m1)$varcor[[1]]
    m2Level2Cov <- summary(m2)$varcor[[1]]
    
    # Save the level 1 variance components from m1 and m2 
    m1Level1Var <- attr(summary(m1)$varcor, "sc")
    m2Level1Var <- attr(summary(m2)$varcor, "sc")
    
    # Define a fomrula object to be passed to lme4 during bootstrap
    m1FormulaChar <- 
      Reduce(paste, deparse(m1@call$formula))
    
    m1Dv <- substr(m1FormulaChar, 1, 1)
    
    #m1Formula <- m1@call$formula
    
    m2FormulaChar <- 
      Reduce(paste, deparse(m2@call$formula))
    
    m2Dv <- substr(m2FormulaChar, 1, 1)

    #m2Formula <- m2@call$formula
    
    # Begin bootstrap for-loop here
    for(i in 1:nsims) {
      
      # Generate L2 and L1 residuals from the variance components 
      m1L2Residuals <- 
        mvtnorm::rmvnorm(
          n = length(clusterId),
          mean = rep(0, nrow(m1Level2Cov)),
          sigma = m1Level2Cov
        )
      
      m1L1Residuals <- 
        rnorm(
          n = nrow(data),
          mean = 0, 
          sd = sqrt(m1Level1Var)
        )
      
      m2L2Residuals <- 
        mvtnorm::rmvnorm(
          n = length(clusterId),
          mean = rep(0, nrow(m2Level2Cov)),
          sigma = m2Level2Cov
        )
      
      m2L1Residuals <- 
        rnorm(
          n = nrow(data),
          mean = 0,
          sd = sqrt(m2Level1Var)
        )
      
      # Create a residual data frame 
      l2ResidualData <- 
        data.frame(
          M1_L2_RESID = m1L2Residuals,
          M2_L2_RESID = m2L2Residuals,
          CLUSTER = clusterId
        )
      
      names(l2ResidualData)[3] <- cluster
      
      l1ResidualData <- 
        data.frame(
          M1_L1_RESID = m1L1Residuals,
          M2_L1_RESID = m2L1Residuals
        )
      
      bootData1 <- 
        bootData %>%
        dplyr::left_join(
          l2ResidualData,
          by = cluster
        ) %>%
        dplyr::bind_cols(
          l1ResidualData
        )
      
      bootData1$m1Dv <- 
        bootData1$M1_PRED_MEAN + 
        bootData1$M1_L2_RESID + 
        bootData1$M1_L1_RESID
      
      bootData1$m2Dv <- 
        bootData1$M2_PRED_MEAN + 
        bootData1$M2_L2_RESID + 
        bootData1$M2_L1_RESID

      names(bootData1)[which(names(bootData1) == "m1Dv")] <- m1Dv
      names(bootData1)[which(names(bootData1) == "m2Dv")] <- m2Dv
      
      bootResult <- 
        func(
          m1Formula = m1Formula,
          m2Formula = m2Formula,
          data = bootData1,
          iv = "X_CWC",
          med = "M_CWC"
        )
      
      if(i == 1) {
        bootVec <- bootResult
      } else {
        bootVec <- c(bootVec, bootResult)
      }
      
    }
  }
  
  return(bootVec)
}
