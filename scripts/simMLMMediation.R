#######################################################
# This script contains the code needed to create the  #
# simMLMMediation function.                           #  
#######################################################

simMLMMediation <- function(
  mlm.mediation.model = "2-1-1", 
  total.sample.size = 1000, 
  number.clusters = 200,
  model.one.fixed.effect = c(M.X = .40), 
  model.two.fixed.effect = c(Y.M_CLUSTER_MEAN = .14, Y.M_CWC = -.39),
  conditional.icc = c(.17, .17)
) {
  
  # Check if the user specified one of the three correct mlm.mediation.model
  # types.
  if(!(mlm.mediation.model %in% c("2-1-1",
                                  "1-1-1",
                                  "2-2-1"))) {
    stop("User needs to specify mlm.mediation.model.")
  }
  
  # Check if total.sample.size is divisble by number.clusters 
  if((total.sample.size %% number.clusters) != 0) {
    stop("Total sample size needs to be divisble by the number of clusters.")
  }
  
  ##### Define objects that are used across all mediation models ##### 
  # Turn model parameter vectors into matrices
  userBx <- as.matrix(c(INT = 0, model.one.fixed.effect))
  userBm <- as.matrix(c(INT = 0, model.two.fixed.effect))
  
  ##### Create random-effects design matrix for models 1 and 2: Z
  clusters <- rep(
    1:number.clusters,
    each = total.sample.size / number.clusters 
  )
  
  Z <- model.matrix(~ as.factor(clusters) - 1)
  
  ###### Begin 2-1-1 Mediation Model simulation #####
  if(mlm.mediation.model == "2-1-1") {
    
    ##### Create Beta matrices #####  
    
    # Check to see if correct effects have been provided 
    betaXCheck <- sum(
      c("M.X", "M.X_MEAN") %in% names(model.one.fixed.effect)
    )
    
    betaMCheck <- sum(
      c("Y.M_CLUSTER_MEAN", "Y.M_CWC") %in% names(model.two.fixed.effect)
    )
    
    if(betaXCheck == 0) {
      stop("User has not provided the correct model one parameter names.")
    }
    
    if(betaMCheck < 2) {
      stop("User has not provided the correct model two parameter names.")
    }
    
    # Create X and M fixed-effect (beta) matrices
    Bx <- as.matrix(rep(0, nrow(userBx)),
                    nrow = nrow(userBx),
                    ncol = 1
    )
    
    Bm <- as.matrix(rep(0, nrow(userBm)),
                    nrow = nrow(userBm),
                    ncol = 1
    )
    
    # Organize the Bx and Bm matrices so we can identify the fixed-effects
    # of interest: Effect of X on Mediator
    Bx[1, 1] <- userBx["INT", 1]
    Bx[2, 1] <- userBx["M.X", 1]
    
    Bm[1, 1] <- userBm["INT", 1]
    Bm[2, 1] <- userBm["Y.M_CLUSTER_MEAN", 1]
    Bm[3, 1] <- userBm["Y.M_CWC", 1]
    
    ##### Generate model 1 M ~ X #####
    # Randomly generate X variables from a multivariate normal dist. 
    X <- rnorm(
      n = number.clusters,
      mean = 0,
      sd = 1
    )
    
    # Repeat the X matrix so that each unit in a cluster recieves the same 
    # X value
    X <- rep(X, each = (total.sample.size / number.clusters))
    
    # Add column vector of 1s for the intercept term
    X <- cbind(1, X)
    
    # Create fixed f(x) 
    xPred <- X%*%Bx
    
    # Using the variance of mPred and the user-provided conditional ICCs, 
    # determine the L1 and L2 error variances so that the total variance
    # of M = 1
    
    # Variance in the mediator not explained by XB
    mUnexpVar <- 1 - var(xPred)
    
    if(mUnexpVar <= 0) {
      stop("The fixed-effects of X are too large.")
    }
    
    # Conditional ICC for mediator -- ICC after controlling for X
    iccM <- conditional.icc[1]
    
    # Variable replacement to make the equations easier to work with
    z <- (1 - iccM) / iccM
    
    # Mediator Level 2 Error
    mL2ErrorVar <- (mUnexpVar) / (1 + z)
    
    # Mediator Level 1 Error
    mL1ErrorVar <- (mUnexpVar * z) / (1 + z)
    
    # Generate L1 and L2 (cluster) residuals 
    mL1Error <- rnorm(
      total.sample.size,
      mean = 0,
      sd = sqrt(mL1ErrorVar)
    ) %>%
      as.matrix()
    
    mL2Error <- rnorm(
      number.clusters,
      mean = 0,
      sd = sqrt(mL2ErrorVar)
    ) %>%
      as.matrix()
    
    # Generate the mediator variable 
    m <- X%*%Bx + Z%*%mL2Error + mL1Error 
    
    ##### Generate model 2 Y ~ M #####
    
    # Generate the cluster means for the mediator 
    mClusterMean <- tapply(
      m, 
      clusters,
      mean
    ) %>%
      rep(
        each = (total.sample.size / number.clusters)
      )
    
    # Center mediator using cluster means 
    mCWC <- m - mClusterMean
    
    # Create a design matrix that includes the mediator cluster mean and 
    # the CWC mediator 
    Xm <- cbind(1, mClusterMean, mCWC)
    
    # Create fixed f(x) 
    mPred <- Xm%*%Bm
    
    # Using the variance of mPred and the user-provided conditional ICCs, 
    # determine the L1 and L2 error variances so that the total variance
    # of Y = 1
    
    # Variance in the mediator not explained by XB
    yUnexpVar <- 1 - var(mPred)
    
    if(yUnexpVar <= 0) {
      stop("The fixed-effects for the mediator variables are too large.")
    }
    
    # Conditional ICC for y -- ICC after controlling for mediator effects
    iccY <- conditional.icc[2]
    
    # Variable replacement to make the equations easier to work with
    z <- (1 - iccY) / iccY
    
    # Mediator Level 2 Error
    yL2ErrorVar <- (yUnexpVar) / (1 + z)
    
    # Mediator Level 1 Error
    yL1ErrorVar <- (yUnexpVar * z) / (1 + z)
    
    # Generate L1 and L2 (cluster) residuals 
    yL1Error <- rnorm(
      total.sample.size,
      mean = 0,
      sd = sqrt(yL1ErrorVar)
    ) %>%
      as.matrix()
    
    yL2Error <- rnorm(
      number.clusters,
      mean = 0,
      sd = sqrt(yL2ErrorVar)
    ) %>%
      as.matrix()
    
    # Generate the mediator variable 
    y <- Xm%*%Bm + Z%*%yL2Error + yL1Error 
    
    # Collect the simulated data into one data frame: data
    data <- 
      data.frame(
        Y = y, 
        M = m,
        M_CLUST_MEAN = mClusterMean,
        M_CWC = mCWC, 
        X = X[, 2],
        CLUSTER = clusters
      ) 
  }
  
  ###### Begin 2-2-1 Mediation Model simulation #####
  if(mlm.mediation.model == "2-2-1") {
    
    ##### Create Beta matrices #####  
    
    # Check to see if correct effects have been provided 
    betaXCheck <- sum(
      c("M.X") %in% names(model.one.fixed.effect)
    )
    
    betaMCheck <- sum(
      c("Y.M") %in% names(model.two.fixed.effect)
    )
    
    if(betaXCheck != 1) {
      stop("User has not provided the correct model one parameters.")
    }
    
    if(betaMCheck != 1) {
      stop("User has not provided the correct model two parameter.")
    }
    
    # Create X and M fixed-effect (beta) matrices
    Bx <- as.matrix(rep(0, nrow(userBx)),
                    nrow = nrow(userBx),
                    ncol = 1
    )
    
    Bm <- as.matrix(rep(0, nrow(userBm)),
                    nrow = nrow(userBm),
                    ncol = 1
    )
    
    # Organize the Bx and Bm matrices so we can identify the fixed-effects
    # of interest: Effect of X on Mediator
    Bx[1, 1] <- userBx["INT", 1]
    Bx[2, 1] <- userBx["M.X", 1]
    
    Bm[1, 1] <- userBm["INT", 1]
    Bm[2, 1] <- userBm["Y.M", 1]
    
    ##### Generate model 1 M ~ X (Linear Regression Model) #####
    # Randomly generate X variables from a multivariate normal dist. 
    X <- rnorm(
      n = number.clusters,
      mean = 0,
      sd = 1
    )
    
    # Add column vector of 1s for the intercept term
    X <- cbind(1, X)
    
    # Create fixed f(x) 
    xPred <- X%*%Bx
  
    # Variance in the mediator not explained by XB.
    # Unconditional mediator variance is scaled to equal 1. 
    # This allows us to interpret user-provided parameters as 
    # standardized effect sizes (beta coefficients).
    mErrorVar <- 1 - var(xPred)
    
    if(mErrorVar <= 0) {
      stop("The fixed-effects of X are too large.")
    }
    
    # Generate mediator residuals
    mError <- rnorm(
      number.clusters,
      mean = 0,
      sd = sqrt(mErrorVar)
    ) %>%
      as.matrix()
    
    # Generate the mediator variable 
    m <- X%*%Bx + mError 
    
    ##### Generate model 2 Y ~ M #####
    
    # Create a design matrix Xm for Y ~ M
    m <- rep(m, each = (total.sample.size / number.clusters))
    Xm <- cbind(1, m)
    
    # Create fixed f(x) 
    mPred <- Xm%*%Bm
    
    # Using the variance of mPred and the user-provided conditional ICCs, 
    # determine the L1 and L2 error variances so that the total variance
    # of Y = 1
    
    # Variance in the mediator not explained by XB
    yUnexpVar <- 1 - var(mPred)
    
    if(yUnexpVar <= 0) {
      stop("The fixed-effects for the mediator variables are too large.")
    }
    
    # Conditional ICC for y -- ICC after controlling for mediator effects
    iccY <- conditional.icc
    
    # Variable replacement to make the equations easier to work with
    z <- (1 - iccY) / iccY
    
    # Mediator Level 2 Error
    yL2ErrorVar <- (yUnexpVar) / (1 + z)
    
    # Mediator Level 1 Error
    yL1ErrorVar <- (yUnexpVar * z) / (1 + z)
    
    # Generate L1 and L2 (cluster) residuals 
    yL1Error <- rnorm(
      total.sample.size,
      mean = 0,
      sd = sqrt(yL1ErrorVar)
    ) %>%
      as.matrix()
    
    yL2Error <- rnorm(
      number.clusters,
      mean = 0,
      sd = sqrt(yL2ErrorVar)
    ) %>%
      as.matrix()
    
    # Generate the mediator variable 
    y <- Xm%*%Bm + Z%*%yL2Error + yL1Error 
    
    # Collect the simulated data into one data frame: data
    data <- 
      data.frame(
        Y = y, 
        M = m,
        X = rep(X[, 2], each = (total.sample.size / number.clusters)),
        CLUSTER = clusters
      ) 
  }
  row.names(data) <- NULL
  return(data)
}

