
# Fetch relevant model and participant data from noisy_perception model
source('samples_model-fxns_basic.R')


### README ###
#' This file includes supporting functions for fitting lines to output from the samples_model.
#' These functions require first running the model using functions in `samples_model-fxns_basic.R`.
#' Follow the example in `samples_model-analysis_drift.R`.
#' 
#' Note many of the supporting functions are copied or adapted from:
#' `numberline/numexpt/R/fit.block-models.2014-06-02.R`
#' `numberline/numexpt/R/fit.models.2014-06-02.R`
#'  



###############
### GLOBALS ###
###############

BLOCKSIZE = 30 # number of trials per block


#################
### FUNCTIONS ###
#################

# bi-linear power-law mapping
map.bipower = function(x, a, b) {
  crit = a
  slope = 10 ^ b
  lx = log10(x)
  ly = ((lx > crit) * (crit + (lx - crit) * slope) + (lx <= crit) * lx);
  return(10 ^ ly)
}

# a trivially modified version of the above where `b` is not log transformed,
# allowing for computing CIs with `a` set to +/- `b`
map.bipower.ci = function(x, a, b) {
  crit = a
  slope = b
  lx = log10(x)
  ly = ((lx > crit) * (crit + (lx - crit) * slope) + (lx <= crit) * lx);
  return(10 ^ ly)
}

# general log likelihood function (with robustness)
loglik = function(x, y, map.fx, a, b, s) {
  sum(
    pmax(-6, dnorm(log10(y) - log10(map.fx(x, a, b)), 0, s, log = T))
  )
}

# Compute best fitting params for estimate data
brutefit = function(tmp) {
  nLL = function(a, b, s) {
    -loglik(tmp$num_dots, tmp$answer, map.bipower, a, b, 10 ^ s) + PRIORS[[1]](a) + PRIORS[[2]](b) + PRIORS[[3]](s)
  }
  
  iter = 0
  fits = NULL
  fit = NULL
  while (is.null(fits)) {
    try(fit <- summary(mle(nLL, # NB: this needs to be assigned with `<-` for some reason
                           start = list(a = runif(1, PARAMS["ma"], PARAMS["sa"]), 
                                      b = runif(1, PARAMS["mb"], PARAMS["sb"]), 
                                      s = rnorm(1, PARAMS["ms"], PARAMS["ss"])))), TRUE)
    iter = iter + 1
    
    if (!is.null(fit)) {
      # fits = c(tmp$subject[1], -0.5 * fit@m2logL, length(tmp$num_dots), fit@coef[,"Estimate"])
      fits = c(tmp$subject[1], -0.5 * fit@m2logL, length(tmp$num_dots), fit@coef[,"Estimate"], fit@coef[,"Std. Error"])
    } else {
      if (iter > 50) {
        # fits = c(tmp$subject[1], -9999, 0, 0, 0, 0)
        fits = c(tmp$subject[1], -9999, 0, 0, 0, 0, 0, 0, 0)
      }      
    }
  }
  # names(fits) = c("subject", "logL", "n", "a", "b", "s")
  names(fits) = c("subject", "logL", "n", "a", "b", "s", "se.a", "se.b", "se.s")
  
  return(fits)
}

# Util function used when fitting slopes: splits trial data into blocks
splitBlock = function(trial, n, total.trials) {
  floor((trial - 1) / (total.trials / n))
}

# Util function used when calculating correlation matrix
namedSlopes = function(x) {
  z = data.frame(x$b)
  rownames(z) = x$subject
  return(z)
}

# Util function used when calculating correlation matrix
cbind.fill = function(...) {
  nm = list(...) 
  rnames = unique(unlist(lapply(nm, rownames)))
  nm = lapply(nm, function(x) {
    newrows = rnames[!rnames %in% rownames(x)]
    newentries = matrix(nrow = length(newrows), ncol = ncol(x))
    rownames(newentries) = newrows
    colnames(newentries) = colnames(x)
    x = rbind(x, newentries)
    return(x)
  })
  nm = lapply(nm, function(x) {
    y = data.frame(x[order(as.numeric(rownames(x))),])
    rownames(y) = as.character(sort(as.numeric(rownames(x))))
    colnames(y) = colnames(x)
    return(y)
  })
  return(do.call(cbind, nm))
}

# Compute best fitting slopes for each participant x trial block of size `blocksizes` in `data`
fit.slopes = function(blocksizes, data) {
  n = blocksizes[1]
  fitsBlock = list()
  for (k in 1:n) {
    tmp = subset(data, splitBlock(data$trial, n, TRIALS) == (k - 1))
    fitsBlock[[k]] = data.frame(do.call(rbind, by(tmp, tmp$subject, brutefit)))
  }

  return(fitsBlock)
}

# Get m blocks by m blocks matrix of fitted slope correlations between each block
get.cor.matrix = function(slopes) {
  # Calculate slope correlations
  block.slopes = do.call(cbind.fill, lapply(slopes, namedSlopes)) # m subjects by n blocks slope values
  cor.matrix = cor(block.slopes, block.slopes, use = "pairwise.complete.obs") # n blocks by n blocks slope correlation matrix
  rownames(cor.matrix) = c() # NB: clearing out rownames and colnames is necessary for the processing below
  colnames(cor.matrix) = c()
  
  return(cor.matrix)
}

# Convert m blocks by m blocks matrix of fitted slope correlations to data frame
# Has column for each block and correlation between those blocks
# Initially has m^2 rows for each pair of blocks, then prunes out lower half
get.cor.df = function(cor.matrix) {
  slope.cor.df = reshape::melt(cor.matrix) # data frame with columns for block x, block y, and slope correlation b/n those blocks, n blocks x n blocks rows
  names(slope.cor.df) = c("block1", "block2", "slope.corr")
  slope.cor.df = slope.cor.df[slope.cor.df$block1 <= slope.cor.df$block2,] # remove redundant lower half of matrix
  slope.cor.df$slope.corr[slope.cor.df$block1 == slope.cor.df$block2] = NA # set correlation to NA in identical blocks
  
  return(slope.cor.df)
}

# Format data frame of correlations by block pair to show mean, se of correlations by block distance across blocks
# NB: keep dplyr:: function designations to ensure that summarize is based on grouping
# See https://stackoverflow.com/questions/26923862/why-are-my-dplyr-group-by-summarize-not-working-properly-name-collision-with
get.distance.cors = function(cor.df, blocksizes) {
  dist.df = cor.df %>%
    dplyr::mutate(block.dist = block2 - block1,
           trial.dist = (TRIALS / blocksizes[1]) * block.dist) %>%
    dplyr::group_by(trial.dist) %>%
    dplyr::summarize(mean.cor = mean(slope.corr, na.rm = TRUE),
              se.cor = sd(slope.corr, na.rm = TRUE) / sqrt(length(slope.corr)))

  return(dist.df[-1,]) # remove first row, will have NaN as mean.cor and NA as se.cor
}

