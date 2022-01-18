
#'
#' Experiment analysis for estimation drift manuscript
#'


setwd(paste0(here::here(), "/vullab/estimation_drift/"))
# rm(list=ls())

library(cocor)
library(Matrix)
library(stats4)
library(tidyverse)




### GLOBALS ====================================================================

# Globals for experiment setup and data
DATA_PATH = "/data/"
MIN_ESTIMATE = 1 # lowest number in the range of dots
MAX_ESTIMATE = 1000 # highest number in the range of dots
TRIALS = 300 # number of trials in the experiment
BLOCKSIZE = 30 # number of blocks in split-nth fitting



### ANALYSIS FUNCTIONS =========================================================

# Data reading
read_data = function(filepath, ntrials) {
  # read data
  fp = paste0(getwd(), filepath)
  files = list.files(fp)
  dat = data.frame()
  subject = 1
  for (f in files) {
    q = read.csv2(paste(fp, f, sep = ""), sep = ",", header = T, colClasses = "character")
    q$subject = subject
    dat = rbind(dat, q)
    subject = subject + 1
  }

  # format existing columns
  to_num = function(x) {as.numeric(as.character(x))} # supporting util
  dat$run = to_num(dat$run)
  dat$index = to_num(dat$index)
  dat$num_dots = to_num(dat$num_dots)
  dat$answer1 = to_num(dat$answer1)
  dat$answer2 = to_num(dat$answer2)
  dat$points1 = to_num(dat$points1)
  dat$points2 = to_num(dat$points2)
  dat$time = to_num(dat$time)

  # add relevant new columns
  dat$avg_answer = 10 ^ (log10(pmax(1, dat$answer1)) / 2 + log10(pmax(1, dat$answer2)) / 2) # blended average of participant answers for this array
  dat$answer = dat$answer1 # choose default response to use for analyses
  dat$trial = 0
  for (s in unique(dat$subject)) {
    dat$trial[dat$subject == s] = 1:ntrials
  }
  return(dat)
}

# BIC function
get_BIC = function(num_params, num_obs, ll) {
  return((num_params * log(num_obs)) - (2 * ll))
}

# Model fitting functions
# simple power-law mapping
map_power = function(x, a, b) {
  10^(a + 10^b * log10(x))
}

# bi-linear mapping
map_bipower = function(x, a, b) {
  crit = a
  slope = 10^b
  lx = log10(x)
  ly = ((lx > crit) * (crit + (lx - crit) * slope) + (lx <= crit) * lx);
  return(10^ly)
}

## general log likelihood function (with robustness)
loglik = function(x, y, map_fx, a, b, s) {
  sum(
    pmax(-6, dnorm(log10(y) - log10(map_fx(x, a, b)), 0, s, log = T))
  )
}

# Log likelihood with variance following a linear function
# Instead of simply fitting variance `s`, we fit intercept `s.a` and linear slope `s.b`
# using the function sigma(x) = exp(log(x)*slope + intercept)
loglik_var = function(x, y, map_fx, a, b, s_a, s_b) {
  sum(
    pmax(-6, dnorm(log10(y) - log10(map_fx(x, a, b)), 0, 10^(s_a + s_b * log10(x)), log = T))
  )
}


# fit function
# `tmp` is individual subject's data, all at once or in trial blocks
brutefit = function(tmp) {
  nLL = function(a, b, s) {
    -loglik(tmp$num_dots, tmp$answer, usefx, a, b, 10^s) +
      priors[[1]](a) + priors[[2]](b) + priors[[3]](s)
  }

  iter = 0
  fits = NULL
  fit = NULL
  while (is.null(fits)) {
    try(fit <- summary(mle(nLL,
                          start = list(a = runif(1, ps["ma"], ps["sa"]),
                                       b = runif(1, ps["mb"], ps["sb"]),
                                       s = rnorm(1, ps["ms"], ps["ss"])))), TRUE)
    iter = iter + 1

    if (!is.null(fit)) {
      fits = c(tmp$subject[1], -0.5*fit@m2logL, length(tmp$num_dots), fit@coef[,"Estimate"])
    } else {
      if (iter > 50) {
        fits = c(tmp$subject[1], -9999, 0, 0, 0, 0)
      }
    }
  }
  names(fits) = c("subject", "logL", "n", "a", "b", "s")
  return(fits)
}



# fit function that assumes linear variance slope
# `tmp` is individual subject's data
# Instead of simply fitting variance `s`, we fit intercept and slope `s.a`, `s.b`
brutefit_var = function(tmp) {
  nLL = function(a, b, s_a, s_b) {
    -loglik_var(tmp$num_dots, tmp$answer, usefx, a, b, s_a, s_b) +
      priors[[1]](a) + priors[[2]](b) + priors[[3]](s_a) + priors[[4]](s_b)
  }

  iter = 0
  fits = NULL
  fit = NULL
  while (is.null(fits)) {
    try(fit <- summary(mle(nLL,
                           start = list(a = runif(1, ps["ma"], ps["sa"]),
                                        b = runif(1, ps["mb"], ps["sb"]),
                                        s_a = rnorm(1, ps["ms_a"], ps["ss_a"]),
                                        s_b = rnorm(1, ps["ms_b"], ps["ss_b"])))), TRUE)
    iter = iter + 1

    if (!is.null(fit)) {
      fits = c(tmp$subject[1], -0.5*fit@m2logL, length(tmp$num_dots), fit@coef[,"Estimate"])
    } else {
      if (iter > 500) {
        fits = c(tmp$subject[1], -9999, 0, 0, 0, 0, 0)
      }
    }
  }
  names(fits) = c("subject", "logL", "n", "a", "b", "s_a", "s_b")
  return(fits)
}


# Util function used when fitting slopes: splits trial data into sequential blocks
split_block = function(trial, n, total_trials) {
  floor((trial - 1) / (total_trials / n))
}

# Util function used when fitting slopes: splits trial data into modular blocks
split_mod = function(trial, n){
  trial %% n
}

# Util function used when calculating correlation matrix
named_slopes = function(x) {
  z = data.frame(x$b)
  rownames(z) = x$subject
  return(z)
}

# Util function used when calculating correlation matrix
cbind_fill = function(...) {
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
fit_block_slopes = function(blocksizes, dat, total_trials) {
  n = blocksizes[1]
  fits_block = list()
  fits_mod = list()
  for (k in 1:n) {
    tmp = subset(dat, split_block(dat$trial, n, total_trials) == (k - 1))
    fits_block[[k]] = data.frame(do.call(rbind, by(tmp, tmp$subject, brutefit)))

    tmp = subset(dat, split_mod(dat$trial, n) == (k - 1))
    fits_mod[[k]] = data.frame(do.call(rbind, by(tmp, tmp$subject, brutefit)))
  }

  ret = list()
  ret[['block']] = fits_block
  ret[['mod']] = fits_mod
  return(ret)
}

# Get m blocks by m blocks matrix of fitted slope correlations between each block
get_cor_matrix = function(slopes) {
  # Calculate slope correlations
  block_slopes = do.call(cbind_fill, lapply(slopes, named_slopes)) # m subjects by n blocks slope values
  cor_matrix = cor(block_slopes, block_slopes, use = "pairwise.complete.obs") # n blocks by n blocks slope correlation matrix
  rownames(cor_matrix) = c() # NB: clearing out rownames and colnames is necessary for subsequent processing
  colnames(cor_matrix) = c()

  return(cor_matrix)
}

# Convert m blocks by m blocks matrix of fitted slope correlations to data frame
# Has column for each block and correlation between those blocks
# Initially has m^2 rows for each pair of blocks, then prunes out lower half
get_cor_df = function(cor_matrix) {
  slope_cor_df = reshape::melt(cor_matrix) # data frame with columns for block x, block y, and slope correlation b/n those blocks, n blocks x n blocks rows
  names(slope_cor_df) = c("block1", "block2", "slope_corr")
  slope_cor_df = slope_cor_df[slope_cor_df$block1 <= slope_cor_df$block2,] # remove redundant lower half of matrix
  slope_cor_df$slope_corr[slope_cor_df$block1 == slope_cor_df$block2] = NA # set correlation to NA in identical blocks

  return(slope_cor_df)
}

# Format data frame of correlations by block pair to show mean, se of correlations by block distance across blocks
# NB: keep dplyr:: function designations to ensure that summarize is based on grouping
# See https://stackoverflow.com/questions/26923862/why-are-my-dplyr-group-by-summarize-not-working-properly-name-collision-with
get_distance_cors = function(cor_df, total_trials, split_blocks) {
  dist_df = cor_df %>%
    dplyr::mutate(block_dist = block2 - block1,
                  trial_dist = (total_trials / split_blocks) * block_dist) %>%
    dplyr::group_by(trial_dist) %>%
    dplyr::summarize(mean_cor = mean(slope_corr, na.rm = TRUE),
                     se_cor = sd(slope_corr, na.rm = TRUE) / sqrt(length(slope_corr)))

  return(dist_df)
}

# Shuffle function for shuffling trial number
shuffle_data = function(dat) {
  dat %>%
    group_by(subject) %>%
    mutate(trial = sample(trial, length(trial), replace = F))
}


### GRAPHING FUNCTIONS =========================================================

# Plotting functions
my_log_breaks = function(lims) {
  majors = seq(floor(log10(lims[1])), ceiling(log10(lims[2])), by = 1)
  minors = log10(unlist(lapply(majors[-1], function(x){seq(10^(x - 1), 9 * 10^(x - 1), by = 10^(x - 1))})))
  return(list(majors, minors))
}

mylogx = function(lims) {
  breaks = my_log_breaks(lims)
  scale_x_log10(limits = lims,
                breaks = 10^breaks[[1]],
                minor_breaks = breaks[[2]])
}

mylogy = function(lims) {
  breaks = my_log_breaks(lims)
  scale_y_log10(limits = lims,
                breaks = 10^breaks[[1]],
                minor_breaks = breaks[[2]])
}


# theme for plots of individual data
individ_plot_theme = theme(
  # titles
  plot.title = element_text(face = "bold", size = 28),
  axis.title.y = element_text(face = "bold", size = 28),
  axis.title.x = element_text(face = "bold", size = 28),
  legend.title = element_text(face = "bold", size = 16),
  # axis text
  axis.text.y = element_text(size = 16),
  axis.text.x = element_text(size = 14, hjust = 0, vjust = 0),
  # legend text
  legend.text = element_text(size = 24),
  # facet text
  strip.text = element_text(face = "bold", size = 12),
  # backgrounds, lines
  panel.background = element_blank(),
  strip.background = element_blank(),

  panel.grid = element_line(color = "gray"),
  axis.line = element_line(color = "black"),
  # positioning
  legend.position = "bottom",
  legend.key = element_rect(color = "transparent", fill = "transparent")
)



### ANALYSIS: Bilinear estimate calibration ====================================

data = read_data(DATA_PATH, TRIALS)

# Check fits for stimuli < 100
# data = data %>%
#   filter(num_dots <= 100) # removes approx. 7200 rows!

# Check fits only above subitizing range
# data = data %>%
  # filter(num_dots > 4) # removes approx. 500 rows

usefx = map_power
ps = c(0.2, 0.4, -0.3, 0.3, -0.7, 0.2)
names(ps) = c("ma", "sa", "mb", "sb", "ms", "ss")
priors = list()
priors[[1]] = function(x){0}
priors[[2]] = function(x){-dnorm(x, 0, 0.2, log = T)}
priors[[3]] = function(x){0}
power_fits = data.frame(do.call(rbind, by(data, data$subject, brutefit)))
print(paste("Failed power fits:", sum(power_fits$logL == -9999)))

usefx = map_bipower
ps = c(0.7, 1.5, -0.5, 0.2, -0.7, 0.2)
names(ps) = c("ma", "sa", "mb", "sb", "ms", "ss")
priors = list()
priors[[1]] = function(x){-dnorm(x, 2, 3, log = T)} #
priors[[2]] = function(x){-dnorm(x, 0, 0.5, log = T)} #
priors[[3]] = function(x){-dnorm(x, -1, 0.25, log = T)} #
bipower_fits = data.frame(do.call(rbind, by(data, data$subject, brutefit)))
print(paste("Failed bipower fits:", sum(bipower_fits$logL == -9999)))

predictions = data.frame()
for (s in unique(data$subject)) {
  stims = seq(1, max(data$num_dots), by = 1)
  powparams = power_fits[power_fits$subject == s,]
  powpred = (map_power(stims, powparams$a, powparams$b))
  biparams = bipower_fits[bipower_fits$subject == s,]
  bipred = (map_bipower(stims, biparams$a, biparams$b))
  predictions = rbind(predictions,
                       data.frame(subject = s,
                                  num_dots = stims,
                                  powpred = powpred,
                                  bipred = bipred))
}


# Compare model fits

prediction_eval = data %>%
  left_join(predictions, by = c("subject", "num_dots"))

cor_output = data.frame(
  subject = numeric(),
  bipred_cor = numeric(),
  powpred_cor = numeric()
)

for (s in unique(prediction_eval$subject)) {
  cor_output = rbind(cor_output,
                     data.frame(
                       subject = s,
                       bipred_cor = cor.test(prediction_eval$answer[prediction_eval$subject == s],
                                             prediction_eval$bipred[prediction_eval$subject == s])$estimate,
                       powpred_cor = cor.test(prediction_eval$answer[prediction_eval$subject == s],
                                              prediction_eval$powpred[prediction_eval$subject == s])$estimate))
}


# Calculate and compare r-squared for bilinear and power law functions
cor_output = cor_output %>%
  mutate(bipred_rsq = bipred_cor^2,
         powpred_rsq = powpred_cor^2)

mean(cor_output$bipred_rsq)
mean(cor_output$powpred_rsq)

binom.test(x = sum(cor_output$bipred_rsq > cor_output$powpred_rsq),
           n = length(cor_output$subject),
           p = 0.5)

# Calculate and compare BIC for bilinear and power law functions

power_bic = get_BIC(2, power_fits$n, power_fits$logL)
bipower_bic = get_BIC(2, bipower_fits$n, bipower_fits$logL)

power_fits = power_fits %>%
  rowwise() %>%
  mutate(BIC = get_BIC(2, n, logL))

bipower_fits = bipower_fits %>%
  rowwise() %>%
  mutate(BIC = get_BIC(2, n, logL))


mean(bipower_fits$BIC)
mean(power_fits$BIC)

sum(bipower_fits$BIC < power_fits$BIC)


# Summary stats about cutoff

bipower_fits = bipower_fits %>%
  mutate(cutoff_linear = 10^a,
         slope_linear = 10^b)

mean(bipower_fits$cutoff_linear)
mean(bipower_fits$a)
sd(bipower_fits$a)
10^mean(bipower_fits$a)

mean(bipower_fits$b)
sd(bipower_fits$b)
10^mean(bipower_fits$b)


### PLOTS: Bilinear estimate calibration =======================================

# FIGURE: individual data
data_subj = data %>%
  filter(subject %in% c(18)) # sample subject


ggplot(data_subj, aes(x = num_dots, y = answer)) +
  geom_point(alpha = 0.75, color = "red", size = 2) +
  geom_abline() +
  mylogx(c(1, 200)) +
  mylogy(c(1, 200)) +
  xlab("Number presented") +
  ylab("Number reported") +
  individ_plot_theme


# well, this is hacky AF...
strip_labels = c("1" = "Subject 1", "2" = "Subject 2", "3" = "Subject 3", "4" = "Subject 4", "5" = "Subject 5", "6" = "Subject 6",
                 "7" = "Subject 7", "8" = "Subject 8", "9" = "Subject 9", "10" = "Subject 10", "11" = "Subject 11", "12" = "Subject 12",
                 "13" = "Subject 13", "14" = "Subject 14", "15" = "Subject 15", "16" = "Subject 16", "17" = "Subject 17", "18" = "Subject 18",
                 "19" = "Subject 19", "20" = "Subject 20", "21" = "Subject 21", "22" = "Subject 22", "23" = "Subject 23", "24" = "Subject 24")

ggplot(data, aes(x = num_dots, y = answer)) +
  geom_point(alpha = 0.25, color = "red") +
  geom_line(data = predictions, aes(x = num_dots, y = bipred), color = rgb(0, 0.6, 0), size = 1) +
  geom_line(data = predictions, aes(x = num_dots, y = powpred), color = "blue", size = 1) +
  geom_abline(position = "identity") +
  mylogx(c(1, max(data$num_dots))) +
  # mylogy(c(1, max(data$answer))) +
  mylogy(c(1, 300)) +
  xlab("Number presented") +
  ylab("Number reported") +
  individ_plot_theme +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  facet_wrap(~subject, ncol = 6,
             labeller = labeller(subject = strip_labels))


### ANALYSIS: Stable individual differences ====================================

# Add modular and blocked split half column
split_half_data = data %>%
  mutate(modular_split_half = as.numeric(trial %% 2 == 0),
         blocked_split_half = as.numeric(trial <= TRIALS / 2))
# sanity check
table(split_half_data$subject, split_half_data$modular_split_half)
table(split_half_data$subject, split_half_data$blocked_split_half)

# Fit slopes to each modular split half
# NB: this uses more restricted priors than above (notably on cutoff)
# this is intentional so we can evaluate the slopes while holding fitted cutoffs fairly constant
usefx = map_bipower
ps = c(0.7, 1.5, -0.5, 0.2, -0.7, 0.2)
names(ps) = c("ma", "sa", "mb", "sb", "ms", "ss")
priors = list()
priors[[1]] = function(x){-dnorm(x, 1.14, 0.1, log = T)} #
priors[[2]] = function(x){-dnorm(x, -0.1, 0.25, log = T)} #
priors[[3]] = function(x){-dnorm(x, -1, 0.05, log = T)} #

d0 = split_half_data %>%
  filter(modular_split_half == 0)
d1 = split_half_data %>%
  filter(modular_split_half == 1)

bipower_fits_d0 = data.frame(do.call(rbind, by(d0, d0$subject, brutefit)))
print(paste("Failed bipower fits:", sum(bipower_fits_d0$logL == -9999)))
bipower_fits_d1 = data.frame(do.call(rbind, by(d1, d1$subject, brutefit)))
print(paste("Failed bipower fits:", sum(bipower_fits_d1$logL == -9999)))

cor_mod = cor.test(bipower_fits_d0$b, bipower_fits_d1$b)
cor_mod


### ANALYSIS: Calibration drift ================================================

# First, compare modular split half above to blocked split half
# NB: this uses a fairly restricted prior on cutoff for same reason as cited above
usefx = map_bipower
ps = c(0.7, 1.5, -0.5, 0.2, -0.7, 0.2)
names(ps) = c("ma", "sa", "mb", "sb", "ms", "ss")
priors = list()
priors[[1]] = function(x){-dnorm(x, 1.14, 0.1, log = T)} #
priors[[2]] = function(x){-dnorm(x, -0.1, 0.25, log = T)} #
priors[[3]] = function(x){-dnorm(x, -1, 0.05, log = T)} #

h1 = split_half_data %>%
  filter(blocked_split_half == 1)
h2 = split_half_data %>%
  filter(blocked_split_half == 0)

bipower_fits_h1 = data.frame(do.call(rbind, by(h1, h1$subject, brutefit)))
print(paste("Failed bipower fits:", sum(bipower_fits_h1$logL == -9999)))
bipower_fits_h2 = data.frame(do.call(rbind, by(h2, h2$subject, brutefit)))
print(paste("Failed bipower fits:", sum(bipower_fits_h2$logL == -9999)))

cor_block = cor.test(bipower_fits_h1$b, bipower_fits_h2$b)
cor_block

# Compare block and modular split-halfs
cocor.indep.groups(r1.jk = cor_block$estimate, r2.hm = cor_mod$estimate, n1 = 24, n2 = 24)
# z sig. < 0 implies r1 is sig. < r2



# Next, look at drift for split-Nths
# NB: this uses a fairly restricted prior on cutoff for same reason as cited above
usefx = map_bipower
ps = c(0.7, 1.5, -0.5, 0.2, -0.7, 0.2)
names(ps) = c("ma", "sa", "mb", "sb", "ms", "ss")
priors = list()
priors[[1]] = function(x){-dnorm(x, 1.14, 0.1, log = T)} #
priors[[2]] = function(x){-dnorm(x, -0.1, 0.25, log = T)} #
priors[[3]] = function(x){-dnorm(x, -1, 0.05, log = T)} #

# Fit slopes
slopes = fit_block_slopes(c(BLOCKSIZE), data, TRIALS) # NB: this can take ~10s
fits_block_subj = slopes[["block"]]
# Get matrix of fitted slope correlations
cor_matrix_subj = get_cor_matrix(fits_block_subj)
# Format correlation matrix as data frame to plot slope correlations by trial block in analysis section below
slope_cor_df_subj = get_cor_df(cor_matrix_subj)
# Process slope correlations by trial distance to get mean, se across participants
cor_means_df_blocks_subj = get_distance_cors(slope_cor_df_subj, TRIALS, BLOCKSIZE)
cor_means_df_blocks_subj = cor_means_df_blocks_subj %>% filter(trial_dist > 0)

# Fit line to drift
dist_df_subj = slope_cor_df_subj %>%
  dplyr::mutate(block_dist = block2 - block1,
                trial_dist = (TRIALS / BLOCKSIZE) * block_dist) %>%
  dplyr::group_by(trial_dist)

model_drift = lm(data = dist_df_subj, slope_corr ~ trial_dist)
summary(model_drift)
confint.lm(model_drift, level = 0.95)

preds = data.frame(trial_dist = seq(1, TRIALS)) %>%
  mutate(corr_pred = model_drift$coefficients[1] + trial_dist * model_drift$coefficients[2])


fits_mod_subj = slopes[["mod"]]
# Get matrix of fitted slope correlations
mod_cor_matrix_subj = get_cor_matrix(fits_mod_subj)
# Format correlation matrix as data frame to plot slope correlations by trial block in analysis section below
mod_slope_cor_df_subj = get_cor_df(mod_cor_matrix_subj)
# Process slope correlations by trial distance to get mean, se across participants
mod_cor_means_df_blocks_subj = get_distance_cors(mod_slope_cor_df_subj, TRIALS, BLOCKSIZE)
mod_cor_means_df_blocks_subj = mod_cor_means_df_blocks_subj %>% filter(trial_dist > 0)

# Shuffled fit
data_shuffle = shuffle_data(data)
# Fit slopes
fits_block_shuffle = fit_block_slopes(c(BLOCKSIZE), data_shuffle, TRIALS)[["block"]] # NB: this can take ~10s
# Get matrix of fitted slope correlations
cor_matrix_shuffle = get_cor_matrix(fits_block_shuffle)
# Format correlation matrix as data frame to plot slope correlations by trial block in analysis section below
slope_cor_df_shuffle = get_cor_df(cor_matrix_shuffle)
# Process slope correlations by trial distance to get mean, se across participants
cor_means_df_blocks_shuffle = get_distance_cors(slope_cor_df_shuffle, TRIALS, BLOCKSIZE)
cor_means_df_blocks_shuffle = cor_means_df_blocks_shuffle %>% filter(trial_dist > 0)

# Fit line to shuffled drift
dist_df_shuffle = slope_cor_df_shuffle %>%
  dplyr::mutate(block_dist = block2 - block1,
                trial_dist = (TRIALS / BLOCKSIZE) * block_dist) %>%
  dplyr::group_by(trial_dist)

model_drift_shuffle = lm(data = dist_df_shuffle, slope_corr ~ trial_dist)
summary(model_drift_shuffle)
confint.lm(model_drift_shuffle, level = 0.95)



### PLOTS: Calibration drift ===================================================

# FIGURE: slope correlation matrix
slope_cor_df_subj %>%
  ggplot(aes(x = as.factor(block1), y = as.factor(block2), fill = slope_corr)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = 0.5, limits = c(0, 1)) + # midpoint = 0.25, limits = c(0, 1)
  xlab("Trial block") + ylab("Trial block") +
  ggtitle("Trial block slope correlations") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 16),
        axis.text.y = element_text(angle = 90),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"),
        title = element_text(size = 18, face = "bold"),
        panel.grid = element_blank())


# FIGURE: slope correlations by trial distance
ggplot() +
  # block split correlations
  geom_point(data = cor_means_df_blocks_subj,
             aes(x = trial_dist,
                 y = mean_cor,
                 color = "subjects")) +
  geom_errorbar(data = cor_means_df_blocks_subj,
                aes(x = trial_dist,
                    ymin = mean_cor - se_cor,
                    ymax = mean_cor + se_cor,
                    color = "subjects"),
                width = 5) +
  # shuffle correlations
  geom_point(data = cor_means_df_blocks_shuffle,
             aes(x = trial_dist,
                 y = mean_cor,
                 color = "shuffle")) +
  geom_errorbar(data = cor_means_df_blocks_shuffle,
                aes(x = trial_dist,
                    ymin = mean_cor - se_cor,
                    ymax = mean_cor + se_cor,
                    color = "shuffle"),
                width = 5) +
  #geom_line(data = preds, aes(x = trial.dist, y = corr.pred)) +
  ylim(0.45, 0.95) + xlim(0, 300) +
  labs(x = "Distance (trials)", y = "Correlation") +
  ggtitle("Drift in slope correlation") +
  scale_color_manual(name = "",
                     labels = c("subjects" = "block split", "shuffle" = "shuffled trials"),
                     values = c("subjects" = "red", "shuffle" = "seagreen")) +
  individ_plot_theme




### ANALYSIS: Increasing CoV ===================================================

# TODO break this out into analysis and plots, as above

# CoV fits: linear function for variance, evaluate fitted (variance) slopes
# NB: this uses a fairly restricted prior on cutoff for same reason as cited above
# But a loose prior on the variables of interest for this analysis, s.a and s.b
usefx = map_bipower
# Starting vals assume variance intercept ~10^-1 and slope ~10^0
ps = c(0.7, 1.5, -0.5, 0.2, -1, 0.1, 0, 0.1)
names(ps) = c("ma", "sa", "mb", "sb", "ms_a", "ss_a", "ms_b", "ss_b")
priors = list()
priors[[1]] = function(x){-dnorm(x, 1.14, 0.1, log = T)} #
priors[[2]] = function(x){-dnorm(x, -0.1, 0.25, log = T)} #
priors[[3]] = function(x){-dnorm(x, -1, 0.25, log = T)} #
priors[[4]] = function(x){-dnorm(x, 0, 0.25, log = T)} # s.b prior for slope of 1

# First test with a single subject
single_subj = data %>% filter(subject == 1)
bipower_fits_var = data.frame(do.call(rbind, by(single_subj, single_subj$subject, brutefit_var)))
bipower_fits_var

# All subjects
bipower_fits_var = data.frame(do.call(rbind, by(data, data$subject, brutefit_var)))
print(paste("Failed bipower + linear var fits:", sum(bipower_fits_var$logL == -9999)))

# get CI
var_slopes_test = t.test(bipower_fits_var$s_b, mu = 0)
var_slopes_test$statistic
var_slopes_test$parameter
var_slopes_test$p.value
ci = var_slopes_test$conf.int
ci

mean_slope = mean(bipower_fits_var$s_b)
bipower_fits_var %>%
  ggplot(aes(x = s_b)) +
  geom_histogram(binwidth = 0.1, fill = I("grey"), col = I("black"), alpha = 0.75) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 2, color = "red") +
  geom_vline(aes(xintercept = mean_slope), linetype = "solid", size = 1, color = "black") +
  geom_vline(aes(xintercept = ci[1]), linetype = "dashed", size = 1, color = "black") +
  geom_vline(aes(xintercept = ci[2]), linetype = "dashed", size = 1, color = "black") +
  labs(x = "Log slope of variance") +
  scale_x_continuous(breaks = seq(-0.1, 0.8, 0.1)) +
  individ_plot_theme +
  theme(axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 0.2))






### APPENDIX: Posterior predictive sampling of log slope of variance ===========

# Posterior predictive sampling of variance data
xrange = 2:1000
nsamples = 10000
usefx = map_bipower

mean_cutoff = mean(bipower_fits_var$a)
mean_slope = mean(bipower_fits_var$b)
mean_var_int = mean(bipower_fits_var$s_a)
mean_var_slope = mean(bipower_fits_var$s_b)
mean_var_static = mean(bipower_fits$s)

post_data = data.frame(
  magnitude = numeric(),
  linear_variance_estimate = numeric(),
  static_variance_estimate = numeric()
)

for(mag in xrange) {
  est_linear_var = 10^rnorm(1,
                            log10(usefx(mag, mean_cutoff, mean_slope)),
                            10^(mean_var_int + mean_var_slope * log10(mag)))
  est_static_var = 10^rnorm(1,
                            log10(usefx(mag, mean_cutoff, mean_slope)),
                            10^mean_var_static)
  post_data = rbind(post_data,
                    data.frame(magnitude = mag,
                               linear_variance_estimate = est_linear_var,
                               static_variance_estimate = est_static_var))
}

post_data %>%
  ggplot(aes(x = magnitude)) +
  geom_point(aes(y = linear_variance_estimate),
             alpha = 0.75, color = "red", size = 2) +
  # Can exclude this
  geom_point(aes(y = static_variance_estimate),
             alpha = 0.75, color = "blue", size = 2) +
  geom_abline() +
  # mylogx(c(1, 200)) +
  # mylogy(c(1, 200)) +
  mylogx(c(1, 1000)) +
  mylogy(c(1, 1000)) +
  xlab("Number presented") +
  ylab("Number reported") +
  individ_plot_theme


# Calculate CoV of resampled data
magnitude_binsize = 25
post_data = post_data %>%
  mutate(
    estimate_bin = case_when(
      magnitude < 5 ~ 0,
      magnitude >= 5 & magnitude < 10 ~ 1,
      magnitude >= 10 & magnitude < 20 ~ 2,
      TRUE ~ floor(magnitude / magnitude_binsize) + 2),
    estimate_bin_value = case_when(
      estimate_bin == 0 ~ 3,
      estimate_bin == 1 ~ 7,
      estimate_bin == 2 ~ 14.5,
      TRUE ~ (estimate_bin - 2) * magnitude_binsize + 9.5)
  )

post_data_summary = post_data %>%
  group_by(estimate_bin, estimate_bin_value) %>%
  summarize(cov_linear = sd(linear_variance_estimate) / mean(linear_variance_estimate),
            cov_static = sd(static_variance_estimate) / mean(static_variance_estimate))

cov_mod = lm(data = post_data_summary,
             cov_linear ~ estimate_bin_value)
summary(cov_mod)
cov_mod$coefficients['(Intercept)']
cov_mod$coefficients['estimate_bin_value']

post_data_summary = post_data_summary %>%
  mutate(cov_linear_predictions = cov_mod$coefficients['(Intercept)'] + (estimate_bin_value * cov_mod$coefficients['estimate_bin_value']))


post_data_summary %>%
  ggplot(aes(x = estimate_bin_value)) +
  geom_point(aes(y = cov_linear), color = "red", alpha = 0.75, size = 2) +
  # geom_point(aes(y = cov_static), color = "blue", alpha = 0.75, size = 2) +
  geom_line(aes(y = cov_linear_predictions)) +
  labs(x = "Magnitude", y = "CoV estimate") +
  individ_plot_theme



### APPENDIX: Calibration accuracy, precision ==================================

# Analysis: make sure people aren't getting better at the task
# Computes average fitted bilinear slopes by trial block
df_slopes = data.frame()
for (x in 1:length(fits_block_subj)) {
  mean_slope = mean(10^fits_block_subj[[x]]$b)
  se_slopes = sd(10^fits_block_subj[[x]]$b) / sqrt(length(fits_block_subj[[x]]$b))
  df_slopes = rbind(df_slopes, data.frame(x, mean_slope, se_slopes))
}

# Plot
ggplot(data = df_slopes, aes(x = x, y = mean_slope)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  geom_ribbon(aes(ymin = mean_slope - se_slopes, ymax = mean_slope + se_slopes),
              color = "gray", alpha = 0.3) +
  labs(x = "trial block", y = "mean slope estimate") +
  individ_plot_theme


# Analysis: make sure people aren't getting better at the task (v2)
# Computes average subject mean squared error by trial block
data = data %>%
  mutate(
    block = split_block(trial, BLOCKSIZE[1], TRIALS) + 1,
    trial_resid = (answer - num_dots)^2)

rmse_subj_summary = data %>%
  group_by(subject, block) %>%
  summarize(mean_subj_resid = mean(trial_resid),
            trials = n(),
            se_resid = sd(trial_resid) / sqrt(trials))

rmse_block_summary = rmse_subj_summary %>%
  group_by(block) %>%
  summarize(mean_block_resid = mean(mean_subj_resid),
            subjects = n(),
            se_block_resid = sd(mean_subj_resid) / sqrt(subjects))

# Plot
rmse_block_summary %>%
  ggplot(aes(x = block, y = mean_block_resid)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  geom_ribbon(aes(ymin = mean_block_resid - se_block_resid,
                  ymax = mean_block_resid + se_block_resid),
              color = "gray", alpha = 0.3) +
  labs(x = "trial block", y = "residuals^2 against number") +
  individ_plot_theme


# Analysis: make sure people's estimates aren't getting more "fine tuned"
# Computes average subject mean squared error *against fitted prediction* by block
data = data %>%
  # need to group by variables that give us one row at a time in `data` for list indexing in fits_block_subj
  group_by(subject, trial) %>%
  mutate(block_cutoff = fits_block_subj[[block]]$a[fits_block_subj[[block]]$subject == subject],
         block_slope = fits_block_subj[[block]]$b[fits_block_subj[[block]]$subject == subject],
         block_pred = map_bipower(num_dots, block_cutoff, block_slope),
         block_residsq = (block_pred - answer)^2)
  # filter(block_residsq < 100000) # is this necessary?

data_block_preds_summary = data %>%
  group_by(block) %>%
  summarize(mean_residsq = mean(block_residsq),
            num_subj = length(unique(subject)),
            se_residsq = sd(block_residsq) / sqrt(num_subj))

# Plot
data_block_preds_summary %>%
  ggplot(aes(x = block, y = mean_residsq)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  geom_ribbon(aes(ymin = mean_residsq - se_residsq, ymax = mean_residsq + se_residsq),
              color = "gray", alpha = 0.3) +
  labs(x = "trial block", y = "residuals^2 against predictions") +
  individ_plot_theme


