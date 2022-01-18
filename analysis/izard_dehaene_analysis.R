
#'
#' Script for analyzing data from Izard & Dehaene (2008)
#'

setwd("/Users/erikbrockbank/web/vullab/estimation_drift/")
# rm(list = ls())

library(cocor)
library(Matrix)
library(stats4)
library(tidyverse)

# NB: this takes approx. 1min.
# source('analysis/exp1-analysis_drift.R')



# Functions ====
# NB: these are copied from `exp1-analysis_drift.R` so we can modify as needed

# BIC function
get_BIC = function(num_params, num_obs, ll) {
  return((num_params * log(num_obs)) - (2 * ll))
}

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

# general log likelihood function (with robustness)
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
      if (iter > 500) {
        print(paste("Failed to fit data for subject: ", tmp$subject[1], " n = ", length(tmp$num_dots)))
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




# Initialize data ====

e1_data = read_tsv("analysis/id2008/experiment1.txt")
e2_data = read_tsv("analysis/id2008/data_experiment2.txt")

# Rename columns to match experiment data so functions can be re-used
e1_data = e1_data %>%
  filter(!is.nan(resp)) %>% # Some NaNs that prevent  model fitting
  rename(subject = Subject,
         num_dots = stim,
         answer = resp)

e2_data = e2_data %>%
  filter(!is.nan(resp)) %>% # Some NaNs that prevent model fitting
  rename(subject = Subject,
         num_dots = stim,
         answer = resp)

#' E1 data: conditions
#' condition == 0:
#' dotsize varies
#' area varies but seems it's meant to be constant (84 values, 576-6400)
#' Hull sort of constant (4 values, 320-340)
#' density varies (99 values, 0.000010-0.000865)

#' condition == 1:
#' dotsize constant
#' area varies (100 values, 100-10000)
#' Hull constant (340)
#' density varies (100 values, 0.000009-0.000865)

#' condition == 5:
#' dotsize constant
#' area varies (100 values, 100-10000)
#' Hull varies (13 values, 100-340)
#' density constant (76 values, but all hovering very close to .0001)



e1_data %>%
  group_by(subject) %>%
  summarize(estimates = n())

e1_data %>%
  filter(answer == 30) %>%
  group_by(subject) %>%
  summarize(mean(num_dots))

#' Based on the above:
#' 'ML' = subject 5
#' 'PQ' = subject 6
#' 'AL' = subject 1
#' 'DC' = subject 2
#' 'BF' = subject 3
#' subject 4 = ????

# Add subject identifiers based on the above
subject_lookup = c("1" = "AL", "2" = "DC", "3" = "BF", "5" = "ML", "6" = "PQ")

e1_data = e1_data %>%
  filter(subject != 4) %>%
  mutate(subject_inits = subject_lookup[as.character(subject)])
e1_data$subject_inits = factor(e1_data$subject_inits,
                               levels = c("ML", "PQ", "AL", "DC", "BF"))


e2_data %>%
  group_by(subject) %>%
  summarize(estimates = n())




# PLOTS: Data overview ====

e1_data %>%
  ggplot(aes(x = num_dots, y = answer)) +
  geom_point(alpha = 0.5) +
  geom_abline(linetype = "dashed", size = 1) +
  mylogx(c(1, max(e1_data$num_dots))) +
  mylogy(c(1, max(e1_data$answer))) +
  facet_wrap(~subject_inits,
             ncol = 5) +
  labs(x = "Number presented", y = "Number reported") +
  individ_plot_theme +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


e2_data %>%
  filter(inducer == 30) %>% # Alternatives: 25, 39
  ggplot(aes(x = num_dots, y = answer)) +
  geom_point(alpha = 0.25, color = "red") +
  geom_abline(linetype = "dashed", size = 1) +
  mylogx(c(1, max(e2_data$num_dots))) +
  mylogy(c(1, max(e2_data$answer))) +
  facet_wrap(~subject,
             ncol = 6) +
  labs(x = "Number presented", y = "Number reported") +
  individ_plot_theme +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))



# ANALYSIS: Fit bilinear models ====

# Use same starting values and priors as original analysis
usefx = map_power
ps = c(0.2, 0.4, -0.3, 0.3, -0.7, 0.2)
names(ps) = c("ma", "sa", "mb", "sb", "ms", "ss")
priors = list()

# Empty priors
priors[[1]] = function(x){0}
priors[[2]] = function(x){0}
priors[[3]] = function(x){0}
# Informative priors
# priors[[1]] = function(x){0}
# priors[[2]] = function(x){-dnorm(x, 0, 0.2, log = T)}
# priors[[3]] = function(x){0}

power_fits_id1 = data.frame(do.call(rbind, by(e1_data, e1_data$subject, brutefit)))
print(paste("Failed power fits:", sum(power_fits_id1$logL == -9999)))

usefx = map_bipower
ps = c(0.7, 1.5, -0.5, 0.2, -0.7, 0.2)
names(ps) = c("ma", "sa", "mb", "sb", "ms", "ss")
priors = list()
# Empty priors (same result)
priors[[1]] = function(x){0}
priors[[2]] = function(x){0}
priors[[3]] = function(x){0}
# Informative priors
# priors[[1]] = function(x){-dnorm(x, 2, 3, log = T)} #
# priors[[2]] = function(x){-dnorm(x, 0, 0.5, log = T)} #
# priors[[3]] = function(x){-dnorm(x, -1, 0.25, log = T)} #
bipower_fits_id1 = data.frame(do.call(rbind, by(e1_data, e1_data$subject, brutefit)))
print(paste("Failed bipower fits:", sum(bipower_fits_id1$logL == -9999)))



# ANALYSIS: Compare model fits ====

# Convert params to linear space for easier interpretation
power_fits_id1 = power_fits_id1 %>%
  mutate(intercept = 10^a,
         slope = 10^b)
bipower_fits_id1 = bipower_fits_id1 %>%
  mutate(cutoff = 10^a,
         slope = 10^b)


# Get predictions from each model
predictions = data.frame()
max_est = max(e1_data$num_dots)
for (s in unique(power_fits_id1$subject)) { # NB: this should be e1_data$subject with all fits
  stims = seq(1, max_est, by = 1)
  powparams = power_fits_id1[power_fits_id1$subject == s,]
  powpred = map_power(stims, powparams$a, powparams$b)
  biparams = bipower_fits_id1[bipower_fits_id1$subject == s,]
  bipred = map_bipower(stims, biparams$a, biparams$b)
  predictions = rbind(predictions,
                      data.frame(subject = s,
                                 num_dots = stims,
                                 powpred = powpred,
                                 bipred = bipred))
}


# Correlate predictions with empirical estimates

prediction_eval = e1_data %>%
  inner_join(predictions, by = c("subject", "num_dots"))

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


# Evaluate BIC for each model

bipower_fits_id1 = bipower_fits_id1 %>%
  rowwise() %>%
  mutate(BIC = get_BIC(2, n, logL))

power_fits_id1 = power_fits_id1 %>%
  rowwise() %>%
  mutate(BIC = get_BIC(2, n, logL))


mean(bipower_fits_id1$BIC)
mean(power_fits_id1$BIC)

sum(bipower_fits_id1$BIC < power_fits_id1$BIC)


# PLOTS: Predictions for each model alongside true estimates ====

predictions = predictions %>%
  mutate(subject_inits = subject_lookup[as.character(subject)])

e1_data %>%
  ggplot(aes(x = num_dots, y = answer)) +
  geom_point(alpha = 0.25) +
  geom_abline(linetype = "dashed", size = 1) +
  geom_line(data = predictions, aes(x = num_dots, y = powpred),
            size = 1, color = "blue") +
  geom_line(data = predictions, aes(x = num_dots, y = bipred),
            size = 1, color = rgb(0, 0.6, 0)) +
  mylogx(c(1, max(e1_data$num_dots))) +
  mylogy(c(1, max(e1_data$answer))) +
  facet_wrap(~subject_inits,
             ncol = 5) +
  labs(x = "Number presented", y = "Number reported") +
  individ_plot_theme +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))




# ANALYSIS: Increasing CoV model for E1 and E2 ====


# Fit E1 data
usefx = map_bipower
# Starting vals assume variance intercept ~10^-1 and slope ~10^0
ps = c(0.7, 1.5, -0.5, 0.2, -1, 0.1, 0, 0.1)
names(ps) = c("ma", "sa", "mb", "sb", "ms_a", "ss_a", "ms_b", "ss_b")
priors = list()
priors[[1]] = function(x){-dnorm(x, 1, 0.1, log = T)} #
priors[[2]] = function(x){-dnorm(x, -0.1, 0.25, log = T)} #
priors[[3]] = function(x){-dnorm(x, -1, 0.25, log = T)} #
priors[[4]] = function(x){-dnorm(x, 0, 0.25, log = T)} # s.b prior for slope of 1


e1_bipower_fits_var = data.frame(do.call(rbind, by(e1_data, e1_data$subject, brutefit_var)))
print(paste("Failed bipower + linear var fits:", sum(e1_bipower_fits_var$logL == -9999)))
e1_bipower_fits_var


# get CI
e1_var_slopes_test = t.test(e1_bipower_fits_var$s_b, mu = 0)
e1_var_slopes_test$statistic
e1_var_slopes_test$parameter
e1_var_slopes_test$p.value
e1_ci = e1_var_slopes_test$conf.int
e1_ci


# Fit E2 data
e2_calibrated = e2_data %>%
  filter(inducer == 30)

usefx = map_bipower
# Starting vals assume variance intercept ~10^-1 and slope ~10^0
ps = c(0.7, 1.5, -0.5, 0.2, -1, 0.1, 0, 0.1)
names(ps) = c("ma", "sa", "mb", "sb", "ms_a", "ss_a", "ms_b", "ss_b")
priors = list()
priors[[1]] = function(x){-dnorm(x, 1, 0.1, log = T)} #
priors[[2]] = function(x){-dnorm(x, -0.1, 0.25, log = T)} #
priors[[3]] = function(x){-dnorm(x, -1, 0.5, log = T)} #
priors[[4]] = function(x){-dnorm(x, 0, 0.25, log = T)} # s.b prior for slope of 1


e2_bipower_fits_var = data.frame(do.call(rbind, by(e2_calibrated, e2_calibrated$subject, brutefit_var)))
print(paste("Failed bipower + linear var fits:", sum(e2_bipower_fits_var$logL == -9999)))
e2_bipower_fits_var


# get CI
e2_var_slopes_test = t.test(e2_bipower_fits_var$s_b, mu = 0)
e2_var_slopes_test$statistic
e2_var_slopes_test$parameter
e2_var_slopes_test$p.value
e2_ci = e2_var_slopes_test$conf.int
e2_ci



# PLOTS: Distributions of linear variance slopes ====

mean_slope_e1 = mean(e1_bipower_fits_var$s_b)
e1_bipower_fits_var %>%
  ggplot(aes(x = s_b)) +
  geom_histogram(binwidth = 0.1, fill = I("grey"), col = I("black"), alpha = 0.75) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 2, color = "red") +
  geom_vline(aes(xintercept = mean_slope_e1), linetype = "solid", size = 1, color = "black") +
  geom_vline(aes(xintercept = e1_ci[1]), linetype = "dashed", size = 1, color = "black") +
  geom_vline(aes(xintercept = e1_ci[2]), linetype = "dashed", size = 1, color = "black") +
  labs(x = "Log slope of variance") +
  # scale_x_continuous(breaks = seq(-0.1, 0.8, 0.1)) +
  individ_plot_theme +
  theme(axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 0.2))


mean_slope_e2 = mean(e2_bipower_fits_var$s_b)
e2_bipower_fits_var %>%
  ggplot(aes(x = s_b)) +
  geom_histogram(binwidth = 0.1, fill = I("grey"), col = I("black"), alpha = 0.75) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 2, color = "red") +
  geom_vline(aes(xintercept = mean_slope_e2), linetype = "solid", size = 1, color = "black") +
  geom_vline(aes(xintercept = e2_ci[1]), linetype = "dashed", size = 1, color = "black") +
  geom_vline(aes(xintercept = e2_ci[2]), linetype = "dashed", size = 1, color = "black") +
  labs(x = "Log slope of variance") +
  # scale_x_continuous(breaks = seq(-0.1, 0.8, 0.1)) +
  individ_plot_theme +
  theme(axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 0.2))





