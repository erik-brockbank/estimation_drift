
### README ###
#' This file does an analysis of the *bilinear fit* across subject estimates
#' and model estimates for different versions of the model.
#'
#' It takes the output from the model in `samples_model-fxns_basic` and
#' then uses more advanced functions in `samples_model-fxns_drift` to
#' fit lines to the model and human estimates, for easy analysis of bilinear parameters
#' and comparison between model and participant bilinear fit.
#'


setwd(paste0(here::here(), "/vullab/estimation_drift/"))
rm(list=ls())

library(Rmisc) # needed for call to `multiplot`
# Fetch relevant model functions from samples_model
source('analysis/estimation_model-fxns_basic.R')
# Fetch relevant functions for fitting lines to model data
source('analysis/estimation_model-fxns_drift.R')


### GLOBALS ====================================================================

OUTPUT_FILE = 'estimation_model-bilinear.RData'
# load(paste('analysis/', OUTPUT_FILE, sep = ""))

### ANALYSIS FUNCTIONS =========================================================

# NB: this over-rides equivalent functions from fxns_drift.R
# Compute best fitting params for estimate data
brutefit = function(tmp) {
  nLL = function(a, b, s) {
    -loglik(tmp$num_dots, tmp$answer, usefx, a, b, 10 ^ s) + PRIORS[[1]](a) + PRIORS[[2]](b) + PRIORS[[3]](s)
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

# simple power-law mapping
map.power = function(x, a, b) {
  10^(a + 10^b * log10(x))
}


### GRAPHING FUNCTIONS =========================================================

# Graphing functions
my.log.breaks = function(lims){
  majors = seq(floor(log10(lims[1])), ceiling(log10(lims[2])), by = 1)
  minors = log10(unlist(lapply(
    majors[-1],
    function(x){seq(10 ^ (x - 1), 9 * 10 ^ (x - 1), by = 10 ^ (x - 1))})))
  return(list(majors, minors))
}

mylogx = function(lims){
  breaks = my.log.breaks(lims)
  scale_x_log10(limits = lims,
                breaks = 10 ^ breaks[[1]],
                minor_breaks = breaks[[2]])
}

mylogy = function(lims){
  breaks = my.log.breaks(lims)
  scale_y_log10(limits = lims,
                breaks = 10 ^ breaks[[1]],
                minor_breaks = breaks[[2]])
}

individ_plot_theme = theme(
  # titles
  plot.title = element_text(face = "bold", size = 28),
  axis.title.y = element_text(face = "bold", size = 24),
  axis.title.x = element_text(face = "bold", size = 24),
  legend.title = element_text(face = "bold", size = 16),
  # axis text
  axis.text.y = element_text(size = 14),
  axis.text.x = element_text(size = 14, angle = 90, hjust = 0, vjust = 0),
  # legend text
  legend.text = element_text(size = 14),
  # facet text
  strip.text = element_text(face = "bold", size = 20),
  # backgrounds, lines
  panel.background = element_blank(),
  strip.background = element_blank(),

  panel.grid = element_line(color = "gray"),
  axis.line = element_line(color = "black"),
  # positioning
  legend.position = "bottom"
)

# NB: this is the same as "drift_theme" in `estimate_model-analysis_drift.R`
large_text_theme = theme(
  # titles
  plot.title = element_text(face = "bold", size = 32),
  axis.title.y = element_text(face = "bold", size = 32),
  axis.title.x = element_text(face = "bold", size = 32),
  legend.title = element_text(face = "bold", size = 16),
  # axis text
  axis.text.y = element_text(size = 20),
  axis.text.x = element_text(size = 20, hjust = 1), #, angle = 60, hjust = 0, vjust = 0.1
  # legend text
  legend.text = element_text(size = 24),
  # facet text
  strip.text = element_text(face = "bold", size = 28),
  # backgrounds, lines
  panel.background = element_blank(),
  #strip.background = element_blank(),

  panel.grid = element_line(color = "gray"),
  axis.line = element_line(color = "black"),
  # positioning
  legend.position = "bottom"
)

plot.fits.all.subjects = function(subj.data, model.data, predictions) {
  ggplot() +
    geom_point(data = subj.data,
               aes(x = num_dots, y = answer, color = "subject"), alpha = 0.25) +
    geom_point(data = model.data,
               aes(x = num_dots, y = answer, color = "model"), alpha = 0.25) +
    geom_line(data = predictions,
              aes(x = num_dots, y = bipred.subj, color = "subject")) +
    geom_line(data = predictions,
              aes(x = num_dots, y = bipred.mod, color = "model")) +
    geom_abline() +
    mylogx(c(1, 300)) +
    mylogy(c(1, 300)) +
    ggtitle("Sample estimates") +
    xlab("Number presented") +
    ylab("Number reported") +
    scale_color_manual(name = "Estimates",
                       values = c("subject" = "blue", "model" = "red")) +
    individ_plot_theme + # TODO take in theme as a param?
    facet_wrap(~subject, ncol = 6)
}

plot.fits.sample.subjects = function(subj.data, model.data, predictions, subj.labels) {
  ggplot() +
    geom_point(data = subj.data,
               aes(x = num_dots, y = answer, color = "subject"), alpha = 0.25, size = 2) +
    geom_point(data = model.data,
               aes(x = num_dots, y = answer, color = "model"), alpha = 0.25, size = 2) +
    geom_line(data = predictions,
              aes(x = num_dots, y = bipred.subj, color = "subject"), size = 2) +
    geom_line(data = predictions,
              aes(x = num_dots, y = bipred.mod, color = "model"), size = 2) +
    geom_abline() +
    mylogx(c(1, 300)) +
    mylogy(c(1, 1000)) +
    ggtitle("Sample estimates") +
    xlab("Number presented") +
    ylab("Number reported") +
    scale_color_manual(name = "Estimates",
                       values = c("subject" = "blue", "model" = "red")) +
    individ_plot_theme +
    facet_wrap(~subject, ncol = 3,
               labeller = labeller(subject = subj.labels))
}

scattterplot.fits.comparison = function(subj.fits, model.fits.summary) {
  ggplot() +
    geom_point(data = subj.fits, aes(x = cutoff.trans, y = slope.trans, color = "subjects"),
               size = 3) +
    geom_errorbar(data = subj.fits,
                  aes(x = cutoff.trans, color = "subjects",
                                                ymin = slope.trans - slope.error.trans,
                                                ymax = slope.trans + slope.error.trans),
                  width = 0) +
    geom_point(data = model.fits.summary,
               aes(x = mean.cutoff, y = mean.slope, color = "model"),
               size = 3) +
    geom_errorbar(data = model.fits.summary,
                  aes(x = mean.cutoff, color = "model",
                                                       ymin = mean.slope - se.slope,
                                                       ymax = mean.slope + se.slope),
                  width = 0.5) +
    geom_errorbarh(data = model.fits.summary,
                   aes(y = mean.slope, color = "model",
                                                        xmin = mean.cutoff - se.cutoff,
                                                        xmax = mean.cutoff + se.cutoff),
                   height = 0.05) +
    geom_hline(yintercept = 1.0, linetype = "dashed") +
    scale_color_manual(name = "Data source",
                       values = c("subjects" = "red", "model" = "blue"),
                       labels = c("subjects" = "subjects", "model" = "model")) +
    #ylim(0, 1.5) +
    #xlim(0, 100) +
    #mylogx(c(1, 100)) +
    ggtitle("Model bilinear fit is consistent with subjects") +
    labs(x = "Fitted bilinear cutoff", y = "Fitted slope above cutoff") +
    individ_plot_theme
}

histogram.slopes.comparison = function(subj.fits, model.fits.summary) {
  subj.fits %>%
    ggplot(aes(x = slope.trans)) +
    geom_histogram(binwidth = 0.05, fill = I("grey"), col = I("black")) +
    geom_vline(data = model.fits.summary,
               aes(xintercept = mean.slope), color = "red", linetype = "dashed", size = 2) +
    labs(x = "Fitted slopes (model in red)") +
    ggtitle("Model underestimation is similar to subjects") +
    individ_plot_theme +
    theme(axis.title = element_text(face = "bold", size = 16),
          plot.title = element_text(face = "bold", size = 18))

}

histogram.cutoff.comparison = function(subj.fits, model.fits.summary) {
  subj.fits %>%
    ggplot(aes(x = cutoff.trans)) +
    geom_histogram(binwidth = 3, fill = I("grey"), col = I("black")) +
    geom_vline(data = model.fits.summary,
               aes(xintercept = mean.cutoff), color = "red", linetype = "dashed", size = 2) +
    labs(x = "Fitted cutoffs (model in red)") +
    ggtitle("Model cutoff is similar to subjects") +
    individ_plot_theme +
    theme(axis.title = element_text(face = "bold", size = 16),
          plot.title = element_text(face = "bold", size = 18))
}


plot.fits.marginalized = function(subj.data, model.data, predictions) {
  ggplot() +
    geom_point(data = model.data,
               aes(x = num_dots, y = answer, color = "model"), size = 2, alpha = 0.05) +
    geom_point(data = subj.data,
               aes(x = num_dots, y = answer, color = "subjects"), size = 2, alpha = 0.05) +
    geom_line(data = predictions,
              aes(x = num_dots, y = prediction.subj, color = "subjects"), size = 2, alpha = 0.85) +
    geom_line(data = predictions,
              aes(x = num_dots, y = prediction.mod, color = "model"), size = 2, alpha = 0.85) +
    geom_ribbon(data = predictions,
                mapping = aes(x = num_dots, ymin = prediction.ll.mod, ymax = prediction.ul.mod), inherit.aes = FALSE, alpha = 0.2) +
    geom_ribbon(data = predictions,
                mapping = aes(x = num_dots, ymin = prediction.ll.subj, ymax = prediction.ul.subj), inherit.aes = FALSE, alpha = 0.2) +
    geom_abline(position = "identity") +
    mylogx(c(1, MAX_ESTIMATE)) +
    mylogy(c(1, MAX_ESTIMATE)) +
    xlab("Number presented") +
    ylab("Number reported") +
    ggtitle("Comparing model and human bilinear fits") +
    scale_color_manual(name = element_blank(),
                       values = c("subjects" = "red", "model" = "blue"),
                       labels = c("subjects" = "subjects", "model" = "model")) +
    large_text_theme +
    theme(legend.position = c(0.8, 0.26),
          legend.background = element_rect(color = "gray80", size = 0.5, linetype = "solid"))

}



### RUN MODELS =================================================================

data = run.model.baseline()

# Select relevant data for subjects, base model, individual variability model
subj.data = data %>%
  select(subject, trial, num_dots, answer)

model.data = data %>%
  mutate(answer = model.answer) %>% # align column names to match participant data
  select(subject, trial, num_dots, answer)



### ANALYSIS ===================================================================

# Fit slopes to each subject's data (in aggregate, not by blocks)
usefx = map.bipower
PARAMS = c(0.7, 1.5, -0.5, 0.2, -0.7, 0.2)
names(PARAMS) = c("ma", "sa", "mb", "sb", "ms", "ss")
PRIORS = list()
PRIORS[[1]] = function(x){-dnorm(x, 1.5, 0.1, log = T)} #
PRIORS[[2]] = function(x){-dnorm(x, -0.2, 0.1, log = T)} #
PRIORS[[3]] = function(x){-dnorm(x, -1, 0.1, log = T)} #

# Fit subject data
bipower.fits.subj = data.frame(do.call(rbind, by(subj.data, subj.data$subject, brutefit)))
print(paste("Failed bipower fits:", sum(bipower.fits.subj$logL == -9999)))
bipower.fits.subj

# Fit model data
bipower.fits.mod = data.frame(do.call(rbind, by(model.data, model.data$subject, brutefit)))
print(paste("Failed bipower fits:", sum(bipower.fits.mod$logL == -9999)))
bipower.fits.mod

# Fit power law functions to model estimates for comparison with bilinear fits
usefx = map.power
PARAMS = c(0.2, 0.4, -0.3, 0.3, -0.7, 0.2)
names(PARAMS) = c("ma", "sa", "mb", "sb", "ms", "ss")
PRIORS = list()
PRIORS[[1]] = function(x){0}
PRIORS[[2]] = function(x){-dnorm(x, 0, 0.2, log = T)}
priors[[3]] = function(x){0}
power.fits.mod = data.frame(do.call(rbind, by(model.data, model.data$subject, brutefit)))
print(paste("Failed power fits:", sum(power.fits.mod$logL == -9999)))
power.fits.mod


# Add columns transforming fits to non-log space
bipower.fits.subj = bipower.fits.subj %>%
  mutate(cutoff.trans = 10^a,
         slope.trans = 10^b,
         cutoff.error.trans = se.a,
         slope.error.trans = se.b)

bipower.fits.mod = bipower.fits.mod %>%
  mutate(cutoff.trans = 10^a,
         slope.trans = 10^b,
         cutoff.error.trans = se.a,
         slope.error.trans = se.b)

# Summary data structure used in subject/model comparison scatterplot
bipower.fits.mod.summary = bipower.fits.mod %>%
  summarize(mean.cutoff = mean(cutoff.trans),
            mean.slope = mean(slope.trans),
            se.cutoff = sd(cutoff.trans) / sqrt(length(cutoff.trans)),
            se.slope = sd(slope.trans) / sqrt(length(slope.trans)))


# Get predictions based on model fits
predictions = data.frame()
for (s in unique(subj.data$subject)){
  stims = seq(1, 300, by = 1)
  # subject params, static
  biparams.subj = bipower.fits.subj[bipower.fits.subj$subject == s,]
  bipred.subj = (map.bipower(stims, biparams.subj$a, biparams.subj$b))
  # model params
  biparams.mod = bipower.fits.mod[bipower.fits.mod$subject == s,]
  bipred.mod = (map.bipower(stims, biparams.mod$a, biparams.mod$b))
  # model params: power
  powparams.mod = power.fits.mod[power.fits.mod$subject == s,]
  powpred.mod = (map.power(stims, powparams.mod$a, powparams.mod$b))
  predictions = rbind(predictions,
                      data.frame(subject = s,
                                 num_dots = stims,
                                 bipred.subj = bipred.subj,
                                 bipred.mod = bipred.mod,
                                 powpred.mod = powpred.mod))
}

# Compare model bilinear fits to power law fits
prediction_eval = model.data %>%
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
                                             prediction_eval$bipred.mod[prediction_eval$subject == s])$estimate,
                       powpred_cor = cor.test(prediction_eval$answer[prediction_eval$subject == s],
                                              prediction_eval$powpred.mod[prediction_eval$subject == s])$estimate))
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





# Fit slopes to estimates across all subjects, marginalizing over subject
model.data.agg = model.data
model.data.agg$subject = 1 # set one subject for all data points
subj.data.agg = subj.data
subj.data.agg$subject = 1 # set one subject for all data points


# Fit slopes for model estimates in aggregate
fits.agg.mod = data.frame(do.call(rbind, by(model.data.agg, model.data.agg$subject, brutefit)))
print(paste("Failed bipower fits:", sum(fits.agg.mod$logL == -9999)))
fits.agg.mod

# Fit slopes for subject estimates in aggregate
fits.agg.subj = data.frame(do.call(rbind, by(subj.data.agg, subj.data.agg$subject, brutefit)))
print(paste("Failed bipower fits:", sum(fits.agg.subj$logL == -9999)))
fits.agg.subj


true_vals = 1:MAX_ESTIMATE
predictions.agg = data.frame('num_dots' = true_vals,
                         'prediction.mod' = map.bipower.ci(true_vals, fits.agg.mod$a, 10^fits.agg.mod$b),
                         'prediction.subj' = map.bipower.ci(true_vals, fits.agg.subj$a, 10^fits.agg.subj$b),
                         'prediction.ul.mod' = map.bipower.ci(true_vals, fits.agg.mod$a, 10^(fits.agg.mod$b + fits.agg.mod$se.b)),
                         'prediction.ll.mod' = map.bipower.ci(true_vals, fits.agg.mod$a, 10^(fits.agg.mod$b - fits.agg.mod$se.b)),
                         'prediction.ul.subj' = map.bipower.ci(true_vals, fits.agg.subj$a, 10^(fits.agg.subj$b + fits.agg.subj$se.b)),
                         'prediction.ll.subj' = map.bipower.ci(true_vals, fits.agg.subj$a, 10^(fits.agg.subj$b - fits.agg.subj$se.b)))


### PLOT RESULTS ===============================================================

### 1. Sample subjects plots: model and subject underestimation side by side ###
# NB: these plots not used in manuscript

# Look at all participants matched to model
fits.comparison.overall = plot.fits.all.subjects(subj.data, model.data, predictions)
fits.comparison.overall

# Choose sample subjects from graph above
sample.subjects = c(19, 20, 21)
subj.labels = c("19" = "Subject 19", "20" = "Subject 20", "21" = "Subject 21")
sdat = subset(subj.data, subj.data$subject %in% sample.subjects)
mod.sdat = subset(model.data, model.data$subject %in% sample.subjects)
spredictions = subset(predictions, predictions$subject %in% sample.subjects)

fits.comparison.sample = plot.fits.sample.subjects(sdat, mod.sdat, spredictions, subj.labels)
fits.comparison.sample



### 2. Plot bilinear fit scatter plots and histograms ###

# Bilinear fit comparison scatterplot (subjects and static model)
fits.comparison.scatter = scattterplot.fits.comparison(bipower.fits.subj, bipower.fits.mod.summary)
fits.comparison.scatter


# Bilinear fit comparison histograms

slopes.plot = histogram.slopes.comparison(bipower.fits.subj, bipower.fits.mod.summary)
slopes.plot

cutoff.plot = histogram.cutoff.comparison(bipower.fits.subj, bipower.fits.mod.summary)
cutoff.plot

multiplot(cutoff.plot, slopes.plot, cols = 2)



### 3. Plot overall model underestimation ###

# Plot model and subject data aggregate fits, marginalized over individual subjects
marginal.comparison.plot = plot.fits.marginalized(subj.data.agg, model.data.agg, predictions.agg)
marginal.comparison.plot


save(subj.data, model.data,
     predictions,
     bipower.fits.subj, bipower.fits.mod.summary,
     subj.data.agg, model.data.agg,
     predictions.agg, file = paste("analysis/", OUTPUT_FILE, sep = ""))

