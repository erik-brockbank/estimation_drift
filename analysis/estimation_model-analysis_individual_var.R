
### README ###
#' This file does an analysis of the *individual variability* across subject estimates
#' and model estimates for different versions of the model.
#'
#' It takes the output from the model in `samples_model-fxns_basic` and
#' then uses more advanced functions in `samples_model-fxns_drift` to
#' fit lines to the model and human estimates, for easy analysis variability
#' in bilinear fit parameters across subject and model participants
#'


setwd("/Users/erikbrockbank/web/vullab/numberline/erikb-2018/")
rm(list=ls())

library(viridis)

# Fetch relevant model functions from samples_model
source('samples_model-fxns_basic.R')
# Fetch relevant functions for fitting lines to model data
source('samples_model-fxns_drift.R')

DATA_FILE = "samples_model_individual_var.RData" # global here to keep reads/saves consistent
# load(DATA_FILE)

##########################
### ANALYSIS FUNCTIONS ###
##########################

get.split.half.cor = function(data) {
  data = data %>%
    group_by(subject) %>%
    # modular split half below
    # mutate(split.half = as.numeric(trial %% 2 == 0)) %>% # this is just a stand-in, we shuffle these below
    # mutate(split.half = replace(split.half, values = sample(split.half, length(split.half), replace = F)))
    # blocked split half below
    mutate(split.half = trial <= max(data$trial) / 2)

  # Fit subject split half estimates
  d0 = data %>%
    filter(split.half == 0)
  d1 = data %>%
    filter(split.half == 1)

  split.half.fits.0 = data.frame(do.call(rbind, by(d0, d0$subject, brutefit)))
  print(paste("Failed split half fits:", sum(split.half.fits.0$logL == -9999)))

  split.half.fits.1 = data.frame(do.call(rbind, by(d1, d1$subject, brutefit)))
  print(paste("Failed split half fits:", sum(split.half.fits.1$logL == -9999)))

  split.half.corrs = cor.test(split.half.fits.0$b, split.half.fits.1$b)

  return(split.half.corrs)
}


# Graphing functions
my.log.breaks = function(lims){
  majors = seq(floor(log10(lims[1])), ceiling(log10(lims[2])), by = 1)
  minors = log10(unlist(lapply(majors[-1], function(x){seq(10 ^ (x - 1), 9 * 10 ^ (x - 1), by = 10 ^ (x - 1))})))
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

plot.estimates.all = function(data, graph.title) {
  data %>%
    ggplot(aes(x = num_dots, y = answer)) +
    geom_point(alpha = 0.25, color = "red", size = 0.75) +
    geom_abline() +
    mylogx(c(MIN_ESTIMATE, MAX_ESTIMATE)) +
    mylogy(c(MIN_ESTIMATE, MAX_ESTIMATE + 200)) +
    xlab("Number presented") +
    ylab("Number reported") +
    ggtitle(graph.title) +
    theme(axis.title = element_text(size = 16, face = "bold"),
          title = element_text(size = 18, face = "bold")) +
    facet_wrap(~subject, ncol = 6)
}

plot.model.comparison.histogram = function(model.slopes.low, model.slopes.high, subj.slopes) {
  ggplot() +
    # subject slopes
    geom_histogram(data = model.slopes.low, aes(x = slope.trans, color = "model.low", fill = "model.low"),
                   alpha = 0.5,
                   binwidth = 0.05) +
    # high var slopes
    geom_histogram(data = model.slopes.high, aes(x = slope.trans, color = "model.high", fill = "model.high"),
                   alpha = 0.5,
                   binwidth = 0.05) +
    geom_vline(xintercept = mean(subj.slopes$slope.trans), color = "black", linetype = "dashed", size = 2) +
    scale_color_manual(name = element_blank(),
                       labels = c("model.low" = "low variability model",
                                  "model.high" = "high variability model"),
                       values = c("model.low" = "red", "model.high" = "blue")) +
    scale_fill_manual(name = element_blank(),
                      labels = c("model.low" = "low variability model",
                                 "model.high" = "high variability model"),
                      values = c("model.low" = "red", "model.high" = "blue")) +
    # scale_color_viridis(name = element_blank(),
    #                    labels = c("model.low" = "low variability model",
    #                               "model.high" = "high variability model"),
    #                    discrete = T,
    #                    begin = 0, end = 0.5) +
    # scale_fill_viridis(name = element_blank(),
    #                   labels = c("model.low" = "low variability model",
    #                              "model.high" = "high variability model"),
    #                   discrete = T,
    #                   begin = 0, end = 0.5) +
    ggtitle("Comparing model slope variability") +
    labs(x = "Fitted bilinear slope") +
    individ_plot_theme
}

plot.model.human.comparison.histogram = function(model.slopes, subj.slopes, plot.title) {
  ggplot() +
    # subject slopes
    geom_histogram(data = subj.slopes, aes(x = slope.trans, color = "subjects", fill = "subjects"),
                   alpha = 0.5,
                   binwidth = 0.05) +
    # high var slopes
    geom_histogram(data = model.slopes, aes(x = slope.trans, color = "model", fill = "model"),
                   alpha = 0.5,
                   binwidth = 0.05) +
    scale_color_manual(name = element_blank(),
                       labels = c("subjects" = "subjects",
                                  "model" = "model"),
                       values = c("subjects" = "red", "model" = "blue")) +
    scale_fill_manual(name = element_blank(),
                      labels = c("subjects" = "subjects",
                                 "model" = "model"),
                      values = c("subjects" = "red", "model" = "blue")) +
    ggtitle(plot.title) +
    labs(x = "Fitted bilinear slope") +
    individ_plot_theme
}

plot.model.comparison.scatter = function(model.slopes.low, model.slopes.high, subj.slopes) {
  ggplot() +
    geom_point(data = model.slopes.low, aes(x = cutoff.trans, y = slope.trans, color = "model.low"),
               size = 3) +
    geom_errorbar(data = model.slopes.low, aes(x = cutoff.trans,
                                                   ymin = slope.trans - slope.error.trans,
                                                   ymax = slope.trans + slope.error.trans,
                                                   color = "model.low"),
                  width = 0) +
    geom_point(data = model.slopes.high, aes(x = cutoff.trans, y = slope.trans, color = "model.high"),
               size = 3) +
    geom_errorbar(data = model.slopes.high, aes(x = cutoff.trans,
                                                    ymin = slope.trans - slope.error.trans,
                                                    ymax = slope.trans + slope.error.trans,
                                                    color = "model.high"),
                  width = 0) +
    geom_hline(yintercept = 1.0, linetype = "dashed") +
    scale_color_manual(name = "Data source",
                       values = c("model.low" = "red", "model.high" = "blue"),
                       labels = c("model.low" = "low variability model", "model.high" = "high variability model")) +
    # ylim(0, 1.5) +
    # xlim(0, 30) +
    ggtitle("Comparing model fits") +
    labs(x = "Fitted bilinear cutoff", y = "Fitted slope above cutoff") +
    individ_plot_theme
}

plot.model.human.comparison.scatter = function(model.slopes, subj.slopes, plot.title) {
  ggplot() +
    geom_point(data = subj.slopes, aes(x = cutoff.trans, y = slope.trans, color = "subjects"),
               size = 3) +
    geom_errorbar(data = subj.slopes, aes(x = cutoff.trans, color = "subjects",
                                                ymin = slope.trans - slope.error.trans,
                                                ymax = slope.trans + slope.error.trans),
                  width = 0) +
    geom_point(data = model.slopes, aes(x = cutoff.trans, y = slope.trans, color = "model"),
               size = 3) +
    geom_errorbar(data = model.slopes, aes(x = cutoff.trans,
                                                   ymin = slope.trans - slope.error.trans,
                                                   ymax = slope.trans + slope.error.trans,
                                                   color = "model"),
                  width = 0) +
    geom_hline(yintercept = 1.0, linetype = "dashed") +
    scale_color_manual(name = "Data source",
                       values = c("subjects" = "red", "model" = "blue"),
                       labels = c("subjects" = "subjects", "model" = "model")) +
    ylim(0, 1.5) +
    xlim(0, 30) +
    ggtitle(plot.title) +
    labs(x = "Fitted bilinear cutoff", y = "Fitted slope above cutoff") +
    individ_plot_theme
}

plot.model.human.split.half = function(split.half.corrs) {
  split.half.corrs %>%
    # ggplot(aes(x = fct_reorder(source, mean.corr, .desc = TRUE), y = mean.corr, fill = source)) +
    ggplot(aes(x = source, y = mean.corr, fill = source)) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_errorbar(aes(ymin = lower.corr, ymax = upper.corr), width = 0.25) +
    labs(x = "", y = "Fitted slope correlation") +
    ggtitle("Split half model comparison") +
    individ_plot_theme +
    scale_fill_viridis(discrete = T,
                       labels = c("subjects" = "subjects",
                                  "model.high.var" = "high variability model",
                                  "model.low.var" = "low variability model",
                                  "model.baseline" = "baseline model")) +
    theme(axis.text.x = element_blank(),
          legend.title = element_blank())
}



###############
### GLOBALS ###
###############

PARAMS = c(0.7, 1.5, -0.5, 0.2, -0.7, 0.2)
names(PARAMS) = c("ma", "sa", "mb", "sb", "ms", "ss")

PRIORS = list()
PRIORS[[1]] = function(x){-dnorm(x, 1.5, 0.1, log = T)} #
PRIORS[[2]] = function(x){-dnorm(x, -0.2, 0.1, log = T)} #
PRIORS[[3]] = function(x){-dnorm(x, -1, 0.1, log = T)} #

MODEL_RUNS = 10 # NB: takes about 90s per model run



###########################
### ANALYSIS: SLOPE FIT ###
###########################

# Run models
subj.data = run.model.baseline()
subj.data = subj.data %>%
  select(subject, trial, num_dots, answer)

model.data.high.var = model.wrapper(nruns = 1,
                                    model.fxn = run.model.individ.memories,
                                    n.memories = 10)
model.data.low.var = model.wrapper(nruns = 1,
                                   model.fxn = run.model.individ.memories,
                                   n.memories = 1000)

# Fit cutoff, slope params
bipower.fits.subj = data.frame(do.call(rbind, by(subj.data, subj.data$subject, brutefit)))
bipower.fits.mod.high = data.frame(do.call(rbind, by(model.data.high.var, model.data.high.var$subject, brutefit)))
bipower.fits.mod.low = data.frame(do.call(rbind, by(model.data.low.var, model.data.low.var$subject, brutefit)))

# Add columns transforming fits to non-log space
bipower.fits.subj = bipower.fits.subj %>%
  mutate(cutoff.trans = 10^a,
         slope.trans = 10^b)

bipower.fits.mod.high = bipower.fits.mod.high %>%
  mutate(cutoff.trans = 10^a,
         slope.trans = 10^b)

bipower.fits.mod.low = bipower.fits.mod.low %>%
  mutate(cutoff.trans = 10^a,
         slope.trans = 10^b)


# Histogram of model fits compared against each other
plot.model.comparison.histogram(bipower.fits.mod.low, bipower.fits.mod.high, bipower.fits.subj)



######################################
### ANALYSIS: SPLIT HALF SLOPE FIT ###
######################################

# Overall structure for storing correlation values
model.split.half.corrs = data.frame(source = character(),
                              model.run = numeric(),
                              corr = numeric(),
                              lower.ci = numeric(),
                              upper.ci = numeric())

for (x in seq(1:MODEL_RUNS)) {
  print(paste("######## CYCLE: ", x, " ########")) # approx. 5 mins / cycle
  # Run models
  model.data.high.var = model.wrapper(nruns = 1,
                                      model.fxn = run.model.individ.memories,
                                      n.memories = 10)
  model.data.low.var = model.wrapper(nruns = 1,
                                     model.fxn = run.model.individ.memories,
                                     n.memories = 1000)
  model.data.baseline = model.wrapper(nruns = 1,
                                      model.fxn = run.model.baseline)

  # Get split half correlations
  model.data.high.var.cors = get.split.half.cor(model.data.high.var)
  model.data.low.var.cors = get.split.half.cor(model.data.low.var)
  model.baseline.cors = get.split.half.cor(model.data.baseline)

  # Update aggregate data frame of split half corrs
  model.split.half.corrs = rbind(model.split.half.corrs, data.frame(source = "model.high.var",
                                                                    model.run = x,
                                                                    corr = model.data.high.var.cors$estimate,
                                                                    lower.ci = model.data.high.var.cors$conf.int[1],
                                                                    upper.ci = model.data.high.var.cors$conf.int[2]))
  model.split.half.corrs = rbind(model.split.half.corrs, data.frame(source = "model.low.var",
                                                                    model.run = x,
                                                                    corr = model.data.low.var.cors$estimate,
                                                                    lower.ci = model.data.low.var.cors$conf.int[1],
                                                                    upper.ci = model.data.low.var.cors$conf.int[2]))
  model.split.half.corrs = rbind(model.split.half.corrs, data.frame(source = "model.baseline",
                                                                    model.run = x,
                                                                    corr = model.baseline.cors$estimate,
                                                                    lower.ci = model.baseline.cors$conf.int[1],
                                                                    upper.ci = model.baseline.cors$conf.int[2]))
}

# Summarize many model runs above
detach(package:Rmisc)
detach(package:plyr)
model.split.half.summary = model.split.half.corrs %>%
  group_by(source) %>%
  summarize(mean.corr = mean(corr),
            n.obs = n(),
            corr.sd = sd(corr),
            corr.se = sd(corr) / sqrt(n.obs),
            lower.corr = mean.corr - corr.se,
            upper.corr = mean.corr + corr.se)

# Generate single run of subject data
subj.data = run.model.baseline()
subj.data = subj.data %>%
  select(subject, trial, num_dots, answer)
subj.split.half.corrs = get.split.half.cor(subj.data)

# Add subject data to summary
model.split.half.summary = as.data.frame(model.split.half.summary)
split.half.summary = rbind(data.frame(source = "subjects",
                                                      mean.corr = subj.split.half.corrs$estimate,
                                                      n.obs = NA,
                                                      corr.sd = NA,
                                                      corr.se = NA,
                                                      lower.corr = subj.split.half.corrs$conf.int[1],
                                                      upper.corr = subj.split.half.corrs$conf.int[2]),
                           model.split.half.summary)



# Bar plot of correlations over multiple runs above
plot.model.human.split.half(split.half.summary)



# Save data
save(model.split.half.corrs, subj.split.half.corrs, split.half.summary, bipower.fits.mod.low, bipower.fits.mod.high, bipower.fits.subj, file = DATA_FILE)

