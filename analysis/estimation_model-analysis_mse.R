
#'
#' Modeling for estimation drift manuscript
#'

setwd(paste0(here::here(), "/vullab/estimation_drift/"))
rm(list=ls())

source('analysis/estimation_model-fxns_basic.R')


### GLOBALS ====================================================================

SAVE_MSE_DATA = TRUE
MSE_DATA_FILE = 'estimation_model-mse.RData'
# This will load mse.df into the environment with saved data
# load(paste('analysis/', MSE_DATA_FILE, sep = ""))

### ANALYSIS FUNCTIONS =========================================================

# General function to calculate MSE of participant and model estimates
calculate.mse = function(data) {
  data = data %>%
    group_by(subject) %>%
    mutate(subj.sq.error = (answer - num_dots) ^ 2,
           model.sq.error = (model.answer - num_dots) ^ 2)
  return(data)
}


### GRAPHING FUNCTIONS =========================================================

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
  plot.title = element_text(face = "bold", size = 32),
  axis.title.y = element_text(face = "bold", size = 32),
  axis.title.x = element_text(face = "bold", size = 32),
  legend.title = element_text(face = "bold", size = 16),
  # axis text
  axis.text.y = element_text(size = 16),
  axis.text.x = element_text(size = 14, angle = 90, hjust = 0, vjust = 0),
  # legend text
  legend.text = element_text(size = 24),
  # facet text
  strip.text = element_text(face = "bold", size = 28),
  strip.background = element_blank(),
  # backgrounds, lines
  panel.background = element_blank(),
  #strip.background = element_blank(),

  panel.grid = element_line(color = "gray"),
  axis.line = element_line(color = "black"),
  # positioning
  legend.position = "bottom"
)

individ_plot_theme_mse = theme(
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


plot.sample.subjects = function(data, sample.subjects = c(), label) {
  sample.subject.plot = data %>%
    filter(subject %in% sample.subjects,
           trial > N_SAMPLES) %>%
    ggplot(aes(x = num_dots)) +
    geom_point(aes(y = answer, color = "subject"), alpha = 0.5) +
    geom_point(aes(y = model.answer, color = "model"), alpha = 0.5) +
    geom_abline() +
    #ggtitle(paste("Model estimates, calibration = ", N_SAMPLES, " samples, p = ", P_BUMPER)) +
    ggtitle(paste("Sample estimates")) +
    labs(x = "Number presented", y = "Number reported") +
    scale_color_manual(name = "Estimates",
                       values = c("subject" = "blue", "model" = "red")) +
    mylogx(c(MIN_ESTIMATE, MAX_ESTIMATE)) +
    mylogy(c(MIN_ESTIMATE, MAX_ESTIMATE)) +
    individ_plot_theme + # TODO pass this in to the parent function
    facet_wrap(~subject, ncol = 3, labeller = label)

  return(sample.subject.plot)
}


plot.mse = function(data) {
  data %>% ggplot(aes(x = N)) +
    geom_line(aes(y = mse.subj, color = "subjects")) +
    geom_line(aes(y = mse.mod, color = "model")) +
    geom_point(aes(y = mse.mod, color = "model"), size = 2) +
    scale_color_manual(name = element_blank(),
                       values = c("subjects" = "blue", "model" = "red")) +
    scale_y_continuous(limits = c(0, 4000)) +
    ggtitle("Model performance with increasing samples") +
    labs(x = "Samples", y = "Estimate error (MSE)") +
    individ_plot_theme_mse + # TODO consider making this a param in the function
    theme(legend.position = c(0.2, 0.2),
          legend.background = element_rect(color = "gray80", size = 0.5, linetype = "solid"))

}



### ANALYSIS: baseline model ===================================================
p.bumper_baseline = 1.0
# Run baseline model
data.base = run.model.baseline(p.bumper = p.bumper_baseline)

# Graph three sample subjects and model estimates for a close look
sample.subjects = c(3, 7, 9)
label = labeller(subject = c("3" = "Subject 3", "7" = "Subject 7", "9" = "Subject 9"))
sample.subj.plot.base = plot.sample.subjects(data.base, sample.subjects, label)
sample.subj.plot.base


# Calculate MSE for baseline model over different numbers of samples
# NB: this can take a while!! (approx. 60 seconds per run, but this increases as N_SAMPLES increases)
mse.df = data.frame('N' = numeric(), 'p' = numeric(), 'mse.subj' = numeric(), 'mse.mod' = numeric())
for (n in seq(1, 25)) {
  print(paste("Running model with N_SAMPLES = ", n, ", P_BUMPER = ", p.bumper_baseline))
  data = run.model.baseline(n.samples = n,
                            p.bumper = p.bumper_baseline)
  data = calculate.mse(data)
  mse.df = rbind(mse.df, data.frame('N' = n, 'p' = p.bumper_baseline,
                                    'mse.subj' = mean(data$subj.sq.error),
                                    'mse.mod' = mean(data$model.sq.error[data$trial > n])))
}

# Graph MSE results
mse.plot = plot.mse(mse.df)
mse.plot
# Display sample subject and MSE graphs together
sample.subj.plot.base / mse.plot

# optional: save output data
if (SAVE_MSE_DATA) {
  save(mse.df, file = paste('analysis/', MSE_DATA_FILE, sep = ""))
}


