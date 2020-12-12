
### README ###
#' This file contains supporting functions for the samples model.
#'
#' The model attempts to simulate a process of estimating quantity whereby the participant
#' has a noisy estimate of the magnitude of the number of dots shown and chooses an exact numerical estimate
#' based on the probability that this magnitude is greater or less than the magnitudes on previously seen trials.
#' In addition, it has a parameter for K 'bumpers' of rough magnitude -> number mappings which are assumed to align
#' with people's existing knowledge of certain magnitudes and prevent the model from wildly over or underestimating when
#' it sees a number which is higher or lower than any previous trials.
#'

library(Matrix)
library(patchwork)
library(stats4)
library(tidyverse)



###############
### GLOBALS ###
###############

# Globals for experiment setup and data
DATA = "/data/"
MIN_ESTIMATE = 1 # lowest number in the range of dots
MAX_ESTIMATE = 1000 # highest number in the range of dots
TRIALS = 300 # number of trials in the experiment (easier to set this as a global than compute it with e.g. max(trial))


############################
### MODEL DEFAULT PARAMS ###
############################

PERCEIVED_DOTS_NOISE_SD = 0.1 # log noise of perceived magnitudes

PRIOR_EXP = -1 # exponential slope parameter for prior

# params for sample from posterior
POSTERIOR_SAMPLE_EXP = 4 # exponent `l` to use when sampling from the posterior with probability p^l

# params for sample of previous trials
TRIALS_EXP = -1 # explonential slope parameter for sampling previous trials

# params for sample of real world "bumpers"
BUMPER_EXP = -4 # exponential slope parameter for sampling "bumper" values

# hyperparams: number of samples and probability that samples are derived from real world bumpers
N_SAMPLES = 20 # total number of samples to take from previous trials and "bumper" distribution
P_BUMPER = 0.5 # probability p that a given sample comes from previous trials or real world "bumper" distribution

N_MEMORIES = 30 # number of memories to use for individual variability model



#################
### FUNCTIONS ###
#################

# Data reading
read.data = function(data_dir, n.trials) {
  # read data
  fp = paste0(getwd(), data_dir)
  files = list.files(fp)
  data = data.frame()
  subject = 1
  for (f in files) {
    q = read.csv2(paste(fp, f, sep = ""), sep = ",", header = T, colClasses = "character")
    q$subject = subject
    data = rbind(data, q)
    subject = subject + 1
  }

  # format existing columns
  to.num = function(x) {as.numeric(as.character(x))} # supporting util
  data$run = to.num(data$run)
  data$index = to.num(data$index)
  data$num_dots = to.num(data$num_dots)
  data$answer1 = to.num(data$answer1)
  data$answer2 = to.num(data$answer2)
  data$points1 = to.num(data$points1)
  data$points2 = to.num(data$points2)
  data$time = to.num(data$time)

  # add relevant new columns
  data$avg_answer = 10 ^ (log10(pmax(1, data$answer1)) / 2 + log10(pmax(1, data$answer2)) / 2) # blended average of participant answers for this array
  data$answer = data$avg_answer # choose default response to use for analyses
  data$trial = 0
  for (s in unique(data$subject)) {
    data$trial[data$subject == s] = 1:n.trials
  }
  data$model.answer = 0
  data$bumper.trial = FALSE

  return(data)
}

# General function to calculate Coefficient of Variation based on standard deviation of sampled magnitude values
calculate.cov = function(weber_noise) {
  x = 1 # arbitrary
  samples = 10000 # arbitrary
  magnitude = 10 ^ rnorm(samples, log10(x), sd = weber_noise)
  return(sd(magnitude) / x) # approx. CoV
}


# Retrieve previous trials and estimate "bumpers"
# When mapping.vals is defined, uses an idiosyncratic set of individual memories,
# otherwise samples them for each trial (baseline model)
get.samples = function(trial.i, subj, single.subj.data, n.trials, n.samples, p.bumper, min.est, max.est,
                       perceived.mag.sd, trials.exp, bumper.exp, mapping.vals = NULL) {
  # get set of previous trials and real world bumpers to sample
  sample.bool = rbinom(n.samples, 1, p.bumper) # n.samples length vector with 1 for bumper samples, 0 for previous trials
  num.bumpers = sum(sample.bool)
  num.trials = n.samples - num.bumpers
  if (num.trials >= trial.i) { # more previous trial samples needed than available
    # use up to number of available trials for previous trial sampling, then fill in rest with bumpers
    num.trials = trial.i - 1
    num.bumpers = n.samples - num.trials
  }
  # fetch previous trials for magnitude comparison
  prev.trials = sample.previous.trials(trial.i, single.subj.data, num.trials, trials.exp)
  # add in "bumpers" as if they were previous trials (note we're providing a fixed set from which to sample here)
  prev.trials = sample.bumper.mappings(subj, prev.trials, num.bumpers, min.est, max.est, bumper.exp, perceived.mag.sd, mapping.vals)

  return(prev.trials)
}


# Sample previous trials based on exponential decay function with slope trials.exp
sample.previous.trials = function(trial.i, single.subject.data, n.trials, trials.exp) {
  if (n.trials > 0) {
    # set of number mappings from which to sample bumpers
    trial_set = data.frame("sub_index" = seq(1, trial.i - 1)) # number of trials *back* from current trial to sample
    trial_set$sub_index.scaled = trial_set$sub_index ^ trials.exp
    trial_set$ptrial_sub = trial_set$sub_index.scaled / sum(trial_set$sub_index.scaled)
    prev.trial.sub.indices = sort(with(trial_set, sample(sub_index, n.trials, prob = ptrial_sub, replace = F)))

    # get only trials m indices prior to trial.i for m in prev.trial.sub.indices
    prev.trials = single.subject.data %>%
      filter(trial %in% (trial.i - prev.trial.sub.indices))
  } else {
    prev.trials = single.subject.data[0,] # return formatting of single subject data but with no previous trials
  }

  return(prev.trials)
}

# Sample real world "bumper" mappings; return magnitude, value mapping and associated probability
get.bumper.set = function(bumper.min, bumper.max, bumper.exp, samples) {
  world_vals = data.frame("num" = seq(bumper.min, bumper.max))
  world_vals$num.scaled = world_vals$num ^ bumper.exp
  world_vals$pnum = world_vals$num.scaled / sum(world_vals$num.scaled)

  bumper_set = sample_n(tbl = world_vals, size = samples, replace = F, weight = world_vals$pnum)
  return(bumper_set)
}

# Sample "bumper" mappings by adding known magnitude -> number mappings to previous trials
# that are being used for likelihood calculation
sample.bumper.mappings = function(subj, prev.trials, k_bumpers, bumper.min, bumper.max, bumper.exp, mag.sd, bumper.set.fixed = NULL) {
  if (k_bumpers > 0) {
    bumper_samples = data.frame(num_dots = numeric(), estimate = numeric()) # sampled value of bumpers for mapping
    if (is.null(bumper.set.fixed)) {
      # fetch set of sampled bumper number mappings: get sample of true numbers, create noisy mapping
      bumper_set = get.bumper.set(bumper.min, bumper.max, bumper.exp, k_bumpers)
      for (bumper in bumper_set$num) {
        bumper_samples = rbind(bumper_samples, data.frame(num_dots = bumper,
                                                          estimate = round(10 ^ rnorm(1, log10(bumper), mag.sd), digits = 0)))
      }
    } else {
      if (length(unique(bumper.set.fixed$num_dots)) <= k_bumpers) {
        bumper_samples = bumper.set.fixed
      } else {
        bumper_samples = sample_n(bumper.set.fixed, k_bumpers, replace = F, weight = bumper.set.fixed$prob)
      }
    }

    prev.trials = prev.trials %>%
      add_row(subject = subj, # this gets passed in so we don't rely on previous trials if there are none
              num.dots = bumper_samples$num_dots,
              model.answer = bumper_samples$estimate,
              bumper.trial = TRUE)
  }

  return(prev.trials)
}


# calculate log likelihood using magnitudes and number estimates from previous trials
calculate.loglik = function(prev.trials, n.samples, min.est, max.est) {
  if (nrow(prev.trials) < n.samples) {
    n.samples = nrow(prev.trials) # If we have fewer previous trials than samples, watch out!
  }
  # matrix of number estimates initialized for use in calculation of likelihoods in loop below
  # NB: this matrix A doesn't need to be created every time we call this function but it doesn't slow it down much
  # so I'm keeping it here for compactness
  A = matrix(rep(min.est:max.est, each = n.samples), nrow = n.samples)

  # discretize possible answer estimates, set log likelihood at each value to be either p.mag.greater or 1 - p.mag.greater
  # use matrix with a row for each of the previous trials and a column for each candidate estimate value.
  # each row, col element in X will be the log likelihood that the magnitude on the current trial maps to the column value, given the magnitude
  #   and corresponding estimate in that row's trial
  X = matrix(0, nrow = n.samples, ncol = max.est) # Initialize matrix for storing log likelihoods
  # matrix storing proability in each col that magnitude on this trial is < magnitude in trial indicated by that row
  P = matrix(rep(1 - prev.trials$p.mag.greater[1:(n.samples)], max.est), nrow = n.samples)
  # matrix storing numerical estimate generated for each previous trial
  B = matrix(rep(floor(prev.trials$model.answer[1:(n.samples)]), max.est), nrow = n.samples)
  X[A <= B] = log10(P[A <= B]) # NB: if we switch this to < and the next line to >=, it makes a big diff, not totally clear why
  X[A > B] = log10(1 - P[A > B])
  # across each previous trial that we've calculated a log likelihood vector for above, calculate aggregate likelihood for estimates
  loglik = data.frame('number.est' = seq(min.est:max.est), 'sum.loglik' = colSums(X))

  # process log likelihood calculated above
  loglik = loglik %>%
    mutate(
      sum.loglik = sum.loglik - max(sum.loglik), #add a constant C to each value in loglike to make it more reasonable
      loglik.probability.raw = 10 ^ sum.loglik # reconvert out of log space
    )

  return(loglik)
}

# calculate posterior by taking in log likelihood data, assuming power law prior
# NB: this can be modified to make different assumptions about the prior
calculate.posterior = function(loglik, prior.exp) {
  posterior = loglik %>%
    mutate(prior = number.est ^ prior.exp, # assume power law probability of seeing a particular number
           prior.norm = prior / sum(prior), # normalize
           posterior.raw = loglik.probability.raw * prior.norm, # multiply loglik by the prior
           posterior.norm = posterior.raw / sum(posterior.raw)) # normalize

  return(posterior)
}


# calculate MAP for a given posterior distribution over possible number estimates
posterior.map = function(posterior) {
  map = posterior$number.est[posterior$posterior.norm == max(posterior$posterior.norm)]
  return(map)
}

# generate a sample from the posterior distribution over possible number estimates
posterior.sample = function(posterior) {
  post.sample = with(posterior, sample(number.est, 1, prob = posterior.norm))
  return(post.sample)
}

# calculate mean of the posterior distribution over possible number estimates
posterior.mean = function(posterior) {
  post.mean = sum(posterior$number.est * posterior$posterior.norm)
  return(post.mean)
}

# calculate median of the posterior distribution over possible number estimates
# NB: median of the posterior is used for similar problems in Griffiths, Tenenbaum 2006
posterior.median = function(posterior) {
  estimate.med = posterior %>%
    mutate(post_sum = cumsum(posterior.norm)) %>%
    filter(post_sum >= 0.5) %>%
    summarize(est = min(number.est)) %>%
    pull(est)

  return(estimate.med)
}

# this is a function for sampling according to p^l for some l
# posterior.sample is a special case of this function where l=1
# higher l will sample closer to MAP, lower l will be closer to mean
posterior.sample.pl = function(posterior, sample.exp) {
  post.sample = with(posterior, sample(number.est, 1, prob = (posterior.raw ^ sample.exp)))
  return(post.sample)
}


# Baseline model assumes no individual differences over separate experimental runs
run.model.baseline = function(data.path = DATA,
                              n.trials = TRIALS,
                              n.samples = N_SAMPLES,
                              p.bumper = P_BUMPER,
                              min.est = MIN_ESTIMATE,
                              max.est = MAX_ESTIMATE,
                              perceived.mag.sd = PERCEIVED_DOTS_NOISE_SD,
                              trials.exp = TRIALS_EXP,
                              bumper.exp = BUMPER_EXP,
                              prior.exp = PRIOR_EXP,
                              post.exp = POSTERIOR_SAMPLE_EXP,
                              est.fxn = posterior.sample.pl) {

  data = read.data(data.path, n.trials)
  subjects = unique(data$subject)

  for (subj in subjects) {
    print(paste("Running with participant: ", subj))
    print(paste("Bumper probability: ", p.bumper))
    print(paste("Samples: ", n.samples))

    # reformat data to be consistent with util functions above
    single.subj.data = data %>%
      mutate(participant.answer = answer,
             num.dots = num_dots) %>%
      select(subject, trial, num.dots, participant.answer, model.answer, bumper.trial) %>%
      filter(subject == subj)

    for (trial.i in seq(1, max(data$trial[data$subject == subj]))) {
      prev.trials = get.samples(trial.i, subj, single.subj.data, n.trials, n.samples, p.bumper,
                                min.est, max.est, perceived.mag.sd, trials.exp, bumper.exp)
      # compare magnitude of trial.i to previous trial magnitudes
      # -> assumes both have a mean around the true number of dots, sd = PERCEIVED_DOTS_NOISE_SD
      curr.trial.mag.mean = single.subj.data$num.dots[single.subj.data$trial == trial.i]
      prev.trials = prev.trials %>%
        mutate(p.mag.greater = 1 - pnorm(0, mean = log10(curr.trial.mag.mean) - log10(num.dots),
                                         sd = sqrt(2 * perceived.mag.sd)))

      # calculate initial log likelihood for this trial based on p.mag.greater from previous trials
      loglik = calculate.loglik(prev.trials, n.samples, min.est, max.est)
      # calculate posterior using log likelihood calculated above
      posterior = calculate.posterior(loglik, prior.exp)

      # Generate an estimate from the posterior
      # -> (functions above for mean, median, probability matched sample, or MAP)
      estimate = est.fxn(posterior, post.exp) # all estimate functions above take in posterior

      # add estimate to single subject results
      single.subj.data$model.answer[single.subj.data$trial == trial.i] = estimate
    }
    # add estimates for this subject to original data
    data$model.answer[data$subject == subj] = single.subj.data$model.answer

  }

  return(data)
}


# Individual variability model samples memories for each individual run of the model
run.model.individ.memories = function(data.path = DATA,
                                      n.trials = TRIALS,
                                      n.samples = N_SAMPLES,
                                      p.bumper = P_BUMPER,
                                      min.est = MIN_ESTIMATE,
                                      max.est = MAX_ESTIMATE,
                                      perceived.mag.sd = PERCEIVED_DOTS_NOISE_SD,
                                      trials.exp = TRIALS_EXP,
                                      bumper.exp = BUMPER_EXP,
                                      prior.exp = PRIOR_EXP,
                                      post.exp = POSTERIOR_SAMPLE_EXP,
                                      n.memories = N_MEMORIES,
                                      est.fxn = posterior.sample.pl) {

  data = read.data(data.path, n.trials)
  subjects = unique(data$subject)

  for (subj in subjects) {
    print(paste("Running with participant: ", subj))
    print(paste("Bumper probability: ", p.bumper))
    print(paste("Samples: ", n.samples))
    print(paste("Memories: ", n.memories))

    # reformat data to be consistent with util functions above
    single.subj.data = data %>%
      mutate(participant.answer = answer,
             num.dots = num_dots) %>%
      select(subject, trial, num.dots, participant.answer, model.answer, bumper.trial) %>%
      filter(subject == subj)

    # before running trials, determine participant's mapping "memories" to sample from
    mapping.set = get.bumper.set(MIN_ESTIMATE, MAX_ESTIMATE, BUMPER_EXP, n.memories)
    mapping.vals = data.frame(num_dots = numeric(), estimate = numeric(), prob = numeric())
    for (bumper in mapping.set$num) {
      mapping.vals = rbind(mapping.vals, data.frame(num_dots = bumper,
                                                    estimate = round(10 ^ rnorm(1, log10(bumper), PERCEIVED_DOTS_NOISE_SD), digits = 0),
                                                    prob = mapping.set$pnum[mapping.set$num == bumper]))
    }
    mapping.vals$prob = mapping.vals$prob / sum(mapping.vals$prob) # re-normalize probabilities

    for (trial.i in seq(1, max(data$trial[data$subject == subj]))) {
      # get set of previous trials and real world bumpers to sample
      prev.trials = get.samples(trial.i, subj, single.subj.data, n.trials, n.samples, p.bumper,
                                min.est, max.est, perceived.mag.sd, trials.exp, bumper.exp, mapping.vals)
      # compare magnitude of trial.i to previous trial magnitudes
      # -> assumes both have a mean around the true number of dots, sd = PERCEIVED_DOTS_NOISE_SD
      curr.trial.mag.mean = single.subj.data$num.dots[single.subj.data$trial == trial.i]
      prev.trials = prev.trials %>%
        mutate(p.mag.greater = 1 - pnorm(0, mean = log10(curr.trial.mag.mean) - log10(num.dots),
                                         sd = sqrt(2 * perceived.mag.sd)))

      # calculate initial log likelihood for this trial based on p.mag.greater from previous trials
      loglik = calculate.loglik(prev.trials, n.samples, min.est, max.est)
      # calculate posterior using log likelihood calculated above
      posterior = calculate.posterior(loglik, prior.exp)

      # Generate an estimate from the posterior
      # -> (functions above for mean, median, probability matched sample, or MAP)
      estimate = est.fxn(posterior, post.exp) # all estimate functions above take in posterior

      # add estimate to single subject results
      single.subj.data$model.answer[single.subj.data$trial == trial.i] = estimate
    }
    # add estimates for this subject to original data
    data$model.answer[data$subject == subj] = single.subj.data$model.answer

  }

  return(data)
}


model.wrapper = function(nruns, model.fxn, ...) {
  model.data = data.frame(
    model.run = numeric(),
    subject = character(),
    trial = numeric(),
    num_dots = numeric(),
    answer = numeric(),
    model.answer = numeric()
  )
  for (x in seq(1:nruns)) {
    print(paste("######## Model run: ", x, " ########"))
    data = model.fxn(...)
    data = data %>%
      mutate(model.run = x,
             subj.repeated = subject + ((x - 1) * max(data$subject)))
    model.data = rbind(model.data,
                       data.frame(
                         model.run = rep(x, nrow(data)),
                         subject = data$subj.repeated,
                         trial = data$trial,
                         num_dots = data$num_dots,
                         answer = data$model.answer
                         ))
  }

  return(model.data)
}
