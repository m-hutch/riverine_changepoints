#--------------------------#
#### MCMC Moving Window ####
#--------------------------#
# A numerical solution to Bayesian Online Changepoint Detection
# by sampling a log Pearson III distribution across a moving window
# this is a simplified moving window update where the posterior mean/SD
# of the parameters become the prior parameters for the next step.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(cmdstanr
               , posterior
               , dplyr
               , logr)
source("../utils.R")

#### Set up logging ####
log.fp <- file.path(getwd(), paste0("mcmc_moving_window.log"))
lf <- log_open(log.fp)

#### Load the Stan model ####
file <- file.path("logpearson3.stan")
mod <- cmdstan_model(file)

#### Load Data ####
#data <- readRDS('stream_gage_values.RData')
data <- data.pull(datestring="2025-06-01T00:00:00Z/2025-08-31T00:00:00Z")
data <- data.resolveNA(data)
saveRDS(data, 'example_data.RDS')

# prior year
#data <- readRDS('stream_gage_values_2024.RData')
data.prev <- data.pull(datestring="2024-06-01T00:00:00Z/2024-08-31T00:00:00Z")
data.prev <- data.resolveNA(data.prev)
saveRDS(data.prev, 'example_prior_data.RDS')

# use last year's data to develop a prior belief 
plot.density(data.prev)

#### Split data by gage point ####
gage.data <- split(data$discharge_max, data$monitoring_location_id)
gage.data.prev <- split(data.prev$discharge_max, data.prev$monitoring_location_id)

#### Input Data Setup ####
#gage <- names(gage.data)[1]
gage.results <-list()
for(gage in names(gage.data)){
gage_values <- gage.data[[gage]]

N_total <- length(gage_values)
window_size <- 10
results_list <- list()

#### Initial Prior Setup ####
lmoms <- est.lp3(gage.data.prev[[gage]])[['para']]
current_priors <- list(
  prior_shape = 4/lmoms[['gamma']]^2
  , prior_scale = 2*lmoms[['sigma']]/lmoms['gamma']
  , prior_location = lmoms[['mu']]-(4/lmoms[['gamma']]^2)*(2*lmoms[['sigma']]/lmoms['gamma'])
)

#### Changepoint Setup ####
# run_length_probs <- matrix(0, nrow = N_total+2, ncol = N_total+2)
# run_length_probs[10,10] <- 1
# changepoint_probs <- numeric(N_total)
# hazard_lambda <-  5

run_mcmc_segment <- function(segment_data, t, shape, scale, location) {
  segment_data <- pmax(as.numeric(segment_data), 1e-8) # ensure positive
  fit <-mod$sample(
                  data = list(N = length(segment_data), 
                              y = as.array(segment_data),
                              shape_prior=as.numeric(shape), 
                              scale_prior=as.numeric(scale), 
                              location_prior=as.numeric(location))
                  , chains = 1, parallel_chains = 1)
  log_print(fit)
  return(fit$draws())
}


#### Sequential Bayesian Update Loop (Moving Window) ####
for (t in window_size:N_total) {
  # Define the current window data
  start_index <- t - window_size + 1
  current_window_data <- log10(gage_values[start_index:t])
  
  new_rl_probs <- numeric(window_size+1)
  segments <- lapply(0:window_size, function(rl) gage_values[(t-rl):t])
  
  # Prepare data list for Stan
  stan_data <- list(
    N = window_size
    , y = current_window_data
    , shape_prior = current_priors$prior_shape
    , scale_prior = current_priors$prior_scale
    , location_prior = current_priors$prior_location
  )
  
  log_print(paste("Processing time step", t, "..."))
  
  # posteriors <- parallel::mclapply(segments, run_mcmc_segment,
  #                        t=t, mc.cores = 12,
  #                        shape=current_priors$prior_shape, 
  #                        scale=current_priors$prior_scale, 
  #                        location=current_priors$prior_location)
  # 
  # for (rl in 1:(window_size)) {
  #   post <- colMeans(posteriors[[rl]])
  #   pred_prob <- mean(FAdist::dgamma3(current_window_data, 
  #                                     shape = post[2], 
  #                                     scale = post[3]+1e-6,
  #                                     thres = post[4])) # post$location
  #   new_rl_probs[rl] <- run_length_probs[rl, t] * pred_prob * (1 - 1/hazard_lambda)
  # }
  # 
  # # changepoint probability
  # cp_prob <- sum(run_length_probs[(t-window_size):t, t] * (1/hazard_lambda))
  # new_rl_probs[1] <- cp_prob
  # log_print(paste("Change Point Probability", cp_prob))
  # #normalize
  # new_rl_probs <- new_rl_probs / sum(new_rl_probs)
  # run_length_probs[(t-window_size+1):(t+1), t+1] <- new_rl_probs
  # changepoint_probs[t] <- new_rl_probs[1]

  # Fit the model using cmdstanr
  fit <- mod$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    refresh = 0, # Suppress intermediate output
    show_messages = FALSE
  )
  
  log_print(fit)
  # Extract posterior samples
  draws <- fit$draws(variables = c("shape", "scale", "location")
                     , format = "df")
  
  # Store results (e.g., mean of posteriors)
  results_list[[t]] <- colMeans(draws)
  
  # Update priors for the next iteration (Sequential Updating)
  current_priors$shape_prior <- mean(draws$shape)
  current_priors$scale_prior <- mean(draws$scale)
  current_priors$location_prior <- mean(draws$location)
  
  }

# Analyze Results
results_df <- bind_rows(results_list)
results_df <- results_df %>%
  mutate(time = window_size:N_total)

# Calculate a simple "changepoint score" by looking at the first difference of mu
# (The magnitude of the jump in the parameter mean)
results_df <- results_df %>%
  mutate(
    # rolling Z-score
    shape_zscore = (shape - zoo::rollapply(shape, 5, mean, fill = NA, align = "right")) /
      zoo::rollapply(shape, 5, sd, fill = NA, align = "right")
    , scale_zscore = (scale - zoo::rollapply(scale, 5, mean, fill = NA, align = "right")) /
      zoo::rollapply(scale, 5, sd, fill = NA, align = "right")
    , location_zscore = (location - zoo::rollapply(location, 5, mean, fill = NA, align = "right")) /
      zoo::rollapply(location, 5, sd, fill = NA, align = "right")
  )

gage.results[[gage]] <- results_df
}


# plot results 
par(mfrow = c(3,1), mar = c(4,4,4,4))
for(gage in names(gage.data)){
  results_df<-gage.results[[gage]]
# Plot the Z-score to see where changes happen
plot(results_df$time, abs(results_df$shape_zscore), type = 'l', col = 'blue',
     main = "Changepoint Indicator (Shape Z-score)"
     , xlab = "Time Step", ylab = "Z-Score"
     , ylim=c(0, 2.5))
abline(h = 1.5, col = 'red', lty = 2) # Threshold at 2 standard deviations

plot(results_df$time, abs(results_df$scale_zscore), type = 'l', col = 'blue',
     main = "Changepoint Indicator (Scale Z-score)"
     , xlab = "Time Step", ylab = "Z-Score"
     , ylim=c(0, 2.5))
abline(h = 1.5, col = 'red', lty = 2) # Threshold at 2 standard deviations

plot(results_df$time, abs(results_df$location_zscore), type = 'l', col = 'blue',
     main = "Changepoint Indicator (Location Z-score)"
     , xlab = "Time Step", ylab = "Z-Score"
     , ylim=c(0, 2.5))
abline(h = 1.5, col = 'red', lty = 2) # Threshold at 2 standard deviations
}

par(mfrow = c(4,1), mar = c(2,2,2,2))
for(gage in names(gage.results)){
  results_df <- gage.results[[gage]]

location_cps <- results_df %>%
  filter(abs(location_zscore) > 1.7)
scale_cps <- results_df %>%
  filter(abs(scale_zscore) > 1.7)
shape_cps <- results_df %>%
  filter(abs(shape_zscore) > 1.7)

plot(log10(gage_values), type='l'
     , main = paste('Gage',gage)
     , xlab='timestep', ylab='discharge')
abline(v = location_cps$time, col='red', lty=2)
abline(v = shape_cps$time, col='blue', lty=3)
abline(v = scale_cps$time, col='green', lty=4)
}
plot.new()
legend("center", legend=c('location change'
                          ,'shape chage'
                          ,'scale change')
       ,col = c('red', 'blue', 'green')
       ,lty = c(2,3,4))

saveRDS(gage.results, "results.RDS")





