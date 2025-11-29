#-----------------------------#
#### Gamma Conjugate Prior ####
#-----------------------------#
# An analytical solution to Bayesian Online Changepoint Detection
# using a gamma-gamma conjugate prior
if (!require("pacman")) install.packages("pacman")
pacman::p_load(Rcpp)


#### Define Functions ####
source("../utils.R")
sourceCpp(file='bocpd_gamma.cpp')

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

# use 2024 data to develop a prior belief 
plot.density(data.prev)

#### split data by gage point ####
gage.data <- split(data$discharge_max, data$monitoring_location_id)
gage.data.prev <- split(data.prev$discharge_max, data.prev$monitoring_location_id)

# stream flows from gage 1 (northernmost) to gage 2 to gage 3 (southernmost)
# gage1 <- gage_data[["USGS-08166140"]]
# gage2 <- gage_data[["USGS-08166200"]]
# gage3 <- gage_data[["USGS-08166250"]]

results <-  list()
for(gage in names(gage.data)){
  # gage values to assess in "real-time"
  gage_values <- gage.data[[gage]]
  
  # prior belief only based on past observed values
  prior <- est.lp3(gage.data.prev[[gage]])[['para']]
  
  # translate L-moments into log Pearson III parameters
  scale <- 2*prior[['sigma']]/prior[['gamma']]
  shape <- 4/prior[['gamma']]^2
  location <- prior[['mu']] - scale*shape
  
  results[[gage]] <- bocpd_pearsonIII(data = log10(gage_values)
                    , prior_shape = shape
                    , prior_scale = scale
                    , location = location
                    , hazard_lambda = 2
                    , adaptive_hazard = TRUE)
  
}

saveRDS(results, 'bocpd_conjugate_results.RData')

# plot results
image_filename <- "bocpd_cojugate_cpp_withData.png"
plot.changepoints(gage.data, results, image_filename)
