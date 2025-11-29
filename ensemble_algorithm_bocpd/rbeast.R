#----------------------#
#### Rbeast Example ####
#----------------------#
# example using the Rbeast package
if (!require("pacman")) install.packages("pacman")
pacman::p_load(Rbeast,logr)
source("../utils.R")

# get data
data <- data.pull()
data <- data.resolveNA(data)
saveRDS(data, 'example_data.RDS')
gage.data <- split(data$discharge_max, data$monitoring_location_id)

# set up logging
log.fp <- file.path(getwd(), paste0("rbeast.log"))
lf <- log_open(log.fp)

# run Rbeast
results <- list()
for(gage in names(gage.data)){
  o = beast(gage.data[[gage]], season='none')
  image_filename <- paste0("rbeast_", gage, ".png")
  png(image_filename, width = 800, height = 1200, res = 150)
  plot(o)
  dev.off()
  log_print(o)
  results[[gage]] <- o
}

image_filename <- paste0("changepoints.png")
png(image_filename, width = 800, height = 1200, res = 150)
par(mfrow = c(3,1), mar = c(4,4,4,4))
for(gage in names(results)){
  o <- results[[gage]]
  gage_values <- log10(gage.data[[gage]])
  
  cp <- as.numeric(na.omit(o$trend$cp))
  
  plot(log10(gage_values), type='l'
       , main = paste('Gage',gage)
       , xlab='timestep', ylab='discharge')
  abline(v = cp, col='red', lty=1)
}
dev.off()