#-------------------------#
#### Utility Functions ####
#-------------------------#
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dataRetrieval
               , dplyr
               , tidyr
               , lmom
               , lmomco
               , kdensity
               )

#' Pull and Format Data
#'
#' Uses USGS dataRetrieval API to pull daily summary values from stream gages.
#' 
#' @param datestring A string indicating the start and end dates
#' @param selected_sites A list of USGS monitoring site IDs
#' @param parameters A list of USGS parameter codes
#' @return data frame of reported gage values
data.pull <- function(
    datestring="2025-06-01T00:00:00Z/2025-08-31T00:00:00Z",
    selected_sites = c('USGS-08166250','USGS-08166200', 'USGS-08166140'),
    parameters = c("00060", #  Discharge (cubic feet per second)
                   "00065") # Gage height (feet)
){
  data <- dataRetrieval::read_waterdata_daily(
    monitoring_location_id =  selected_sites,
    parameter_code = parameters ,
    time=datestring
  )
  
  formatted_data <- data %>% 
    mutate(parameter = case_match(parameter_code, 
                                  "00060"~"discharge", 
                                  "00065"~"gage_height"),
           statistic = case_match(statistic_id,
                                  "00001" ~ "max",
                                  "00002" ~ "min",
                                  "00003" ~ "mean")) %>%
    pivot_wider(id_cols = c(monitoring_location_id, time, geometry),
                names_from = c(parameter, statistic), values_from=value)
  
  return(formatted_data)
}

#' Handling Missing Data
#' 
#' consistently handle NAs the same way every time.
#' 
#' @param df data frame from the data.pull() function
#' @return data frame with no NA values
data.resolveNA<- function(df){
  numeric_cols <- names(df)[sapply(df, is.numeric)]
  
  # Loop through each numeric column
  for (col_name in numeric_cols) {
    
    # Get the column as a vector
    x <- df[[col_name]]
    n <- length(x)
    i <- 1
    
    # Loop through the rows of the current column
    while (i <= n) {
      if (is.na(x[i])) {
        # Find the starting index of the NA sequence
        start_na <- i
        
        # Find the index of the next non-NA value
        end_na <- i
        while (end_na <= n && is.na(x[end_na])) {
          end_na <- end_na + 1
        }
        
        # Get the value of the previous non-NA row
        prev_val <- if (start_na > 1) x[start_na - 1] else NA
        
        # Get the value of the next non-NA row
        next_val <- if (end_na <= n) x[end_na] else NA
        
        # Replace the sequential NAs with the average
        if (!is.na(prev_val) && !is.na(next_val)) {
          fill_value <- (prev_val + next_val) / 2
          x[start_na:(end_na - 1)] <- fill_value
        }
        
        # Move the index past the filled NA sequence
        i <- end_na
      } else {
        # Move to the next element if it's not NA
        i <- i + 1
      }
    }
    
    # Assign the modified vector back to the data frame
    df[[col_name]] <- x
  }
  return(df)
}

#' Estimate Log-PearsonIII parameters
#' 
#' uses L-moments to estimate the mean, variance, and skewness
#'
#' @param flow_values vector of reported flow discharge values
#' @return An L-moment object with mu, sigma, and gamma 
est.lp3 <- function(flow_values){
  log_data <- log10(flow_values)
  sample_lmoms <- lmom::samlmu(log_data)
  slmom_vec <- lmomco::vec2lmom(sample_lmoms)
  lp3_parameters <- lmomco::parpe3(slmom_vec)
  return(lp3_parameters)
}

#' Plot gage value densities
#' 
#' Plot empirical density with estimated log-PearsonIII theoretical distribution
#'
#' @param data data frame from data.pull()
#' @param variable column name of variable to plot density of
#' @param kernel name of kernel type for plotting empirical density
#' @param adjust numerical value of smoothing bandwidth for empirical density
#' @param logx if TRUE, uses the log-density
plot.density <- function(data, variable='discharge_max',
                         kernel='gamma', adjust = 1, 
                         logx=TRUE){
  sensor_ids <- unique(data$monitoring_location_id)
  par(mfrow=c(length(sensor_ids),1))
  for(sensor in sensor_ids){
    sensor_data <- subset(data, 
                          monitoring_location_id == sensor)
    x <- sensor_data[[variable]]
    x <- x[!is.na(x)]
    x_vals <- seq(min(x), max(x), length.out = length(x))
    if(logx){
      x[x==0] <- 0.1
      x <- log10(x)
      x <- x[x>=0]
    }
    # Plot empirical density
    dens <- kdensity(x, bw='SJ', adjust=adjust, kernel=kernel)
    plot(dens, main = paste("Sensor:", sensor)
         , xlab = variable
         , ylab = "Density"
         , col = "black", lwd = 2)
    
    flow_data <- sensor_data[[variable]]
    lp3_parameters <- est.lp3(flow_data)
    log_bounds <- quape3(c(0.001, 0.999), lp3_parameters)
    log_x_values <- seq(log_bounds[1], log_bounds[2], length.out = 100)
    
    # Convert log values back to the original scale
    x_values <- exp(log_x_values)
    # Calculate density on log scale
    density_log_scale <- pdfpe3(log_x_values, lp3_parameters)
    
    # Transform back to the original scale density
    lp3_density <- density_log_scale
    lines(log_x_values, lp3_density, col = 'red', lwd = 2, lty = 2)
    
  }
}


plot.changepoints <- function(gage.data, results, image_filename=""){
  if(length(image_filename)>1){
    png(image_filename
        , width = 800
        , height = 400*length(gage.data)
        , res = 150)
  }
  par(mfrow = c(length(gage.data),1), mar = c(4, 4, 4, 4))
  
  for(gage in names(gage.data)){
    log_data <- log10(gage.data[[gage]])
    result <- results[[gage]]
    scale_factor <- max(log_data) / max(result$changepoint_probs)
    
    # Plot log data
    plot(log_data, type = "l", col = "blue", lwd = 2,
         ylab = "Log-Transformed Data", xlab = "Time Step",
         main = paste("Sensor", gage))
    
    # Add changepoint probabilities on secondary axis
    par(new = TRUE)
    z<-result$changepoint_probs * scale_factor
    labels <- pretty(result$changepoint_probs, n=6)
    plot(result$changepoint_probs * scale_factor, type = "l", col = "orange", lwd = 2, axes = FALSE, xlab = "", ylab = "")
    axis(side = 4, at = seq(min(z), max(z), length.out = length(labels)), 
         labels = labels)
    mtext("Changepoint Probability", side = 4, line = 3)
    
    # Add legend
    legend("bottomright", legend = c("Log Data", "Changepoint Probability"),
           col = c("blue", "orange"), lwd = 2)
  }
  if(length(image_filename)>1){dev.off()}
  par(mfrow = c(1,1), mar = c(1, 1, 1, 1)) # reset plot settings
}
