#------------------------#
#### Data Exploration ####
#------------------------#
# Exploration of the example River Gage data from Kerr County Texas

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dataRetrieval
               ,nhdplusTools
               ,leaflet
               ,dplyr
               ,tidyr
               ,sf)


par(mfrow=c(1,1))

#### USGS Gauge Sites ####
sites <- read_waterdata_monitoring_location(
  state_code = "48", # Texas
  county_code = "265", # Kerr County
  site_type_code ='ST')  # Stream

# map
leaflet(data = sites) %>%
  addTiles() %>%
  addMarkers(data=sites, popup = ~monitoring_location_id) %>%
  addCircleMarkers(data=sites,
                   radius = 5, 
                   color = "blue", 
                   stroke = FALSE, 
                   fillOpacity = 0.7) %>%
  setView(lng = -99.14, lat = 30.05, zoom = 10)

#### Flowlines ####
site_ids<- unique(sites$monitoring_location_id)
pour_points <- data.frame(featureSource = rep("nwissite", length(site_ids)), 
                    featureID = site_ids) 
upstream_network<-apply(pour_points, 1, navigate_nldi, mode="upstreamMain")
flowlines <- lapply(upstream_network, function(x) x$UM_flowlines)
plot(sites$geometry, main='Flowlines')
for(i in seq(1,length(upstream_network))){
  if(length(upstream_network[[i]])==2){
    plot(upstream_network[[i]]$UM_flowlines, 
          color = upstream_network[[i]]$nhdplus_comid, lwd = 2, add=T)
  }
}

#### Discharge and Gage Height ####
# parameter codes:
# 00060 = Discharge (cubic feet per second)
# 00065 = Gage height (feet)
parameters <- c("00060", "00065")
date_range <- "2025-06-01T00:00:00Z/2025-08-31T00:00:00Z"

data <- read_waterdata_daily(
  monitoring_location_id =  sites$monitoring_location_id,
  parameter_code = parameters ,
  time=date_range
)

# Split data by parameter
discharge_data <- data %>% 
  subset(parameter_code == "00060") %>%
  subset(statistic_id=='00003') #mean value
gage_height_data <- data %>% 
  subset(parameter_code == "00065") %>%
  subset(statistic_id=='00003')

# Plot Discharge
plot(NULL, 
     xlim = range(discharge_data$time, na.rm=T), 
     ylim = range(discharge_data$value, na.rm=T),
     xlab = "Date", ylab = "Discharge (cfs)", 
     main = "Discharge Over Time by Site")
discharge_sites <- unique(discharge_data$monitoring_location_id)
colors <- rainbow(length(discharge_sites))
for (i in seq_along(discharge_sites)) {
  site_data <- subset(discharge_data, 
                      monitoring_location_id == discharge_sites[i])
  lines(site_data$time, site_data$value, col = colors[i], lwd = 2)
}
legend("right", legend = discharge_sites, col = colors, lwd = 2
       , cex=0.8, text.width = 10)

# Plot Gage Height
plot(NULL, 
     xlim = range(gage_height_data$time, na.rm=T), 
     ylim = range(gage_height_data$value, na.rm=T),
     xlab = "Date", ylab = "Gage Height (ft)", 
     main = "Gage Height Over Time by Site")
gh_sites <- unique(gage_height_data$monitoring_location_id)
colors <- rainbow(length(gh_sites))
for (i in seq_along(gh_sites)) {
  site_data <- subset(gage_height_data, 
                      monitoring_location_id == gh_sites[i])
  lines(site_data$time, site_data$value, col = colors[i], lwd = 2 )
}
legend("topright", legend = gh_sites, col = colors, lwd = 2
       , cex=0.8, text.width = 10)

discharge_sites==gh_sites

flow_sites <- subset(sites, monitoring_location_id %in% gh_sites)

plot(flow_sites$geometry, col=colors, main="Gage locations")
legend("topright", legend = flow_sites$monitoring_location_id, 
       col = colors, lwd = 2, cex=0.5)

site_ids<- unique(flow_sites$monitoring_location_id)
pour_points <- data.frame(featureSource = rep("nwissite", length(site_ids)), 
                          featureID = site_ids) 
upstream_network<-apply(pour_points, 1, navigate_nldi, mode="upstreamMain")
flowlines <- upstream_network$UM_flowlines
plot(flow_sites$geometry,  col=colors, pch=16, main="Gage locations with flowlines")
legend("topright", legend = flow_sites$monitoring_location_id, 
       col = colors, lwd = 2, cex=0.5)
for(i in seq(1,length(upstream_network))){
  if(length(upstream_network[[i]])==2){
    plot(upstream_network[[i]]$UM_flowlines, 
         col = "cornflowerblue", lwd = 2, add=T)
  }
}

#### Data Selection & Formatting ####
selected_sites <- c('USGS-08166250','USGS-08166200','USGS-08166140')
clean_data <- data %>% 
  subset(monitoring_location_id %in% selected_sites) %>%
  mutate(parameter = case_match(parameter_code, 
                                "00060"~"discharge", 
                                "00065"~"gage_height"),
         statistic = case_match(statistic_id,
                                "00001" ~ "max",
                                "00002" ~ "min",
                                "00003" ~ "mean")) %>%
  pivot_wider(id_cols = c(monitoring_location_id, time, geometry),
              names_from = c(parameter, statistic), values_from=value)

#write.csv(clean_data, "stream_gage_values.csv")
#saveRDS(clean_data, "stream_gage_values.RData")
# USGS-08166140 feeds into -> USGS-08166200 feeds into -> 08166250
#save.image(file="USGS_dataset.RData")
