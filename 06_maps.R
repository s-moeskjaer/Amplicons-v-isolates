library(ggplot2)
library(dplyr)
require(maps)

setwd('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/')

#https://www.datanovia.com/en/blog/how-to-create-a-map-using-ggplot2/
# Scale bar: http://oswaldosantos.github.io/ggsn/
## Sampling site points 
sample_sites <- read.csv('zz_sampling_sites.csv', sep = ';')


world_map <- map_data("world")
ggplot(world_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill="lightgray", colour = "white")

######################
# ALL SAMPLING SITES #
######################
# EU Contries to be included 
some.eu.countries <- c(
  "France", "Switzerland", "Germany",
  "Austria", "Belgium", "UK", "Netherlands",
  "Denmark", "Luxembourg"
)
# Retrievethe map data
some.eu.maps <- map_data("world", region = some.eu.countries)

# Compute the centroid as the mean longitude and lattitude
# Used as label coordinate for country's names
region.lab.data <- some.eu.maps %>%
  group_by(region) %>%
  summarise(long = mean(long), lat = mean(lat))

p = ggplot(some.eu.maps, aes(x = long, y = lat)) +
  geom_polygon(fill = "lightgray", colour = "white", aes( group = group))+
  geom_text(aes(label = region), data = region.lab.data,  size = 3, hjust = 0.5)+
  #scale_fill_viridis_d()+
  theme_void()+
  theme(legend.position = "none")+
  geom_point(data = sample_sites, aes(x=long, y=X...lat, colour="red"), alpha = 1/2) + 
  scalebar(some.eu.maps, dist = 100, dist_unit = "km",
           transform = TRUE, model = "WGS84")

ggsave('map_scale.pdf', plot = p, width = 15, height = 20, unit = 'cm', useDingbats = FALSE)


######################
     # DK ONLY #
######################

# Some EU Contries
some.eu.countries <- c(
  "Denmark"
)
# Retrievethe map data
some.eu.maps <- map_data("world", region = some.eu.countries)

# Compute the centroid as the mean longitude and lattitude
# Used as label coordinate for country's names
region.lab.data <- some.eu.maps %>%
  group_by(region) %>%
  summarise(long = mean(long), lat = mean(lat))

#remove bottom three rows from sample sites 
sample_sitesDK <- sample_sites[1:50,]

p = ggplot(some.eu.maps, aes(x = long, y = lat)) +
  geom_polygon(fill = "lightgray", colour = "white", aes( group = group))+
  geom_text(aes(label = region), data = region.lab.data,  size = 3, hjust = 0.5)+
  #scale_fill_viridis_d()+
  theme_void()+
  theme(legend.position = "none")+
  geom_point(data = sample_sitesDK, aes(x=long, y=X...lat, colour="red"), alpha = 1/2) +
  scalebar(some.eu.maps, dist = 100, dist_unit = "km",
           transform = TRUE, model = "WGS84")

ggsave('map_DK_scale.pdf', plot = p, width = 5, height = 5, unit = 'cm', useDingbats = FALSE)

