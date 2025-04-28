######### Study area maú
library(sf) 
library(ggplot2)
library(ggspatial)
library(leaflet)
library(htmlwidgets)


study_area <- st_read("resources/koda_study_area.kml")
transect<- st_read("resources/koda_point_transect.kml")
transect_points <- st_cast(transect, "POINT")


study_area <- st_transform(study_area, 3857) 
transect_points <- st_transform(transect_points, 3857)


buffer_size <- 200  
study_area_buffered <- st_buffer(study_area, dist = buffer_size)


static_map_area <- ggplot() +
  annotation_map_tile(type = "osm", zoom = 15) +
  geom_sf(data = study_area, fill = "blue", alpha = 0.3, color = "darkblue", size = 1) +
  geom_sf(data = transect_points, color = "red", size = 4, shape = 16) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tr", height = unit(1, "cm"), width = unit(1, "cm")) +
  coord_sf(xlim = c(st_bbox(study_area_buffered)[1], st_bbox(study_area_buffered)[3]), 
           ylim = c(st_bbox(study_area_buffered)[2], st_bbox(study_area_buffered)[4]),
           expand = FALSE) +
  labs(caption = "Background: OpenStreetMap") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())



# Print map
print(static_map_area)

# Save map to file
ggsave("map_point_transect.png", static_map_area, width = 10, height = 5, dpi = 300)



# Calculate study area in km2
print(st_crs(study_area))
st_area_projected <- st_transform(study_area, 5514)
print(st_crs(st_area_projected))
study_area_km2 <- st_area(st_area_projected)
print(study_area_km2/1000000)



################# Locations map for Meta ###############################

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

coords <- data.frame(
  Latitude = c(-27.23, -18.87, -12.50, -35.65, -12.47, 30.56, 36.49, 43.94, 
               29.71, 11.09, 23.45, 27.29, 22.43, 5.11, 10.20, 40.85, 49.93, 
               -29.48, -30.49, -36.93, -32.78, -26.40, -41.30, -33.80),
  Longitude = c(152.74, 146.13, 130.95, 143.64, 130.83, -109.75, -121.70, -71.70, 
                -82.46, 76.77, 120.97, 88.68, 114.18, 100.99, 77.50, -3.94, 14.11, 
                146.13, 151.74, 149.33, 116.95, 116.23, 148.12, 150.53)
)


world <- ne_countries(scale = "medium", returnclass = "sf")¨¨

map <- ggplot(data = world) +
  geom_sf(fill = "gray90", color = "black") +  # Draw world map
  geom_jitter(data = coords, aes(x = Longitude, y = Latitude), 
              color = "red", size = 4, width = 1, height = 1) +       # Plot study sites
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "")


print(map)

# Save as a png
ggsave("study_sites_map.png", map, width = 8, height = 6, dpi = 300)




