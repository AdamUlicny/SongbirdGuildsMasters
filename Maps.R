# Load required packages
library(sf) 
library(ggplot2)
library(ggspatial)
library(leaflet)
library(htmlwidgets)

# Read KML files
study_area <- st_read("resources/koda_study_area.kml")
transect<- st_read("resources/koda_point_transect.kml")
transect_points <- st_cast(transect, "POINT")

# Transform data to appropriate CRS for ggspatial
study_area <- st_transform(study_area, 3857)  # Web Mercator
transect_points <- st_transform(transect_points, 3857)

# Create a buffer around the study area for better context
buffer_size <- 200  # meters
study_area_buffered <- st_buffer(study_area, dist = buffer_size)


static_map_area <- ggplot() +
  # Add OpenStreetMap tiles (this uses a different method than before)
  annotation_map_tile(type = "osm", zoom = 15) +
  # Study area polygon
  geom_sf(data = study_area, fill = "blue", alpha = 0.3, color = "darkblue", size = 1) +
  # Point transect points
  geom_sf(data = transect_points, color = "red", size = 4, shape = 16) +
  # Add scale bar
  annotation_scale(location = "bl", width_hint = 0.3) +
  # Add north arrow
  annotation_north_arrow(location = "tr", height = unit(1, "cm"), width = unit(1, "cm")) +
  # Crop to the buffered bounding box
  coord_sf(xlim = c(st_bbox(study_area_buffered)[1], st_bbox(study_area_buffered)[3]), 
           ylim = c(st_bbox(study_area_buffered)[2], st_bbox(study_area_buffered)[4]),
           expand = FALSE) +
  # Add title and theme
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

# Load required packages
library(ggplot2)
library(sf)
library(rnaturalearth)

library(rnaturalearthdata)

# Define study site coordinates
coords <- data.frame(
  Latitude = c(-27.23, -18.87, -12.50, -35.65, -12.47, 30.56, 36.49, 43.94, 
               29.71, 11.09, 23.45, 27.29, 22.43, 5.11, 10.20, 40.85, 49.93, 
               -29.48, -30.49, -36.93, -32.78, -26.40, -41.30, -33.80),
  Longitude = c(152.74, 146.13, 130.95, 143.64, 130.83, -109.75, -121.70, -71.70, 
                -82.46, 76.77, 120.97, 88.68, 114.18, 100.99, 77.50, -3.94, 14.11, 
                146.13, 151.74, 149.33, 116.95, 116.23, 148.12, 150.53)
)


# Load a world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create the plot
map <- ggplot(data = world) +
  geom_sf(fill = "gray90", color = "black") +  # Draw world map
  geom_jitter(data = coords, aes(x = Longitude, y = Latitude), 
              color = "red", size = 4, width = 1, height = 1) +       # Plot study sites
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "")

# Display the plot
print(map)

# Save as an image file
ggsave("study_sites_map.png", map, width = 8, height = 6, dpi = 300)



################## Rayshader 3D Map ###############################
# Load required packages
library(sf)
library(rayshader)
library(elevatr)
library(raster)

# Read KML files
study_area <- st_read("resources/koda_study_area.kml")
transect <- st_read("resources/koda_point_transect.kml")
transect_points <- st_cast(study_area, "POINT")

# Transform data to WGS84 (for elevation)
study_area_wgs84 <- st_transform(study_area, 4326)
transect_points_wgs84 <- st_transform(transect_points, 4326)

# Create a buffered bounding box for elevation data
buffer_size_degrees <- 0.9  # Buffer in degrees
study_area_buffered <- st_buffer(study_area_wgs84, dist = buffer_size_degrees)
bbox <- st_bbox(study_area_buffered)

# Get elevation data
elev_raster <- get_elev_raster(study_area_buffered, z = 12, src = "aws")
elev_matrix <- raster_to_matrix(elev_raster)

# Define extent for mapping coordinates
extent_mat <- c(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"])

# ---- Green Terrain Shading ----
texture_map <- elev_matrix %>%
  sphere_shade(texture = "imhof1") %>%  # "imhof1" is a green terrain texture
  add_shadow(ray_shade(elev_matrix, zscale = 3, lambert = TRUE), 0.5) %>%
  add_shadow(ambient_shade(elev_matrix), 0.5)

# ---- Plot 3D Terrain ----
texture_map %>%
  plot_3d(
    elev_matrix, zscale = 5, fov = 0, theta = 20, phi = 35,
    windowsize = c(1000, 800), zoom = 0.75
  )

# ---- Add Transect Points ----
points_coords <- st_coordinates(transect_points_wgs84)
point_elevations <- raster::extract(elev_raster, points_coords[, 1:2]) + 50  # Small vertical offset

for (i in 1:nrow(points_coords)) {
  render_points(
    extent = extent_mat,
    lat = points_coords[i, 2],
    long = points_coords[i, 1],
    altitude = point_elevations[i],
    zscale = 10,
    size = 6,
    color = "red"
  )
}

# ---- Add Study Area Polygon ----
poly_coords <- st_coordinates(study_area_wgs84)

# Ensure poly_coords has multiple parts before using L1
if ("L1" %in% colnames(poly_coords)) {
  unique_parts <- unique(poly_coords[, "L1"])  # Correct indexing
  
  for (part in unique_parts) {
    part_coords <- poly_coords[poly_coords[, "L1"] == part, , drop = FALSE]  # Drop = FALSE to retain matrix structure
    path_long <- part_coords[, "X"]
    path_lat <- part_coords[, "Y"]
    path_alt <- rep(max(point_elevations, na.rm = TRUE) + 100, length(path_lat))
    
    render_path(
      extent = extent_mat,
      lat = path_lat,
      long = path_long,
      altitude = path_alt,
      zscale = 10,
      color = "blue",
      linewidth = 3
    )
  }
} else {
  # Single-part polygon case
  path_long <- poly_coords[, "X"]
  path_lat <- poly_coords[, "Y"]
  path_alt <- rep(max(point_elevations, na.rm = TRUE) + 100, length(path_lat))
  
  render_path(
    extent = extent_mat,
    lat = path_lat,
    long = path_long,
    altitude = 4690,
    zscale = 50,
    color = "blue",
    linewidth = 3 
  )
}

# ---- Add Compass and Scale Bar ----
render_compass(position = "E")
render_scalebar(position = "S", distance = 1000, unit = "m")

# ---- Save Snapshot ----
render_snapshot("3d_study_area_green.png", clear = TRUE)

# ---- Render High-Quality Image ----
render_highquality(
  filename = "3d_study_area_green_hq.png",
  width = 3000,
  height = 2000,
  samples = 256,
  clear = TRUE
)

