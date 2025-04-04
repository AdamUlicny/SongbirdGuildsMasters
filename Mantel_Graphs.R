# Load necessary library
library(ggplot2)

# --- Configuration ---

# Define Node Positions (Adjust these coordinates to fine-tune spacing if needed)
# These values create a reasonably tight equilateral triangle.
positions <- data.frame(
  name = c("Fylogeneze", "Morfologie", "Gildy"),
  x = c(0, -0.4, 0.4),      # X-coordinates
  y = c(0.4, -0.2, -0.2)   # Y-coordinates (adjusted slightly for visual balance)
  # You could make y = c(0.5, -0.25, -0.25) for a perfect equilateral triangle
)

# --- Data Preparation ---

# !! REPLACE THESE WITH YOUR ACTUAL MANTEL STATISTIC VALUES !!
# Placeholder Mantel test statistics (replace with your actual data)
mantel_Global_P_G1 <- list(statistic = 0.149) # Example value
mantel_Global_P_M  <- list(statistic = 0.110) # Example value
mantel_Global_M_G1 <- list(statistic = 0.091) # Example value

# Define Edges connecting the nodes
edges <- data.frame(
  from = c("Fylogeneze", "Morfologie", "Gildy"),
  to = c("Gildy", "Fylogeneze", "Morfologie"),
  weight = c(mantel_Australia_P_G1[["statistic"]],
             mantel_Australia_P_M[["statistic"]],
             mantel_Australia_M_G1[["statistic"]])
)

# Merge node positions with edge data to get coordinates for geom_segment
df_edges <- merge(edges, positions, by.x = "from", by.y = "name")
df_edges <- merge(df_edges, positions, by.x = "to", by.y = "name", suffixes = c("_from", "_to"))

# Calculate midpoint coordinates for placing edge labels
df_edges$label_x <- (df_edges$x_from + df_edges$x_to) / 2
df_edges$label_y <- (df_edges$y_from + df_edges$y_to) / 2

# --- Plotting ---

# Create the plot
gg_mantel <- ggplot() +
  # Draw edges (segments) between nodes
  # Edge thickness is mapped to 'weight'
  geom_segment(
    data = df_edges,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to, size = weight),
    color = "grey60", # Slightly darker gray
    show.legend = FALSE # Hide legend for edge size
  ) +
  # Define the range of edge thickness
  scale_size_continuous(range = c(1, 4)) + # Adjust range (min, max thickness)
  
  # Draw nodes (points)
  geom_point(
    data = positions,
    aes(x = x, y = y),
    color = c("lightblue", "lightgreen", "lightcoral"), # Node colors
    size = 40  # REDUCED node size
  ) +
  # Add text labels inside nodes
  geom_text(
    data = positions,
    aes(x = x, y = y, label = name),
    color = "black",
    fontface = "bold",
    size = 5 # Node label size (adjust as needed)
  ) +
  # Add text labels for edge weights near the middle of edges
  geom_text(
    data = df_edges,
    aes(x = label_x, y = label_y, label = round(weight, 3)),
    color = "black",
    size = 6, # Edge label size (adjust as needed, often slightly smaller than node labels)
    nudge_y = 0.03,
    nudge_x = 0.07# Optional: slightly move labels off the edge line
  ) +
  # Use a minimal theme
  theme_void() +
  # Set plot limits to ensure all elements fit, with some padding
  # Adjust padding (e.g., 0.1 or 0.15) if elements get cut off
  xlim(min(positions$x) - 0.15, max(positions$x) + 0.15) +
  ylim(min(positions$y) - 0.15, max(positions$y) + 0.15) +
  # Ensure aspect ratio is 1:1 so the triangle looks correct
  coord_fixed(ratio = 1)

# Print the plot to the viewer (optional)
# print(gg_mantel)

# --- Save the Plot ---

# Save the plot as an SVG file with REASONABLE dimensions
svg("mantel_graph_North_America.svg", width = 7, height = 7) # REDUCED dimensions
print(gg_mantel)
dev.off()

# You can also save as PNG:
# ggsave("mantel_graph_global_improved.png", plot = gg_mantel, width = 7, height = 7, dpi = 300)

