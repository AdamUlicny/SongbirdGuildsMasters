library(ggplot2)

positions <- data.frame(
  name = c("Fylogeneze", "Morfologie", "Gildy"),
  x = c(0, -0.4, 0.4),     
  y = c(0.4, -0.2, -0.2)  
)

# Change _Asia_ to whatever else. Need prepared DF's from Meta.R
edges <- data.frame(
  from = c("Fylogeneze", "Morfologie", "Gildy"),
  to = c("Gildy", "Fylogeneze", "Morfologie"),
  weight = c(mantel_Asia_P_G1[["statistic"]],
             mantel_Asia_P_M[["statistic"]],
             mantel_Asia_M_G1[["statistic"]])
)

df_edges <- merge(edges, positions, by.x = "from", by.y = "name")
df_edges <- merge(df_edges, positions, by.x = "to", by.y = "name", suffixes = c("_from", "_to"))


df_edges$label_x <- (df_edges$x_from + df_edges$x_to) / 2
df_edges$label_y <- (df_edges$y_from + df_edges$y_to) / 2

gg_mantel <- ggplot() +
  geom_segment(
    data = df_edges,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to, size = weight),
    color = "grey60", 
    show.legend = FALSE 
  ) +
  scale_size_continuous(range = c(1, 4)) + 
  geom_point(
    data = positions,
    aes(x = x, y = y),
    color = c("lightblue", "lightgreen", "lightcoral"),
    size = 40  
  ) +
  geom_text(
    data = positions,
    aes(x = x, y = y, label = name),
    color = "black",
    fontface = "bold",
    size = 5 
  ) +
  geom_text(
    data = df_edges,
    aes(x = label_x, y = label_y, label = round(weight, 3)),
    color = "black",
    size = 6,
    nudge_y = 0.03,
    nudge_x = 0.07
  ) +
  theme_void() +
  xlim(min(positions$x) - 0.15, max(positions$x) + 0.15) +
  ylim(min(positions$y) - 0.15, max(positions$y) + 0.15) +
  coord_fixed(ratio = 1)


# don't forget to change _Asia to specific continent
svg("mantel_graph_Asia.svg", width = 7, height = 7) 
print(gg_mantel)
dev.off()



