# This code will produce a heatmap/dot plot showing the relative read abundance
# and number of ASVs at the bottom, mid, and surface levels for different animal
# groups. There are two `heatmap_data` options: one for plotting animals and the
# other (commented out) for plotting non-animal taxa.

# Set libraries
library(tidyverse)
library(ggtext)

# Prepare data for heatmap for animals
heatmap_data <- melted_all %>%
  dplyr::filter( 
    !is.na(Phylum),
    !is.na(Class),
    Kingdom == "Animalia", 
    Phylum != "Myzozoa"
    ) %>%
  # Remove unnecessary columns
  dplyr::select(OTU, Depth_class, Site, Abundance, Phylum, Class) %>% 
  # Add column for abundance and number of ASVs (OTUs)
  group_by(Depth_class, Site, Phylum, Class) %>% 
  dplyr::summarise(Abundance = sum(Abundance),
                   n_ASV = n_distinct(OTU),
                   .groups = "drop") %>%
  # Add column for relative abundance
  group_by(Depth_class, Site) %>%
  dplyr::mutate(relative_abundance = Abundance/sum(Abundance)) %>%
  ungroup()

# Prepare data for heatmap for animals for non-animal taxa
# heatmap_data <- melted_all %>%
#   mutate(Kingdom = case_when(Phylum == "Myzozoa" ~ "Chromista",
#                              .default = Kingdom)) %>%
#   dplyr::filter( 
#     !is.na(Kingdom),
#     !is.na(Phylum),
#     Kingdom != "Animalia") %>%
#   dplyr::select(OTU, Depth_class, Site, Abundance, Kingdom, Phylum) %>% 
#   group_by(Depth_class, Site, Kingdom, Phylum) %>% 
#   dplyr::summarise(Abundance = sum(Abundance),
#                    n_ASV = n_distinct(OTU),
#                    .groups = "drop") %>%
#   group_by(Depth_class, Site) %>%
#   dplyr::mutate(relative_abundance = Abundance/sum(Abundance)) %>%
#   ungroup()

# Create labels for phyla
heatmap_data_labels <- heatmap_data %>%
  # Find the (alphabetically) first class in each phylum
  distinct(Phylum, Class) %>% 
  arrange(Phylum, Class) %>%
  group_by(Phylum) %>%
  slice_head(n = 1) %>%
  # Add column with (bolded) phylum name and class name
  mutate(taxa_labels = paste0("<span style = 'color:#0B0405FF;'>**", Phylum, "**</span><br>",Class,"<br>"))

# Combine labels with data 
# This will replace the first class label with the combine phylum/class label
# created above
heatmap_data_with_labels <- heatmap_data %>%
  left_join(heatmap_data_labels, by = join_by(Phylum, Class)) %>%
  mutate(taxa_labels = case_when(is.na(taxa_labels) ~ Class,
                                 .default = taxa_labels))


# Plot data
heatmap_data_with_labels %>% 
  ggplot(aes(Site, forcats::fct_rev(taxa_labels), fill = relative_abundance, size = n_ASV)) +
  # Set point shape and stroke
  geom_point(stroke = 0.2, shape = 21, color = "black") +
  # Set point size range
  scale_size_continuous(range = c(2,7), breaks = c(5, 50, 500)) +
  # Set viridis color palette
  viridis::scale_fill_viridis(option = "mako", 
                              direction = -1, 
                              oob = scales::squish) +
  # Facet by phylum and depth
  facet_grid(Phylum~Depth_class, scales = "free", space = "free", switch = "y") + 
  # Set legend labels
  labs(fill = "Relative read\nabundance\n(per site)", size = "Number of\nASVs") + 
  # Set theming
  theme_classic() +
  theme(text = element_text(size = 10, family = "Arial"),
        axis.line = element_blank(),
        panel.spacing = unit(10, "pt"),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey95"),
        axis.text = element_markdown(size = 10, family = "Arial"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text.y.left = element_blank(),
        strip.background.x = element_rect(fill = "grey95", color = "white", linewidth = 0),
        strip.text.x = element_text(size = 10, family = "Arial", face = "bold"),
        plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
        legend.position = "right",
        legend.justification = "bottom",
        legend.title.position = "top",
        legend.title = element_text(size = 10, hjust = 0.5, vjust = 0.5))

# Save plot
ggsave(here::here("plots/heatmap_all_animals.svg"), width = 8, height = 9, units = "in", dpi = 300) 


