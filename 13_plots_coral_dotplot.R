# This code will produce a dot plot showing coral eDNA and video detections 
# grouped by species.

# Set libraries
library(tidyverse)
library(ggtext)

# Set theme elements
theme_simple <- function() {
  theme_classic() +
    theme(
      text = element_text(family = "Arial", size = 10),
      # lighten axis lines and ticks
      axis.line = element_line(linewidth = 0.25, colour = "grey80"), 
      axis.ticks = element_line(linewidth = 0.25, colour = "grey80"))
}

# Site order for plot
site_order <- c("AA", "LW2", "BF2", "BF3", "SR2", "SR3", "BF5", "PE1", "DS2")

# Import video data
video_data <- read.table(here::here("additional_files/coral_dotplot_video_data.tsv"), sep = "\t", header = TRUE, row.names = NULL) %>%
  dplyr::filter(Abundance > 0) %>%
  # set count size for video
  dplyr::mutate(Count = 3)
  
# Filter data for plot
ps_coral_dotplot_data <- melted_all %>%
  # Remove anything without a genus, abundance below zero and not 28S
  dplyr::filter(!is.na(Genus), Abundance > 0, Marker == "28S") %>%
  # Add a colum for detection type
  dplyr::mutate(Type = "eDNA") %>%
  # Add "sp." designation to genus-level detections
  dplyr::mutate(Species = case_when(!is.na(Genus) & is.na(Species) ~ paste0(Genus, " sp."), .default = Species)) %>%
  # Keep only coral genera
  dplyr::filter(str_detect(Genus, "Antipathes|Ellisella|Thesea|Scleracis|Bebryce|Muricea|Swiftia|Tanacetipathes|Callistephanus|Leptogorgia|Paramuricea|Cirrhipathes|Flabellum")) %>%
  # Add columns for number of ASVs (OTU) and summed abundance for species and sites
  dplyr::group_by(Species, Site) %>%
  dplyr::mutate(Count = n_distinct(OTU), Abundance = sum(Abundance)) %>%
  # Remove duplicate entries
  dplyr::distinct(Species, Site, .keep_all = TRUE) %>%
  # Remove unnecessary columns
  dplyr::select(Site, Species, Abundance, Count, Type) %>%
  ungroup() %>%
  # Fix site naming issues
  dplyr::filter(Site != "PE02", Site != "PE03", Site != "PE04", Site != "PE05") %>%
  # Join with the video data
  dplyr::full_join(video_data, by = join_by(Site, Species, Abundance, Type, Count)) %>%
  # Add column for group (helps with sorting)
  dplyr::mutate(Group = paste0(Species, " ", Type), 
                # Add markdown to italicize species names
                Species = paste0("*",Species,"*"))

# Plot data
ps_coral_dotplot_data %>% 
  # Plot data, with fill color as abundance and shape as detection type
  ggplot(aes(x = factor(Site, levels = site_order), 
             y = forcats::fct_rev(Group), 
             fill = Abundance, 
             shape = Type),
         color = "black") + 
  # Set point and stroke size
  geom_point(size = 4, stroke = 0.2) +
  # Set viridis color palette for fill
  scale_fill_viridis_c(option = "mako", direction = -1) +
  # Set shapes for detection types
  scale_shape_manual(values = c(21, 24)) +
  # Facet by species to group detection types together
  facet_grid(Species ~ ., scales = "free", space = "free", switch = "y") +
  # Make the shapes bigger on the legend
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  # Theming
  theme_simple() +
  theme(axis.text = element_text(family = "Arial", size = 10),
        axis.title = element_text(family = "Arial", size = 10),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.25, linetype = 3, color = "gray80"),
        legend.text = element_text(family = "Arial", size = 10),
        legend.justification = "bottom",
        strip.background = element_blank(),
        strip.text.y.left = element_markdown(family = "Arial", size = 10, angle = 0, hjust = 1),
        legend.title = element_text(hjust = 0.5)) +
  # Nicer legend labels
  labs(fill = "Number of\nreads (eDNA) or\ncounts (video)", 
       shape = "Detection type") +
  # Remove x and y axis labels
  xlab("") +
  ylab("")

# Save plot
ggsave(here::here("plots/coral_dotplot_redux.svg"), plot = last_plot(), width = 6.5, height = 6, units = "in", dpi = 300)
