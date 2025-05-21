# This code will produce a plot showing fish eDNA detections by depth. Note: Fish 
# order names were replaced with fish icons in Adobe Illustrator 2024.

# Set libraries
library(tidyverse)
library(ggtext) # element_markdown

# Import fish depth data
fish <- read.table(here::here("additional_files/fish_depths_NEW_DATABASE.tsv"), header=TRUE, sep="\t")
background <- read.table(here::here("additional_files/fish_depth_background.tsv"), header=TRUE, sep="\t")

# Prepare data for plotting
fish_plot_data <- melted_all %>%
  # Only keep 12S data with at least species-level taxonomy and abundance greater than zero
  dplyr::filter(Marker == "12S", !is.na(Species), Abundance > 0) %>%
  # Add fish depth data
  left_join(fish) %>%
  # Sum abundance per species and sample
  group_by(Species, Sample) %>%
  mutate(Abundance = sum(Abundance)) %>%
  ungroup() %>%
  # Keep only unique values of species, sample and abundance
  distinct(Species, Sample, Abundance, .keep_all = TRUE) %>%
  select(Order, Sample, Species, Abundance, Depth, min, max) %>%
  arrange(Species) %>%
  # Add icons to some species
  # "Seen on video" (camera icon) and "better ID from eDNA" (plus symbol)
  mutate(Species = case_when( 
    # better than video (camera +)
    str_detect(Species, "Rhynchoconger flavus|Myrophis punctatus|Ophichthus gomesii|Selar crumenophthalmus|Seriola dumerili|Seriola rivoliana|Etrumeus teres|Malacanthus plumieri|Diaphus splendidus|Chromis enchrysurus|Chromis insolata|Pronotogrammus martinicensis|Trichopsetta ventralis|Citharichthys cornutus|Syacium gunteri|Symphurus civitatium|Stomias affinis|Sphoeroides dorsalis") ~ 
      paste0("<span style='font-family:\"Font Awesome 6 Free Solid\"'>&#xf083;</span>+ ", Species), 
    # seen on video only (camera)
    str_detect(Species, "Caranx crysos|Chloroscombrus chrysurus|Decapterus punctatus|Selene setapinnis|Trachurus lathami|Leiostomus xanthurus|Micropogonias undulatus|Pareques iwamotoi|Pristipomoides aquilonaris|Serranus atrobranchus|Priacanthus arenatus|Auxis rochei|Auxis thazard|Lagodon rhomboides|Pagrus pagrus|Fistularia petimba|Mullus auratus|Upeneus parvus|Sphoeroides spengleri") ~ 
      paste0("<span style='font-family:\"Font Awesome 6 Free Solid\"'>&#xf083;</span> ", Species),
    str_detect(Species, "Rhomboplites aurorubens|Balistes capriscus|Pomatomus saltatrix|Lutjanus campechanus") ~ 
      paste0("<span style='color:maroon'>**VU** </span><span style='font-family:\"Font Awesome 6 Free Solid\"'>&#xf083;</span>+ ", Species),
    .default = Species))

# Plot data
fish_plot_data %>%
  ggplot(aes(Depth, Species)) +
  # Add background bars to separate fish orders (set manual in separate file)
  geom_rect(data = background, aes(x = NULL, y = NULL, 
                                   xmin=0, xmax=Inf, ymin=ymin, ymax=ymax, fill=odd), 
            alpha = 0.2, show.legend = FALSE) +
  # Add depth ranges (blue bars)
  geom_linerange(data = distinct(fish_plot_data, Species, .keep_all = TRUE), 
                 aes(xmin = min, xmax = max), 
                 linetype = "solid", linewidth = 4, alpha = 0.6, color="#89B7CF") +
  # Add detection points
  geom_point(aes(x = Depth, y = Species, size = Abundance), color = "black", alpha = 0.5) +
  # Group fish species by order
  facet_grid(~ Order, scales = "free", space = "free") +
  # Flip the x/y coordinates so fish species labels are on bottom
  coord_flip() +
  # Reverse the depth scale so deeper depths are at the bottom
  scale_x_continuous(trans = "reverse", expand = c(0,0), limits = c(112, 0)) +
  # Set scale for detection points
  scale_size_continuous(range = c(2,6)) +
  # Set colors for background bars
  scale_fill_manual(values = c("#ffffff","gray80")) +
  # Set theming
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = "Arial", color = "black"),
    axis.text = element_markdown(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank(), 
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    panel.grid.major.y = element_line(linewidth = 0.5, colour = "grey80", linetype = "dotted"),
    strip.background = element_blank(), 
    strip.text.x = element_text(angle = 90),
    # Remove space between facet panels
    panel.spacing = unit(0, "lines"),
    legend.position = "bottom",
    legend.justification = "right",
    legend.title = element_text(size = 10, family = "Arial"),
    plot.margin = unit(c(0.25,0.25,0.25,1), "in"),
  ) +
  # Set x axis label
  xlab("Depth (m)") 

# Save plot
ggsave(
  here::here("plots/fish_depths_NEW_DATABASE.svg"), 
  plot = last_plot(), 
  width = 14, 
  height = 8,
  units = "in",
  dpi = 300
) 
