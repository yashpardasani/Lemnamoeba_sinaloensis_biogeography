# Loading the libraries
library(readr)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

## File 3 Eukbank_18S_V4_samples.tsv.gz: Contextual data table describing the samples
#look for sample here and check the information: then go to project 
# Specify the path to your gzipped TSV file
eukbank_file3 <- "/path/eukbank_18S_V4_samples.tsv.gz"
# Read the gzipped TSV file
file3_data <- read.delim(eukbank_file3, header = TRUE, sep = "\t")

## File 4 ukbank_18S_V4_counts.tsv.gz: Number of read sassociated to each ASV in Eukbank sample
#linking file, check for sample, amplicon are asvs
eukbank_file4 <- "path/eukbank_18S_V4_counts.tsv.gz"
file4_data <- read.delim(eukbank_file4, header = TRUE, sep = "\t")

#Checking the total number of samples
length(unique(file4_data$sample))

####Making a tsv file with only the amplicons of interest, i.e., only the amplicons that were closely related to Lemnamoeba in the 18S tree

#Give a text file with the names of amplicons/asvs
M6MM_in_eukbank <- "~/path/M6MM_clade_in_tree.txt"

amplicons_of_interest <- tolower(readLines(M6MM_in_eukbank))

#just to check how many are there, and we have everything
head(amplicons_of_interest)

#file4 data: the one that we loaded earlier, this step is to filter that and only keep just Lemnamoeba related amplicons
filtered_data <- file4_data[file4_data$amplicon %in% tolower(amplicons_of_interest), ]
### Now I am adding sample information to the above filtered_data and making a bigger tsv file, merging filtered data, i.e., file 4 data with file3 data
merged_metadata <- merge(filtered_data, file3_data, by = "sample", all.x = TRUE)

#checking the number of samples in which lemnamoeba-related ASVs are found
length(unique(merged_metadata$sample))

#Loading the world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Removing anything with no value for latitude and longitude
asv_df <- merged_metadata %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  mutate(amplicon = as.factor(amplicon))

##Adding manual colours, ASV names are manually changed later, check Supplementary table for details
custom_colors <- c(
  # Closely related clade (3 very distinct greens)
  "67979c40272eac9966eb759e6673671b12cf2d27" = "#00441b",  # dark forest
  "f57714ffce3df4c25df0cc4045d71086a82eb5c4" = "#78c679",  # bright lime
  "2ca3e2d13a8691c548b20a86d270d2df9807f9f4" = "#1c9099",  # teal-green
  
  # Distant ASVs
  "38ee549aeb136514faa8066a30dc9eaaa5d53b71" = "#d95f02",  # orange
  "37b8a3088b48540df95c352ee17dc9d87a4373c7" = "#e7298a"   # magenta  
)

#Adding relative abundance
#I made a new row in asv_df: nreads.x is the number of reads of that particular asv in that particular sample and nreads.y is the total number of reads of all asvs in that sample

asv_df <- merged_metadata |>
  filter(!is.na(latitude), !is.na(longitude)) |>
  mutate(
    amplicon = as.factor(amplicon),
    envplot  = as.factor(envplot),
    rel_abundance = nreads.x / nreads.y         # proportion
  )

asv_df <- asv_df |>
  mutate(percent_abundance = 100 * rel_abundance)

###############Plotting according to co-ordinates ################

asv_max <- asv_df %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  group_by(amplicon, longitude, latitude, envplot) %>%
  summarise(
    percent_abundance = max(percent_abundance, na.rm = TRUE),
    .groups = "drop"
  )

#plotting with respect to co-ordinates, only the max 
asv_distribution_plot <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70") +
  geom_point(
    data = asv_max,
    aes(
      x     = longitude,
      y     = latitude,
      size  = percent_abundance,
      color = amplicon,
      shape = envplot
    ),
    alpha = 0.8
  ) +
  scale_size_continuous(
    name   = "Max percent abundance (%)",
    range  = c(1.5, 8),
    trans  = "sqrt",
    breaks = c(0.1, 0.5, 1, 2, 3)
  ) +
  scale_color_manual(values = custom_colors, name = "Amplicon (ASV)") +
  scale_shape_manual(
    values = c("marine_sediment" = 16, "marine_water" = 17),
    name   = "Environment"
  ) +
  coord_sf() +
  theme_minimal() +
  labs(
    x = "longitude",
    y = "latitude",
    title = "Global distribution of ASVs\n(size = max relative abundance per site)"
  )

ggsave(
  filename = "/path/asv_distribution.pdf",
  plot     = asv_distribution_plot,
  device   = pdf,       
  width    = 6,              
  height   = 3,
  dpi      = 300
)

########################### Depth analysis #####################################

asv_depth <- asv_df %>%
  filter(!is.na(depth)) %>%
  mutate(
    amplicon = factor(
      amplicon,
      levels = c(
        "67979c40272eac9966eb759e6673671b12cf2d27",
        "f57714ffce3df4c25df0cc4045d71086a82eb5c4",
        "2ca3e2d13a8691c548b20a86d270d2df9807f9f4",
        "38ee549aeb136514faa8066a30dc9eaaa5d53b71",
        "37b8a3088b48540df95c352ee17dc9d87a4373c7"
      )
    )
  )

## Stretching y co-ordinate
stretch_depth <- function(d) {
  ifelse(
    d <= 500,
    d * 3,                    # shallow depths stretched
    1500 + (d - 500) / 2      # deep depths compressed
  )
}

#
asv_depth_disp <- asv_depth %>%
  mutate(depth_display = stretch_depth(depth))

# Real depths we want tick marks for:
breaks_real    <- c(0, 100, 200, 300, 400, 500, 1000, 1500, 2000)
breaks_display <- stretch_depth(breaks_real)

#plot
depth_plot_v2 <- ggplot(asv_depth_disp, aes(x = amplicon, y = depth_display)) +
  geom_point(
    aes(
      size  = percent_abundance,
      color = amplicon,
      shape = envplot
    ),
    position = position_jitter(width = 0.12, height = 0),
    alpha = 0.85
  ) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(
    values = c(
      "marine_sediment" = 16,
      "marine_water"    = 17
    )
  ) +
  scale_size_continuous(
    name   = "Percent abundance (%)",
    range  = c(1.5, 8),
    trans  = "sqrt",
    breaks = c(0.1, 0.5, 1, 2, 3)
  ) +
  scale_y_reverse(
    limits = range(breaks_display),
    breaks = breaks_display,
    labels = breaks_real,    # show real depths on the axis
    name   = "Depth (m)"
  ) +
  labs(
    x     = "Amplicon (ASV)",
    color = "Amplicon (ASV)",
    shape = "Environment",
    title = "Depth distribution of ASVs\n(size = relative abundance)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right"
  )

ggsave(
  filename = "/path/depth_plot.pdf",
  plot     = depth_plot_v2,
  device   = pdf,       
  width    = 6,            
  height   = 3,
  dpi      = 300
)
