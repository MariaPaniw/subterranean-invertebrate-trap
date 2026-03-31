###Subterranean traps descriptive analysis ###

### COMPLETE QUALITATIVE ANALYSIS ###

# 1. INSTALL AND LOAD REQUIRED PACKAGES -----------------------------------
required_packages <- c(
  "tidyverse", "lubridate", "ggpubr", "RColorBrewer", "ggvenn", "patchwork",
  "networkD3", "flextable", "officer", "ggiraph", "scales", "viridis", "glue"
)

# Install missing packages
new_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
invisible(lapply(required_packages, require, character.only = TRUE))

# 2. DATA PREPARATION -----------------------------------------------------
# Read and clean both datasets
subterranean <- read.csv("TableS1_forAnalysis.csv") %>%
  mutate(
    Timestamp = as.Date(Timestamp),
    Month = month(Timestamp, label = TRUE, abbr = FALSE),
    Year = year(Timestamp),
    Method = "Subterranean"
  ) %>%
  filter(Year >= 2024)  # Filter for 2024-present data

# Read and prepare pitfall data
pitfall <- read.csv("Pitfall_data2.csv") %>%
  mutate(
    Timestamp = as.Date(Timestamp),
    Month = month(Timestamp, label = TRUE, abbr = FALSE),
    Year = year(Timestamp),
    Method = "Pitfall"
  ) %>%
  filter(Year >= 2024)  # Filter for 2024-present data

# Combine and clean data
combined_clean <- bind_rows(subterranean, pitfall) %>%
  mutate(
    Family = ifelse(is.na(Family), "Unidentified", Family),
    Genera = ifelse(is.na(Genera), "Unidentified", Genera),
    TaxonID = paste(Family, Genera, sep = ": "),
    LifeStage = case_when(
      is.na(LifeStage) ~ "Unidentified",
      LifeStage == "" ~ "Unidentified",
      TRUE ~ LifeStage
    )
  ) %>%
  filter(Family != "Unidentified") %>%
  mutate(Presence = ifelse(Number_caught > 0, 1, 0))

# 3. QUALITATIVE STATISTICS TABLE ----------------------------------------
top_taxa <- combined_clean %>%
  group_by(Method, TaxonID) %>%
  summarise(Total = sum(Number_caught), .groups = "drop") %>%
  group_by(TaxonID) %>% 
  slice_max(Total, n = 10)
view(top_taxa)

qual_table <- combined_clean %>%
  filter(TaxonID %in% top_taxa) %>%
  group_by(Method, TaxonID, LifeStage) %>%
  summarise(Count = sum(Number_caught), .groups = "drop") %>%
  group_by(Method) %>%
  mutate(Method_Percent = round(Count/sum(Count)*100, 1)) %>%
  group_by(TaxonID) %>%
  mutate(Overall_Percent = round(Count/sum(Count)*100, 1)) %>%
  arrange(Method, desc(Count))

# 3. Summarise and calculate proportions
taxon_table <- combined_clean %>%
  group_by(
    Method,
    Class,
    Order,
    Family,
    Subfamily,
    Tribe,
    Genera,
    LifeStage
  ) %>%
  summarise(n = sum(Number_caught), .groups = "drop") %>%
  group_by(Method) %>%
  mutate(Proportion = round(n / sum(n), 2)) %>%
  ungroup() %>%
  # Rename columns for final output
  rename(
    Genus = Genera,
    `Life stage` = LifeStage
  ) %>%
  select(Class, Order, Family, Subfamily, Tribe, Genus, `Life stage`, Method, n, Proportion)

library(DT)
# 4. Show as interactive DT table
datatable(
  taxon_table,
  extensions = "Buttons",
  options = list(
    pageLength = 25,
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'excel'),
    scrollX = TRUE
  ),
  rownames = FALSE
)

# 4. VISUALIZATIONS ------------------------------------------------------

## 4a. Enhanced Venn Diagram
# Create presence/absence matrix for Venn diagram
presence_absence <- combined_clean %>%
  group_by(Method, TaxonID) %>%
  summarise(Present = as.numeric(sum(Number_caught) > 0)) %>%
  pivot_wider(names_from = Method, values_from = Present, values_fill = 0)

# Venn diagram of taxa overlap
venn_data <- list(
  Subterranean = presence_absence %>% filter(Subterranean == 1) %>% pull(TaxonID),
  Pitfall = presence_absence %>% filter(Pitfall == 1) %>% pull(TaxonID)
)

venn_plot <- ggvenn(venn_data, 
       fill_color = c("darkblue", "#C45508"),
      stroke_size = 0.5,  
      set_name_size = 6,  
      text_size = 7,
      text_color = "white"
) + 
  ggtitle("Taxonomic overlap between both sampling methods") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

venn_plot

ggsave(filename = "C:/Users/Walt/Nextcloud/anything like this is private/Walt/BUGS/Paper/Images&Tables/venn_plot.png",
       plot = venn_plot,
       width = 8, height = 6, units = "in", dpi = 300)


#Creating a static sankey for Taxon difference
# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
# Install and load ggsankey if not already installed:
install.packages("remotes"); remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)

# ––– 1. Identify top‐8 in each category –––
top8_subterranean <- combined_clean %>%
  filter(Method == "Subterranean") %>%
  count(TaxonID, wt = Number_caught, name = "n") %>%
  slice_max(n, n = 8) %>%
  pull(TaxonID)

top8_pitfall <- combined_clean %>%
  filter(Method == "Pitfall") %>%
  count(TaxonID, wt = Number_caught, name = "n") %>%
  slice_max(n, n = 8) %>%
  pull(TaxonID)

# “Shared” = taxa caught by *both* methods; then pick the top‐6 of those by total captures
shared_candidates <- intersect(
  combined_clean$TaxonID[combined_clean$Method=="Subterranean"],
  combined_clean$TaxonID[combined_clean$Method=="Pitfall"]
)

top6_shared <- combined_clean %>%
  filter(TaxonID %in% shared_candidates) %>%
  count(TaxonID, wt = Number_caught, name = "n") %>%
  slice_max(n, n = 6) %>%
  pull(TaxonID)

# ––– 2. Union and filter –––
top_taxa <- unique(c(top8_subterranean, top8_pitfall, top6_shared))

top_data <- combined_clean %>%
  filter(TaxonID %in% top_taxa)

# 3. Calculate total captures per method (for percentage calculations)
method_totals <- combined_clean %>%
  group_by(Method) %>%
  summarise(total_method = sum(Number_caught), .groups = "drop")
total_pitfall    <- method_totals$total_method[method_totals$Method == "Pitfall"]
total_subterrane <- method_totals$total_method[method_totals$Method == "Subterranean"]

# ––– 3. Then rerun your summary + sankey pipeline on `top_data` –––
taxon_summary <- top_data %>%
  group_by(TaxonID, Method) %>%
  summarise(count = sum(Number_caught), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Method, values_from = count, values_fill = 0)

# 5. Compute percentages for each taxon's contributions and create labels
taxon_summary <- taxon_summary %>%
  mutate(pit_pct = Pitfall / total_pitfall * 100,
         sub_pct = Subterranean / total_subterrane * 100,
         label = case_when(
           Subterranean == 0 ~ paste0(TaxonID, "\n[Pitfall]: n=", Pitfall, ", ", round(pit_pct, 1), "%"),
           Pitfall == 0      ~ paste0(TaxonID, "\n[Subterranean]: n=", Subterranean, ", ", round(sub_pct, 1), "%"),
           TRUE              ~ paste0(TaxonID, "\nPitfall: n=", Pitfall, ", ", round(pit_pct, 1), "% ",
                                      "Subterranean: n=", Subterranean, ", ", round(sub_pct, 1), "%")
         ),
         category = case_when(
           Subterranean == 0 ~ "Pitfall-only",
           Pitfall == 0      ~ "Subterranean-only",
           TRUE              ~ "Shared"
         ))

# 6. Prepare flows data for Sankey (one row per Method→Taxon link with its count)
flows_df <- taxon_summary %>%
  select(TaxonID, Pitfall, Subterranean) %>%
  tidyr::pivot_longer(cols = c(Pitfall, Subterranean), names_to = "Method", values_to = "Count") %>%
  filter(Count > 0)   # remove zero-count links (taxon not caught by that method)

# 7. Transform to long format required by ggsankey (Method as first stage, TaxonID as second)
long_df <- flows_df %>% 
  make_long(Method, TaxonID, value = Count)

# Replace TaxonID nodes with the detailed label text
taxon_label_map <- setNames(taxon_summary$label, taxon_summary$TaxonID)
long_df$node      <- ifelse(long_df$node %in% names(taxon_label_map),
                            taxon_label_map[long_df$node], long_df$node)
long_df$next_node <- ifelse(long_df$next_node %in% names(taxon_label_map),
                            taxon_label_map[long_df$next_node], long_df$next_node)

# 8. Define colors for each node (blue for Subterranean-only, orange for Pitfall-only, green for shared)
color_map <- c(
  # Taxon nodes by category:
  setNames(rep("#FBC687", sum(taxon_summary$category == "Pitfall-only")),      # orange
           taxon_summary$label[taxon_summary$category == "Pitfall-only"]),
  setNames(rep("#A3C1E0", sum(taxon_summary$category == "Subterranean-only")), # blue
           taxon_summary$label[taxon_summary$category == "Subterranean-only"]),
  setNames(rep("#A6D9CE", sum(taxon_summary$category == "Shared")),            # green
           taxon_summary$label[taxon_summary$category == "Shared"]),
  # Method nodes (match their respective colors):
  "Pitfall"      = "#FBC687",  # rust red
  "Subterranean" = "#A3C1E0"   # navy blue
)


# 9. Create the Sankey diagram
sankey_plot <- ggplot(long_df, aes(x = x, next_x = next_x, 
                                   node = node, next_node = next_node,
                                   value = value^0.02, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.7, node.alpha = 0.9, width = 0.22, node.color = "black") +              # Sankey flows and nodes
  geom_sankey_label(aes(fill = node), size = 4.5, label.size = 0, colour = "black") +  # Node labels with white background
  scale_fill_manual(values = color_map) +                                    # Apply manual colors
  theme_sankey(base_size = 14) +                                             # Clean Sankey theme for elegant look
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20)) +
  labs(x = NULL)# No legend (nodes are directly labeled)

# Display the plot
print(sankey_plot)

ggsave(filename = "C:/Users/Walt/Nextcloud/anything like this is private/Walt/BUGS/Paper/Images&Tables/sankey_plot.png",
       plot = sankey_plot,
       width = 13, height = 13, units = "in", dpi = 500)

#Now a sankey using Family and life stage

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
# Install and load ggsankey if not already installed:
install.packages("remotes"); remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)

# Combine and clean data - changed TaxonID 
combined_clean_life <- bind_rows(subterranean, pitfall) %>%
  mutate(
    Family = ifelse(is.na(Family), "Unidentified", Family),
    Genera = ifelse(is.na(Genera), "Unidentified", Genera),
    TaxonID = paste(Family, LifeStage, sep = ": "),
    LifeStage = case_when(
      is.na(LifeStage) ~ "Unidentified",
      LifeStage == "" ~ "Unidentified",
      TRUE ~ LifeStage
    )
  ) %>%
  filter(Family != "Unidentified") %>%
  mutate(Presence = ifelse(Number_caught > 0, 1, 0))

# ––– 1. Identify top‐8 in each category –––
top8_subterranean <- combined_clean_life %>%
  filter(Method == "Subterranean") %>%
  count(TaxonID, wt = Number_caught, name = "n") %>%
  slice_max(n, n = 8) %>%
  pull(TaxonID)

top8_pitfall <- combined_clean_life %>%
  filter(Method == "Pitfall") %>%
  count(TaxonID, wt = Number_caught, name = "n") %>%
  slice_max(n, n = 8) %>%
  pull(TaxonID)

# “Shared” = taxa caught by *both* methods; then pick the top‐6 of those by total captures
shared_candidates <- intersect(
  combined_clean_life$TaxonID[combined_clean_life$Method=="Subterranean"],
  combined_clean_life$TaxonID[combined_clean_life$Method=="Pitfall"]
)

top6_shared <- combined_clean_life %>%
  filter(TaxonID %in% shared_candidates) %>%
  count(TaxonID, wt = Number_caught, name = "n") %>%
  slice_max(n, n = 6) %>%
  pull(TaxonID)

# ––– 2. Union and filter –––
top_taxa <- unique(c(top8_subterranean, top8_pitfall, top6_shared))

top_data <- combined_clean_life %>%
  filter(TaxonID %in% top_taxa)

# ––– 3. Then rerun your summary + sankey pipeline on `top_data` –––
taxon_summary <- top_data %>%
  group_by(TaxonID, Method) %>%
  summarise(count = sum(Number_caught), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Method, values_from = count, values_fill = 0)

# 5. Compute percentages for each taxon's contributions and create labels
taxon_summary <- taxon_summary %>%
  mutate(pit_pct = Pitfall / total_pitfall * 100,
         sub_pct = Subterranean / total_subterrane * 100,
         label = case_when(
           Subterranean == 0 ~ paste0(TaxonID, "\n[Pitfall]: n=", Pitfall, ", ", round(pit_pct, 1), "%"),
           Pitfall == 0      ~ paste0(TaxonID, "\n[Subterranean]: n=", Subterranean, ", ", round(sub_pct, 1), "%"),
           TRUE              ~ paste0(TaxonID, "\nPitfall: n=", Pitfall, ", ", round(pit_pct, 1), "% ",
                                      "Subterranean: n=", Subterranean, ", ", round(sub_pct, 1), "%")
         ),
         category = case_when(
           Subterranean == 0 ~ "Pitfall-only",
           Pitfall == 0      ~ "Subterranean-only",
           TRUE              ~ "Shared"
         ))

# 6. Prepare flows data for Sankey (one row per Method→Taxon link with its count)
flows_df <- taxon_summary %>%
  select(TaxonID, Pitfall, Subterranean) %>%
  tidyr::pivot_longer(cols = c(Pitfall, Subterranean), names_to = "Method", values_to = "Count") %>%
  filter(Count > 0)   # remove zero-count links (taxon not caught by that method)

# 7. Transform to long format required by ggsankey (Method as first stage, TaxonID as second)
long_df <- flows_df %>% 
  make_long(Method, TaxonID, value = Count)

# Replace TaxonID nodes with the detailed label text
taxon_label_map <- setNames(taxon_summary$label, taxon_summary$TaxonID)
long_df$node      <- ifelse(long_df$node %in% names(taxon_label_map),
                            taxon_label_map[long_df$node], long_df$node)
long_df$next_node <- ifelse(long_df$next_node %in% names(taxon_label_map),
                            taxon_label_map[long_df$next_node], long_df$next_node)

# 8. Define colors for each node (blue for Subterranean-only, orange for Pitfall-only, green for shared)
color_map <- c(
  # Taxon nodes by category:
  setNames(rep("#FBC687", sum(taxon_summary$category == "Pitfall-only")),      # orange
           taxon_summary$label[taxon_summary$category == "Pitfall-only"]),
  setNames(rep("#A3C1E0", sum(taxon_summary$category == "Subterranean-only")), # blue
           taxon_summary$label[taxon_summary$category == "Subterranean-only"]),
  setNames(rep("#A6D9CE", sum(taxon_summary$category == "Shared")),            # green
           taxon_summary$label[taxon_summary$category == "Shared"]),
  # Method nodes (match their respective colors):
  "Pitfall"      = "#FBC687",  # rust red
  "Subterranean" = "#A3C1E0"   # navy blue
)

# 9. Create the Sankey diagram
sankeylife_plot <- ggplot(long_df, aes(x = x, next_x = next_x, 
                                   node = node, next_node = next_node,
                                   value = value^0.12, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.7, node.alpha = 0.9, width = 0.22, node.color = "black") +              # Sankey flows and nodes
  geom_sankey_label(aes(fill = node), size = 4.5, label.size = 0, colour = "black") +  # Node labels with white background
  scale_fill_manual(values = color_map) +                                    # Apply manual colors
  theme_sankey(base_size = 14) +                                             # Clean Sankey theme for elegant look
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20)) +
  labs(x = NULL)# No legend (nodes are directly labeled)

# Display the plot
print(sankeylife_plot)

ggsave(filename = "C:/Users/Walt/Nextcloud/anything like this is private/Walt/BUGS/Paper/Images&Tables/sankeylife_plot2.png",
       plot = sankeylife_plot,
       width = 13, height = 10, units = "in", dpi = 400)
