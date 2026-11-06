library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(here)

# 1. Load data
df1 <- read.csv(here("results", "tables", "section_3_4", "taxa_E1_core.csv"))
df4 <- read.csv(here("results", "tables", "section_3_4", "taxa_E4_core.csv"))

# 2. Extract taxonomy and calculate mean log2FC
process_df <- function(df, lfc_cols, type_name) {
  df %>%
    mutate(
      Phylum = str_extract(feature, "p__[^;]+") %>% str_replace("p__", ""),
      Family = str_extract(feature, "f__[^;]+") %>% str_replace("f__", ""),
      Species = str_extract(feature, "s__[^;]+") %>% str_replace("s__", "")
    ) %>%
    mutate(mean_lfc = rowMeans(select(., all_of(lfc_cols)))) %>%
    select(feature, Phylum, Family, Species, mean_lfc) %>%
    mutate(CoreType = type_name)
}

t1 <- process_df(df1, c("log2FoldChange_E1vE2", "log2FoldChange_E1vE3"), "Type 1")
t4 <- process_df(df4, c("log2FoldChange_E4vE2", "log2FoldChange_E4vE3"), "Type 4")

# Identify shared
shared_features <- intersect(t1$feature, t4$feature)

# Combine and assign membership stringently
# For the plot, we need a single bar value for each taxon.
# If it's shared, we can use the average or specifically Type 1's just for layout, 
# but since they are convergent (both upregulated), we'll just plot T1's and mark it Shared.
df_combined <- bind_rows(t1, t4) %>%
  mutate(Membership = case_when(
    feature %in% shared_features ~ "Shared",
    CoreType == "Type 1" ~ "Type 1 Exclusive",
    CoreType == "Type 4" ~ "Type 4 Exclusive"
  )) %>%
  # Remove the duplicate shared feature (keep T1's lfc for plotting)
  distinct(feature, .keep_all = TRUE)

# Clean up Species names
df_combined$Species <- str_replace_all(df_combined$Species, "_", " ")

# Group logically: First by Membership, then Phylum
df_combined <- df_combined %>%
  arrange(Membership, Phylum, mean_lfc)

# 3. Add Empty Bars for separation
# Calculate how many gaps we need. We'll gap by Phylum AND Membership to keep it clean.
df_combined$Group <- paste(df_combined$Membership, df_combined$Phylum, sep="_")

empty_bar <- 2
unique_groups <- unique(df_combined$Group)
to_add <- data.frame(matrix(NA, empty_bar * length(unique_groups), ncol(df_combined)))
colnames(to_add) <- colnames(df_combined)
to_add$Group <- rep(unique_groups, each = empty_bar)

df_plot <- rbind(df_combined, to_add) %>%
  arrange(Group, is.na(mean_lfc), mean_lfc) %>%
  mutate(id = row_number())

# 4. Prepare label angles and positions
label_data <- df_plot %>% filter(!is.na(mean_lfc))
number_of_bar <- nrow(df_plot)
angle <- 90 - 360 * (label_data$id - 0.5) / number_of_bar
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle + 180, angle)

ylim_min <- min(df_plot$mean_lfc, na.rm=TRUE) - 10
ylim_max <- max(df_plot$mean_lfc, na.rm=TRUE) + 8
label_data$y_pos <- pmax(label_data$mean_lfc, 0) + 0.5 

# Add outer marker positions
label_data$marker_y <- ylim_max - 2

# Colors
custom_colors <- c("#003f5c", "#bc5090", "#ffa600", "#58508d", "#ff6361", "#2f4b7c", "#a05195", "#d45087")
# Membership colors aligned with the global enterotype palette
# (Type 1 = red, Type 4 = blue; see scripts/utils/palettes.R). Shared = neutral grey.
membership_colors <- c("Type 1 Exclusive" = "#E64B35", "Type 4 Exclusive" = "#3C5488", "Shared" = "#999999")

# 5. Plotting
p <- ggplot(df_plot, aes(x = as.factor(id), y = mean_lfc)) +
  geom_hline(yintercept = 0, color = "gray30", linewidth=0.5, linetype="dashed") +
  # The actual bars colored by Phylum
  geom_bar(aes(fill = Phylum), stat = "identity", alpha = 0.85, width=0.8) +
  scale_fill_manual(values = custom_colors, name="Phylum", na.translate = FALSE) +
  # Add colored points for Membership on the outer ring
  geom_point(data=label_data, aes(x=id, y=marker_y, color=Membership), size=2, alpha=0.9) +
  scale_color_manual(values = membership_colors, name="Core Status") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size=8),
    legend.title = element_text(size=9, face="bold"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  coord_polar(start = 0) +
  geom_text(data=label_data, 
            aes(x=id, y=y_pos, label=Species, hjust=hjust), 
            color="black", fontface="italic", size=2.5, 
            angle=label_data$angle, inherit.aes=FALSE)
  
# Output goes to supplementary/figures/ as FigureS4A (consumed by assemble_figureS4.R)
ggsave(here("supplementary", "figures", "FigureS4A_Radial_Taxa.pdf"), p, width=12, height=10)
ggsave(here("supplementary", "figures", "FigureS4A_Radial_Taxa.png"), p, width=12, height=10, dpi=300)

cat("Successfully generated Panel B: Radial Bar Chart for Taxa with Membership.\n")
