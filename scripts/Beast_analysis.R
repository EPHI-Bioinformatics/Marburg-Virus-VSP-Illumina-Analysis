# ============================================================================
# Phylogenetic Tree Visualization using BEAST Output
# ============================================================================
# Author: Maysa Niazy
# Version: 2.0
# ============================================================================

# Load required libraries with error checking
required_packages <- c("ape", "ggtree", "tidyverse", "phytools", "ggplot2", 
                       "ggnewscale", "RColorBrewer", "scico", "ggrepel")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

# Install ggtree from Bioconductor if needed
if (!require("ggtree", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("ggtree")
  library(ggtree)
}

# ============================================================================
# 1. CONFIGURATION SECTION 
# ============================================================================

# File paths (update these)
tree_file <- "mrv_mcc.tre"  # BEAST MCC tree file
metadata_file <- "all_seq_metadata.csv"      
output_dir <- "phylogenetic_results"              # Output directory

# Analysis parameters
node_support_threshold <- 0.7                     # Minimum node support to display
tip_label_size <- 2.5                             # Size of tip labels
point_size <- 3                                    # Size of tip points
clade_offset <- 0.05                              # Offset for clade labels
tree_width <- 12                                   # Plot width in inches
tree_height <- 10                                  # Plot height in inches
dpi_resolution <- 300                              # Output resolution

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# ============================================================================
# 2. READ AND VALIDATE DATA
# ============================================================================

# Function to read tree file with error handling
read_tree_file <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("Tree file not found: ", file_path)
  }
  
  # Try different formats
  tree <- tryCatch({
    if (grepl("\\.nexus$|\\.nex$|\\.nxs$", file_path, ignore.case = TRUE)) {
      read.nexus(file_path)
    } else if (grepl("\\.tre$|\\.tree$|\\.nwk$", file_path, ignore.case = TRUE)) {
      read.tree(file_path)
    } else {
      # Try both formats
      tryCatch(read.nexus(file_path), error = function(e) read.tree(file_path))
    }
  }, error = function(e) {
    stop("Failed to read tree file: ", e$message)
  })
  
  cat("Successfully loaded tree with", Ntip(tree), "tips and", Nnode(tree), "nodes\n")
  return(tree)
}

# Read tree and metadata
tree <- read_tree_file(tree_file)
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

# Validate metadata structure
required_columns <- c("tip_label", "country", "host")
if (!all(required_columns %in% colnames(metadata))) {
  stop("Metadata missing required columns: ", 
       paste(setdiff(required_columns, colnames(metadata)), collapse = ", "))
}

# ============================================================================
# 3. CONFIGURE COLOR PALETTES AND SHAPES
# ============================================================================

# Your custom color palette for countries
country_palette <- c(
  "Kenya" = "#FF7F00",
  "DRC" = "#4DAF4A",
  "Sierra Leone" = "#E78AC3",
  "Uganda" = "#377EB8",
  "Angola" = "#984EA3",
  "Netherlands" = "black",
  "Germany" = "#FFFF33",
  "Ethiopia" = "#E41A1C",
  "Guinea" = "cyan2",
  "Ghana" = "#999999",
  "Tanzania" = "maroon2",
  "Rwanda" = "#A6D854",
  "South Africa" = "maroon4"
)

# Define host shapes
host_shapes <- c("Human" = 16, "Bat" = 18)

# Define clade annotations (update node numbers based on your tree)
clade_annotations <- data.frame(
  node = c(150, 155, 160, 165, 170, 175),  # Replace with actual node numbers
  label = c("RAVV.1", "RAVV.2", "MARV.A.1", "MARV.A.3", "MARV.B", "MARV.B.2"),
  color = c("deepskyblue2", "brown4", "deepskyblue4", "darkslategray1", "coral2", "lightsalmon"),
  stringsAsFactors = FALSE
)

# ============================================================================
# 4. PREPARE METADATA FOR PLOTTING
# ============================================================================

# Ensure all tips in tree have metadata
tips_in_tree <- tree$tip.label
metadata_filtered <- metadata %>%
  filter(tip_label %in% tips_in_tree)

# Check for missing metadata
missing_tips <- setdiff(tips_in_tree, metadata_filtered$tip_label)
if (length(missing_tips) > 0) {
  warning("Missing metadata for ", length(missing_tips), " tips. Creating dummy entries.")
  dummy_metadata <- data.frame(
    tip_label = missing_tips,
    country = "Unknown",
    host = "Unknown",
    stringsAsFactors = FALSE
  )
  metadata_filtered <- bind_rows(metadata_filtered, dummy_metadata)
}

# Create tip data frame
tip_data <- data.frame(tip_label = tree$tip.label, stringsAsFactors = FALSE) %>%
  left_join(metadata_filtered, by = "tip_label") %>%
  mutate(
    country = ifelse(is.na(country), "Unknown", country),
    host = ifelse(is.na(host), "Unknown", host)
  )

# ============================================================================
# 5. CREATE BASE TREE PLOT
# ============================================================================

# Create base tree
p <- ggtree(tree, layout = "rectangular", size = 0.5, color = "gray30") +
  # Add scale bar
  geom_treescale(x = 0, y = -2, fontsize = 3, offset = 1) +
  # Add node support values
  geom_nodepoint(
    aes(subset = as.numeric(label) > node_support_threshold & 
          !is.na(as.numeric(label))),
    size = 2, color = "black", alpha = 0.6
  ) +
  geom_nodelab(
    aes(label = ifelse(as.numeric(label) > node_support_threshold, 
                       paste0(round(as.numeric(label), 3)), "")),
    size = 2.5, hjust = -0.2, color = "gray30"
  )

# ============================================================================
# 6. ADD METADATA ANNOTATIONS
# ============================================================================

p <- p %<+% tip_data +
  geom_tippoint(
    aes(color = country, shape = host),
    size = point_size,
    stroke = 0.5
  ) +
  scale_color_manual(
    name = "Geographic Origin",
    values = country_palette,
    na.value = "gray80"
  ) +
  scale_shape_manual(
    name = "Host Source",
    values = host_shapes,
    labels = c("Human" = "Human (Circle)", "Bat" = "Bat (Diamond)")
  )

# Optional: Add tip labels (commented by default for cleaner plots)
# p <- p + geom_tiplab(size = tip_label_size, offset = 0.01, align = FALSE, linesize = 0.2)

# ============================================================================
# 7. ADD CLADE ANNOTATIONS
# ============================================================================

# Function to safely add clade labels
add_clade_labels <- function(plot, clade_df, tree_obj, offset = clade_offset) {
  for (i in 1:nrow(clade_df)) {
    # Validate node exists
    if (clade_df$node[i] <= Ntip(tree_obj) + Nnode(tree_obj)) {
      plot <- plot + geom_cladelabel(
        node = clade_df$node[i],
        label = clade_df$label[i],
        color = clade_df$color[i],
        offset = offset,
        offset.text = offset * 0.4,
        barsize = 1.2,
        fontsize = 3.5,
        angle = 90,
        hjust = 0.5,
        alpha = 0.8
      )
    } else {
      warning("Node ", clade_df$node[i], " not found in tree. Skipping clade label.")
    }
  }
  return(plot)
}

# Add clade annotations if provided
if (exists("clade_annotations") && nrow(clade_annotations) > 0) {
  p <- add_clade_labels(p, clade_annotations, tree)
}

# ============================================================================
# 8. CUSTOMIZE THEME AND LAYOUT
# ============================================================================

p_final <- p +
  theme_tree2() +
  theme(
    # Legend positioning and styling
    legend.position = "right",
    legend.justification = "top",
    legend.box = "vertical",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.6, "cm"),
    legend.spacing.y = unit(0.15, "cm"),
    legend.background = element_rect(fill = "white", color = "gray80", size = 0.3),
    
    # Plot margins
    plot.margin = margin(t = 20, r = 60, b = 20, l = 20, unit = "pt"),
    
    # Background
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Phylogenetic Analysis of Marburg and Ravn Viruses",
    subtitle = "Maximum Clade Credibility Tree from BEAST Analysis",
    x = "Time (years)",
    y = "",
    caption = paste("Node labels show posterior support values (≥", 
                    node_support_threshold, " shown)", sep = "")
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray50", margin = margin(b = 10)),
    plot.caption = element_text(size = 8, hjust = 0, color = "gray50", margin = margin(t = 10))
  )

# Display the plot
print(p_final)

# ============================================================================
# 9. SAVE PLOT IN MULTIPLE FORMATS
# ============================================================================

# Function to save plots in multiple formats
save_plot_formats <- function(plot, base_filename, output_dir, width, height, dpi) {
  # PDF (vector format - best for publications)
  ggsave(
    filename = file.path(output_dir, paste0(base_filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    device = "pdf"
  )
  cat("Saved PDF:", base_filename, ".pdf\n")
  
  # PNG (high resolution raster)
  ggsave(
    filename = file.path(output_dir, paste0(base_filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi * 2,  # Double DPI for PNG
    device = "png"
  )
  cat("Saved PNG:", base_filename, ".png\n")
  
  # TIFF (for journals)
  ggsave(
    filename = file.path(output_dir, paste0(base_filename, ".tiff")),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    compression = "lzw",
    device = "tiff"
  )
  cat("Saved TIFF:", base_filename, ".tiff\n")
  
  # SVG (for web/editing)
  if (require("svglite", quietly = TRUE)) {
    svglite::svglite(
      file.path(output_dir, paste0(base_filename, ".svg")),
      width = width,
      height = height
    )
    print(plot)
    dev.off()
    cat("Saved SVG:", base_filename, ".svg\n")
  }
}

# Save main plot
save_plot_formats(p_final, "MARV_RAVV_phylogenetic_tree", 
                  output_dir, tree_width, tree_height, dpi_resolution)

