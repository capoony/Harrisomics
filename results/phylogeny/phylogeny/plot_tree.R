# Load necessary R libraries
suppressMessages({
    library('ggtree')
    library('ggplot2')
    library('ape')
    if (!require('phytools', quietly = TRUE)) {
        # Use ape functions if phytools is not available
        nodeHeights <- function(tree) {
            return(max(node.depth.edgelength(tree)))
        }
    }
})

# Get tree file and outgroup from command line arguments
args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
outgroup <- args[2]
output_prefix <- args[3]

if (is.na(tree_file) || !file.exists(tree_file)) {
    stop("Tree file not found: ", tree_file)
}

# Load tree file
tree <- read.tree(tree_file)

# Calculate tree height (on x-axis)
Xmax <- max(node.depth.edgelength(tree))

# Root tree with outgroup if specified and present
if (!is.na(outgroup) && outgroup %in% tree$tip.label) {
    tree <- root(tree, outgroup = outgroup)
    cat("Tree rooted with outgroup:", outgroup, "\n")
} else {
    cat("Outgroup not found or not specified, using midpoint rooting\n")
    tree <- midpoint(tree)
}

# Plot tree
PLOT.tree <- ggtree(tree) +
    ggtitle('Partridge and Quail Phylogeny (BUSCO genes)') +
    theme_tree2() +
    theme_bw() +
    xlim(0, Xmax + 0.005) +
    xlab('avg. substitutions/site') +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    geom_tiplab(size = 3)

# Add bootstrap support if available
if ("node.label" %in% names(tree) && !all(is.na(tree$node.label))) {
    PLOT.tree <- PLOT.tree + geom_nodelab(size = 2)
}

# Export tree
ggsave(filename = paste0(output_prefix, '.pdf'), PLOT.tree, width = 12, height = 8)
ggsave(filename = paste0(output_prefix, '.png'), PLOT.tree, width = 12, height = 8, dpi = 300)

cat("Tree plots saved to:", paste0(output_prefix, '.pdf'), "and", paste0(output_prefix, '.png'), "\n")
