
# This is to plot the maximum likelihood biogeographic states on the trees 
# as yielded by BioGeoBEARS analyses.

setwd()

# Load packages ----
libs <- c("treeio", "phytools", "geiger", "tidytree", "dplyr", "tidyverse", "doParallel",
          "RColorBrewer", "ggtree", "scico", "MCMCtreeR", "viridis")
lapply(libs, require, character.only = TRUE)

# Import list of clade trees and remove Lytorhynchus ----
consensus_genera <- readRDS("objects/consensus_genera.rds")


# Import list of tip and node states for each clade and remove Lytorhynchus ----
consensus_states <- readRDS("objects/consensus_states.rds")

# Define clade names 
clade_names <- names(consensus_states)

# Set the colors of each biogeographic state ----
state_colors <- c('F' = "#652581", R = "#efb004",
                  S = "#2b8258", E = "#ff5959", 
                  U = "#8CD9FF", FR = "#ff937c", 
                  RE = "#ff8433", RS = '#7ea233',
                  SE = '#978215', FS = '#00629c', 
                  FE = '#be0073', FU = '#317dd4', 
                  RU = '#4ddca3', SU = '#2fb0b3', 
                  EU = '#d39bf2')
plot(1:15, 1:15, col=state_colors, cex=5, pch=16)


# Set parameters
branch_color <- "black"
branch_size <- 0.4
max_lim <- 30

# Loop for plotting all clades ----
pdf("biogeo_plots_ggtree/biogeo_results_bestmodel.pdf", paper = "a4")
for (i in 1:length(consensus_genera)){
  if (length(consensus_genera[[i]]@phylo$tip.label) >= 3){
    tree <- consensus_genera[[i]]
    states <- consensus_states[[i]]
    states_df <- data.frame(node = c(1:length(states)), biogeo = states)
    
    # Create ggtree object: tree + states data frame
    t <- ggtree(tree, layout="rectangular", color=branch_color, size=branch_size) %<+%
      states_df
    
    tree_plot <- revts(t +
                         geom_point(aes(color=biogeo), size=5, alpha=1, shape=15) +
                         scale_color_manual(values=state_colors,
                                            name = 'biogeographic state',
                                            breaks = c('F', 'R', 'FR', 'E', 'S', 'U', 'RE', 'RS', 'SE', 'FS', 'FE', 'FU'),
                                            labels = c('F (Africa)', 'R (Arabia)', 'FR', 'E (Eastern Mediterranean)', 
                                                       'S (Asia/Oceania)', 'U (Europe)', 'RE', 'RS', 'SE',
                                                       'FS', 'FE', 'FU')) +
                         geom_text(aes(label=biogeo), size=3) +
                         geom_tiplab(size=3, offset=1, fontface="italic") +
                         theme_tree2(legend.position='right') +
                         xlab("Ma") +
                         scale_x_continuous(breaks=seq(-100, 0, 10), labels=abs(seq(-100, 0, 10)),
                                            limits = c(-tree[rootnode(tree@phylo),]$height-10, max_lim)) +
                         ggtitle(names(consensus_genera)[i]) +
                         theme(plot.title = element_text(face = "italic"))
    )
    
    # Print the plot
    print(tree_plot)
    
  }
    
}
dev.off()


# Plot node numbers ----
pdf("biogeo_plots_ggtree/node_numbers.pdf", paper = "a4")
for (i in 1:length(consensus_genera)){
  if (length(consensus_genera[[i]]@phylo$tip.label) >= 3){
    tree <- consensus_genera[[i]]
    states <- consensus_states[[i]]
    states_df <- data.frame(node = c(1:length(states)), biogeo = states)
    
    # Create ggtree object: tree + states data frame
    t <- ggtree(tree, layout="rectangular", color=branch_color, size=branch_size) %<+%
      states_df
    
    tree_plot <- revts(t +
                         geom_point(aes(color=biogeo), size=5, alpha=1, shape=15) +
                         scale_color_manual(values=state_colors,
                                            name = 'biogeographic state',
                                            breaks = c('F', 'R', 'FR', 'E', 'S', 'U', 'RE', 'RS', 'SE', 'FS', 'FE', 'FU'),
                                            labels = c('F (Africa)', 'R (Arabia)', 'FR', 'E (Eastern Mediterranean)', 
                                                       'S (Asia/Oceania)', 'U (Europe)', 'RE', 'RS', 'SE',
                                                       'FS', 'FE', 'FU')) +
                         geom_text(aes(label=node), size=3) +
                         geom_tiplab(size=3, offset=1, fontface="italic") +
                         theme_tree2(legend.position='right') +
                         xlab("Ma") +
                         scale_x_continuous(breaks=seq(-100, 0, 10), labels=abs(seq(-100, 0, 10)),
                                            limits = c(-tree[rootnode(tree@phylo),]$height-10, max_lim)) +
                         ggtitle(names(consensus_genera)[i]) +
                         theme(plot.title = element_text(face = "italic"))
    )
    
    # Print the plot
    print(tree_plot)
    
  }
   
}
dev.off()



