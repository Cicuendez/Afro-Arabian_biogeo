
# DESCRIPTION ----
## Plot Afro-Arabia tree with tip labels indicating area and with 
## groups of interest highlited


# Load packages ----
libs <- c("treeio", "phytools", "geiger", "tidytree", "dplyr", "tidyverse", "doParallel",
          "RColorBrewer", "ggtree", "scico", "viridis")
lapply(libs, require, character.only = TRUE)


# Working directory
setwd()

# Import phylogeny
tree <- read.beast(file="data/AfAr_tree.nex")
tree_phylo <- tree@phylo
tip_labels <- tree_phylo$tip.label

# Import species-distribution information
tiplabels_biogeo <- read.table("data/tiplabels_biogeo.csv", sep=";", header = TRUE, 
                               stringsAsFactors = FALSE)

tiplabels_biogeo_0 <- tiplabels_biogeo
tiplabels_biogeo <- tiplabels_biogeo_0 %>%
  mutate(biogeo = case_when(in_out == 'Outgroup' ~ "Out", T ~ biogeo))
tiplabels_biogeo <- tiplabels_biogeo %>%
  mutate(biogeo = case_when(biogeo == 'Asia/Oceania' ~ 'Asia_Oceania', T ~ biogeo)) %>%
  mutate(biogeo = case_when(biogeo == 'Asia' ~ 'Asia_Oceania', T ~ biogeo))

# Number of ingroups and outgroups ----
outgroups <- tiplabels_biogeo %>%
  filter(in_out == "Outgroup")
nrow(outgroups)

ingroups <- tiplabels_biogeo %>%
  filter(in_out == "Ingroup")
nrow(ingroups)


biogeo_table <- tiplabels_biogeo %>%
  separate(biogeo, c("biogeo1", "biogeo2"), sep = "/", remove = TRUE)

biogeo_table <- biogeo_table %>%
  mutate(biogeo2 = case_when(is.na(biogeo2) ~ 'Out', T ~ biogeo2)) %>%
  mutate(biogeo2 = case_when(biogeo2 == 'Asia' ~ 'Asia_Oceania', T ~ biogeo2)) %>%
  mutate(biogeo1 = as.factor(biogeo1)) %>%
  mutate(biogeo2 = as.factor(biogeo2))

levels(biogeo_table$biogeo1)
levels(biogeo_table$biogeo2)

# Set colors
colors_biogeo <- c(Africa = "#652581", Arabia = "#efb004", 
                   Asia_Oceania ="#2b8258", EasternMed = "#ff5959", 
                   Europe = "#8CD9FF", Out = NA)
plot(c(1:6), c(1:6), col = colors_biogeo, cex = 4, pch = 16)

## Trying ggtree: node and tip number, biogeo tip colors, highlited clade, tiplabels
tree@phylo$tip.label <- gsub("_", " ", tree@phylo$tip.label)
biogeo_table$species <- gsub("_", " ", biogeo_table$species)


angle <- 20
t <- ggtree(tree, layout="fan", open.angle=angle) %<+%
  biogeo_table

# plot node numbers
t + 
  geom_tippoint(aes(color=biogeo1), size=0.5, alpha=1, shape=16) +
  geom_tippoint(aes(x = x+2.2, color=biogeo2), size = 0.5, alpha=1, shape=16) +
  scale_color_manual(values=colors_biogeo) +
  geom_nodepoint(color="lightblue", size=0.5) +
  geom_text2(aes(label=node), size=0.5) +
  geom_tiplab2(size=0.5, offset=2.5) +
  #  scale_color_brewer("Biogeo", palette="Spectral") +
  theme_tree2(legend.position='right') +
  ggsave("plots/tree/ggtree0.pdf", width=8.27, height=11.69)



## Node numbers of the clades ----
node_pristurus <- 1074
node_ptyodactylus <- 960
node_tropiocolotes <- 996
node_stenodactylus <- 984
node_hemidactylus <- 1008
node_chalcides_scincus <- 959
node_chalcides <- 920
node_scincus <- 942
node_mesalina <- 887
node_acanthodactylus <- 851
node_echis <- 738
node_cerastes <- 735
node_bitis <- 721
node_telescopus <- 826
node_malpolon <- 821
node_psammophis <- 794
node_atractaspis <- 752
node_naja <- 766
node_varanus <- 570
node_chamaeleo <- 681
node_uromastyx <- 667
node_pseudotrapelus <- 641

# Add clade information to the tree ----
clade_list <- split(tiplabels_biogeo$species, tiplabels_biogeo$Clade)
names(clade_list)
clades <- c(node_acanthodactylus, node_atractaspis, node_bitis, node_cerastes, 
               node_chalcides, node_chamaeleo, node_echis, node_malpolon, node_mesalina,
               node_naja, node_pristurus, node_psammophis, node_pseudotrapelus, node_ptyodactylus,
               node_scincus, node_stenodactylus, node_telescopus, node_tropiocolotes, 
               node_uromastyx, node_varanus, node_hemidactylus)
names(clades) <- names(clade_list)[-which(names(clade_list) == "No")]

t_clades <- groupClade(tree, sort(clades), group_name="monophyletic")

## Set branch color and branch width ----
branch_color <- "gray40"
branch_size <- 0.3
tree_clades <- ggtree(t_clades, layout="fan", open.angle = angle, color=branch_color, size=branch_size) %<+%
  biogeo_table

## Plot tree with node and tip numbers, biogeo colors, clades ----
t1 <- tree_clades +
  geom_cladelabel(node=clades[1], label=names(clades)[1], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[2], label=names(clades)[2], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[3], label=names(clades)[3], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[4], label=names(clades)[4], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[5], label=names(clades)[5], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[6], label=names(clades)[6], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[7], label=names(clades)[7], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[8], label=names(clades)[8], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[9], label=names(clades)[9], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[10], label=names(clades)[10], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[11], label=names(clades)[11], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[12], label=names(clades)[12], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[13], label=names(clades)[13], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[14], label=names(clades)[14], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[15], label=names(clades)[15], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[16], label=names(clades)[16], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[17], label=names(clades)[17], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[18], label=names(clades)[18], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[19], label=names(clades)[19], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[20], label=names(clades)[20], fontsize=1, offset = 1.5) +
  geom_cladelabel(node=clades[21], label=names(clades)[21], fontsize=1, offset = 1.5)

tree_contodo <- t1 +
  geom_tippoint(aes(color=biogeo1), size=0.8, alpha=1, shape=16) +
  geom_tippoint(aes(x = x+3, color=biogeo2), size = 0.8, alpha=1, shape=16) +
  scale_color_manual(values=colors_biogeo) +
  geom_nodepoint(color="lightblue", size=0.8, shape=16) +
  geom_text2(aes(label=node), size=0.5) +
  theme_tree2(legend.position='right') +
  ggsave("plots/tree/ggtree00.pdf", width=8.27, height=11.69)


# CALIBRATION NODES ----
cal_nodes <- c(node_gekkota = 950,
               node_scincoidea = 918,
               node_lacertoidea = 845,
               node_serpentes = 707,
               node_anguimorpha = 563,
               node_iguania = 634,
               node_sphaerodactylus = 1109,
               node_teratoscincus = 1107,
               node_naultinus_woodworthia = 954,
               node_naultinus_woodworthia_rhacodactylus_oedodera = 953,
               node_podarcis = 904,
               node_naja = 766,
               node_porthidium = 746,
               node_phelsuma = 1071,
               node_oedodera_rhacodactylus = 955,
               node_intellagama = 637,
               node_gallotia_psammodromus = 905,
               node_hodzhhakulia = 846,
               node_balnealacerta = 559,
               node_sphenodon = 556,
               node_leiolep = 638
)

names(cal_nodes)
cal_names <- gsub("node_", "", names(cal_nodes))
names(cal_nodes) <- cal_names

calibrations <- data.frame(name=names(cal_nodes), node=cal_nodes)
calibrations <- calibrations %>% 
  arrange(node) %>%
  mutate(number=1:length(cal_nodes)) %>%
  select(number, name, node)
write.table(calibrations, "calibrations_number.csv", sep=";", quote=F, row.names = F)

tree_cladelabels <- t1

levels(tree_cladelabels$data$biogeo1)

tree_cladelabels_biogeo <- tree_cladelabels +
  geom_tippoint(aes(color=biogeo1), size=0.7, alpha=0.9) +
  geom_tippoint(aes(x = x+3, color=biogeo2), size = 0.7, alpha=0.9) +
  scale_color_manual(values=colors_biogeo) +
  theme_tree2(legend.position='right') +
  ggsave("plots/tree/ggtree_cladelabels_biogeo.pdf", width=8.27, height=11.69)


## Plot tree with clades, biogeo AND CALIBRATIONS ----
tree_cladelabels_biogeo_calibrations <- tree_cladelabels_biogeo %<+%
  calibrations

t_all <- tree_cladelabels_biogeo_calibrations +
  geom_point2(aes(subset=node %in% calibrations$node), color="orchid2", size=1, alpha=1) +
  geom_text2(aes(label=number), size=0.7) +
  ggsave("plots/tree/ggtree_cladelabels_biogeo_calibrations.pdf", width=8.27, height=11.69)

## Add tiplabels (for supplementary) ----
t_all_tiplabels <- t_all +
  geom_tiplab2(size=0.5, offset=3.8) +
  ggsave("plots/tree/ggtree_cladelabels_biogeo_calibrations_tiplabels.pdf", width=8.27, height=11.69)





