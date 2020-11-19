# Description ----
# This script uses the output tree from BEAST2 to extract individual clades 
# and biogeographic data from those clades to reconstruct 
# ancestral areas (with the BioGeoBEARS package).
# It simulates biogeographic histories on those clades.
# It gets information on the specific nodes where the ancestral reconstruction 
# and the simulated biogeographic histories yielded different types of 
# biogeographic events: dispersal, vicariance and range contraction, involvin 
# specific areas (Africa and Arabia).


# Working directory ----
setwd()

# PACKAGES ----
library(dplyr)
library(DescTools)
#library(ips)
#library(cowplot) #package 'cowplot' requires R >= 3.5.0
library(ggplot2)
library(ggtree)
library(tidytree)
library(treeio)
library(phytools)
library(geiger)
#library(treeman)
library(optimx)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(rexpokit)
library(cladoRcpp)
library(devtools)
library(BioGeoBEARS)

calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)# slight speedup hopefully

# Import tree object (treeio format) ----
bigtree <- read.beast(file="data/AfAr_tree.nex")
bigtree_phylo <- as.phylo(bigtree)
bigtree_tibble <- as_tibble(bigtree)

# Define clades based on genera names ----
genera <- c("Acanthodactylus","Atractaspis","Bitis","Cerastes",
            "Chalcides","Chamaeleo","Echis","arid_Hemidactylus",
            "Malpolon","Mesalina","Naja","Pristurus",
            "Psammophis","Pseudotrapelus|Acanthocercus|Xenagama","Ptyodactylus",
            "Scincus|Scincopus|Eumeces",
            "Stenodactylus", "Telescopus", "Tropiocolotes","Uromastyx|Saara",
            "Varanus")

# Get the node for each clade from the big tree ----
nodelist_consensus <- vector(mode="numeric", length(genera))
names(nodelist_consensus) <- genera
for (i in 1:length(genera)){
  nodelist_consensus[i] <- MRCA(bigtree_tibble, grep(genera[i], bigtree_phylo$tip.label))$node
}

# Extract clades from the big tree ----
consensus_genera <- vector("list", length(genera))
names(consensus_genera) <- genera
for (i in 1:length(nodelist_consensus)){
  consensus_genera[[i]] <- tree_subset(bigtree, nodelist_consensus[i], levels_back=0)
}
## In consensus_genera we have the subtrees (treedata class) for each clade of the consensus tree.


# READING LOCATION FILES AND STORING THEM IN A LIST (EACH ELEMENT IS CLASS "tipranges" from BioGeoBEARS) ----
# We have a folder called "locations" in the working directory.
# Within "locations", we have one file per clade with information on species ranges.
location_files <- list.files("locations", full.names=TRUE)
tipranges_list <- vector("list", length(genera))
names(tipranges_list) <- genera

#setwd("locations/") # Not necessary if full.names=TRUE
for(i in 1:length(location_files)){
  tipranges_list[[i]] = getranges_from_LagrangePHYLIP(lgdata_fn=location_files[i])
}
# tipranges_list is a list with the locations for each genus.
# BioGeoBEARS will call each element of that list to reconstruct the biogeography of each genus.

# Set the maximum number of areas any species may occupy; this cannot be larger
# than the number of areas you set up, but it can be smaller.
max_range_size = 2

# Since BioGeoBEARS doesn't work with actual R objects, but with file paths, we'll use the location_files list.
# Furthermore, each time in the loop, the correspondent tree will be exported and then the file path will be used.

# BIOGEOGRAPHIC RECONSTRUCTION (BioGeoBEARS) ----
DECconsensus_output_biogeo <- vector("list", length(genera))
names(DECconsensus_output_biogeo) <- genera
DIVAconsensus_output_biogeo <- vector("list", length(genera))
names(DIVAconsensus_output_biogeo) <- genera
BAYAREAconsensus_output_biogeo <- vector("list", length(genera))
names(BAYAREAconsensus_output_biogeo) <- genera
BESTMODELconsensus <- vector("list", length(genera))
names(BESTMODELconsensus) <- genera
Qmat_consensus <- vector("list", length(genera))
names(Qmat_consensus) <- genera
COOmat_consensus <- vector("list", length(genera))
names(COOmat_consensus) <- genera
rootstate_consensus <- vector("list", length(genera))
names(rootstate_consensus) <- genera
consensus_genera_ape <- vector("list", length(genera))
names(consensus_genera_ape) <- genera
wd <- getwd()
trfn <- np(paste(wd,"/biogeo_tmp/","tree_tmp.nex",sep=""))
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
cores <- 4 # number of cores for BioGeoBEARS analysis

# Functions from BioGeoBEARS to work with a tree in NEXUS format ----
source("functions/bears_optim_run_nexus.R")
source("functions/check_BioGeoBEARS_run_nexus.R")
source("functions/check_trfn_nexus.R")
source("functions/get_Qmat_COOmat_from_BioGeoBEARS_run_object_nexus.R")
source("functions/get_Qmat_COOmat_from_res_nexus.R")
source("functions/readfiles_BioGeoBEARS_run_nexus.R")
source("functions/BioGeoBEARS_extract_Qmat_COOmat_v1.R")

#################################################
##### BIOGEOBEARS LOOP FOR CONSENSUS GENERA #####
#################################################
for (i in 1:length(genera)){
  if (length(as.phylo(consensus_genera[[i]])$tip.label) >= 3){
    write.beast(consensus_genera[[i]], "biogeo_tmp/tree_tmp.nex")
    tr <- read.nexus("biogeo_tmp/tree_tmp.nex")
    consensus_genera_ape[[i]] <- tr
    #trfn <- np(paste(wd,"/biogeo_tmp/","tree_tmp.tre",sep=""))
    print(genera[i])
    geogfn <- location_files[i]
    tipranges <- tipranges_list[[i]]
    # DEC MODEL
    print(paste(genera[i], "DEC"))
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = trfn
    BioGeoBEARS_run_object$geogfn = geogfn
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc. (see Massana et al.)
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = cores
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    check_BioGeoBEARS_run_nexus(BioGeoBEARS_run_object)
    runslow = TRUE
    resfn = "biogeo_tmp/output_DEC.RData"
    if (runslow)
    {
      res = bears_optim_run_nexus(BioGeoBEARS_run_object)
      res
      
      save(res, file=resfn)
      resDEC = res
    } else {
      # Loads to "res"
      load(resfn)
      resDEC = res
    }
    
    #resDEC$condlikes_of_each_state
    # store the output table
    resDEC_list[[i]] <- resDEC
    DECconsensus_output_biogeo[[i]] <- resDEC$ML_marginal_prob_each_state_at_branch_top_AT_node
    #resDEC$outputs@params_table["d","est"] #rate of dispersal (range expansion)
    #resDEC$outputs@params_table["e", "est"] # rate of "extinction" (range contraction or extirpation)
    
    # Plot
    #plot(tr, cex=.5)
    #nodelabels(frame="circle", cex=0.4)
    #tiplabels(frame="circle", cex=0.4)
    analysis_titletxt = paste(genera[i], "DEC")
    results_object = resDEC
    resDEC$inputs
    # States
    pdf(paste("DEC_plots/", genera[i],"_dec.pdf", sep=""))
    res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    # Pie chart
    plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    dev.off()
    
    # DEC+J model
    print(paste(genera[i], "DEC+J"))
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = trfn
    BioGeoBEARS_run_object$geogfn = geogfn
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O???Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DEC+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart = resDEC$outputs@params_table["d","est"]
    estart = resDEC$outputs@params_table["e","est"]
    jstart = 0.0001
    
    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
    
    # Add j as a free parameter
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
    
    check_BioGeoBEARS_run_nexus(BioGeoBEARS_run_object)
    
    resfn = "biogeo_tmp/output_DECJ.RData"
    runslow = TRUE
    if (runslow)
    {
      #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
      
      res = bears_optim_run_nexus(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      
      resDECj = res
    } else {
      # Loads to "res"
      load(resfn)
      resDECj = res
    }
    
    resDECj_list[[i]] <- resDECj
    DECjconsensus_output_biogeo[[i]] <- resDECj$ML_marginal_prob_each_state_at_branch_top_AT_node
    
    #######################################################
    # Plot ancestral states - DECJ
    #######################################################
    analysis_titletxt = paste(genera[i], "DEC+J")
    
    # Setup
    results_object = resDECj
    #scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
    
    pdf(paste("DECJ_plots/", genera[i],"_decj.pdf", sep=""))
    # States
    res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    
    # Pie chart
    plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    
    dev.off()  # Turn off PDF
    
    # DIVA MODEL
    print(paste(genera[i], "DIVA"))
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = trfn
    BioGeoBEARS_run_object$geogfn = geogfn
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc. (see Massana et al.)
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = cores
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
    
    check_BioGeoBEARS_run_nexus(BioGeoBEARS_run_object)
    
    runslow = TRUE
    resfn = "biogeo_tmp/output_DIVA.RData"
    if (runslow)
    {
      res = bears_optim_run_nexus(BioGeoBEARS_run_object)
      res
      
      save(res, file=resfn)
      resDIVALIKE = res
    } else {
      # Loads to "res"
      load(resfn)
      resDIVALIKE = res
    }
    
    # store the output table
    resDIVA_list[[i]] <- resDIVALIKE
    DIVAconsensus_output_biogeo[[i]] <- resDIVALIKE$ML_marginal_prob_each_state_at_branch_top_AT_node
    
    # Plot
    analysis_titletxt = paste(genera[i], "DIVALIKE")
    results_object = resDIVALIKE
    # States
    pdf(paste("DIVA_plots/", genera[i],"_diva.pdf", sep=""), width = 10, height = 15)
    res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    # Pie chart
    plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    dev.off()
    
    # DIVALIKE+J MODEL
    print(paste(genera[i], "DIVA+J"))
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = trfn
    BioGeoBEARS_run_object$geogfn = geogfn
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O???Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DIVALIKE+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart = resDIVALIKE$outputs@params_table["d","est"]
    estart = resDIVALIKE$outputs@params_table["e","est"]
    jstart = 0.0001
    
    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
    
    # Remove subset-sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
    
    # Allow classic, widespread vicariance; all events equiprobable
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
    
    # Add jump dispersal/founder-event speciation
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
    
    # Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999
    
    check_BioGeoBEARS_run_nexus(BioGeoBEARS_run_object)
    
    resfn = "biogeo_tmp/output_DIVAJ.RData"
    runslow = TRUE
    if (runslow)
    {
      #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
      
      res = bears_optim_run_nexus(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      
      resDIVALIKEj = res
    } else {
      # Loads to "res"
      load(resfn)
      resDIVALIKEj = res
    }
    
    # store the output table
    resDIVAj_list[[i]] <- resDIVALIKEj
    DIVAjconsensus_output_biogeo[[i]] <- resDIVALIKEj$ML_marginal_prob_each_state_at_branch_top_AT_node
    
    # Plot
    analysis_titletxt = paste(genera[i], "DIVALIKE+J")
    results_object = resDIVALIKEj
    # States
    pdf(paste("DIVAJ_plots/", genera[i],"_divaj.pdf", sep=""), width = 10, height = 15)
    res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    # Pie chart
    plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    dev.off()
    
    
    # BAYAREA MODEL
    print(paste(genera[i], "BAYAREA"))
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = trfn
    BioGeoBEARS_run_object$geogfn = geogfn
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc. (see Massana et al.)
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = cores
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
    check_BioGeoBEARS_run_nexus(BioGeoBEARS_run_object)
    
    runslow = TRUE
    resfn = "biogeo_tmp/output_BAYAREA.RData"
    if (runslow)
    {
      res = bears_optim_run_nexus(BioGeoBEARS_run_object)
      res
      
      save(res, file=resfn)
      resBAYAREALIKE = res
    } else {
      # Loads to "res"
      load(resfn)
      resBAYAREALIKE = res
    }
    # store the output table
    resBAYAREA_list[[i]] <- resBAYAREALIKE
    BAYAREAconsensus_output_biogeo[[i]] <- resBAYAREALIKE$ML_marginal_prob_each_state_at_branch_top_AT_node
    
    # Plot
    analysis_titletxt = paste(genera[i], "BAYAREALIKE")
    results_object = resBAYAREALIKE
    pdf(paste("BAYAREA_plots/", genera[i],"_bayarea.pdf", sep=""))
    # States
    res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    # Pie chart
    plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    dev.off()
    
    # BAYAREA+J model
    print(paste(genera[i], "BAYAREA+J"))
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = trfn
    BioGeoBEARS_run_object$geogfn = geogfn
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O???Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = "GenSA"
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up BAYAREALIKE+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart = resBAYAREALIKE$outputs@params_table["d","est"]
    estart = resBAYAREALIKE$outputs@params_table["e","est"]
    jstart = 0.0001
    
    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
    
    # No subset sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    # No vicariance
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
    
    # *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
    
    # Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
    
    # Adjust linkage between parameters
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
    
    # Only sympatric/range-copying (y) events allowed, and with 
    # exact copying (both descendants always the same size as the ancestor)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
    
    # NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
    # machines. I can't replicate this on my Mac machines, but it is almost certainly
    # just some precision under-run issue, when optim/optimx tries some parameter value 
    # just below zero.  The "min" and "max" options on each parameter are supposed to
    # prevent this, but apparently optim/optimx sometimes go slightly beyond 
    # these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
    # slightly for each parameter:
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 1e-12
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 1e-12
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 1e-12
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
    
    check_BioGeoBEARS_run_nexus(BioGeoBEARS_run_object)
    
    resfn = "biogeo_tmp/output_BAYAREAJ.RData"
    runslow = TRUE
    if (runslow)
    {
      res = bears_optim_run_nexus(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      
      resBAYAREALIKEj = res
    } else {
      # Loads to "res"
      load(resfn)
      resBAYAREALIKEj = res
    }
    # store the output table
    resBAYAREAj_list[[i]] <- resBAYAREALIKEj
    BAYAREAjconsensus_output_biogeo[[i]] <- resBAYAREALIKEj$ML_marginal_prob_each_state_at_branch_top_AT_node
    
    
    # Plot
    analysis_titletxt = paste(genera[i], "BAYAREALIKE+J")
    results_object = resBAYAREALIKEj
    # States
    pdf(paste("BAYAREAJ_plots/", genera[i],"_bayareaj.pdf", sep=""), width = 10, height = 15)
    res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    # Pie chart
    plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    dev.off()
    
    
    ##### ... Model selection ####
    # Once we have run all the models, we can choose the best fit model according to the AICc
    
    LnL_DEC = get_LnL_from_BioGeoBEARS_results_object(resDEC)
    numparams_DEC = 2
    cols_DEC = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
    
    LnL_DECJ = get_LnL_from_BioGeoBEARS_results_object(resDECj)
    numparams_DECJ = 3
    cols_DECJ = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
    
    LnL_DIVA = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
    numparams_DIVA = 2
    cols_DIVA = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
    
    LnL_DIVAJ = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)
    numparams_DIVAJ = 3
    cols_DIVAJ = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
    
    LnL_BAYAREA = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
    numparams_BAYAREA = 2
    cols_BAYAREA = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
    
    LnL_BAYAREAJ = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)
    numparams_BAYAREAJ = 3
    cols_BAYAREAJ = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
    
    
    # Table with LnL
    restable = rbind(cols_DEC, cols_DECJ, cols_DIVA, cols_DIVAJ, cols_BAYAREA, cols_BAYAREAJ)
    row.names(restable) = c("DEC", "DECJ", "DIVALIKE", "DIVALIKEJ", "BAYAREALIKE", "BAYAREALIKEJ")
    restable
    # Add AIC and AIC weight
    AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
    restable = cbind(restable, AICtable)
    #restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
    #restable_AIC_rellike
    
    if (length(as.phylo(consensus_genera[[i]])$tip.label) >= 4){
      # Add AICc and AICc weight
      samplesize = length(tr$tip.label)
      AICctable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
      restable = cbind(restable, AICctable)
      #restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
      #restable_AICc_rellike
      
      restable_list[[i]] <- restable
      
      # Model selection (min AICc): 1 is DEC, 2 is DEC+J, 3 is DIVA, 4 is DIVA+J,
      # 5 is BAYAREA, 6 is BAYAREA+J
      # and store the Q matrix of the best model for posterior simulations
      if (restable$AICc[1] == min(restable$AICc)){
        BESTMODELconsensus[[i]] <- DECconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DEC"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDEC)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDEC)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDEC$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AICc[2] == min(restable$AICc)){
        BESTMODELconsensus[[i]] <- DECjconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DECj"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDECj)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDECj)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDECj$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AICc[3] == min(restable$AICc)){
        BESTMODELconsensus[[i]] <- DIVAconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DIVA"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKE)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKE)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDIVALIKE$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AICc[4] == min(restable$AICc)){
        BESTMODELconsensus[[i]] <- DIVAjconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DIVAj"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKEj)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKEj)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDIVALIKEj$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AICc[5] == min(restable$AICc)){
        BESTMODELconsensus[[i]] <- BAYAREAconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "BAYAREA"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKE)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKE)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resBAYAREALIKE$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AICc[6] == min(restable$AICc)){
        BESTMODELconsensus[[i]] <- BAYAREAjconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "BAYAREAj"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKEj)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKEj)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resBAYAREALIKEj$relative_probs_of_each_state_at_bottom_of_root_branch)
      }
      
      # For groups with less than 4 species, AICc calculation returns an error. 
      # So we use AIC instead.
    } else {
      restable_list[[i]] <- restable
      
      # Model selection (min AIC): 1 is DEC, 2 is DEC+J, 3 is DIVA, 4 is DIVA+J,
      # 5 is BAYAREA, 6 is BAYAREA+J
      # and store the Q matrix of the best model for posterior simulations
      if (restable$AIC[1] == min(restable$AIC)){
        BESTMODELconsensus[[i]] <- DECconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DEC"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDEC)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDEC)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDEC$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AIC[2] == min(restable$AIC)){
        BESTMODELconsensus[[i]] <- DECjconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DECj"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDECj)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDECj)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDECj$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AIC[3] == min(restable$AIC)){
        BESTMODELconsensus[[i]] <- DIVAconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DIVA"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKE)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKE)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDIVALIKE$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AIC[4] == min(restable$AIC)){
        BESTMODELconsensus[[i]] <- DIVAjconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DIVAj"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKEj)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKEj)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDIVALIKEj$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AIC[5] == min(restable$AIC)){
        BESTMODELconsensus[[i]] <- BAYAREAconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "BAYAREA"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKE)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKE)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resBAYAREALIKE$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AIC[6] == min(restable$AIC)){
        BESTMODELconsensus[[i]] <- BAYAREAjconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "BAYAREAj"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKEj)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKEj)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resBAYAREALIKEj$relative_probs_of_each_state_at_bottom_of_root_branch)
      }
      
      
      
    }
    
  } else {
    print(paste(genera[i], " has less than 3 species! =(", sep=""))
  }
  
}
##### END OF THE CONSENSUS BIOGEOBEARS LOOP ####

saveRDS(resDEC_list, "objects/resDEC_list.rds")
saveRDS(resDECj_list, "objects/resDECj_list.rds")
saveRDS(resDIVA_list, "objects/resDIVA_list.rds")
saveRDS(resDIVAj_list, "objects/resDIVAj_list.rds")
saveRDS(resBAYAREA_list, "objects/resBAYAREA_list.rds")
saveRDS(resBAYAREAj_list, "objects/resBAYAREAj_list.rds")

resDEC_list <- readRDS("objects/resDEC_list.rds")
resDIVA_list <- readRDS("objects/resDIVA_list.rds")
resBAYAREA_list <- readRDS("objects/resBAYAREA_list.rds")

resDECj_list <- readRDS("objects/resDECj_list.rds")
resDIVAj_list <- readRDS("objects/resDIVAj_list.rds")
resBAYAREAj_list <- readRDS("objects/resBAYAREAj_list.rds")

consensus_genera <- readRDS("objects/consensus_genera.rds")

source("functions/modsel.R")
source("functions/get_Qmat_COOmat_from_res_nexus.R")
source("functions/get_Qmat_COOmat_from_BioGeoBEARS_run_object_nexus.R")
source("functions/check_trfn_nexus.R")
source("functions/readfiles_BioGeoBEARS_run_nexus.R")

source("functions/bears_optim_run_nexus.R")
source("functions/check_BioGeoBEARS_run_nexus.R")
source("functions/BioGeoBEARS_extract_Qmat_COOmat_v1.R")

# MODEL SELECTION ----
res_bgb_list_all <- list(resDEC_list, resDECj_list, resDIVA_list, resDIVAj_list,
                         resBAYAREA_list, resBAYAREAj_list)
names(res_bgb_list_all) <- c('DEC', 'DECJ', 'DIVA', 'DIVAJ', 'BAYAREA', 'BAYAREAJ')

# Choose the models to include in model selection:
res_bgb_list <- list(resDEC_list, resDIVA_list, resBAYAREA_list)
names(res_bgb_list) <- c('DEC', 'DIVA', 'BAYAREA')

# Run the model selection function
a <- modsel(results_list = res_bgb_list, treedata_list = consensus_genera2)
length(res_bgb_list$DEC)
names(res_bgb_list$DEC)
names(consensus_genera)
length(consensus_genera)

# Save output:
BESTMODELconsensus <- a$BESTMODELconsensus
restable_list <- a$restable_list
table_bestmodels <- a$table_bestmodels
Qmat_consensus <- a$Qmat_consensus
COOmat_consensus <- a$COOmat_consensus
rootstate_consensus <- a$rootstate_consensus


saveRDS(BESTMODELconsensus, "objects/BESTMODELconsensus.rds")
saveRDS(restable_list, "objects/restable_list.rds")
saveRDS(table_bestmodels, "objects/table_bestmodels.rds")
saveRDS(Qmat_consensus, "objects/Qmat_consensus.rds")
saveRDS(COOmat_consensus, "objects/COOmat_consensus.rds")
saveRDS(rootstate_consensus, "objects/rootstate_consensus.rds")


table_bestmodels <- readRDS("objects/table_bestmodels.rds")
write.table(table_bestmodels, "objects/table_bestmodels.csv", sep=";", 
            row.names = FALSE, quote = FALSE)
BESTMODELconsensus <- readRDS("objects/BESTMODELconsensus.rds")

consensus_genera <- readRDS("objects/consensus_genera.rds")

# In BESTMODELconsensus there are matrices with the probabilities of each biogeographic state
# in each tip and node.


# TRANSFORM PROBABILITY MATRICES INTO VECTORS OF STATES ----
##### TRANSFORM PROBABILITIES MATRICES INTO VECTORS OF STATES #####
for(i in 1:length(location_files)){
  tipranges_list[[i]] = getranges_from_LagrangePHYLIP(lgdata_fn=location_files[i])
}
which(names(tipranges_list) == "Lytorhynchus")
tipranges_list <- tipranges_list[-9]

# In consensus genera
consensus_states <- vector("list", length(genera))
names(consensus_states) <- genera
for (i in 1:length(genera)){
  if (length(as.phylo(consensus_genera[[i]])$tip.label) >= 3){
    areas = getareas_from_tipranges_object(tipranges_list[[i]])
    statenames = areas_list_to_states_list_new(areas, maxareas = max_range_size,
                                               include_null_range = TRUE,
                                               split_ABC = FALSE)
    statenames = unlist(statenames)
    relprobs_matrix <- BESTMODELconsensus[[i]]
    MLprobs = get_ML_probs(relprobs_matrix)
    MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames,
                                           returnwhat = "states", if_ties = "takefirst")
    consensus_states[[i]] <- MLstates
  }
}
# consensus_states is a list with vectors of states for each genera (consensus)
saveRDS(consensus_states, "objects/consensus_states.rds")
consensus_states <- readRDS("objects/consensus_states.rds")

#DEFINE BIOGEOGRAPHIC SCENARIOS ----
# Africa: "F"
# Arabia: "R"
# Africa+Arabia: "FR"
# Establish the different scenarios. Each scenario will be a list of vectors of two/three elements:
# (state1, state2, state3), where state1 = state in parental node, state2 = state in descendendant 1,
# state3 = state in descendant 2.
# For dispersal and extirpation (anagenetic processes), we only need two states (parental and descendant); 
# for vicariance, we need three states: parental and the two descendent nodes. 
dispersal_Af2Ar <- list(c('F', 'R'), c('F', 'FR'), 
                        c('FE', 'FR'), c('FE', 'RE'), 
                        c('FS', 'FR'), c('FS', 'RS'))

dispersal_Ar2Af <- list(c('R', 'F'), c('R', 'FR'), 
                        c('RE', 'FR'), c('RE', 'FE'),
                        c('RS', 'FR'), c('RS', 'FS'))

vicariance <- list(c("FR","F","R"), c("FR","R","F"), 
                   c("FR", "RE", "F"), c("FR", "F", "RE"),
                   c("FR", "FE", "R"), c("FR", "R", "FE"),
                   c("FR", "F", "RS"), c("FR", "RS", "F"),
                   c("FR", "FS", "R"), c("FR", "R", "FS"),
                   c("FR", "FE", "RE"), c("FR", "RE", "FE"),
                   c("FR", "FS", "RS"), c("FR", "RS", "FS"),
                   c("FR", "FE", "RS"), c("FR", "RS", "FE"),
                   c("FR", "FS", "RE"), c("FR", "RE", "FS"))

extirpation_from_Ar <- list(c("FR", "F"), c("FR", "FE"), c("FR", "FS"))

extirpation_from_Af <- list(c("FR", "R"), c("FR", "RE"), c("FR", "RS"))

# Now that we have established the scenarios, we can go node by node categorizing all of them
# in their corresponding scenarios.

##### .-.-.-.-. COUNTING AND CATEGORIZING BIOGEOGRAPHIC EVENTS IN CONSENSUS #####
Af2Ar_nodes_list_cons <- vector("list", length(genera))
names(Af2Ar_nodes_list_cons) <- genera
Af2Ar_nodes_list_cons[1:length(genera)] <- 0

Ar2Af_nodes_list_cons <- vector("list", length(genera))
names(Ar2Af_nodes_list_cons) <- genera
Ar2Af_nodes_list_cons[1:length(genera)] <- 0

vicariance_nodes_list_cons <- vector("list", length(genera))
names(vicariance_nodes_list_cons) <- genera
vicariance_nodes_list_cons[1:length(genera)] <- 0

extirp_from_Af_nodes_list_cons <- vector("list", length(genera))
names(extirp_from_Af_nodes_list_cons) <- genera
extirp_from_Af_nodes_list_cons[1:length(genera)] <- 0

extirp_from_Ar_nodes_list_cons <- vector("list", length(genera))
names(extirp_from_Ar_nodes_list_cons) <- genera
extirp_from_Ar_nodes_list_cons[1:length(genera)] <- 0



# Get vicariance nodes (the parental node)
for (i in 1:length(genera)){
  if (length(as.phylo(consensus_genera[[i]])$tip.label) >= 3){
    for (j in min(as.phylo(consensus_genera[[i]])$edge[,1]):max(as.phylo(consensus_genera[[i]])$edge[,1])){
      v <- c(consensus_states[[i]][j],
             consensus_states[[i]][as.phylo(consensus_genera[[i]])$edge[,2][as.phylo(consensus_genera[[i]])$edge[,1]==j]][1],
             consensus_states[[i]][as.phylo(consensus_genera[[i]])$edge[,2][as.phylo(consensus_genera[[i]])$edge[,1]==j]][2])
      if (list(v) %in% vicariance){
        vicariance_nodes_list_cons[[i]] <- append(vicariance_nodes_list_cons[[i]], j)
      }
    }
  }
}

# Get dispersal and extirpation nodes (the descendent node)
'%!in%' <- function(x,y)!('%in%'(x,y))
for (i in 1:length(genera)){
  if (length(as.phylo(consensus_genera[[i]])$tip.label) >= 3){
    for (j in sort(as.phylo(consensus_genera[[i]])$edge[,2])){
      parent <- consensus_states[[i]][as.phylo(consensus_genera[[i]])$edge[,1][as.phylo(consensus_genera[[i]])$edge[,2]==j]]
      son <- consensus_states[[i]][j]
      d <- c(parent, son)
      if (list(d) %in% dispersal_Af2Ar){
        Af2Ar_nodes_list_cons[[i]] <- append(Af2Ar_nodes_list_cons[[i]], j)
      }
      if (list(d) %in% dispersal_Ar2Af){
        Ar2Af_nodes_list_cons[[i]] <- append(Ar2Af_nodes_list_cons[[i]], j)
      }
      
      # extirpation (if the parent node is not a vicariant node)
      if (as.phylo(consensus_genera[[i]])$edge[,1][as.phylo(consensus_genera[[i]])$edge[,2]==j] %!in% vicariance_nodes_list_cons[[i]]){
        if (list(d) %in% extirpation_from_Ar){
          extirp_from_Ar_nodes_list_cons[[i]] <- append(extirp_from_Ar_nodes_list_cons[[i]], j)
        }
        if (list(d) %in% extirpation_from_Af){
          extirp_from_Af_nodes_list_cons[[i]] <- append(extirp_from_Af_nodes_list_cons[[i]], j)
        }
      }
    }
  }
}

#max(unlist(lapply(Af2Ar_nodes_list_cons, length)))
#max(unlist(lapply(Ar2Af_nodes_list_cons, length)))
#max(unlist(lapply(vicariance_nodes_list_cons, length)))
#max(unlist(lapply(extirp_from_Af_nodes_list_cons, length)))
#max(unlist(lapply(extirp_from_Ar_nodes_list_cons, length)))

#which(unlist(lapply(Af2Ar_nodes_list_cons, length)) == 2)

# Now we have 5 lists, one for each scenario of dispersal, vicariance and extirpation.
# Each of those lists contains vectors of simulated transition nodes for each genus.
# First, we're going to remove the zeros from the vectors with transition nodes.

for (i in 1:length(genera)){
  if (length(Af2Ar_nodes_list_cons[[i]]) != 1){
    Af2Ar_nodes_list_cons[[i]] <- Af2Ar_nodes_list_cons[[i]][-1]
  }
  if (length(Ar2Af_nodes_list_cons[[i]]) != 1){
    Ar2Af_nodes_list_cons[[i]] <- Ar2Af_nodes_list_cons[[i]][-1]
  }
  if (length(vicariance_nodes_list_cons[[i]]) != 1){
    vicariance_nodes_list_cons[[i]] <- vicariance_nodes_list_cons[[i]][-1]
  }
  if (length(extirp_from_Af_nodes_list_cons[[i]]) != 1){
    extirp_from_Af_nodes_list_cons[[i]] <- extirp_from_Af_nodes_list_cons[[i]][-1]
  }
  if (length(extirp_from_Ar_nodes_list_cons[[i]]) != 1){
    extirp_from_Ar_nodes_list_cons[[i]] <- extirp_from_Ar_nodes_list_cons[[i]][-1]
  }
}

# Now we have already removed the zeros from the vectors where there are transition nodes.
# We can make a data frame with them. We're going to make 5 lists of *number of genera* dataframes (some of them
# will not have nodes) with the transition nodes.
Af2Ar_nodes_dfs_cons <- vector("list", length(genera))
names(Af2Ar_nodes_dfs_cons) <- genera
Ar2Af_nodes_dfs_cons <- vector("list", length(genera))
names(Ar2Af_nodes_dfs_cons) <- genera
vicariance_nodes_dfs_cons <- vector("list", length(genera))
names(vicariance_nodes_dfs_cons) <- genera
extirp_from_Af_nodes_dfs_cons <- vector("list", length(genera))
names(extirp_from_Af_nodes_dfs_cons) <- genera
extirp_from_Ar_nodes_dfs_cons <- vector("list", length(genera))
names(extirp_from_Ar_nodes_dfs_cons) <- genera


for (i in 1:length(genera)){
  Af2Ar_nodes_dfs_cons[[i]] <- data.frame(node=Af2Ar_nodes_list_cons[[i]], genus=genera[i],
                                          height=0, min=0, max=0, Af2Ar=1, Ar2Af=0, Vic=0, ExtAf=0, ExtAr=0)
  Ar2Af_nodes_dfs_cons[[i]] <- data.frame(node=Ar2Af_nodes_list_cons[[i]], genus=genera[i],
                                          height=0, min=0, max=0, Af2Ar=0, Ar2Af=1, Vic=0, ExtAf=0, ExtAr=0)
  vicariance_nodes_dfs_cons[[i]] <- data.frame(node=vicariance_nodes_list_cons[[i]], genus=genera[i],
                                               height=0, min=0, max=0, Af2Ar=0, Ar2Af=0, Vic=1, ExtAf=0, ExtAr=0)
  extirp_from_Af_nodes_dfs_cons[[i]] <- data.frame(node=extirp_from_Af_nodes_list_cons[[i]], genus=genera[i],
                                                   height=0, min=0, max=0, Af2Ar=0, Ar2Af=0, Vic=0, ExtAf=1, ExtAr=0)
  extirp_from_Ar_nodes_dfs_cons[[i]] <- data.frame(node=extirp_from_Ar_nodes_list_cons[[i]], genus=genera[i],
                                                   height=0, min=0, max=0, Af2Ar=0, Ar2Af=0, Vic=0, ExtAf=0, ExtAr=1)
  
}

# Now we have 5 lists, one for each scenario, and in each list there are *number of genera* data frames,
# one per genus, with the nodes where there has been a biogeographic event.

# We want a list with one data frame per genus.
#rbind(Af2Ar_nodes_dfs_cons[[1]], Ar2Af_nodes_dfs_cons[[1]])
#Af2Ar_nodes_dfs_cons[[1]]$node == 0
#extirp_from_Af_nodes_dfs_cons[[15]]$node[1] == 0
#rbind(extirp_from_Af_nodes_dfs_cons[[15]], Af2Ar_nodes_dfs_cons[[15]])

event_list_cons <- vector("list", length(genera))
names(event_list_cons) <- genera

for (i in 1:length(genera)){
  event_list_cons[[i]] <- rbind(Af2Ar_nodes_dfs_cons[[i]], Ar2Af_nodes_dfs_cons[[i]],
                                vicariance_nodes_dfs_cons[[i]], extirp_from_Af_nodes_dfs_cons[[i]],
                                extirp_from_Ar_nodes_dfs_cons[[i]])
}


# Now we already have one data frame per genus with the transition nodes. We can remove the zeros:
for (i in 1:length(genera)){
  event_list_cons[[i]] <- event_list_cons[[i]][event_list_cons[[i]]$node!=0,]
}

# Now we need to include information about the divergence times: the height of the nodes with transitions.

##### INCLUDE INFO ON DIVERGENCE TIMES FOR CONSENSUS #####
#consensus_genera[[1]]
#event_list_cons
#as_tibble(consensus_genera[[1]])$height_0.95_HPD[[60]]
#as_tibble(consensus_genera[[1]])$height_0.95_HPD[[60]][1]
#as_tibble(consensus_genera[[1]])$height_0.95_HPD[[60]][2]
#as_tibble(consensus_genera[[1]])$height_0.95_HPD[[event_list_cons[[1]]$node[2]]][2]
#as_tibble(consensus_genera[[1]])$height_median[[event_list_cons[[1]]$node[2]]]
#nrow(event_list_cons[[1]])

# Divergence times and 95% HPD (min and max) for each vicariance node 
# and branch length for dispersal and extirpation nodes.
for (i in 1:length(genera)){
  if (nrow(event_list_cons[[i]]) > 0){
    for (j in 1:length(event_list_cons[[i]]$node)){
      
      # Vicariance nodes
      if (event_list_cons[[i]][j, "Vic"] == 1){
        event_list_cons[[i]]$min[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[event_list_cons[[i]]$node[j]]][1]
        event_list_cons[[i]]$max[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[event_list_cons[[i]]$node[j]]][2]
        event_list_cons[[i]]$height[j] <- as_tibble(consensus_genera[[i]])$height_median[[event_list_cons[[i]]$node[j]]]
      }
      
      # Dispersal and extirpation nodes
      if (event_list_cons[[i]][j, "Af2Ar"] == 1 | event_list_cons[[i]][j, "Ar2Af"] == 1 | event_list_cons[[i]][j, "ExtAf"] == 1 | event_list_cons[[i]][j, "ExtAr"] == 1){
        dat <- as_tibble(consensus_genera[[i]])
        # Min is the age of the descendent node
        desc_age <- dat$height[dat$node == event_list_cons[[i]]$node[j]]
        event_list_cons[[i]]$min[j] <- desc_age
        
        # Max is the age of the parental node
        parent_age <- dat$height[dat$parent[dat$node == event_list_cons[[i]]$node[j]]]
        event_list_cons[[i]]$max[j] <- parent_age
        
        # let's set the height as the midpoint between the parental and the descendent nodes
        event_list_cons[[i]]$height[j] <- mean(c(desc_age, parent_age))
        
      }
    }
  }
}

# With this object we can represent the age and type of each biogeographic event for each clade.
saveRDS(event_list_cons, "objects/event_list_cons.rds")

##### EVENTS PER MILLION YEAR FOR EACH GENUS (CONSENSUS) #####
myr <- 60
ma <- myr
event_myr_list_cons <- vector("list", length(genera))
names(event_myr_list_cons) <- genera
for (i in 1:length(genera)){
  event_myr_list_cons[[i]] <- data.frame(Ma=c(1:60), N=0)
}
event_myr_list_cons[[1]][1,]$N

# We also can create a list of vectors representing intervals of 1 Ma.
ma <- vector("list", ma)
for (i in 1:length(ma)){
  ma[[i]] <- c(i-1, i)
}

# Counting the maximum number of transition in each million year for each genus.
for (i in 1:length(genera)){
  if (nrow(event_list_cons[[i]]) > 0){
    for (j in 1:nrow(event_list_cons[[i]])){
      vi <- c(event_list_cons[[i]]$min[j], event_list_cons[[i]]$max[j])
      for (k in 1:length(ma)){
        if (Overlap(na.omit(vi), ma[[k]]) != 0){
          event_myr_list_cons[[i]][k,]$N <- event_myr_list_cons[[i]][k,]$N + 1
        }
      }
    }
  }
}

# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
event_all_cons <- data.frame(Ma=c(1:myr), N=0)
for (i in 1:length(genera)){
  event_all_cons$N <- event_all_cons$N + event_myr_list_cons[[i]]$N
}

##### AFRICA TO ARABIA PER MILLION YEAR FOR EACH GENUS (CONSENSUS) #####
Af2Ar_myr_list_cons <- vector("list", length(genera))
names(Af2Ar_myr_list_cons) <- genera
for (i in 1:length(genera)){
  Af2Ar_myr_list_cons[[i]] <- data.frame(Ma=c(1:60), N=0)
}
Af2Ar_myr_list_cons[[1]][1,]$N

# Counting the maximum number of Af -> Ar dispersals in each million year for each genus.
for (i in 1:length(genera)){
  if (nrow(event_list_cons[[i]]) > 0){
    for (j in 1:nrow(event_list_cons[[i]])){
      if (event_list_cons[[i]]$Af2Ar[j] > 0){
        vi <- c(event_list_cons[[i]]$min[j], event_list_cons[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            Af2Ar_myr_list_cons[[i]][k,]$N <- Af2Ar_myr_list_cons[[i]][k,]$N + 1
          }
        }
      }
    }
  }
}
Af2Ar_myr_list_cons



# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
Af2Ar_all_cons <- data.frame(Ma=c(1:myr), N=0)
for (i in 1:length(genera)){
  Af2Ar_all_cons$N <- Af2Ar_all_cons$N + Af2Ar_myr_list_cons[[i]]$N
}

##### ARABIA TO AFRICA PER MILLION YEAR FOR EACH GENUS (CONSENSUS) #####
Ar2Af_myr_list_cons <- vector("list", length(genera))
names(Ar2Af_myr_list_cons) <- genera
for (i in 1:length(genera)){
  Ar2Af_myr_list_cons[[i]] <- data.frame(Ma=c(1:60), N=0)
}
Ar2Af_myr_list_cons[[1]][1,]$N

# Counting the maximum number of Ar -> Af dispersals in each million year for each genus.
for (i in 1:length(genera)){
  if (nrow(event_list_cons[[i]]) > 0){
    for (j in 1:nrow(event_list_cons[[i]])){
      if (event_list_cons[[i]]$Ar2Af[j] > 0){
        vi <- c(event_list_cons[[i]]$min[j], event_list_cons[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            Ar2Af_myr_list_cons[[i]][k,]$N <- Ar2Af_myr_list_cons[[i]][k,]$N + 1
          }
        }
      }
    }
  }
}
Ar2Af_myr_list_cons



# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
Ar2Af_all_cons <- data.frame(Ma=c(1:myr), N=0)
for (i in 1:length(genera)){
  Ar2Af_all_cons$N <- Ar2Af_all_cons$N + Ar2Af_myr_list_cons[[i]]$N
}

##### VICARIANCE PER MILLION YEAR FOR EACH GENUS (CONSENSUS) #####
vicariance_myr_list_cons <- vector("list", length(genera))
names(vicariance_myr_list_cons) <- genera
for (i in 1:length(genera)){
  vicariance_myr_list_cons[[i]] <- data.frame(Ma=c(1:60), N=0)
}
vicariance_myr_list_cons[[1]][1,]$N

# Counting the maximum number of vicariance events in each million year for each genus.
for (i in 1:length(genera)){
  if (nrow(event_list_cons[[i]]) > 0){
    for (j in 1:nrow(event_list_cons[[i]])){
      if (event_list_cons[[i]]$Vic[j] > 0){
        vi <- c(event_list_cons[[i]]$min[j], event_list_cons[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            vicariance_myr_list_cons[[i]][k,]$N <- vicariance_myr_list_cons[[i]][k,]$N + 1
          }
        }
      }
    }
  }
}
vicariance_myr_list_cons

# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
vicariance_all_cons <- data.frame(Ma=c(1:myr), N=0)
for (i in 1:length(genera)){
  vicariance_all_cons$N <- vicariance_all_cons$N + vicariance_myr_list_cons[[i]]$N
}

##### EXTIRPATION FROM AFRICA PER MILLION YEAR FOR EACH GENUS (CONSENSUS) #####
extirp_from_Af_myr_list_cons <- vector("list", length(genera))
names(extirp_from_Af_myr_list_cons) <- genera
for (i in 1:length(genera)){
  extirp_from_Af_myr_list_cons[[i]] <- data.frame(Ma=c(1:60), N=0)
}
extirp_from_Af_myr_list_cons[[1]][1,]$N

# Counting the maximum number of vicariance events in each million year for each genus.
for (i in 1:length(genera)){
  if (nrow(event_list_cons[[i]]) > 0){
    for (j in 1:nrow(event_list_cons[[i]])){
      if (event_list_cons[[i]]$ExtAf[j] > 0){
        vi <- c(event_list_cons[[i]]$min[j], event_list_cons[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            extirp_from_Af_myr_list_cons[[i]][k,]$N <- extirp_from_Af_myr_list_cons[[i]][k,]$N + 1
          }
        }
      }
    }
  }
}
extirp_from_Af_myr_list_cons

# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
extirp_from_Af_all_cons <- data.frame(Ma=c(1:myr), N=0)
for (i in 1:length(genera)){
  extirp_from_Af_all_cons$N <- extirp_from_Af_all_cons$N + extirp_from_Af_myr_list_cons[[i]]$N
}

##### EXTIRPATION FROM ARABIA PER MILLION YEAR FOR EACH GENUS (CONSENSUS) #####
extirp_from_Ar_myr_list_cons <- vector("list", length(genera))
names(extirp_from_Ar_myr_list_cons) <- genera
for (i in 1:length(genera)){
  extirp_from_Ar_myr_list_cons[[i]] <- data.frame(Ma=c(1:60), N=0)
}
extirp_from_Ar_myr_list_cons[[1]][1,]$N

# Counting the maximum number of vicariance events in each million year for each genus.
for (i in 1:length(genera)){
  if (nrow(event_list_cons[[i]]) > 0){
    for (j in 1:nrow(event_list_cons[[i]])){
      if (event_list_cons[[i]]$ExtAr[j] > 0){
        vi <- c(event_list_cons[[i]]$min[j], event_list_cons[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            extirp_from_Ar_myr_list_cons[[i]][k,]$N <- extirp_from_Ar_myr_list_cons[[i]][k,]$N + 1
          }
        }
      }
    }
  }
}
extirp_from_Ar_myr_list_cons

# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
extirp_from_Ar_all_cons <- data.frame(Ma=c(1:myr), N=0)
for (i in 1:length(genera)){
  extirp_from_Ar_all_cons$N <- extirp_from_Ar_all_cons$N + extirp_from_Ar_myr_list_cons[[i]]$N
}


#######################
##### SIMULATIONS #####
#######################
file.sources = list.files("functions/", pattern="*.R", full.names=T)
file.sources <- file.sources[-6]
source("functions/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
for (i in 1:length(file.sources)){
  source(file.sources[i])
}
Qmat_consensus <- readRDS("objects/Qmat_consensus.rds")
COOmat_consensus <- readRDS("objects/COOmat_consensus.rds")
rootstate_consensus <- readRDS("objects/rootstate_consensus.rds")


##### Simulations with consensus tree #####
sims_consensus <- vector("list", length(genera))
names(sims_consensus) <- genera
nsim <- 1000
for (i in 1:length(genera)){
  if (Ntip(consensus_genera[[i]]) > 2){
    for (s in 1:nsim){
      sims_consensus[[i]][[s]] <- simulate_biogeog_history(phy=as.phylo(consensus_genera[[i]]), Qmat=Qmat_consensus[[i]],
                                                           COO_probs_columnar=COOmat_consensus[[i]],
                                                           index_Qmat_0based_of_starting_state=rootstate_consensus[[i]]-1)
    }
  }
}

# Add 1 to all the states because the function doesn't count the null state (state 1).
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (!is.null(sims_consensus[[i]])){
      sims_consensus[[i]][[s]] <- sims_consensus[[i]][[s]] + 1
    }
  }
}

# Names of the states for each genus
consensus_states_names <- vector("list", length(genera))
names(consensus_states) <- genera
for (i in 1:length(genera)){
  if (length(as.phylo(consensus_genera[[i]])$tip.label) >= 3){
    areas = getareas_from_tipranges_object(tipranges_list[[i]])
    consensus_states_names[[i]] = areas_list_to_states_list_new(areas, maxareas = max_range_size,
                                                                include_null_range = TRUE,
                                                                split_ABC = FALSE)
    consensus_states_names[[i]] = unlist(consensus_states_names[[i]])
    
  }
}

# Now we can transform the numeric states of sims_consensus into the names (F, R, etc)
sims_consensus_letters <- vector("list", length(genera))
names(sims_consensus_letters) <- genera
for (i in 1:length(genera)){
  sims_consensus_letters[[i]] <- vector("list", nsim)
}

for (i in 1:length(genera)){
  for (s in 1:nsim){
    sims_consensus_letters[[i]][[s]] <- vector(mode="character")
  }
}

for (i in 1:length(genera)){
  if (Ntip(consensus_genera[[i]]) > 2){
    for (s in 1:nsim){
      for (j in 1:length(sims_consensus[[i]][[s]])){
        sims_consensus_letters[[i]][[s]][j] <- consensus_states_names[[i]][sims_consensus[[i]][[s]][j]]
      }
    }
  }
}


##### .-.-.-.-. COUNTING AND CATEGORIZING BIOGEOGRAPHIC EVENTS IN SIMULATIONS (CONSENSUS) #####
Af2Ar_SIM_nodes_list_cons <- vector("list", length(genera))
names(Af2Ar_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){
  Af2Ar_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    Af2Ar_SIM_nodes_list_cons[[i]][[s]] <- 0
  }
}

Ar2Af_SIM_nodes_list_cons <- vector("list", length(genera))
names(Ar2Af_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){
  Ar2Af_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    Ar2Af_SIM_nodes_list_cons[[i]][[s]] <- 0
  }
}

vicariance_SIM_nodes_list_cons <- vector("list", length(genera))
names(vicariance_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){
  vicariance_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    vicariance_SIM_nodes_list_cons[[i]][[s]] <- 0
  }
}

extirpfromAf_SIM_nodes_list_cons <- vector("list", length(genera))
names(extirpfromAf_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){
  extirpfromAf_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    extirpfromAf_SIM_nodes_list_cons[[i]][[s]] <- 0
  }
}

extirpfromAr_SIM_nodes_list_cons <- vector("list", length(genera))
names(extirpfromAr_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){
  extirpfromAr_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    extirpfromAr_SIM_nodes_list_cons[[i]][[s]] <- 0
  }
}

# Get vicariance nodes (the parental node)
for (i in 1:length(genera)){
  if (length(as.phylo(consensus_genera[[i]])$tip.label) >= 3){
    for (s in 1:nsim){
      for (j in min(as.phylo(consensus_genera[[i]])$edge[,1]):max(as.phylo(consensus_genera[[i]])$edge[,1])){
        v <- c(sims_consensus_letters[[i]][[s]][j],
               sims_consensus_letters[[i]][[s]][as.phylo(consensus_genera[[i]])$edge[,2][as.phylo(consensus_genera[[i]])$edge[,1]==j]][1],
               sims_consensus_letters[[i]][[j]][as.phylo(consensus_genera[[i]])$edge[,2][as.phylo(consensus_genera[[i]])$edge[,1]==j]][2])
        if (list(v) %in% vicariance){
          vicariance_SIM_nodes_list_cons[[i]][[s]] <- append(vicariance_SIM_nodes_list_cons[[i]][[s]], j)
        }
      }
    }
  }
}

# Get dispersal and extirpation nodes (the descendent node)
'%!in%' <- function(x,y)!('%in%'(x,y))
for (i in 1:length(genera)){
  if (length(as.phylo(consensus_genera[[i]])$tip.label) >= 3){
    for (s in 1:nsim){
      for (j in sort(as.phylo(consensus_genera[[i]])$edge[,2])){
        parent <- sims_consensus_letters[[i]][[s]][as.phylo(consensus_genera[[i]])$edge[,1][as.phylo(consensus_genera[[i]])$edge[,2]==j]]
        son <- sims_consensus_letters[[i]][[s]][j]
        d <- c(parent, son)
        if (list(d) %in% dispersal_Af2Ar){
          Af2Ar_SIM_nodes_list_cons[[i]][[s]] <- append(Af2Ar_SIM_nodes_list_cons[[i]][[s]], j)
          
        }
        if (list(d) %in% dispersal_Ar2Af){
          Ar2Af_SIM_nodes_list_cons[[i]][[s]] <- append(Ar2Af_SIM_nodes_list_cons[[i]][[s]], j)
        }
        
        # extirpation (if the parent node is not a vicariant node)
        if (as.phylo(consensus_genera[[i]])$edge[,1][as.phylo(consensus_genera[[i]])$edge[,2]==j] %!in% vicariance_SIM_nodes_list_cons[[i]][[s]]){
          if (list(d) %in% extirpation_from_Ar){
            extirpfromAr_SIM_nodes_list_cons[[i]][[s]] <- append(extirpfromAr_SIM_nodes_list_cons[[i]][[s]], j)
          }
          if (list(d) %in% extirpation_from_Af){
            extirpfromAf_SIM_nodes_list_cons[[i]][[s]] <- append(extirpfromAf_SIM_nodes_list_cons[[i]][[s]], j)
          }
        }  
      }
    }
  }
}


# Remove zeros.

for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (length(Af2Ar_SIM_nodes_list_cons[[i]][[s]]) != 1){
      Af2Ar_SIM_nodes_list_cons[[i]][[s]] <- Af2Ar_SIM_nodes_list_cons[[i]][[s]][-1]
    }
    if (length(Ar2Af_SIM_nodes_list_cons[[i]][[s]]) != 1){
      Ar2Af_SIM_nodes_list_cons[[i]][[s]] <- Ar2Af_SIM_nodes_list_cons[[i]][[s]][-1]
    }
    if (length(vicariance_SIM_nodes_list_cons[[i]][[s]]) != 1){
      vicariance_SIM_nodes_list_cons[[i]][[s]] <- vicariance_SIM_nodes_list_cons[[i]][[s]][-1]
    }
    if (length(extirpfromAf_SIM_nodes_list_cons[[i]][[s]]) != 1){
      extirpfromAf_SIM_nodes_list_cons[[i]][[s]] <- extirpfromAf_SIM_nodes_list_cons[[i]][[s]][-1]
    }
    if (length(extirpfromAr_SIM_nodes_list_cons[[i]][[s]]) != 1){
      extirpfromAr_SIM_nodes_list_cons[[i]][[s]] <- extirpfromAr_SIM_nodes_list_cons[[i]][[s]][-1]
    }
  }
}

# Make data frames
Af2Ar_SIM_nodes_dfs_cons <- vector("list", length(genera))
names(Af2Ar_SIM_nodes_dfs_cons) <- genera
for (i in 1:length(genera)){
  Af2Ar_SIM_nodes_dfs_cons[[i]] <- vector("list", nsim)
}

Ar2Af_SIM_nodes_dfs_cons <- vector("list", length(genera))
names(Ar2Af_SIM_nodes_dfs_cons) <- genera
for (i in 1:length(genera)){
  Ar2Af_SIM_nodes_dfs_cons[[i]] <- vector("list", nsim)
}

vicariance_SIM_nodes_dfs_cons <- vector("list", length(genera))
names(vicariance_SIM_nodes_dfs_cons) <- genera
for (i in 1:length(genera)){
  vicariance_SIM_nodes_dfs_cons[[i]] <- vector("list", nsim)
}

extirpfromAf_SIM_nodes_dfs_cons <- vector("list", length(genera))
names(extirpfromAf_SIM_nodes_dfs_cons) <- genera
for (i in 1:length(genera)){
  extirpfromAf_SIM_nodes_dfs_cons[[i]] <- vector("list", nsim)
}

extirpfromAr_SIM_nodes_dfs_cons <- vector("list", length(genera))
names(extirpfromAr_SIM_nodes_dfs_cons) <- genera
for (i in 1:length(genera)){
  extirpfromAr_SIM_nodes_dfs_cons[[i]] <- vector("list", nsim)
}


for (i in 1:length(genera)){
  for (s in 1:nsim){
    Af2Ar_SIM_nodes_dfs_cons[[i]][[s]] <- data.frame(node=Af2Ar_SIM_nodes_list_cons[[i]][[s]], genus=genera[i],
                                                     height=0, min=0, max=0, Af2Ar=1, Ar2Af=0, Vic=0, ExtAf=0, ExtAr=0)
    Ar2Af_SIM_nodes_dfs_cons[[i]][[s]] <- data.frame(node=Ar2Af_SIM_nodes_list_cons[[i]][[s]], genus=genera[i],
                                                     height=0, min=0, max=0, Af2Ar=0, Ar2Af=1, Vic=0, ExtAf=0, ExtAr=0)
    vicariance_SIM_nodes_dfs_cons[[i]][[s]] <- data.frame(node=vicariance_SIM_nodes_list_cons[[i]][[s]], genus=genera[i],
                                                          height=0, min=0, max=0, Af2Ar=0, Ar2Af=0, Vic=1, ExtAf=0, ExtAr=0)
    extirpfromAf_SIM_nodes_dfs_cons[[i]][[s]] <- data.frame(node=extirpfromAf_SIM_nodes_list_cons[[i]][[s]], genus=genera[i],
                                                            height=0, min=0, max=0, Af2Ar=0, Ar2Af=0, Vic=0, ExtAf=1, ExtAr=0)
    extirpfromAr_SIM_nodes_dfs_cons[[i]][[s]] <- data.frame(node=extirpfromAr_SIM_nodes_list_cons[[i]][[s]], genus=genera[i],
                                                            height=0, min=0, max=0, Af2Ar=0, Ar2Af=0, Vic=0, ExtAf=0, ExtAr=1)
  }
}

# We want a list with one data frame per genus per simulation, including all types of events.
sim_event_consensus <- vector("list", length(genera))
names(sim_event_consensus) <- genera
for (i in 1:length(genera)){
  sim_event_consensus[[i]] <- vector("list", nsim)
}

for (i in 1:length(genera)){
  for (s in 1:nsim){
    sim_event_consensus[[i]][[s]] <- rbind(Af2Ar_SIM_nodes_dfs_cons[[i]][[s]], Ar2Af_SIM_nodes_dfs_cons[[i]][[s]],
                                           vicariance_SIM_nodes_dfs_cons[[i]][[s]], extirpfromAf_SIM_nodes_dfs_cons[[i]][[s]],
                                           extirpfromAr_SIM_nodes_dfs_cons[[i]][[s]])
    
  }
}

# Now we already have one data frame per genus per sim with the transition nodes. We can remove the zeros:
for (i in 1:length(genera)){
  for (s in 1:nsim){
    sim_event_consensus[[i]][[s]] <- sim_event_consensus[[i]][[s]][sim_event_consensus[[i]][[s]]$node!=0,]
  }
}

# Now we need to include information about the divergence times: the height of the nodes with transitions.
##### INCLUDE INFO ON DIVERGENCE TIMES FOR SIMS #####

# Divergence times and 95% HPD (min and max) for each node
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (nrow(sim_event_consensus[[i]][[s]]) > 0){
      for (j in 1:length(sim_event_consensus[[i]][[s]]$node)){
        
        # Vicariance nodes:
        if (sim_event_consensus[[i]][[s]][j, "Vic"] == 1){
          sim_event_consensus[[i]][[s]]$min[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[sim_event_consensus[[i]][[s]]$node[j]]][1]
          sim_event_consensus[[i]][[s]]$max[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[sim_event_consensus[[i]][[s]]$node[j]]][2]
          sim_event_consensus[[i]][[s]]$height[j] <- as_tibble(consensus_genera[[i]])$height_median[[sim_event_consensus[[i]][[s]]$node[j]]]
        }
        
        # Dispersal and extirpation nodes:
        if (sim_event_consensus[[i]][[s]][j, "Af2Ar"] == 1 | sim_event_consensus[[i]][[s]][j, "Ar2Af"] == 1 | sim_event_consensus[[i]][[s]][j, "ExtAf"] == 1 | sim_event_consensus[[i]][[s]][j, "ExtAr"] == 1){
          dat <- as_tibble(consensus_genera[[i]])
          # Min is the age of the descendent node
          desc_age <- dat$height[dat$node == sim_event_consensus[[i]][[s]]$node[j]]
          sim_event_consensus[[i]][[s]]$min[j] <- desc_age
          
          # Max is the age of the parental node
          parent_age <- dat$height[dat$parent[dat$node == sim_event_consensus[[i]][[s]]$node[j]]]
          sim_event_consensus[[i]][[s]]$max[j] <- parent_age
          
          # let's set the height as the midpoint between the parental and the descendent nodes
          sim_event_consensus[[i]][[s]]$height[j] <- mean(c(desc_age, parent_age))
        }
        
      }
    }
  }
}
saveRDS(sim_event_consensus, "objects/sim_event_cons.rds")

##### EVENTS PER MILLION YEAR FOR EACH GENUS AND SIMULATION (CONSENSUS) #####
myr <- 60
sim_event_myr_list_cons <- vector("list", length(genera))
names(sim_event_myr_list_cons) <- genera
for (i in 1:length(genera)){
  sim_event_myr_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    sim_event_myr_list_cons[[i]][[s]] <- data.frame(Ma=c(1:60), N=0)
  }
}

# We also can create a list of vectors representing intervals of 1 Ma.
ma <- vector("list", myr)
for(i in 1:myr){
  ma[[i]] <- c(i-1, i)
}

# Counting the maximum number of transition in each million year for each genus and each simulation.
sim_event_myr_list_cons0 <- sim_event_myr_list_cons
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (nrow(sim_event_consensus[[i]][[s]]) > 0){
      for (j in 1:nrow(sim_event_consensus[[i]][[s]])){
        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            sim_event_myr_list_cons[[i]][[s]][k,]$N <- sim_event_myr_list_cons[[i]][[s]][k,]$N + 1
          }
        }
      }
    }
  }
}


# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
sim_event_all_cons <- vector("list", nsim)
for (s in 1:nsim){
  sim_event_all_cons[[s]] <- data.frame(Ma=c(1:myr), N=0)
}

for (s in 1:nsim){
  for (i in 1:length(genera)){
    sim_event_all_cons[[s]]$N <- sim_event_all_cons[[s]]$N + sim_event_myr_list_cons[[i]][[s]]$N
  }
}

##### AFRICA TO ARABIA PER MILLION YEAR FOR EACH GENUS AND SIMULATION (CONSENSUS) #####
Af2Ar_sim_myr_list_cons <- vector("list", length(genera))
names(Af2Ar_sim_myr_list_cons) <- genera
for (i in 1:length(genera)){
  Af2Ar_sim_myr_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    Af2Ar_sim_myr_list_cons[[i]][[s]] <- data.frame(Ma=c(1:60), N=0)
  }
}


# Counting the maximum number of Af -> Ar dispersals in each million year for each genus.
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (nrow(sim_event_consensus[[i]][[s]]) > 0){
      for (j in 1:nrow(sim_event_consensus[[i]][[s]])){
        if (sim_event_consensus[[i]][[s]]$Af2Ar[j] > 0){
          vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
          for (k in 1:length(ma)){
            if (Overlap(na.omit(vi), ma[[k]]) != 0){
              Af2Ar_sim_myr_list_cons[[i]][[s]][k,]$N <- Af2Ar_sim_myr_list_cons[[i]][[s]][k,]$N + 1
            }
          }
        }
      }
    }
  }
}

# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
Af2Ar_sim_all_cons <- vector("list", nsim)
for (s in 1:nsim){
  Af2Ar_sim_all_cons[[s]] <- data.frame(Ma=c(1:myr), N=0)
}

for (s in 1:nsim){
  for (i in 1:length(genera)){
    Af2Ar_sim_all_cons[[s]]$N <- Af2Ar_sim_all_cons[[s]]$N + Af2Ar_sim_myr_list_cons[[i]][[s]]$N
  }
}

##### ARABIA TO AFRICA PER MILLION YEAR FOR EACH GENUS AND SIMULATION (CONSENSUS) #####
Ar2Af_sim_myr_list_cons <- vector("list", length(genera))
names(Ar2Af_sim_myr_list_cons) <- genera
for (i in 1:length(genera)){
  Ar2Af_sim_myr_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    Ar2Af_sim_myr_list_cons[[i]][[s]] <- data.frame(Ma=c(1:60), N=0)
  }
}


# Counting the maximum number of Ar -> Af dispersals in each million year for each genus.
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (nrow(sim_event_consensus[[i]][[s]]) > 0){
      for (j in 1:nrow(sim_event_consensus[[i]][[s]])){
        if (sim_event_consensus[[i]][[s]]$Ar2Af[j] > 0){
          vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
          for (k in 1:length(ma)){
            if (Overlap(na.omit(vi), ma[[k]]) != 0){
              Ar2Af_sim_myr_list_cons[[i]][[s]][k,]$N <- Ar2Af_sim_myr_list_cons[[i]][[s]][k,]$N + 1
            }
          }
        }
      }
    }
  }
}

# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
Ar2Af_sim_all_cons <- vector("list", nsim)
for (s in 1:nsim){
  Ar2Af_sim_all_cons[[s]] <- data.frame(Ma=c(1:myr), N=0)
}

for (s in 1:nsim){
  for (i in 1:length(genera)){
    Ar2Af_sim_all_cons[[s]]$N <- Ar2Af_sim_all_cons[[s]]$N + Ar2Af_sim_myr_list_cons[[i]][[s]]$N
  }
}

##### VICARIANCE PER MILLION YEAR FOR EACH GENUS AND SIMULATION (CONSENSUS) #####
vicariance_sim_myr_list_cons <- vector("list", length(genera))
names(vicariance_sim_myr_list_cons) <- genera
for (i in 1:length(genera)){
  vicariance_sim_myr_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    vicariance_sim_myr_list_cons[[i]][[s]] <- data.frame(Ma=c(1:60), N=0)
  }
}


# Counting the maximum number of vicariance in each million year for each genus.
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (nrow(sim_event_consensus[[i]][[s]]) > 0){
      for (j in 1:nrow(sim_event_consensus[[i]][[s]])){
        if (sim_event_consensus[[i]][[s]]$Vic[j] > 0){
          vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
          for (k in 1:length(ma)){
            if (Overlap(na.omit(vi), ma[[k]]) != 0){
              vicariance_sim_myr_list_cons[[i]][[s]][k,]$N <- vicariance_sim_myr_list_cons[[i]][[s]][k,]$N + 1
            }
          }
        }
      }
    }
  }
}

# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
vicariance_sim_all_cons <- vector("list", nsim)
for (s in 1:nsim){
  vicariance_sim_all_cons[[s]] <- data.frame(Ma=c(1:myr), N=0)
}

for (s in 1:nsim){
  for (i in 1:length(genera)){
    vicariance_sim_all_cons[[s]]$N <- vicariance_sim_all_cons[[s]]$N + vicariance_sim_myr_list_cons[[i]][[s]]$N
  }
}

##### EXTIRPATION FROM AFRICA PER MILLION YEAR FOR EACH GENUS AND SIMULATION (CONSENSUS) #####
extirp_from_Af_sim_myr_list_cons <- vector("list", length(genera))
names(extirp_from_Af_sim_myr_list_cons) <- genera
for (i in 1:length(genera)){
  extirp_from_Af_sim_myr_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    extirp_from_Af_sim_myr_list_cons[[i]][[s]] <- data.frame(Ma=c(1:60), N=0)
  }
}


# Counting the maximum number of extirpation from Africa in each million year for each genus.
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (nrow(sim_event_consensus[[i]][[s]]) > 0){
      for (j in 1:nrow(sim_event_consensus[[i]][[s]])){
        if (sim_event_consensus[[i]][[s]]$ExtAf[j] > 0){
          vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
          for (k in 1:length(ma)){
            if (Overlap(na.omit(vi), ma[[k]]) != 0){
              extirp_from_Af_sim_myr_list_cons[[i]][[s]][k,]$N <- extirp_from_Af_sim_myr_list_cons[[i]][[s]][k,]$N + 1
            }
          }
        }
      }
    }
  }
}

# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
extirp_from_Af_sim_all_cons <- vector("list", nsim)
for (s in 1:nsim){
  extirp_from_Af_sim_all_cons[[s]] <- data.frame(Ma=c(1:myr), N=0)
}

for (s in 1:nsim){
  for (i in 1:length(genera)){
    extirp_from_Af_sim_all_cons[[s]]$N <- extirp_from_Af_sim_all_cons[[s]]$N + extirp_from_Af_sim_myr_list_cons[[i]][[s]]$N
  }
}

##### EXTIRPATION FROM ARABIA PER MILLION YEAR FOR EACH GENUS AND SIMULATION (CONSENSUS) #####
extirp_from_Ar_sim_myr_list_cons <- vector("list", length(genera))
names(extirp_from_Ar_sim_myr_list_cons) <- genera
for (i in 1:length(genera)){
  extirp_from_Ar_sim_myr_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    extirp_from_Ar_sim_myr_list_cons[[i]][[s]] <- data.frame(Ma=c(1:60), N=0)
  }
}


# Counting the maximum number of Extirpation from Arabia in each million year for each genus.
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (nrow(sim_event_consensus[[i]][[s]]) > 0){
      for (j in 1:nrow(sim_event_consensus[[i]][[s]])){
        if (sim_event_consensus[[i]][[s]]$ExtAr[j] > 0){
          vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
          for (k in 1:length(ma)){
            if (Overlap(na.omit(vi), ma[[k]]) != 0){
              extirp_from_Ar_sim_myr_list_cons[[i]][[s]][k,]$N <- extirp_from_Ar_sim_myr_list_cons[[i]][[s]][k,]$N + 1
            }
          }
        }
      }
    }
  }
}

# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
extirp_from_Ar_sim_all_cons <- vector("list", nsim)
for (s in 1:nsim){
  extirp_from_Ar_sim_all_cons[[s]] <- data.frame(Ma=c(1:myr), N=0)
}

for (s in 1:nsim){
  for (i in 1:length(genera)){
    extirp_from_Ar_sim_all_cons[[s]]$N <- extirp_from_Ar_sim_all_cons[[s]]$N + extirp_from_Ar_sim_myr_list_cons[[i]][[s]]$N
  }
}


##### Save event objects #####
saveRDS(sim_event_all_cons, "objects/sim_event_all_cons.rds")
saveRDS(event_all_cons, "objects/event_all_cons.rds")
saveRDS(Af2Ar_sim_all_cons, "objects/Af2Ar_sim_all_cons.rds")
saveRDS(Af2Ar_all_cons, "objects/Af2Ar_all_cons.rds")
saveRDS(Ar2Af_sim_all_cons, "objects/Ar2Af_sim_all_cons.rds")
saveRDS(Ar2Af_all_cons, "objects/Ar2Af_all_cons.rds")
saveRDS(vicariance_sim_all_cons, "objects/vicariance_sim_all_cons.rds")
saveRDS(vicariance_all_cons, "objects/vicariance_all_cons.rds")
saveRDS(extirp_from_Af_sim_all_cons, "objects/extirp_from_Af_sim_all_cons.rds")
saveRDS(extirp_from_Af_all_cons, "objects/extirp_from_Af_all_cons.rds")
saveRDS(extirp_from_Ar_sim_all_cons, "objects/extirp_from_Ar_sim_all_cons.rds")
saveRDS(extirp_from_Ar_all_cons, "objects/extirp_from_Ar_all_cons.rds")

# THE END ----

