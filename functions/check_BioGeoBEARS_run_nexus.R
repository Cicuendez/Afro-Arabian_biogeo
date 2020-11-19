check_BioGeoBEARS_run_nexus <- function (inputs, allow_huge_ranges = FALSE, allow_null_range_tips = NULL) 
{
  runjunk = "\n\tinputs = BioGeoBEARS_run_object\n\tallow_huge_ranges=FALSE\n\tallow_null_range_tips=FALSE\n\t"
  if (is.null(inputs$allow_null_tips)) {
    inputs$allow_null_tips = FALSE
  }
  if (is.null(allow_null_range_tips) == FALSE) {
    if (allow_null_range_tips == FALSE) {
      inputs$allow_null_tips = FALSE
    }
    if (allow_null_range_tips == TRUE) {
      inputs$allow_null_tips = TRUE
    }
  }
  if (inputs$use_detection_model == TRUE) {
    if ((is.null(inputs$detects_fn) == TRUE) || (is.null(inputs$controls_fn) == 
                                                 TRUE)) {
      stoptxt = paste0("STOP ERROR in check_BioGeoBEARS_run(): You set inputs$use_detection_model == TRUE, but inputs$detects_fn and/or inputs$controls_fn is NULL. Using the detection model requires a detections counts file (inputs$detects_fn) and a taphonomic control counts file (inputs$controls_fn).\n\nBackground on finding files in R:\n\nCheck your inputs and your working directory (wd). Use commands like:\n - 'getwd()' to get your current R working directory (wd)\n - '?setwd' to see how to set your R working directory\n - 'list.files()' to list the files in your current R working directory\n - and use 'system('open .')' to open your file browsing program from the R command line (this works on macs, at least).")
      cat("\n\n")
      cat(stoptxt)
      cat("\n\n")
      stop(stoptxt)
    }
  }
  tmptr = check_trfn_nexus(trfn = inputs$trfn)
  if (exists("tmptr") == FALSE) {
    stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: no readable Newick tree at:\n", 
                    inputs$trfn, "\n", sep = "")
    cat(stoptxt)
    stop(stoptxt)
  }
  tipnames = tmptr$tip.label
  uniq_tipnames = unique(tipnames)
  if (length(uniq_tipnames) != length(tipnames)) {
    stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: your tree has non-unique tipnames. Make all tipnames unique in both your tree file and geography file.  Current tipnames are listed below:\n\n", 
                    sep = "")
    for (i in 1:length(uniq_tipnames)) {
      TF = uniq_tipnames[i] %in% tipnames
      if (sum(TF) > 1) {
        cat(uniq_tipnames[i])
        cat("\n")
      }
    }
    stop(stoptxt)
  }
  TF = grepl(" ", tipnames)
  if (sum(TF) > 0) {
    stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: your tree has tipnames with spaces. Take these out of your Newick file and re-run. Current tipnames are listed below.\n\n", 
                    sep = "")
    print(tipnames[TF])
    stop(stoptxt)
  }
  TF = grepl("'", tipnames)
  if (sum(TF) > 0) {
    stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: your tree has tipnames with apostrophes ('). Take these out of your Newick file and re-run. Current tipnames are listed below.\n\n", 
                    sep = "")
    print(tipnames[TF])
    stop(stoptxt)
  }
  blren_equal_below_0_TF = tmptr$edge.length <= 0
  if (sum(blren_equal_below_0_TF) > 0) {
    tmptr_table = prt(tmptr, printflag = FALSE)
    rows_w_BL0_TF = tmptr_table$edge.length <= 0
    rows_w_BL0_TF[is.na(rows_w_BL0_TF)] = FALSE
    nodenums = tmptr_table$node[rows_w_BL0_TF]
    branchlengths = tmptr_table$edge.length[rows_w_BL0_TF]
    edge_nums = tmptr_table$parent_br[rows_w_BL0_TF]
    tmptable = cbind(nodenums, branchlengths, edge_nums)
    tmptable = as.data.frame(tmptable, stringsAsFactors = FALSE)
    tmptxt = paste0(nodenums, collapse = ",")
    tmptxt2 = paste0(edge_nums, collapse = ",")
    stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in check_BioGeoBEARS_run(): the input tree has branchlengths <= 0, at these nodes:\n\n", 
                    tmptxt, "\n\n...And, at these edge numbers:\n\n", 
                    tmptxt2, "\n\n")
    cat(stoptxt)
    print(tmptable)
    cat("\nThis can sometimes happen in e.g. MCC (majority clade consensus) trees output by BEAST's TreeAnnotator.\nYou must fix the Newick file. See ?check_BioGeoBEARS_run, and PhyloWiki, for comments. One option is to try impose_min_brlen() and then write the tree to a file.\n", 
        sep = "")
    cat(stoptxt)
    stop(stoptxt)
  }
  if (is.binary(tmptr) == FALSE) {
    stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: Your tree not bifurcating, i.e. is.binary(tmptr) returns FALSE.\n", 
                    "\nYou must fix the Newick file. APE's multi2di() function is an option.  See ?check_BioGeoBEARS_run for comments.\n", 
                    sep = "")
    cat(stoptxt)
    stop(stoptxt)
  }
  if (is.rooted(tmptr) == FALSE) {
    stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: Your tree is not rooted, i.e. is.rooted(tmptr) returns FALSE.\n", 
                    "\nYou must fix the Newick file. APE's root() function is an option, or e.g. re-rooting by hand in FigTree.  However, you will want to make sure that all your tips still come up to the present (assuming you have a typical molecular tree, i.e. no fossils). \n", 
                    sep = "")
    cat(stoptxt)
    stop(stoptxt)
  }
  if (inputs$use_detection_model == FALSE) {
    tipranges = getranges_from_LagrangePHYLIP(inputs$geogfn)
  }
  if (inputs$use_detection_model == TRUE) {
    if (is.character(inputs$detects_fn)) {
      inputs$detects_df = read_detections(inputs$detects_fn, 
                                          OTUnames = NULL, areanames = NULL, tmpskip = 0, 
                                          phy = tmptr)
      tmp_blah = as.matrix(inputs$detects_df)
      tmp_blah[isblank_TF(tmp_blah)] = 0
      tmp_blah[tmp_blah > 0] = 1
      tmp_input = adf2(tmp_blah)
      tipranges_object = define_tipranges_object(tmpdf = tmp_input)
      tipranges_object@df = adf2(tipranges_object@df)
      rownames(tipranges_object@df) = rownames(tmp_blah)
      tipranges = tipranges_object
    }
    else {
      txt = paste0("STOP ERROR in bears_optim_run(): you set use_detection_model=TRUE, so an input file is required for input 'detects_fn'. This is required to set up the tipranges internally.")
      cat("\n\n")
      cat(txt)
      cat("\n\n")
      stop(txt)
    }
  }
  if (exists("tipranges") == FALSE) {
    stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: no readable tipranges text file at inputs$geogfn = '", 
                    inputs$geogfn, "', nor at inputs$detects_fn = '", 
                    inputs$detects_fn, "'. Check your inputs and your working directory (wd). Use commands like:\n - 'getwd()' to get your current R working directory (wd)\n - '?setwd' to see how to set your R working directory\n - 'list.files()' to list the files in your current R working directory\n - and use 'system('open .')' to open your file browsing program from the R command line (this works on macs, at least).\n", 
                    sep = "")
    cat(stoptxt)
    stop(stoptxt)
  }
  tipranges_colnames_TF = is.na(colnames(tipranges@df))
  if (sum(tipranges_colnames_TF) > 0) {
    catstr = paste(colnames(tipranges@df), collapse = " ", 
                   sep = "")
    stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: tipranges area names (columns) have NAs:\n", 
                    catstr, "\nThis probably means your input tipranges file is missing ", 
                    sum(tipranges_colnames_TF), " areaname(s).\n", sep = "")
    cat(stoptxt)
    moref(inputs$geogfn)
    stop(stoptxt)
  }
  tipnames = sort(tmptr$tip.label)
  geogtaxa = sort(rownames(tipranges@df))
  match_mismatch_table_df = compare_two_name_lists(names1 = geogtaxa, 
                                                   names2 = tipnames, listdesc1 = "geography file", listdesc2 = "Newick tree file", 
                                                   list_or_file_txt = "files")
  tipnames_in_geogfile_TF = tipnames %in% geogtaxa
  if ((sum(tipnames_in_geogfile_TF) == length(tipnames_in_geogfile_TF)) == 
      FALSE) {
    tipnames_NOT_in_geogfile_TF = tipnames_in_geogfile_TF == 
      FALSE
    numtips_not_in_geogfn = sum(tipnames_NOT_in_geogfile_TF)
    tips_not_in_geogfile = tipnames[tipnames_NOT_in_geogfile_TF]
    tips_not_in_geogfile_txt = paste(tips_not_in_geogfile, 
                                     collapse = ", ", sep = "")
    if (length(geogtaxa) == length(tipnames)) {
      match_TF = geogtaxa == tipnames
      match_TF_txt = paste(match_TF, collapse = " ", sep = "")
    }
    else {
      match_TF_txt = paste("Cannot display TRUE/FALSE, as length(tipnames)=", 
                           length(tipnames), " & length(geogtaxa)=", length(geogtaxa), 
                           ".\n", sep = "")
    }
    stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in input files: Your geography file has a list of species, and your newick file has a list of species.  These lists have to match *exactly*.  This error message is saying that you have one or more species names that are missing, mis-spelled, differing due to underscores vs. spaces, contain inappropriate characters like apostrophes, etc. Error message: ", 
                    numtips_not_in_geogfn, " tree tip(s) are not in the geographic ranges file.  These are listed.\n", 
                    tips_not_in_geogfile_txt, "\n", "TRUE/FALSE between sort(tmptr$tip.label)==sort(rownames(tipranges@df)):\n", 
                    match_TF_txt, "\n", sep = "")
    cat(stoptxt)
    cat("\n\nList of matches and mis-matches between the geography file and Newick tree file:\n\n")
    print(match_mismatch_table_df)
    cat("\n\nNOTE: Computers are very literal. Things like space vs. '_' vs. '.' matter, as do tabs versus spaces, and multiple spaces or tabs versus single ones. Also, you should NOT have ' or similar characters in either your Newick tree file or your geography file.\n\nTo edit these files, view them in a REAL plain-text editor, not in Word or whatever.\n\nInformation about plain-text editors: http://phylo.wikidot.com/biogeobears#texteditors .\n\nInformation about BioGeoBEARS file formats, with examples: http://phylo.wikidot.com/biogeobears#links_to_files .\n\n")
    stop(stoptxt)
  }
  tipnames = sort(tmptr$tip.label)
  geogtaxa = sort(rownames(tipranges@df))
  geogtaxa_in_treetips_TF = geogtaxa %in% tipnames
  if ((sum(geogtaxa_in_treetips_TF) == length(geogtaxa_in_treetips_TF)) == 
      FALSE) {
    geogtaxa_NOT_in_treetips_TF = geogtaxa_in_treetips_TF == 
      FALSE
    num_geogtaxa_not_in_treetips = sum(geogtaxa_NOT_in_treetips_TF)
    geogtaxa_not_in_treetips = geogtaxa[geogtaxa_NOT_in_treetips_TF]
    geogtaxa_not_in_treetips_txt = paste(geogtaxa_not_in_treetips, 
                                         collapse = ", ", sep = "")
    if (length(geogtaxa) == length(tipnames)) {
      match_TF = geogtaxa == tipnames
      match_TF_txt = paste(match_TF, collapse = " ", sep = "")
    }
    else {
      match_TF_txt = paste("Cannot display TRUE/FALSE, as length(tipnames)=", 
                           length(tipnames), " & length(geogtaxa)=", length(geogtaxa), 
                           ".\n", sep = "")
    }
    stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: ", 
                    num_geogtaxa_not_in_treetips, " taxa in the geography file are not in the tree tips.  These are listed.  This error message is saying you have one or more species names that are missing, mis-spelled, differing due to underscores vs. spaces, contain inappropriate characters like apostrophes, etc.\n", 
                    geogtaxa_not_in_treetips_txt, "\n", "TRUE/FALSE between sort(tmptr$tip.label)==sort(rownames(tipranges@df)):\n", 
                    match_TF_txt, "\n", sep = "")
    cat("\n\nList of matches and mis-matches between the geography file and Newick tree file:\n\n")
    print(match_mismatch_table_df)
    cat(stoptxt)
    cat("\n\nNOTE: Computers are very literal. Things like space vs. '_' vs. '.' matter, as do tabs versus spaces, and multiple spaces or tabs versus single ones. Also, you should NOT have ' or similar characters in either your Newick tree file or your geography file.\n\nTo edit these files, view them in a REAL plain-text editor, not in Word or whatever.\n\nInformation about plain-text editors: http://phylo.wikidot.com/biogeobears#texteditors .\n\nInformation about BioGeoBEARS file formats, with examples: http://phylo.wikidot.com/biogeobears#links_to_files .\n\n")
    stop(stoptxt)
  }
  tmp_tipranges = tipranges@df
  tmp_tipranges[tmp_tipranges == "?"] = 0
  tmp_tipranges = dfnums_to_numeric(tmp_tipranges)
  if (!is.na(inputs$max_range_size)) {
    max_tipsize = max(rowSums(tmp_tipranges))
    if (max_tipsize > inputs$max_range_size) {
      tips_too_big_TF = rowSums(tmp_tipranges) > inputs$max_range_size
      tipranges_too_big = tipranges@df[tips_too_big_TF, 
                                       ]
      stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: max_tipsize=", 
                      max_tipsize, " > inputs$max_range_size=", inputs$max_range_size, 
                      ". Examples:\n", sep = "")
      cat(stoptxt)
      print(tipranges_too_big)
      cat("\n")
      stop(stoptxt)
    }
  }
  rangesize_is_ZERO_TF = rowSums(tmp_tipranges) == 0
  ranges_have_NAs_TF = is.na(unlist(tipranges@df))
  if ((inputs$allow_null_tips == FALSE) && ((sum(rangesize_is_ZERO_TF) > 
                                             0) || (sum(ranges_have_NAs_TF) > 0))) {
    stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: your tipranges file has NAs and/or tips with rangesize 0!\n See tipranges printed below:\n\n", 
                    sep = "")
    cat(stoptxt)
    print(tipranges@df)
    stoptxt2 = paste("\nAnd your problematic tipranges rows are:\n", 
                     sep = "")
    cat(stoptxt2)
    if (sum(rangesize_is_ZERO_TF) > 0) {
      TF = rangesize_is_ZERO_TF
      print(tipranges@df[TF, ])
    }
    if (sum(ranges_have_NAs_TF) > 0) {
      TF = rep(FALSE, nrow(tipranges@df))
      for (i in 1:nrow(tipranges@df)) {
        tmprow = tipranges@df[i, ]
        if (sum(is.na(tmprow)) > 0) {
          TF[i] = TRUE
        }
      }
      print(tipranges@df[TF, ])
    }
    stop(stoptxt)
  }
  tmp_numareas = ncol(tipranges@df)
  if (!is.na(inputs$max_range_size)) {
    max_tipsize = max(rowSums(tmp_tipranges))
  }
  else {
    max_tipsize = tmp_numareas
  }
  tmp_numstates1 = length(inputs$states_list)
  if (length(tmp_numstates1) == 0) {
    tmp_numstates1 = numstates_from_numareas(numareas = tmp_numareas, 
                                             maxareas = max_tipsize, include_null_range = TRUE)
  }
  if (tmp_numstates1 > 2500) {
    stoptxt = paste("\ncheck_BioGeoBEARS_run() says: Your setup has ", 
                    tmp_numstates1, " states (# states = # combinations of geographic ranges). This will be veerry slow. \n", 
                    "In check_BioGeoBEARS_run(), set allow_huge_ranges=TRUE to proceed, but you probably shouldn't bother.  See e.g. ?numstates_from_numareas.\n", 
                    sep = "")
    cat(stoptxt)
    if (allow_huge_ranges == TRUE) {
      pass = 1
    }
    else {
      stop(stoptxt)
    }
  }
  numrows = nrow(inputs$BioGeoBEARS_model_object@params_table)
  list_of_is = NULL
  for (i in 1:numrows) {
    if (inputs$BioGeoBEARS_model_object@params_table[i, "type"] == 
        "free") {
      TF1 = inputs$BioGeoBEARS_model_object@params_table[i, 
                                                         "init"] >= inputs$BioGeoBEARS_model_object@params_table[i, 
                                                                                                                 "min"]
      TF2 = inputs$BioGeoBEARS_model_object@params_table[i, 
                                                         "init"] <= inputs$BioGeoBEARS_model_object@params_table[i, 
                                                                                                                 "max"]
      if ((TF1 + TF2) == 2) {
        (next)()
      }
      else {
        list_of_is = c(list_of_is, i)
      }
    }
    if (length(list_of_is) > 0) {
      stoptxt = paste0("check_BioGeoBEARS_run() says: STOP ERROR: ", 
                       length(list_of_is), " row(s) of the BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table are set to be 'free', but they have starting values ('init') outside of the specified min/max.\n\nFix manually, or run fix_BioGeoBEARS_params_minmax(). Printing these rows to screen....\n\n")
      cat("\n\n")
      cat(stoptxt)
      print(inputs$BioGeoBEARS_model_object@params_table[list_of_is, 
                                                         ])
      stop(stoptxt)
    }
  }
  TF1 = inputs$BioGeoBEARS_model_object@params_table["x", "type"] == 
    "free"
  TF2 = inputs$BioGeoBEARS_model_object@params_table["x", "init"] != 
    0
  TF3 = inputs$BioGeoBEARS_model_object@params_table["x", "est"] != 
    0
  if (TF1 || TF2 || TF3) {
    if (is.character(inputs$distsfn) == FALSE) {
      stoptxt = paste("\nFATAL ERROR: Your 'x' parameter is free or nonzero, but you have input\n", 
                      "no distances file.\n\n", sep = "")
      cat(stoptxt)
      stop(stoptxt)
    }
  }
  TF1 = inputs$BioGeoBEARS_model_object@params_table["w", "type"] == 
    "free"
  TF2 = inputs$BioGeoBEARS_model_object@params_table["w", "init"] != 
    1
  TF3 = inputs$BioGeoBEARS_model_object@params_table["w", "est"] != 
    1
  if (TF1 || TF2 || TF3) {
    if (is.character(inputs$dispersal_multipliers_fn) == 
        FALSE) {
      stoptxt = paste("\nFATAL ERROR: Your 'w' parameter is not set to '1', or is free, but you have input\n", 
                      "no manual dispersal multipliers file.\n\n", 
                      sep = "")
      cat(stoptxt)
      stop(stoptxt)
    }
  }
  TF1 = inputs$BioGeoBEARS_model_object@params_table["n", "type"] == 
    "free"
  TF2 = inputs$BioGeoBEARS_model_object@params_table["n", "init"] != 
    0
  TF3 = inputs$BioGeoBEARS_model_object@params_table["n", "est"] != 
    0
  if (TF1 || TF2 || TF3) {
    if (is.character(inputs$envdistsfn) == FALSE) {
      stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: Your 'w' parameter is not set to '1', or is free, but you have input\n", 
                      "no environmental distances file.\n\n", sep = "")
      cat(stoptxt)
      stop(stoptxt)
    }
  }
  if (is.character(inputs$distsfn)) {
    dims = dim(inputs$list_of_distances_mats[[1]])
    if (dims[1] != dims[2]) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the distance matrix is not square! Instead, yours has ", 
                       dims[1], " rows x ", dims[2], " columns. Check your input file. Printing the matrix below.\n")
      cat(stoptxt)
      print(inputs$list_of_distances_mats[[1]])
      cat("\n\n")
      stop(stoptxt)
    }
    areanames_from_distsfn = colnames(inputs$list_of_distances_mats[[1]])
    areanames = names(tipranges@df)
    if (length(areanames) != length(areanames_from_distsfn)) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas list in your geography file and in your distances file are not the same length!\n")
      cat(stoptxt)
      cat("\nPrinting the two lists below.\n\n")
      cat("Areas in geography file (names(tipranges@df)):\n")
      print(areanames)
      cat("Areas in distances file (colnames(inputs$list_of_distances_mats[[1]]) ):\n")
      print(areanames_from_distsfn)
      cat("\n\n")
      stop(stoptxt)
    }
    TFs = areanames == areanames_from_distsfn
    if (all(TFs) == FALSE) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the area names in your geography file and in your distances file do not match! Either they are in a different order, or they are different. They must be the same, and be in the same order. In these matrices, the areas must be in the same order in BOTH rows and columns, but only the columns get area labels. For example files, see: http://phylo.wikidot.com/biogeobears#links_to_files . For advice on using a *PLAIN-TEXT* editor (not e.g. Word) to edit text files, see: http://phylo.wikidot.com/biogeobears#texteditors \n")
      cat(stoptxt)
      cat("\nPrinting the two lists below.\n\n")
      cat("Areas in geography file (names(tipranges@df)):\n")
      print(areanames)
      cat("Areas in distances file (colnames(inputs$list_of_distances_mats[[1]]) ):\n")
      print(areanames_from_distsfn)
      cat("\n\n")
      stop(stoptxt)
    }
  }
  if (is.character(inputs$envdistsfn)) {
    dims = dim(inputs$list_of_envdistances_mats[[1]])
    if (dims[1] != dims[2]) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the envdistances matrix is not square! Instead, yours has ", 
                       dims[1], " rows x ", dims[2], " columns. Check your input file. Printing the matrix below.\n")
      cat(stoptxt)
      print(inputs$list_of_envdistances_mats[[1]])
      cat("\n\n")
      stop(stoptxt)
    }
    areanames_from_envdistsfn = colnames(inputs$list_of_envdistances_mats[[1]])
    areanames = names(tipranges@df)
    if (length(areanames) != length(areanames_from_envdistsfn)) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas list in your geography file and in your environmental distances file are not the same length!\n")
      cat(stoptxt)
      cat("\nPrinting the two lists below.\n\n")
      cat("Areas in geography file (names(tipranges@df)):\n")
      print(areanames)
      cat("Areas in environmental distances file (colnames(inputs$list_of_envdistances_mats[[1]]) ):\n")
      print(areanames_from_envdistsfn)
      cat("\n\n")
      stop(stoptxt)
    }
    TFs = areanames == areanames_from_envdistsfn
    if (all(TFs) == FALSE) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the area names in your geography file and in your environmental distances file do not match! Either they are in a different order, or they are different. They must be the same, and be in the same order. In these matrices, the areas must be in the same order in BOTH rows and columns, but only the columns get area labels. For example files, see: http://phylo.wikidot.com/biogeobears#links_to_files . For advice on using a *PLAIN-TEXT* editor (not e.g. Word) to edit text files, see: http://phylo.wikidot.com/biogeobears#texteditors \n")
      cat(stoptxt)
      cat("\nPrinting the two lists below.\n\n")
      cat("Areas in geography file (names(tipranges@df)):\n")
      print(areanames)
      cat("Areas in environmental distances file (colnames(inputs$list_of_envdistances_mats[[1]]) ):\n")
      print(areanames_from_envdistsfn)
      cat("\n\n")
      stop(stoptxt)
    }
  }
  if (is.character(inputs$dispersal_multipliers_fn)) {
    dims = dim(inputs$list_of_dispersal_multipliers_mats[[1]])
    if (dims[1] != dims[2]) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the dispersal multipliers matrix is not square! Instead, yours has ", 
                       dims[1], " rows x ", dims[2], " columns. Check your input file. Printing the matrix below.\n")
      cat(stoptxt)
      print(inputs$list_of_dispersal_multipliers_mats[[1]])
      cat("\n\n")
      stop(stoptxt)
    }
    areanames_from_dispersal_multipliers_fn = colnames(inputs$list_of_dispersal_multipliers_mats[[1]])
    areanames = names(tipranges@df)
    if (length(areanames) != length(areanames_from_dispersal_multipliers_fn)) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas list in your geography file and in your manual dispersal multipliers file are not the same length!\n")
      cat(stoptxt)
      cat("\nPrinting the two lists below.\n\n")
      cat("Areas in geography file (names(tipranges@df)):\n")
      print(areanames)
      cat("Areas in manual dispersal multipliers file (colnames(inputs$list_of_dispersal_multipliers_mats[[1]]) ):\n")
      print(areanames_from_dispersal_multipliers_fn)
      cat("\n\n")
      stop(stoptxt)
    }
    TFs = areanames == areanames_from_dispersal_multipliers_fn
    if (all(TFs) == FALSE) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the area names in your geography file and in your manual dispersal multipliers file do not match! Either they are in a different order, or they are different. They must be the same, and be in the same order. In these matrices, the areas must be in the same order in BOTH rows and columns, but only the columns get area labels. For example files, see: http://phylo.wikidot.com/biogeobears#links_to_files . For advice on using a *PLAIN-TEXT* editor (not e.g. Word) to edit text files, see: http://phylo.wikidot.com/biogeobears#texteditors \n")
      cat(stoptxt)
      cat("\nPrinting the two lists below.\n\n")
      cat("Areas in geography file (names(tipranges@df)):\n")
      print(areanames)
      cat("Areas in manual dispersal multipliers file (colnames(inputs$list_of_dispersal_multipliers_mats[[1]]) ):\n")
      print(areanames_from_dispersal_multipliers_fn)
      cat("\n\n")
      stop(stoptxt)
    }
  }
  if (is.character(inputs$area_of_areas_fn)) {
    areanames_from_area_of_areas_fn = colnames(inputs$list_of_area_of_areas[[1]])
    areanames = names(tipranges@df)
    if (length(areanames) != length(areanames_from_area_of_areas_fn)) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas list in your geography file and in your area of areas file are not the same length!\n")
      cat(stoptxt)
      cat("\nPrinting the two lists below.\n\n")
      cat("Areas in geography file (names(tipranges@df)):\n")
      print(areanames)
      cat("Areas in area of areas file (colnames(inputs$list_of_area_of_areas[[1]]) ):\n")
      print(areanames_from_area_of_areas_fn)
      cat("\n\n")
      stop(stoptxt)
    }
    TFs = areanames == areanames_from_area_of_areas_fn
    if (all(TFs) == FALSE) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the area names in your geography file and in your area of areas file do not match! Either they are in a different order, or they are different. They must be the same, and be in the same order. In these matrices, the areas must be in the same order in BOTH rows and columns, but only the columns get area labels. For example files, see: http://phylo.wikidot.com/biogeobears#links_to_files . For advice on using a *PLAIN-TEXT* editor (not e.g. Word) to edit text files, see: http://phylo.wikidot.com/biogeobears#texteditors \n")
      cat(stoptxt)
      cat("\nPrinting the two lists below.\n\n")
      cat("Areas in geography file (names(tipranges@df)):\n")
      print(areanames)
      cat("Areas in area of areas file (colnames(inputs$list_of_area_of_areas[[1]]) ):\n")
      print(areanames_from_area_of_areas_fn)
      cat("\n\n")
      stop(stoptxt)
    }
  }
  if (is.character(inputs$areas_allowed_fn)) {
    dims = dim(inputs$list_of_areas_allowed_mats[[1]])
    if (dims[1] != dims[2]) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas allowed matrix is not square! Instead, yours has ", 
                       dims[1], " rows x ", dims[2], " columns. Check your input file. Printing the matrix below.\n")
      cat(stoptxt)
      print(inputs$list_of_areas_allowed_mats[[1]])
      cat("\n\n")
      stop(stoptxt)
    }
    areanames_from_areas_allowed_fn = colnames(inputs$list_of_areas_allowed_mats[[1]])
    areanames = names(tipranges@df)
    if (length(areanames) != length(areanames_from_areas_allowed_fn)) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas list in your geography file and in your areas allowed file are not the same length!\n")
      cat(stoptxt)
      cat("\nPrinting the two lists below.\n\n")
      cat("Areas in geography file (names(tipranges@df)):\n")
      print(areanames)
      cat("Areas in areas allowed file (colnames(inputs$list_of_areas_allowed_mats[[1]]) ):\n")
      print(areanames_from_areas_allowed_fn)
      cat("\n\n")
      stop(stoptxt)
    }
    TFs = areanames == areanames_from_areas_allowed_fn
    if (all(TFs) == FALSE) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the area names in your geography file and in your areas allowed file do not match! Either they are in a different order, or they are different. They must be the same, and be in the same order. In these matrices, the areas must be in the same order in BOTH rows and columns, but only the columns get area labels. For example files, see: http://phylo.wikidot.com/biogeobears#links_to_files . For advice on using a *PLAIN-TEXT* editor (not e.g. Word) to edit text files, see: http://phylo.wikidot.com/biogeobears#texteditors \n")
      cat(stoptxt)
      cat("\nPrinting the two lists below.\n\n")
      cat("Areas in geography file (names(tipranges@df)):\n")
      print(areanames)
      cat("Areas in areas allowed file (colnames(inputs$list_of_areas_allowed_mats[[1]]) ):\n")
      print(areanames_from_areas_allowed_fn)
      cat("\n\n")
      stop(stoptxt)
    }
  }
  if (is.character(inputs$areas_adjacency_fn)) {
    dims = dim(inputs$list_of_areas_adjacency_mats[[1]])
    if (dims[1] != dims[2]) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas adjacency matrix is not square! Instead, yours has ", 
                       dims[1], " rows x ", dims[2], " columns. Check your input file. Printing the matrix below.\n")
      cat(stoptxt)
      print(inputs$list_of_areas_adjacency_mats[[1]])
      cat("\n\n")
      stop(stoptxt)
    }
    areanames_from_areas_adjacency_fn = colnames(inputs$list_of_areas_adjacency_mats[[1]])
    areanames = names(tipranges@df)
    if (length(areanames) != length(areanames_from_areas_adjacency_fn)) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas list in your geography file and in your areas adjacency file are not the same length!\n")
      cat(stoptxt)
      cat("\nPrinting the two lists below.\n\n")
      cat("Areas in geography file (names(tipranges@df)):\n")
      print(areanames)
      cat("Areas in areas adjacency file (colnames(inputs$list_of_areas_adjacency_mats[[1]]) ):\n")
      print(areanames_from_areas_adjacency_fn)
      cat("\n\n")
      stop(stoptxt)
    }
    TFs = areanames == areanames_from_areas_adjacency_fn
    if (all(TFs) == FALSE) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the area names in your geography file and in your areas adjacency file do not match! Either they are in a different order, or they are different. They must be the same, and be in the same order. In these matrices, the areas must be in the same order in BOTH rows and columns, but only the columns get area labels. For example files, see: http://phylo.wikidot.com/biogeobears#links_to_files . For advice on using a *PLAIN-TEXT* editor (not e.g. Word) to edit text files, see: http://phylo.wikidot.com/biogeobears#texteditors \n")
      cat(stoptxt)
      cat("\nPrinting the two lists below.\n\n")
      cat("Areas in geography file (names(tipranges@df)):\n")
      print(areanames)
      cat("Areas in areas adjacency file (colnames(inputs$list_of_areas_adjacency_mats[[1]]) ):\n")
      print(areanames_from_areas_adjacency_fn)
      cat("\n\n")
      stop(stoptxt)
    }
  }
  if (is.character(inputs$timesfn)) {
    timevals = inputs$timeperiods
    if (identical(timevals, sort(timevals)) == FALSE) {
      stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: Your timeperiods are not in order from youngest to oldest. They need to be.\n")
      cat(stoptxt)
      stop(stoptxt)
    }
    trtable = prt(tmptr, printflag = FALSE)
    root_age = max(trtable$time_bp)
    oldest_time = timevals[length(timevals)]
    if (root_age >= oldest_time) {
      stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: The root of your tree is age=", 
                       root_age, ". But your oldest time in your times file is only time=", 
                       oldest_time, ". The oldest time in the times file must be older than the root age.\n\n(Times in the times file represent time-bin bottoms. For example, if your first time is 10 Ma, this means the first time bin stretches from 0-10 Ma. (Ma = mega-annum = millions of years ago.)\n")
      cat(stoptxt)
      stop(stoptxt)
    }
    if (is.character(inputs$distsfn)) {
      if (length(inputs$list_of_distances_mats) < length(inputs$timeperiods)) {
        stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: fewer distances matrices than timeperiods.\n", 
                        "length(inputs$timeperiods)=", length(inputs$timeperiods), 
                        "\n", "length(inputs$list_of_distances_mats)=", 
                        length(inputs$list_of_distances_mats), "\n", 
                        sep = "")
        cat(stoptxt)
        stop(stoptxt)
      }
    }
    if (is.character(inputs$envdistsfn)) {
      if (length(inputs$list_of_envdistances_mats) < length(inputs$timeperiods)) {
        stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: fewer environmental distances matrices than timeperiods.\n", 
                        "length(inputs$timeperiods)=", length(inputs$timeperiods), 
                        "\n", "length(inputs$list_of_envdistances_mats)=", 
                        length(inputs$list_of_envdistances_mats), "\n", 
                        sep = "")
        cat(stoptxt)
        stop(stoptxt)
      }
    }
    if (is.character(inputs$dispersal_multipliers_fn)) {
      if (length(inputs$list_of_dispersal_multipliers_mats) < 
          length(inputs$timeperiods)) {
        stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: fewer manual dispersal multipliers matrices than timeperiods.\n", 
                        "length(inputs$timeperiods)=", length(inputs$timeperiods), 
                        "\n", "length(inputs$list_of_dispersal_multipliers_mats)=", 
                        length(inputs$list_of_dispersal_multipliers_mats), 
                        "\n", sep = "")
        cat(stoptxt)
        stop(stoptxt)
      }
    }
    if (is.character(inputs$area_of_areas_fn)) {
      if (length(inputs$list_of_area_of_areas) < length(inputs$timeperiods)) {
        stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: fewer area-of-areas vectors than timeperiods.\n", 
                        "length(inputs$timeperiods)=", length(inputs$timeperiods), 
                        "\n", "length(inputs$list_of_area_of_areas)=", 
                        length(inputs$list_of_area_of_areas), "\n", 
                        sep = "")
        cat(stoptxt)
        stop(stoptxt)
      }
    }
    if (is.character(inputs$areas_allowed_fn)) {
      if (length(inputs$list_of_areas_allowed_mats) != 
          length(inputs$timeperiods)) {
        stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: different number of area-allowed matrices than timeperiods.\n", 
                        "length(inputs$timeperiods)=", length(inputs$timeperiods), 
                        "\n", "length(inputs$list_of_areas_allowed_mats)=", 
                        length(inputs$list_of_areas_allowed_mats), 
                        "\n", sep = "")
        cat(stoptxt)
        stop(stoptxt)
      }
      for (i in 1:length(inputs$list_of_areas_allowed_mats)) {
        tmp_areas_allow_mat = inputs$list_of_areas_allowed_mats[[i]]
        if (nrow(tmp_areas_allow_mat) != ncol(tmp_areas_allow_mat)) {
          errortxt = paste0("\n\ncheck_BioGeoBEARS_run() says: STOP ERROR:\n\nAreas-allowed matrices should be square, but your areas-allowed matrix #", 
                            i, " is not square:\n\nncol=", ncol(tmp_areas_allow_mat), 
                            "\nnrow=", nrow(tmp_areas_allow_mat), "\n\n")
          cat(errortxt)
          cat("Printing the offending areas-allowed matrix:\n\n")
          print(tmp_areas_allow_mat)
          cat("\n\n")
          stop(errortxt)
        }
      }
    }
    if (is.character(inputs$areas_adjacency_fn)) {
      if (length(inputs$list_of_areas_adjacency_mats) != 
          length(inputs$timeperiods)) {
        stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: different number of areas-adjacency matrices than timeperiods.\n", 
                        "length(inputs$timeperiods)=", length(inputs$timeperiods), 
                        "\n", "length(inputs$list_of_areas_adjacency_mats)=", 
                        length(inputs$list_of_areas_adjacency_mats), 
                        "\n", sep = "")
        cat(stoptxt)
        stop(stoptxt)
      }
      for (i in 1:length(inputs$list_of_areas_adjacency_mats)) {
        tmp_areas_adjacency_mat = inputs$list_of_areas_adjacency_mats[[i]]
        if (nrow(tmp_areas_adjacency_mat) != ncol(tmp_areas_adjacency_mat)) {
          errortxt = paste0("\n\ncheck_BioGeoBEARS_run() says: STOP ERROR:\n\nAreas-adjacency matrices should be square, but your areas-adjacency matrix #", 
                            i, " is not square:\n\nncol=", ncol(tmp_areas_adjacency_mat), 
                            "\nnrow=", nrow(tmp_areas_adjacency_mat), 
                            "\n\n")
          cat(errortxt)
          cat("Printing the offending areas-adjacency matrix:\n\n")
          print(tmp_areas_adjacency_mat)
          cat("\n\n")
          stop(errortxt)
        }
      }
    }
    if ("tree_sections_list" %in% names(inputs) == FALSE) {
      stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: You have time slices, but you do not have 'inputs$tree_sections_list'.\n", 
                      "Run 'section_the_tree()' to add tree sections to your input BioGeoBEARS_run_object.\n\n(e.g., uncomment the line(s) starting '#section_the_tree()' in the example script!!) (You will have to do this for *each* model, e.g. 6 times for the 6 models in the PhyloWiki example script.)\n\n", 
                      sep = "")
      cat(stoptxt)
      stop(stoptxt)
    }
  }
  if (inputs$use_detection_model == TRUE) {
    if (is.character(inputs$detects_fn) == FALSE) {
      stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR: You have $use_detection_model set to TRUE, but you have no \n", 
                      "detections text file given in '$detects_fn'!\n\n", 
                      sep = "")
      cat(stoptxt)
      stop(stoptxt)
    }
    else {
      if (class(inputs$detects_df) != "data.frame") {
        stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR: You have referenced a detections text file given in '$detects_fn' at:\n\n", 
                        inputs$detects_fn, "\n", "\n...but 'inputs$detects_df' is not a data.frame, perhaps empty!\n\n", 
                        "You should use 'readfiles_BioGeoBEARS_run()' to load these files.\n\n", 
                        "Printing inputs$detects_df below:\n\n", "inputs$detects_df = \n\n", 
                        sep = "")
        cat(stoptxt)
        print(inputs$detects_df)
        stop(stoptxt)
      }
      tr = check_trfn_nexus(trfn = inputs$trfn)
      tipnames = tr$tip.label
      table_rownames = row.names(inputs$detects_df)
      TF = table_rownames == tipnames
      if (sum(TF) != length(TF)) {
        stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the rownames in inputs$detects_df do not match the tip.labels in tr$tip.labels!!\n\n", 
                        "Printing both below:\n\n", sep = "")
        cat(stoptxt)
        print("tipnames:")
        print(tipnames)
        print("table_rownames:")
        print(table_rownames)
        print("match TF:")
        print(TF)
        stop(stoptxt)
      }
    }
  }
  if (inputs$use_detection_model == TRUE) {
    if (is.character(inputs$controls_fn) == FALSE) {
      stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR: You have $use_detection_model set to TRUE, but you have no \n", 
                      "taphonomic controls text file given in '$controls_fn'!\n\n", 
                      sep = "")
      cat(stoptxt)
      stop(stoptxt)
    }
    else {
      if (class(inputs$controls_df) != "data.frame") {
        stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR: You have referenced a taphonomic controls text file given in '$controls_fn' at:\n\n", 
                        inputs$controls_fn, "\n", "\n...but 'inputs$controls_df' is not a data.frame, perhaps empty!\n\n", 
                        "You should use 'readfiles_BioGeoBEARS_run()' to load these files.\n\n", 
                        "Printing inputs$controls_df below:\n\n", "inputs$controls_df = \n\n", 
                        sep = "")
        cat(stoptxt)
        print(inputs$controls_df)
        stop(stoptxt)
      }
      tr = check_trfn_nexus(trfn = inputs$trfn)
      tipnames = tr$tip.label
      table_rownames = row.names(inputs$controls_df)
      TF = table_rownames == tipnames
      if (sum(TF) != length(TF)) {
        stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the rownames in inputs$controls_df do not match the tip.labels in tr$tip.labels!!\n\n", 
                        "Printing both below:\n\n", sep = "")
        cat(stoptxt)
        print("tipnames:")
        print(tipnames)
        print("table_rownames:")
        print(table_rownames)
        print("match TF:")
        print(TF)
        stop(stoptxt)
      }
    }
  }
  return(TRUE)
}
#<bytecode: 0x10b7757e0>
 # <environment: namespace:BioGeoBEARS>