bears_optim_run_nexus <- function (BioGeoBEARS_run_object = define_BioGeoBEARS_run(), 
          skip_optim = FALSE, skip_optim_option = "return_loglike") 
{
  defaults = "\t\n\tBioGeoBEARS_run_object = define_BioGeoBEARS_run()\n\tBioGeoBEARS_run_object\n\tskip_optim=FALSE\n\t\n\tskip_optim=TRUE\n\tskip_optim_option=\"return_all\"\n\t"
  assign("last.warning", NULL, envir = baseenv())
  if (is.null(BioGeoBEARS_run_object$include_null_range)) {
    BioGeoBEARS_run_object$include_null_range = TRUE
  }
  include_null_range = BioGeoBEARS_run_object$include_null_range
  if (is.null(BioGeoBEARS_run_object$allow_null_tips)) {
    BioGeoBEARS_run_object$allow_null_tips = FALSE
  }
  if (exists("tr") == FALSE) {
    tr = read.nexus(BioGeoBEARS_run_object$trfn)
  }
  traitTF = is.null(BioGeoBEARS_run_object$trait) == FALSE
  m = NULL
  if (is.null(BioGeoBEARS_run_object$min_branchlength) == FALSE) {
    min_branchlength = BioGeoBEARS_run_object$min_branchlength
  }
  else {
    min_branchlength = 1e-06
  }
  if (traitTF) {
    if (is.na(BioGeoBEARS_run_object$timesfn) == FALSE) {
      txt = "WARNING: you have loaded a BioGeoBEARS_run_object$trait, but you have a timesfn, indicate a time-stratified analysis. Traits-based dispersal has now been implemented for time-stratified analyses, but is still experimental."
      cat("\n\n")
      cat(txt)
      cat("\n\n")
      warning(txt)
    }
    if (is.null(BioGeoBEARS_run_object$trait_Pmat_txt) == 
        TRUE) {
      txt = "STOP ERROR: you have loaded a BioGeoBEARS_run_object$trait, but you are missing a BioGeoBEARS_run_object$trait_Pmat_txt"
      cat("\n\n")
      cat(txt)
      cat("\n\n")
      stop(txt)
    }
    trait_transition_rates_TF = grepl(pattern = "trait transition rate", 
                                      x = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$desc)
    if (sum(trait_transition_rates_TF) < 1) {
      txt = "STOP ERROR: you have loaded a BioGeoBEARS_run_object$trait, but you need one or more 'trait transition rates' in  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table"
      cat("\n\n")
      cat(txt)
      cat("\n\n")
      stop(txt)
    }
    numtraitstates = ncol(BioGeoBEARS_run_object$trait@df)
    traitbased_dispersal_Ms_TF = grepl(pattern = "trait-based dispersal rate multiplier", 
                                       x = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$desc)
    if (sum(traitbased_dispersal_Ms_TF) != numtraitstates) {
      txt = paste0("STOP ERROR: you have loaded a BioGeoBEARS_run_object$trait, and it has ", 
                   numtraitstates, " states, so you need to have ", 
                   numtraitstates, " multipliers ('m1', 'm2', etc.) with 'desc' field 'trait-based dispersal rate multipliers...' in  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table. Instead, you have only this many: ", 
                   sum(traitbased_dispersal_Ms_TF))
      cat("\n\n")
      cat(txt)
      cat("\n\n")
      stop(txt)
    }
  }
  if (traitTF == TRUE) {
    trait = BioGeoBEARS_run_object$trait
    trait_as_tip_condlikes = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges = trait, 
                                                                              phy = tr, states_list = NULL, maxareas = 1, include_null_range = FALSE, 
                                                                              useAmbiguities = BioGeoBEARS_run_object$useAmbiguities, 
                                                                              trait_as_tip_condlikes = NULL)
    ntrait_states = ncol(trait_as_tip_condlikes)
    BGB_trait_model_params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[traitbased_dispersal_Ms_TF, 
                                                                                                ]
    BGB_trait_Pmat_params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[trait_transition_rates_TF, 
                                                                                               ]
    trait_Pmat_txt = BioGeoBEARS_run_object$trait_Pmat_txt
  }
  else {
    trait_as_tip_condlikes = NULL
    ntrait_states = NULL
    BGB_trait_model_params_table = NULL
    BGB_trait_Pmat_params_table = NULL
  }
  BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
  print_optim = BioGeoBEARS_run_object$print_optim
  if (BioGeoBEARS_run_object$use_detection_model == FALSE) {
    tipranges = getranges_from_LagrangePHYLIP(lgdata_fn = np(BioGeoBEARS_run_object$geogfn))
  }
  if (BioGeoBEARS_run_object$use_detection_model == TRUE) {
    if (BioGeoBEARS_run_object$use_detection_model == TRUE) {
      tipranges = tipranges_from_detects_fn(detects_fn = BioGeoBEARS_run_object$detects_fn)
    }
  }
  speedup = BioGeoBEARS_run_object$speedup
  areas = getareas_from_tipranges_object(tipranges)
  areas_list = seq(0, length(areas) - 1, 1)
  if (is.na(BioGeoBEARS_run_object$max_range_size)) {
    if (is.null(BioGeoBEARS_run_object$states_list)) {
      max_range_size = length(areas)
    }
    else {
      max_range_size = max(sapply(X = BioGeoBEARS_run_object$states_list, 
                                  FUN = length), na.rm = TRUE)
    }
  }
  else {
    max_range_size = BioGeoBEARS_run_object$max_range_size
  }
  max_numareas = max_range_size
  tipranges_df_tmp = tipranges@df
  names(tipranges_df_tmp) = paste0("col", names(tipranges_df_tmp))
  tipranges_df_tmp[tipranges_df_tmp == "?"] = 0
  TF = (rowSums(dfnums_to_numeric(tipranges_df_tmp))) > max_range_size
  if (sum(TF, na.rm = TRUE) > 0) {
    cat("\n\nERROR: Tips with ranges too big:\n", sep = "")
    print(dfnums_to_numeric(tipranges_df_tmp)[TF, ])
    cat("\n\nCheck your input geography file!\n", sep = "")
    txt = paste("ERROR: Some tips (listed above) have range sizes larger than ", 
                max_range_size, sep = "")
    stop(txt)
  }
  if (is.null(BioGeoBEARS_run_object$states_list)) {
    states_list = rcpp_areas_list_to_states_list(areas = areas, 
                                                 maxareas = max_range_size, include_null_range = BioGeoBEARS_run_object$include_null_range)
    states_list
  }
  else {
    states_list = BioGeoBEARS_run_object$states_list
  }
  all_states_list = states_list
  all_geog_states_list_usually_inferred_from_areas_maxareas = all_states_list
  BioGeoBEARS_run_object$all_geog_states_list_usually_inferred_from_areas_maxareas = all_states_list
  if (is.numeric(BioGeoBEARS_run_object$timeperiods) == FALSE) {
    if ((is.null(BioGeoBEARS_run_object$list_of_areas_allowed_mats) == 
         FALSE)) {
      areas_allowed_mat = BioGeoBEARS_run_object$list_of_areas_allowed_mats[[1]]
      original_states_list = states_list
      states_list = prune_states_list(states_list_0based_index = states_list, 
                                      areas_allowed_mat = areas_allowed_mat)
      BioGeoBEARS_run_object$states_list = states_list
      print("Limiting original_states_list using an areas_allowed matrix")
      print("original_states_list")
      print(original_states_list)
      cat("\nlength(original_states_list) = ", length(original_states_list), 
          " states/ranges.\n")
      cat("\n")
      print("states_list")
      print(states_list)
      cat("\nlength(original_states_list) = ", length(original_states_list), 
          " states/ranges.")
      cat("\nlength(states_list) = ", length(states_list), 
          " states/ranges.\n")
    }
    else {
      pass = 1
    }
    if ((is.null(BioGeoBEARS_run_object$list_of_areas_adjacency_mats) == 
         FALSE)) {
      areas_adjacency_mat = BioGeoBEARS_run_object$list_of_areas_adjacency_mats[[1]]
      original_states_list = states_list
      states_list = prune_states_list_by_adjacency(states_list_0based_index = states_list, 
                                                   areas_adjacency_mat = areas_adjacency_mat)
      BioGeoBEARS_run_object$states_list = states_list
      print("Limiting original_states_list using an area adjacency matrix")
      print("original_states_list")
      print(original_states_list)
      print(length(original_states_list))
      cat("\n")
      print("states_list")
      print(states_list)
      print("length(states_list)")
      print(length(states_list))
    }
    else {
      pass = 1
    }
  }
  if (traitTF == TRUE) {
    states_list_ORIG = states_list
  }
  TF1 = (is.null(BioGeoBEARS_run_object$list_of_areas_allowed_mats) == 
           FALSE)
  TF2 = (is.null(BioGeoBEARS_run_object$list_of_areas_adjacency_mats) == 
           FALSE)
  state_space_changing_TF = (TF1 + TF2) > 0
  need_to_print_list_of_states_list = TRUE
  master_states_list = states_list
  if ((is.numeric(BioGeoBEARS_run_object$timeperiods) == TRUE) && 
      (state_space_changing_TF == TRUE) && (is.null(BioGeoBEARS_run_object$lists_of_states_lists_0based) == 
                                            TRUE)) {
    all_geog_states_list_usually_inferred_from_areas_maxareas
    need_to_print_list_of_states_list = FALSE
    ntimes = length(BioGeoBEARS_run_object$timeperiods)
    lists_of_states_lists_0based = list()
    for (ti in 1:ntimes) {
      states_list_for_this_stratum = states_list
      if (TF1 == TRUE) {
        areas_allowed_mat = BioGeoBEARS_run_object$list_of_areas_allowed_mats[[ti]]
        states_list_for_this_stratum = prune_states_list(states_list_0based_index = states_list_for_this_stratum, 
                                                         areas_allowed_mat = areas_allowed_mat)
      }
      else {
        pass = 1
      }
      timeslice_num = ti
      if (timeslice_num == 1) {
        toptime = 0
      }
      else {
        toptime = BioGeoBEARS_run_object$timeperiods[ti - 
                                                       1]
      }
      if (timeslice_num == ntimes) {
        bottime = BioGeoBEARS_run_object$timeperiods[ti]
        catend = "\n\n"
      }
      else {
        bottime = BioGeoBEARS_run_object$timeperiods[ti]
        catend = ""
      }
      txt = paste0("bears_optim_run() note: overall states_list has ", 
                   length(master_states_list), " states/ranges. In stratum #", 
                   ti, " (", toptime, "-", bottime, " mya), states_list_for_this_stratum has ", 
                   length(states_list_for_this_stratum), " states/ranges, due to areas_allowed and/or areas_adjacency matrices. See BioGeoBEARS_run_object$lists_of_states_lists_0based.")
      cat("\n")
      cat(txt)
      cat(catend)
      if (TF2 == TRUE) {
        areas_adjacency_mat = BioGeoBEARS_run_object$list_of_areas_adjacency_mats[[ti]]
        states_list_for_this_stratum = prune_states_list_by_adjacency(states_list_0based_index = states_list_for_this_stratum, 
                                                                      areas_adjacency_mat = areas_adjacency_mat)
      }
      else {
        pass = 1
      }
      lists_of_states_lists_0based[[ti]] = states_list_for_this_stratum
    }
    BioGeoBEARS_run_object$lists_of_states_lists_0based = lists_of_states_lists_0based
  }
  if (is.null(BioGeoBEARS_run_object$lists_of_states_lists_0based) == 
      FALSE) {
    ntimes = length(BioGeoBEARS_run_object$timeperiods)
    states_allowed_TF1 = rep(TRUE, times = length(all_states_list))
    states_allowed_TF2 = rep(TRUE, times = length(all_states_list))
    states_allowed_TF3 = rep(TRUE, times = length(all_states_list))
    for (ntimes_i in 1:ntimes) {
      if ((is.null(BioGeoBEARS_run_object$list_of_areas_allowed_mats) == 
           FALSE)) {
        areas_allowed_mat = BioGeoBEARS_run_object$list_of_areas_allowed_mats[[ntimes_i]]
        states_allowed_TF1 = sapply(X = all_states_list, 
                                    FUN = check_if_state_is_allowed, areas_allowed_mat)
        if (include_null_range == TRUE) {
          states_allowed_TF1[1] = TRUE
        }
      }
      if ((is.null(BioGeoBEARS_run_object$list_of_areas_adjacency_mats) == 
           FALSE)) {
        areas_adjacency_mat = BioGeoBEARS_run_object$list_of_areas_adjacency_mats[[ntimes_i]]
        states_allowed_TF2 = sapply(X = all_states_list, 
                                    FUN = check_if_state_is_allowed_by_adjacency, 
                                    areas_adjacency_mat)
        if (include_null_range == TRUE) {
          states_allowed_TF2[1] = TRUE
        }
      }
      if ((is.null(BioGeoBEARS_run_object$lists_of_states_lists_0based) == 
           FALSE)) {
        states_allowed_TF3 = all_states_list %in% BioGeoBEARS_run_object$lists_of_states_lists_0based[[ntimes_i]]
        if (include_null_range == TRUE) {
          states_allowed_TF3[1] = TRUE
        }
      }
      states_allowed_TF = ((states_allowed_TF1 + states_allowed_TF2 + 
                              states_allowed_TF3) == 3)
      BioGeoBEARS_run_object$lists_of_states_lists_0based[[ntimes_i]] = all_states_list[states_allowed_TF]
    }
    txt = paste0("bears_optim_run() note: BioGeoBEARS_run_object$lists_of_states_lists_0based has been specified. This means there is a different state space in each timebin / stratum / epoch.")
    cat("\n")
    cat(txt)
    cat("\n")
    number_of_lists_of_states = length(BioGeoBEARS_run_object$lists_of_states_lists_0based)
    if (ntimes == number_of_lists_of_states) {
      txt = paste0("bears_optim_run() note: BioGeoBEARS_run_object has ", 
                   ntimes, " timebins and ", number_of_lists_of_states, 
                   " lists of states ranges. Check passed.")
      cat("\n")
      cat(txt)
      cat("\n")
    }
    else {
      txt = paste0("bears_optim_run() STOP ERROR: BioGeoBEARS_run_object has ", 
                   ntimes, " timebins and ", number_of_lists_of_states, 
                   " lists of states ranges. Check FAILED.")
      cat("\n")
      cat(txt)
      cat("\n")
      stop(txt)
    }
    if (need_to_print_list_of_states_list == TRUE) {
      for (ti in 1:ntimes) {
        states_list_for_this_stratum = BioGeoBEARS_run_object$lists_of_states_lists_0based[[ti]]
        timeslice_num = ti
        if (timeslice_num == 1) {
          toptime = 0
        }
        else {
          toptime = BioGeoBEARS_run_object$timeperiods[ti - 
                                                         1]
        }
        if (timeslice_num == ntimes) {
          bottime = BioGeoBEARS_run_object$timeperiods[ti]
          catend = "\n\n"
        }
        else {
          bottime = BioGeoBEARS_run_object$timeperiods[ti]
          catend = ""
        }
        txt = paste0("bears_optim_run() note: overall states_list has ", 
                     length(master_states_list), " states/ranges. In stratum #", 
                     ti, " (", toptime, "-", bottime, " mya), states_list_for_this_stratum has ", 
                     length(states_list_for_this_stratum), " states/ranges, due to user-specified states_lists. See BioGeoBEARS_run_object$lists_of_states_lists_0based.")
        cat("\n")
        cat(txt)
        cat(catend)
      }
    }
  }
  if (is.na(BioGeoBEARS_run_object$force_sparse)) {
    if (length(states_list) > 128) {
      force_sparse = TRUE
      cat("\nNote: force_sparse being set to TRUE, as length(states_list) > 128\n", 
          sep = "")
    }
    else {
      force_sparse = FALSE
    }
  }
  else {
    force_sparse = BioGeoBEARS_run_object$force_sparse
  }
  trfn = np(BioGeoBEARS_run_object$trfn)
  phy = check_trfn_nexus(trfn = trfn)
  if (is.null(BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state) == 
      TRUE) {
    if (BioGeoBEARS_run_object$use_detection_model == FALSE) {
      tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, 
                                                                                             phy, states_list = states_list, maxareas = max_numareas, 
                                                                                             include_null_range = BioGeoBEARS_run_object$include_null_range, 
                                                                                             useAmbiguities = BioGeoBEARS_run_object$useAmbiguities, 
                                                                                             trait_as_tip_condlikes = trait_as_tip_condlikes, 
                                                                                             allow_null_tips = BioGeoBEARS_run_object$allow_null_tips)
    }
    else {
      numareas = length(areas)
      detects_df = BioGeoBEARS_run_object$detects_df
      controls_df = BioGeoBEARS_run_object$controls_df
      mean_frequency = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mf", 
                                                                                    "init"]
      dp = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["dp", 
                                                                        "init"]
      fdp = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["fdp", 
                                                                         "init"]
      tip_condlikes_of_data_on_each_state = tiplikes_wDetectionModel(states_list_0based_index = states_list, 
                                                                     phy = phy, numareas = numareas, detects_df = detects_df, 
                                                                     controls_df = controls_df, mean_frequency = mean_frequency, 
                                                                     dp = dp, fdp = fdp, null_range_gets_0_like = TRUE, 
                                                                     return_LnLs = TRUE, relative_LnLs = TRUE, exp_LnLs = TRUE, 
                                                                     error_check = TRUE)
      if (is.null(BioGeoBEARS_run_object$prior_by_range_size) == 
          FALSE) {
        cat("\n\nNOTE: BioGeoBEARS will multiply the initial tip conditional likelihoods ('tip_condlikes_of_data_on_each_state') by the user-specified 'BioGeoBEARS_run_object$prior_by_range_size'.\n")
        for (iii in 1:nrow(tip_condlikes_of_data_on_each_state)) {
          tip_condlikes_of_data_on_each_state[iii, ] = tip_condlikes_of_data_on_each_state[iii, 
                                                                                           ] * BioGeoBEARS_run_object$prior_by_range_size
        }
        cat("...done.\n\n")
      }
    }
  }
  else {
    tip_condlikes_of_data_on_each_state = BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state
  }
  numstates = ncol(tip_condlikes_of_data_on_each_state)
  if (is.null(BioGeoBEARS_run_object$printlevel)) {
    BioGeoBEARS_run_object$printlevel = 0
  }
  printlevel = BioGeoBEARS_run_object$printlevel
  check_BioGeoBEARS_run_nexus(BioGeoBEARS_run_object)
  if (BioGeoBEARS_run_object$rescale_params == TRUE) {
    BioGeoBEARS_model_object@params_table = scale_BGB_params(orig_params_table = BioGeoBEARS_model_object@params_table, 
                                                             add_smin = 0, add_smax = 1)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
  }
  params = BioGeoBEARS_model_object_to_init_params(BioGeoBEARS_model_object)
  minj = 1e-05
  lower = BioGeoBEARS_model_object_to_params_lower(BioGeoBEARS_model_object)
  upper = BioGeoBEARS_model_object_to_params_upper(BioGeoBEARS_model_object)
  num_cores_to_use = BioGeoBEARS_run_object$num_cores_to_use
  cluster_already_open = BioGeoBEARS_run_object$cluster_already_open
  cluster_was_open = FALSE
  if (.Platform$GUI != "AQUA" && ((is.na(num_cores_to_use) == 
                                   TRUE) || ((is.na(num_cores_to_use) == FALSE) && (num_cores_to_use > 
                                                                                    1)))) {
    num_cores_computer_has = detectCores()
    txt = paste0("Your computer has ", num_cores_computer_has, 
                 " cores.")
    cat("\n")
    cat(txt)
    cat("\n")
    if (is.null(num_cores_to_use)) {
      num_cores_to_use = num_cores_computer_has
    }
    if (num_cores_to_use > num_cores_computer_has) {
      txt = paste0("WARNING from bears_optim_run(): You specified num_cores_to_use=", 
                   num_cores_to_use, " cores, but your computer only has ", 
                   num_cores_computer_has, ". Resetting to ", num_cores_computer_has, 
                   ".")
      cat("\n")
      cat(txt)
      cat("\n")
      warning(txt)
      num_cores_to_use = num_cores_computer_has
    }
    cat("\nYour computer has ", num_cores_computer_has, " cores. You have chosen to use:\nnum_cores_to_use = ", 
        num_cores_to_use, " cores for the matrix exponentiations in the likelihood calculations.\n", 
        sep = "")
    if (is.logical(cluster_already_open) == TRUE) {
      if (cluster_already_open == FALSE) {
        cluster_already_open = makeCluster(rep("localhost", 
                                               num_cores_to_use), type = "SOCK")
        cat("Started cluster with ", num_cores_to_use, 
            " cores.\n\n", sep = "")
        cluster_open = TRUE
        cluster_was_open = FALSE
      }
    }
    else {
      cluster_was_open = TRUE
      cat("Cluster with ", num_cores_to_use, " cores already open.\n\n", 
          sep = "")
    }
  }
  else {
    num_cores_computer_has = detectCores()
    txt = paste0("Your computer has ", num_cores_computer_has, 
                 " cores.")
    cat("\n")
    cat(txt)
    cat("\n")
    if (num_cores_to_use > 1) {
      txt = paste0("WARNING from bears_optim_run(): You specified num_cores_to_use=", 
                   num_cores_to_use, " cores, but in R.app, multicore functionality doesn't work. Resetting num_cores_to_use=1.")
      cat("\n")
      cat(txt)
      cat("\n")
      warning(txt)
      num_cores_to_use = num_cores_computer_has
    }
    cluster_already_open = NULL
    cluster_was_open = FALSE
  }
  if (force_sparse == TRUE) {
    cat("\nNote: force_sparse is set to TRUE; length(states_list)=", 
        length(states_list), "\n", sep = "")
    txt = paste0("\n\nNote: sparse matrix exponentiation is being used. When on exponentiation on a branch is completed, 'L',  'R', 'S', or 'U' will print to screen  (for left and right branches, S segments in time-stratified analyses, U for uppass on a segment/branch). This will help you judge the time this analysis will take.  An ML search takes (at least) 100+ downpass calculations of the log-likelihood (lnL) of the tip data, given on left branch, given the tree, model, and parameters. Each downpass requires a matrix exponentiation on each branch of the tree. Your tree has ", 
                 length(tr$tip.label), " tips, thus ", length(tr$tip.label) + 
                   length(tr$tip.label) - 1, " branches. The transition matrix has ", 
                 numstates, " states (states=possible geographic ranges), so it would be of size ", 
                 numstates, "x", numstates, " if it were dense, but you are using the sparse matrix routine to speed up calculations. Starting now...\n")
    cat(txt)
  }
  if (is.numeric(BioGeoBEARS_run_object$timeperiods)) {
    allareas = areas_list
    all_states_list = states_list
    all_geog_states_list_usually_inferred_from_areas_maxareas
    use_optimx = BioGeoBEARS_run_object$use_optimx
    if ((use_optimx == FALSE) || (use_optimx == "optim")) {
      cat("\n\nNOTE: Before running optim(), here is a test calculation of the data likelihood\nusing calc_loglike_for_optim_stratified() on initial parameter values, with printlevel=2...\nif this crashes, the error messages are more helpful\nthan those from inside optim().\n", 
          sep = "")
      inputs = BioGeoBEARS_run_object
      inputs$printlevel = 2
      loglike = calc_loglike_for_optim_stratified(params = params, 
                                                  BioGeoBEARS_run_object = inputs, phy = phy, tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                                                  print_optim = TRUE, areas_list = areas_list, 
                                                  states_list = states_list, force_sparse = force_sparse, 
                                                  cluster_already_open = cluster_already_open)
      if ((loglike < -1e+10) || (is.finite(loglike) == 
                                 FALSE)) {
        txt = paste0("STOP ERROR #1 in bears_optim_run(). Test calculation of the likelihood with calc_loglike_for_optim_stratified() returned LnL=", 
                     loglike, ", which is not a valid starting likelihood. Probably, you have an overly constrained analysis and have thus made your data impossible. For example, if your tips had ranges A and B, but you disallowed the state AB, then your data would be impossible under DEC, because AB is a required intermediate state. Try removing some of the areas allowed/area adjacency/manual states list constraints.  You can also try changing manual dispersal multipliers from 0 to some small value (e.g. 0.00001) -- note that this can works on dispersal multiplers, but NOT on area constraints, which have to be 0 or 1.  Have a CAREFUL THINK about what you are doing and why you think the list of ranges should be so constrained - do you actually have a good argument for your constraints?")
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        stop(txt)
      }
      cat("\ncalc_loglike_for_optim_stratified() on initial parameters loglike=", 
          loglike, "\n\n\n\nCalculation of likelihood on initial parameters: successful.\n\nNow starting Maximum Likelihood (ML) parameter optimization with optim()...\n\n", 
          sep = "")
      if (skip_optim == TRUE) {
        if (skip_optim_option == "return_loglike") {
          cat("Skipping ML search as skip_optim==TRUE. Returning only the log-likelihood.\n\n", 
              sep = "")
          return(loglike)
        }
        if (skip_optim_option == "return_all") {
          inputs = BioGeoBEARS_run_object
          cat("Skipping ML search as skip_optim==TRUE. Returning everything (ancestral probabilities etc.) conditional on 'est' parameters given in 'BioGeoBEARS_model_object' .\n\n", 
              sep = "")
          optim_result2 = put_params_into_optim_or_optimx_result(BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object, 
                                                                 total_loglikelihood = loglike, use_optimx = BioGeoBEARS_run_object$use_optimx)
        }
      }
      else {
        inputs = BioGeoBEARS_run_object
        optim_result2 = optim(par = params, fn = calc_loglike_for_optim_stratified, 
                              BioGeoBEARS_run_object = inputs, phy = phy, 
                              tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                              print_optim = print_optim, areas_list = areas_list, 
                              states_list = states_list, force_sparse = force_sparse, 
                              cluster_already_open = cluster_already_open, 
                              method = "L-BFGS-B", lower = lower, upper = upper, 
                              control = list(fnscale = -1, trace = 2, maxit = 500))
      }
    }
    if ((use_optimx == TRUE) || (use_optimx == "optimx")) {
      print_optim = BioGeoBEARS_run_object$print_optim
      if (speedup) {
        control_list = list(all.methods = FALSE, maximize = TRUE, 
                            save.failures = TRUE)
        num_free_params = sum(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$type == 
                                "free")
        num_free_params
        itnmax = 50 * num_free_params
      }
      else {
        control_list = list(all.methods = FALSE, maximize = TRUE, 
                            save.failures = TRUE)
        itnmax = NULL
      }
      cat("\n\nNOTE: Before running optimx(), here is a test calculation of the data likelihood\nusing calc_loglike_for_optim_stratified() on initial parameter values, with printlevel=2...\nif this crashes, the error messages are more helpful\nthan those from inside optimx().\n\n", 
          sep = "")
      inputs = BioGeoBEARS_run_object
      inputs$printlevel = 2
      loglike = calc_loglike_for_optim_stratified(params = params, 
                                                  BioGeoBEARS_run_object = inputs, phy = phy, tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                                                  print_optim = TRUE, areas_list = areas_list, 
                                                  states_list = states_list, force_sparse = force_sparse, 
                                                  cluster_already_open = cluster_already_open)
      if ((loglike < -1e+10) || (is.finite(loglike) == 
                                 FALSE)) {
        txt = paste0("STOP ERROR #2 in bears_optim_run(). Test calculation of the likelihood with calc_loglike_for_optim_stratified() returned LnL=", 
                     loglike, ", which is not a valid starting likelihood. Probably, you have an overly constrained analysis and have thus made your data impossible. For example, if your tips had ranges A and B, but you disallowed the state AB, then your data would be impossible under DEC, because AB is a required intermediate state. Try removing some of the areas allowed/area adjacency/manual states list constraints.  You can also try changing manual dispersal multipliers from 0 to some small value (e.g. 0.00001) -- note that this can works on dispersal multiplers, but NOT on area constraints, which have to be 0 or 1.  Have a CAREFUL THINK about what you are doing and why you think the list of ranges should be so constrained - do you actually have a good argument for your constraints?")
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        stop(txt)
      }
      cat("\ncalc_loglike_for_optim_stratified() on initial parameters loglike=", 
          loglike, "\n\n\n\nCalculation of likelihood on initial parameters: successful.\n\nNow starting Maximum Likelihood (ML) parameter optimization with optimx()...\n\n", 
          sep = "")
      cat("\n\nPrinting any warnings() that occurred during calc_loglike_for_optim_stratified():\n\n")
      print(warnings())
      if (skip_optim == TRUE) {
        if (skip_optim_option == "return_loglike") {
          cat("Skipping ML search as skip_optim==TRUE. Returning only the log-likelihood.\n\n", 
              sep = "")
          return(loglike)
        }
        if (skip_optim_option == "return_all") {
          inputs = BioGeoBEARS_run_object
          cat("Skipping ML search as skip_optim==TRUE. Returning everything (ancestral probabilities etc.) conditional on 'est' parameters given in 'BioGeoBEARS_model_object' .\n\n", 
              sep = "")
          optim_result2 = put_params_into_optim_or_optimx_result(BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object, 
                                                                 total_loglikelihood = loglike, use_optimx = BioGeoBEARS_run_object$use_optimx)
        }
      }
      else {
        inputs = BioGeoBEARS_run_object
        scalecheck_results = optimx_scalecheck(par = params, 
                                               lower = lower, upper = upper)
        cat("\n\nResults of optimx_scalecheck() below. Note: sometimes rescaling parameters may be helpful for ML searches, when the parameters have much different absolute sizes. This can be attempted by setting BioGeoBEARS_run_object$rescale_params = TRUE.\n\n")
        print(scalecheck_results)
        minqa_TF = is.element("minqa", installed.packages()[, 
                                                            1])
        if (minqa_TF == FALSE) {
          if (packageVersion("optimx") > 2017) {
            txt = "Warning in bears_optim_run(): optimx version 2018.7.10 requires package 'minqa' to do optimx ML optimization with the 'bobyqa' method (optimization with mix/max limits on parameters). However, optimx 2018.7.10 doesn't load 'minqa' automatically, so you may have to do:\n\ninstall.packages('minqa')\n\n...and re-run, to get rid of this warning, and/or the error where optimx returns NA for the parameter inferences after one step, and crashes the resulting uppass calculations."
            cat("\n\n")
            cat(txt)
            cat("\n\n")
            warning(txt)
            requireNamespace("minqa")
          }
        }
        optim_result2 = optimx(par = params, fn = calc_loglike_for_optim_stratified, 
                               lower = lower, upper = upper, itnmax = itnmax, 
                               method = c("bobyqa"), control = control_list, 
                               BioGeoBEARS_run_object = inputs, phy = phy, 
                               tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                               print_optim = print_optim, areas_list = areas_list, 
                               states_list = states_list, force_sparse = force_sparse, 
                               cluster_already_open = cluster_already_open)
      }
    }
    if (use_optimx == "GenSA") {
      cat("\n\nNOTE: You are optimizing with GenSA::GenSA() ('Generalized Simulated Annealing') instead of optimx() or optim(). GenSA seems to be better for more complex problems (4+ parameters, wildly different scalings). However, it will likely be slower, as it does more calculations of the likelihood to search the parameter space.")
      cat("\n\nNOTE: Before running GenSA::GenSA(), here is a test calculation of the data likelihood\nusing calc_loglike_for_optim_stratified() on initial parameter values, with printlevel=2...\nif this crashes, the error messages are more helpful\nthan those from inside GenSA::GenSA().\n", 
          sep = "")
      inputs = BioGeoBEARS_run_object
      inputs$printlevel = 2
      loglike = calc_loglike_for_optim_stratified(params = params, 
                                                  BioGeoBEARS_run_object = inputs, phy = phy, tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                                                  print_optim = print_optim, areas_list = areas_list, 
                                                  states_list = states_list, force_sparse = force_sparse, 
                                                  cluster_already_open = cluster_already_open)
      if ((loglike < -1e+10) || (is.finite(loglike) == 
                                 FALSE)) {
        txt = paste0("STOP ERROR #3 in bears_optim_run(). Test calculation of the likelihood with calc_loglike_for_optim_stratified() returned LnL=", 
                     loglike, ", which is not a valid starting likelihood. Probably, you have an overly constrained analysis and have thus made your data impossible. For example, if your tips had ranges A and B, but you disallowed the state AB, then your data would be impossible under DEC, because AB is a required intermediate state. Try removing some of the areas allowed/area adjacency/manual states list constraints.  You can also try changing manual dispersal multipliers from 0 to some small value (e.g. 0.00001) -- note that this can works on dispersal multiplers, but NOT on area constraints, which have to be 0 or 1.  Have a CAREFUL THINK about what you are doing and why you think the list of ranges should be so constrained - do you actually have a good argument for your constraints?")
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        stop(txt)
      }
      cat("\ncalc_loglike_for_optim_stratified() on initial parameters loglike=", 
          loglike, "\n\n\n\nCalculation of likelihood on initial parameters: successful.\n\nNow starting Maximum Likelihood (ML) parameter optimization with GenSA::GenSA()...\n\n", 
          sep = "")
      if (skip_optim == TRUE) {
        if (skip_optim_option == "return_loglike") {
          cat("Skipping ML search as skip_optim==TRUE. Returning only the log-likelihood.\n\n", 
              sep = "")
          return(loglike)
        }
        if (skip_optim_option == "return_all") {
          inputs = BioGeoBEARS_run_object
          cat("Skipping ML search as skip_optim==TRUE. Returning everything (ancestral probabilities etc.) conditional on 'est' parameters given in 'BioGeoBEARS_model_object' .\n\n", 
              sep = "")
          optim_result2 = put_params_into_optim_or_optimx_result(BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object, 
                                                                 total_loglikelihood = loglike, use_optimx = BioGeoBEARS_run_object$use_optimx)
        }
      }
      else {
        control_list = list(nb.stop.improvement = 50, 
                            simple.function = TRUE, trace.mat = TRUE)
        if (is.null(BioGeoBEARS_run_object$temperature) == 
            FALSE) {
          temperature = BioGeoBEARS_run_object$temperature
          control_list = c(control_list, list(temperature = temperature))
        }
        if (is.null(BioGeoBEARS_run_object$max.call) == 
            FALSE) {
          max.call = BioGeoBEARS_run_object$max.call
          control_list = c(control_list, list(max.call = max.call))
        }
        else {
          max.call = length(params) * 250
          control_list = c(control_list, list(max.call = max.call))
        }
        inputs = BioGeoBEARS_run_object
        optim_result2 = GenSA::GenSA(par = params, fn = calc_loglike_for_optim_stratified_neg, 
                                     BioGeoBEARS_run_object = inputs, phy = phy, 
                                     tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                                     print_optim = print_optim, areas_list = areas_list, 
                                     states_list = states_list, force_sparse = force_sparse, 
                                     cluster_already_open = cluster_already_open, 
                                     lower = lower, upper = upper, control = control_list)
      }
    }
  }
  else {
    use_optimx = BioGeoBEARS_run_object$use_optimx
    if ((use_optimx == FALSE) || (use_optimx == "optim")) {
      cat("\n\nNOTE: Before running optim(), here is a test calculation of the data likelihood\nusing calc_loglike_for_optim() on initial parameter values, with printlevel=2...\nif this crashes, the error messages are more helpful\nthan those from inside optim().\n\n", 
          sep = "")
      inputs = BioGeoBEARS_run_object
      inputs$printlevel = 2
      loglike = calc_loglike_for_optim(params, BioGeoBEARS_run_object = inputs, 
                                       phy = phy, tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                                       print_optim = BioGeoBEARS_run_object$print_optim, 
                                       areas_list = areas_list, states_list = states_list, 
                                       force_sparse = force_sparse, cluster_already_open = cluster_already_open, 
                                       return_what = "loglike", calc_ancprobs = FALSE)
      if ((loglike < -1e+10) || (is.finite(loglike) == 
                                 FALSE)) {
        txt = paste0("STOP ERROR #4 in bears_optim_run(). Test calculation of the likelihood with calc_loglike_for_optim() returned LnL=", 
                     loglike, ", which is not a valid starting likelihood. Probably, you have an overly constrained analysis and have thus made your data impossible. For example, if your tips had ranges A and B, but you disallowed the state AB, then your data would be impossible under DEC, because AB is a required intermediate state. Try removing some of the areas allowed/area adjacency/manual states list constraints.  You can also try changing manual dispersal multipliers from 0 to some small value (e.g. 0.00001) -- note that this can works on dispersal multiplers, but NOT on area constraints, which have to be 0 or 1.  Have a CAREFUL THINK about what you are doing and why you think the list of ranges should be so constrained - do you actually have a good argument for your constraints?")
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        stop(txt)
      }
      cat("\ncalc_loglike_for_optim() on initial parameters loglike=", 
          loglike, "\n\n\n\nCalculation of likelihood on initial parameters: successful.\n\nNow starting Maximum Likelihood (ML) parameter optimization with optim()...\n\n", 
          sep = "")
      if (skip_optim == TRUE) {
        if (skip_optim_option == "return_loglike") {
          cat("Skipping ML search as skip_optim==TRUE. Returning only the log-likelihood.\n\n", 
              sep = "")
          return(loglike)
        }
        if (skip_optim_option == "return_all") {
          inputs = BioGeoBEARS_run_object
          cat("Skipping ML search as skip_optim==TRUE. Returning everything (ancestral probabilities etc.) conditional on 'est' parameters given in 'BioGeoBEARS_model_object' .\n\n", 
              sep = "")
          optim_result2 = put_params_into_optim_or_optimx_result(BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object, 
                                                                 total_loglikelihood = loglike, use_optimx = BioGeoBEARS_run_object$use_optimx)
        }
      }
      else {
        parscale = (upper - lower)/min(upper - lower)
        print("parscale:")
        print(parscale)
        optim_result2 = optim(par = params, fn = calc_loglike_for_optim, 
                              BioGeoBEARS_run_object = BioGeoBEARS_run_object, 
                              phy = phy, tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                              print_optim = print_optim, areas_list = areas_list, 
                              states_list = states_list, force_sparse = force_sparse, 
                              cluster_already_open = cluster_already_open, 
                              return_what = "loglike", calc_ancprobs = FALSE, 
                              method = "L-BFGS-B", lower = lower, upper = upper, 
                              control = list(fnscale = -1, trace = 2, parscale = parscale))
      }
    }
    if ((use_optimx == TRUE) || (use_optimx == "optimx")) {
      if (speedup) {
        parscale = (upper - lower)/min(upper - lower)
        print("parscale:")
        print(parscale)
        control_list = list(all.methods = FALSE, maximize = TRUE, 
                            save.failures = TRUE)
        num_free_params = sum(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$type == 
                                "free")
        num_free_params
        itnmax = 50 * num_free_params
      }
      else {
        parscale = (upper - lower)/min(upper - lower)
        print("parscale:")
        print(parscale)
        control_list = list(all.methods = FALSE, maximize = TRUE, 
                            save.failures = TRUE)
        itnmax = NULL
      }
      cat("\n\nNOTE: Before running optimx(), here is a test calculation of the data likelihood\nusing calc_loglike_for_optim() on initial parameter values, with printlevel=2...\nif this crashes, the error messages are more helpful\nthan those from inside optimx().\n\n", 
          sep = "")
      inputs = BioGeoBEARS_run_object
      inputs$printlevel = 2
      loglike = calc_loglike_for_optim(params, BioGeoBEARS_run_object = inputs, 
                                       phy = phy, tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                                       print_optim = BioGeoBEARS_run_object$print_optim, 
                                       areas_list = areas_list, states_list = states_list, 
                                       force_sparse = force_sparse, cluster_already_open = cluster_already_open, 
                                       return_what = "loglike", calc_ancprobs = FALSE)
      if ((loglike < -1e+10) || (is.finite(loglike) == 
                                 FALSE)) {
        txt = paste0("STOP ERROR #5 in bears_optim_run(). Test calculation of the likelihood with calc_loglike_for_optim() returned LnL=", 
                     loglike, ", which is not a valid starting likelihood. Probably, you have an overly constrained analysis and have thus made your data impossible. For example, if your tips had ranges A and B, but you disallowed the state AB, then your data would be impossible under DEC, because AB is a required intermediate state. Try removing some of the areas allowed/area adjacency/manual states list constraints.  You can also try changing manual dispersal multipliers from 0 to some small value (e.g. 0.00001) -- note that this can works on dispersal multiplers, but NOT on area constraints, which have to be 0 or 1.  Have a CAREFUL THINK about what you are doing and why you think the list of ranges should be so constrained - do you actually have a good argument for your constraints?")
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        stop(txt)
      }
      cat("\ncalc_loglike_for_optim() on initial parameters loglike=", 
          loglike, "\n\n\n\nCalculation of likelihood on initial parameters: successful.\n\nNow starting Maximum Likelihood (ML) parameter optimization with optimx()...\n\n", 
          sep = "")
      cat("\n\nPrinting any warnings() that occurred during calc_loglike_for_optim():\n\n")
      print(warnings())
      if (skip_optim == TRUE) {
        if (skip_optim_option == "return_loglike") {
          cat("Skipping ML search as skip_optim==TRUE. Returning only the log-likelihood.\n\n", 
              sep = "")
          return(loglike)
        }
        if (skip_optim_option == "return_all") {
          inputs = BioGeoBEARS_run_object
          cat("Skipping ML search as skip_optim==TRUE. Returning everything (ancestral probabilities etc.) conditional on 'est' parameters given in 'BioGeoBEARS_model_object' .\n\n", 
              sep = "")
          optim_result2 = put_params_into_optim_or_optimx_result(BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object, 
                                                                 total_loglikelihood = loglike, use_optimx = BioGeoBEARS_run_object$use_optimx)
        }
      }
      else {
        if (packageVersion("optimx") < 2013) {
          optim_result2 = optimx(par = params, fn = calc_loglike_for_optim, 
                                 gr = NULL, hess = NULL, lower = lower, upper = upper, 
                                 method = c("bobyqa"), itnmax = itnmax, hessian = NULL, 
                                 control = control_list, BioGeoBEARS_run_object = BioGeoBEARS_run_object, 
                                 phy = phy, tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                                 print_optim = print_optim, areas_list = areas_list, 
                                 states_list = states_list, force_sparse = force_sparse, 
                                 cluster_already_open = cluster_already_open, 
                                 return_what = "loglike", calc_ancprobs = FALSE)
        }
        else {
          scalecheck_results = optimx_scalecheck(par = params, 
                                                 lower = lower, upper = upper)
          cat("\n\nResults of optimx_scalecheck() below. Note: sometimes rescaling parameters may be helpful for ML searches, when the parameters have much different absolute sizes. This can be attempted by setting BioGeoBEARS_run_object$rescale_params = TRUE.\n\n")
          print(scalecheck_results)
          minqa_TF = is.element("minqa", installed.packages()[, 
                                                              1])
          if (minqa_TF == FALSE) {
            if (packageVersion("optimx") > 2017) {
              txt = "Warning in bears_optim_run(): optimx version 2018.7.10 requires package 'minqa' to do optimx ML optimization with the 'bobyqa' method (optimization with mix/max limits on parameters). However, optimx 2018.7.10 doesn't load 'minqa' automatically, so you may have to do:\n\ninstall.packages('minqa')\n\n...and re-run, to get rid of this warning, and/or the error where optimx returns NA for the parameter inferences after one step, and crashes the resulting uppass calculations."
              cat("\n\n")
              cat(txt)
              cat("\n\n")
              warning(txt)
              requireNamespace("minqa")
            }
          }
          optim_result2 = optimx(par = params, fn = calc_loglike_for_optim, 
                                 gr = NULL, hess = NULL, lower = lower, upper = upper, 
                                 method = c("bobyqa"), itnmax = itnmax, hessian = FALSE, 
                                 control = control_list, BioGeoBEARS_run_object = BioGeoBEARS_run_object, 
                                 phy = phy, tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                                 print_optim = print_optim, areas_list = areas_list, 
                                 states_list = states_list, force_sparse = force_sparse, 
                                 cluster_already_open = cluster_already_open, 
                                 return_what = "loglike", calc_ancprobs = FALSE)
        }
      }
    }
    if (use_optimx == "GenSA") {
      cat("\n\nNOTE: You are optimizing with GenSA::GenSA() ('Generalized Simulated Annealing') instead of optimx() or optim(). GenSA seems to be better for more complex problems (4+ parameters, wildly different scalings). However, it will likely be slower, as it does more calculations of the likelihood to search the parameter space.")
      cat("\n\nNOTE: Before running GenSA::GenSA(), here is a test calculation of the data likelihood\nusing calc_loglike_for_optim() on initial parameter values, with printlevel=2...\nif this crashes, the error messages are more helpful\nthan those from inside GenSA::GenSA().\n\n", 
          sep = "")
      inputs = BioGeoBEARS_run_object
      inputs$printlevel = 2
      loglike = calc_loglike_for_optim(params, BioGeoBEARS_run_object = inputs, 
                                       phy = phy, tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                                       print_optim = BioGeoBEARS_run_object$print_optim, 
                                       areas_list = areas_list, states_list = states_list, 
                                       force_sparse = force_sparse, cluster_already_open = cluster_already_open, 
                                       return_what = "loglike", calc_ancprobs = FALSE)
      if ((loglike < -1e+10) || (is.finite(loglike) == 
                                 FALSE)) {
        txt = paste0("STOP ERROR #6 in bears_optim_run(). Test calculation of the likelihood with calc_loglike_for_optim() returned LnL=", 
                     loglike, ", which is not a valid starting likelihood. Probably, you have an overly constrained analysis and have thus made your data impossible. For example, if your tips had ranges A and B, but you disallowed the state AB, then your data would be impossible under DEC, because AB is a required intermediate state. Try removing some of the areas allowed/area adjacency/manual states list constraints.  You can also try changing manual dispersal multipliers from 0 to some small value (e.g. 0.00001) -- note that this can works on dispersal multiplers, but NOT on area constraints, which have to be 0 or 1.  Have a CAREFUL THINK about what you are doing and why you think the list of ranges should be so constrained - do you actually have a good argument for your constraints?")
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        stop(txt)
      }
      cat("\ncalc_loglike_for_optim() on initial parameters loglike=", 
          loglike, "\n\n\n\nCalculation of likelihood on initial parameters: successful.\n\nNow starting Maximum Likelihood (ML) parameter optimization with GenSA::GenSA()...\n\n", 
          sep = "")
      cat("\n\nPrinting any warnings() that occurred during calc_loglike_for_optim():\n\n")
      print(warnings())
      if (skip_optim == TRUE) {
        if (skip_optim_option == "return_loglike") {
          cat("Skipping ML search as skip_optim==TRUE. Returning only the log-likelihood.\n\n", 
              sep = "")
          return(loglike)
        }
        if (skip_optim_option == "return_all") {
          inputs = BioGeoBEARS_run_object
          cat("Skipping ML search as skip_optim==TRUE. Returning everything (ancestral probabilities etc.) conditional on 'est' parameters given in 'BioGeoBEARS_model_object' .\n\n", 
              sep = "")
          optim_result2 = put_params_into_optim_or_optimx_result(BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object, 
                                                                 total_loglikelihood = loglike, use_optimx = BioGeoBEARS_run_object$use_optimx)
        }
      }
      else {
        control_list = list(nb.stop.improvement = 50, 
                            simple.function = TRUE, trace.mat = TRUE)
        if (is.null(BioGeoBEARS_run_object$temperature) == 
            FALSE) {
          temperature = BioGeoBEARS_run_object$temperature
          control_list = c(control_list, list(temperature = temperature))
        }
        if (is.null(BioGeoBEARS_run_object$max.call) == 
            FALSE) {
          max.call = BioGeoBEARS_run_object$max.call
          control_list = c(control_list, list(max.call = max.call))
        }
        else {
          max.call = length(params) * 250
          control_list = c(control_list, list(max.call = max.call))
        }
        optim_result2 = GenSA::GenSA(par = params, fn = calc_loglike_for_optim_neg, 
                                     lower = lower, upper = upper, control = control_list, 
                                     BioGeoBEARS_run_object = BioGeoBEARS_run_object, 
                                     phy = phy, tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                                     print_optim = print_optim, areas_list = areas_list, 
                                     states_list = states_list, force_sparse = force_sparse, 
                                     cluster_already_open = cluster_already_open, 
                                     return_what = "loglike", calc_ancprobs = FALSE)
      }
    }
  }
  if ((skip_optim == TRUE) && (skip_optim_option == "return_loglike")) {
    return(loglike)
  }
  optimx_result = optim_result2
  use_optimx = BioGeoBEARS_run_object$use_optimx
  if (printlevel >= 0) {
    cat("\n\nThis is the output from optim, optimx, or GenSA. Check the help on those functions to\ninterpret this output and check for convergence issues:\n\n")
    print(optimx_result)
  }
  if (printlevel >= 1) {
    cat("\n\nReading the optim/optimx/GenSA output into the BioGeoBEARS_model object:\n\nBioGeoBEARS_model_object =\n\n")
  }
  BioGeoBEARS_model_object = update_BioGeoBEARS_model_object_w_optimx_result(BioGeoBEARS_model_object, 
                                                                             optimx_result, use_optimx)
  if (any(is.na(BioGeoBEARS_model_object@params_table$est)) == 
      TRUE) {
    txt = "STOP ERROR in bears_optim_run(). For some reason, your ML optimizer returned one or more NA / NaN values for the estimated parameters. Probably this is a version conflict with an update to one of the optimizer functions/packages (e.g., optim, optimx, minqa, GenSA. Printing BioGeoBEARS_model_object@params_table to screen, below.  Email the BioGeoBEARS Google Group if you cannot figure out the problem."
    cat("\n\n")
    cat(txt)
    cat("\n\n")
    cat("BioGeoBEARS_model_object@params_table:\n\n")
    print(BioGeoBEARS_model_object@params_table)
    cat("\n\n")
    stop(txt)
  }
  if (BioGeoBEARS_run_object$rescale_params == TRUE) {
    cat("\n(Because BioGeoBEARS_run_object$rescale_params == TRUE, using unscale_BGB_params() to return parameter estimates to the original scaling...\n")
    BioGeoBEARS_model_object@params_table = unscale_BGB_params(scaled_params_table = BioGeoBEARS_model_object@params_table)
    if (BioGeoBEARS_run_object$use_optimx == FALSE) {
      optim_result2$par = BioGeoBEARS_model_object@params_table$est[BioGeoBEARS_model_object@params_table$type == 
                                                                      "free"]
    }
    if ((BioGeoBEARS_run_object$use_optimx == TRUE) || (BioGeoBEARS_run_object$use_optimx == 
                                                        "optimx")) {
      if (packageVersion("optimx") >= 2013) {
        param_names = names(optim_result2)
        param_1st_letter = substr(x = param_names, start = 1, 
                                  stop = 1)
        param_TF = param_1st_letter == "p"
        param_names = param_names[param_TF]
        optim_result2[param_names] = BioGeoBEARS_model_object@params_table$est[BioGeoBEARS_model_object@params_table$type == 
                                                                                 "free"]
      }
      if (packageVersion("optimx") < 2013) {
        optim_result2$par[[1]] = BioGeoBEARS_model_object@params_table$est[BioGeoBEARS_model_object@params_table$type == 
                                                                             "free"]
      }
    }
  }
  if (printlevel >= 1) {
    print(BioGeoBEARS_model_object)
  }
  BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
  outputs = BioGeoBEARS_model_object
  if ((is.numeric(BioGeoBEARS_run_object$timeperiods))) {
    return_condlikes_table = BioGeoBEARS_run_object$return_condlikes_table
    calc_ancprobs = BioGeoBEARS_run_object$calc_ancprobs
    fixnode = BioGeoBEARS_run_object$fixnode
    fixlikes = BioGeoBEARS_run_object$fixlikes
    inputs2 = BioGeoBEARS_run_object
    inputs2$BioGeoBEARS_model_object = BioGeoBEARS_model_object
    calc_TTL_loglike_from_condlikes_table = BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table
    model_results = calc_loglike_sp_stratified(tip_condlikes_of_data_on_each_state, 
                                               phy, Qmat = NULL, spPmat = NULL, min_branchlength = min_branchlength, 
                                               return_what = "all", probs_of_states_at_root = NULL, 
                                               rootedge = TRUE, sparse = force_sparse, printlevel = BioGeoBEARS_run_object$printlevel, 
                                               use_cpp = TRUE, input_is_COO = FALSE, spPmat_inputs = NULL, 
                                               cppSpMethod = 3, cluster_already_open = cluster_already_open, 
                                               calc_ancprobs = calc_ancprobs, include_null_range = BioGeoBEARS_run_object$include_null_range, 
                                               fixnode = fixnode, fixlikes = fixlikes, inputs = inputs2, 
                                               allareas = allareas, all_states_list = all_states_list, 
                                               return_condlikes_table = return_condlikes_table, 
                                               calc_TTL_loglike_from_condlikes_table = calc_TTL_loglike_from_condlikes_table)
  }
  else {
    params = BioGeoBEARS_model_object_to_est_params(BioGeoBEARS_model_object)
    calc_ancprobs = BioGeoBEARS_run_object$calc_ancprobs
    model_results = calc_loglike_for_optim(params = params, 
                                           BioGeoBEARS_run_object = BioGeoBEARS_run_object, 
                                           phy = phy, tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                                           print_optim = BioGeoBEARS_run_object$print_optim, 
                                           areas_list = areas_list, states_list = states_list, 
                                           force_sparse = force_sparse, cluster_already_open = cluster_already_open, 
                                           return_what = "all", calc_ancprobs = calc_ancprobs)
  }
  if (cluster_was_open == FALSE) {
    if (exists("cluster_open") && (cluster_open == TRUE)) {
      cat("\n\nStopping cluster with ", num_cores_to_use, 
          " cores.\n\n", sep = "")
      stopCluster(cluster_already_open)
    }
  }
  bears_output = model_results
  bears_output$inputs = BioGeoBEARS_run_object
  bears_output$outputs = outputs
  bears_output$optim_result = optim_result2
  return(bears_output)
}
#<bytecode: 0x115d8fc40>
 # <environment: namespace:BioGeoBEARS>