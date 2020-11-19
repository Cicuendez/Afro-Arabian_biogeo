get_Qmat_COOmat_from_BioGeoBEARS_run_object_nexus <- function (BioGeoBEARS_run_object, BioGeoBEARS_model_object = NULL, 
          numstates = NULL, max_range_size = NULL, include_null_range = TRUE, 
          timeperiod_i = 1) 
{
  setup = "\n\t# Set up a BioGeoBEARS_run_object\n\tBioGeoBEARS_run_object = define_BioGeoBEARS_run()\n\t\n\t# Produces stop error\n\tget_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object, BioGeoBEARS_model_object=NULL, numstates=NULL, include_null_range=TRUE, max_range_size=NULL)\n\t\n\t# Works:\n\t# Get the matrices (trivial case, non-stratified Psychotria example)\n\treturned_mats = get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object, BioGeoBEARS_model_object=NULL, numstates=NULL, max_range_size=4)\n\treturned_mats\n\tnames(returned_mats)\n\t\n\t# Also works:\n\tBioGeoBEARS_run_object$max_range_size = 4\n\treturned_mats = get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object, BioGeoBEARS_model_object=NULL, numstates=NULL, max_range_size=NULL)\n\treturned_mats\n\tnames(returned_mats)\n\t"
  inputs = BioGeoBEARS_run_object
  if (isblank_TF(inputs$max_range_size) == FALSE) {
    max_range_size = inputs$max_range_size
  }
  else {
    if (isblank_TF(max_range_size) == FALSE) {
      inputs$max_range_size = max_range_size
    }
  }
  if (is.null(max_range_size)) {
    txt = "STOP ERROR in get_Qmat_COOmat_from_BioGeoBEARS_run_object(): either max_range_size or BioGeoBEARS_run_object$max_range_size must be specified."
    cat("\n\n")
    cat(txt)
    stop(txt)
    cat("\n\n")
  }
  if (is.null(BioGeoBEARS_model_object) == TRUE) {
    BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
  }
  else {
    BioGeoBEARS_model_object = BioGeoBEARS_model_object
  }
  print_optim = inputs$print_optim
  if (inputs$use_detection_model == FALSE) {
    tipranges = getranges_from_LagrangePHYLIP(lgdata_fn = np(inputs$geogfn))
  }
  if (inputs$use_detection_model == TRUE) {
    tipranges = tipranges_from_detects_fn(detects_fn = inputs$detects_fn)
  }
  speedup = inputs$speedup
  areanames = getareas_from_tipranges_object(tipranges)
  areas = areanames
  areas_list = seq(0, length(areanames) - 1, 1)
  if (is.null(numstates)) {
    numstates = numstates_from_numareas(numareas = length(areanames), 
                                        maxareas = inputs$max_range_size, include_null_range = include_null_range)
  }
  if (is.na(inputs$max_range_size)) {
    if (is.null(inputs$states_list)) {
      max_range_size = length(areas)
    }
    else {
      max_range_size = max(sapply(X = inputs$states_list, 
                                  FUN = length), na.rm = TRUE)
    }
  }
  else {
    max_range_size = inputs$max_range_size
  }
  max_numareas = max_range_size
  TF = (rowSums(dfnums_to_numeric(tipranges@df))) > max_range_size
  if (sum(TF, na.rm = TRUE) > 0) {
    cat("\n\nERROR: Tips with ranges too big:\n", sep = "")
    print(dfnums_to_numeric(tipranges@df)[TF, ])
    cat("\n\nCheck your input geography file!\n", sep = "")
    txt = paste("ERROR: Some tips (listed above) have range sizes larger than ", 
                max_range_size, sep = "")
    stop(txt)
  }
  newstrat = TRUE
  if ((is.null(BioGeoBEARS_run_object$lists_of_states_lists_0based) == 
       FALSE) && (newstrat == TRUE)) {
    area_nums = sort(unique(unlist(inputs$lists_of_states_lists_0based)))
    area_nums
    state_indices_0based_all_timeperiods = unique(unlist(inputs$lists_of_states_lists_0based, 
                                                         recursive = FALSE))
    state_indices_0based_all_timeperiods = sort_list_of_lists_of_numbers(state_indices_0based_all_timeperiods)
    states_list_this_timeperiod = inputs$lists_of_states_lists_0based[[timeperiod_i]]
    states_allowed_this_timeperiod_TF = state_indices_0based_all_timeperiods %in% 
      states_list_this_timeperiod
    states_allowed_this_timeperiod_TF
    states_list = states_list_this_timeperiod
  }
  if ((is.null(BioGeoBEARS_run_object$lists_of_states_lists_0based) == 
       TRUE) || (newstrat == FALSE)) {
    if (is.null(BioGeoBEARS_run_object$states_list)) {
      states_list = rcpp_areas_list_to_states_list(areas = areas, 
                                                   maxareas = max_range_size, include_null_range = include_null_range)
      states_list
    }
    else {
      states_list = BioGeoBEARS_run_object$states_list
    }
  }
  states_list_0based = states_list
  ranges_list = NULL
  for (i in 1:length(states_list_0based)) {
    if ((length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]]))) {
      tmprange = "_"
    }
    else {
      tmprange = paste(areas[states_list_0based[[i]] + 
                               1], collapse = "")
    }
    ranges_list = c(ranges_list, tmprange)
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
  if (force_sparse == TRUE) {
    cat("\nNote: force_sparse is set to TRUE; length(states_list)=", 
        length(states_list), "\n", sep = "")
  }
  trfn = np(inputs$trfn)
  phy = check_trfn_nexus(trfn = trfn)
  inputs = readfiles_BioGeoBEARS_run_nexus(inputs = inputs)
  d = BioGeoBEARS_model_object@params_table["d", "est"]
  e = BioGeoBEARS_model_object@params_table["e", "est"]
  a = BioGeoBEARS_model_object@params_table["a", "est"]
  b = BioGeoBEARS_model_object@params_table["b", "est"]
  phy$edge.length = phy$edge.length^b
  areas = areas_list
  if ((is.null(BioGeoBEARS_run_object$list_of_distances_mats) == 
       FALSE)) {
    distances_mat = BioGeoBEARS_run_object$list_of_distances_mats[[timeperiod_i]]
  }
  else {
    distances_mat = matrix(1, nrow = length(areas), ncol = length(areas))
  }
  x = BioGeoBEARS_model_object@params_table["x", "est"]
  dispersal_multipliers_matrix = distances_mat^x
  if ((is.null(BioGeoBEARS_run_object$list_of_envdistances_mats) == 
       FALSE)) {
    envdistances_mat = BioGeoBEARS_run_object$list_of_envdistances_mats[[timeperiod_i]]
  }
  else {
    envdistances_mat = matrix(1, nrow = length(areas), ncol = length(areas))
  }
  n = BioGeoBEARS_model_object@params_table["n", "est"]
  dispersal_multipliers_matrix = dispersal_multipliers_matrix * 
    envdistances_mat^n
  if ((is.null(BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats) == 
       FALSE)) {
    manual_dispersal_multipliers_matrix = BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats[[timeperiod_i]]
  }
  else {
    manual_dispersal_multipliers_matrix = matrix(1, nrow = length(areas), 
                                                 ncol = length(areas))
  }
  w = BioGeoBEARS_model_object@params_table["w", "est"]
  dispersal_multipliers_matrix = dispersal_multipliers_matrix * 
    manual_dispersal_multipliers_matrix^w
  dmat_times_d = dispersal_multipliers_matrix * matrix(d, nrow = length(areas), 
                                                       ncol = length(areas))
  amat = dispersal_multipliers_matrix * matrix(a, nrow = length(areas), 
                                               ncol = length(areas))
  if ((is.null(BioGeoBEARS_run_object$list_of_area_of_areas) == 
       FALSE)) {
    area_of_areas = BioGeoBEARS_run_object$list_of_area_of_areas[[timeperiod_i]]
  }
  else {
    area_of_areas = rep(1, length(areas))
  }
  u = BioGeoBEARS_model_object@params_table["u", "est"]
  extinction_modifier_list = area_of_areas^(1 * u)
  elist = extinction_modifier_list * rep(e, length(areas))
  Qmat = rcpp_states_list_to_DEmat(areas_list = areas_list, 
                                   states_list = states_list, dmat = dmat_times_d, elist = elist, 
                                   amat = amat, include_null_range = TRUE, normalize_TF = TRUE, 
                                   makeCOO_TF = force_sparse)
  j = BioGeoBEARS_model_object@params_table["j", "est"]
  ysv = BioGeoBEARS_model_object@params_table["ysv", "est"]
  ys = BioGeoBEARS_model_object@params_table["ys", "est"]
  v = BioGeoBEARS_model_object@params_table["v", "est"]
  y = BioGeoBEARS_model_object@params_table["y", "est"]
  s = BioGeoBEARS_model_object@params_table["s", "est"]
  sum_SPweights = y + s + j + v
  maxent_constraint_01 = BioGeoBEARS_model_object@params_table["mx01", 
                                                               "est"]
  maxent_constraint_01v = BioGeoBEARS_model_object@params_table["mx01v", 
                                                                "est"]
  maxent01s_param = BioGeoBEARS_model_object@params_table["mx01s", 
                                                          "est"]
  maxent01v_param = BioGeoBEARS_model_object@params_table["mx01v", 
                                                          "est"]
  maxent01j_param = BioGeoBEARS_model_object@params_table["mx01j", 
                                                          "est"]
  maxent01y_param = BioGeoBEARS_model_object@params_table["mx01y", 
                                                          "est"]
  spPmat_inputs = NULL
  dmat = dispersal_multipliers_matrix
  spPmat_inputs$dmat = dmat
  states_indices = states_list
  if (include_null_range == TRUE) {
    states_indices[1] = NULL
  }
  spPmat_inputs$l = states_indices
  spPmat_inputs$s = s
  spPmat_inputs$v = v
  spPmat_inputs$j = j
  spPmat_inputs$y = y
  spPmat_inputs$maxent01s_param = maxent01s_param
  spPmat_inputs$maxent01v_param = maxent01v_param
  spPmat_inputs$maxent01j_param = maxent01j_param
  spPmat_inputs$maxent01y_param = maxent01y_param
  calc_ancprobs = TRUE
  cppSpMethod = 3
  printmat = FALSE
  if (is.null(spPmat_inputs) == FALSE) {
    spPmat_inputs$l[spPmat_inputs$l == c("_")] = NULL
    spPmat_inputs$l[spPmat_inputs$l == c("-")] = NULL
    spPmat_inputs$l[spPmat_inputs$l == c("-1")] = NULL
  }
  l = spPmat_inputs$l
  s = spPmat_inputs$s
  v = spPmat_inputs$v
  j = spPmat_inputs$j
  y = spPmat_inputs$y
  numareas = max(sapply(X = spPmat_inputs$l, FUN = length), 
                 na.rm = TRUE) + 0
  maxent01s_param = spPmat_inputs$maxent01s_param
  maxent01v_param = spPmat_inputs$maxent01v_param
  maxent01j_param = spPmat_inputs$maxent01j_param
  maxent01y_param = spPmat_inputs$maxent01y_param
  maxent01s = relative_probabilities_of_subsets(max_numareas = numareas, 
                                                maxent_constraint_01 = maxent01s_param, NA_val = 0)
  maxent01v = relative_probabilities_of_vicariants(max_numareas = numareas, 
                                                   maxent_constraint_01v = maxent01v_param, NA_val = 0)
  maxent01j = relative_probabilities_of_subsets(max_numareas = numareas, 
                                                maxent_constraint_01 = maxent01j_param, NA_val = 0)
  maxent01y = relative_probabilities_of_subsets(max_numareas = numareas, 
                                                maxent_constraint_01 = maxent01y_param, NA_val = 0)
  maxprob_as_function_of_ancsize_and_decsize = mapply(FUN = max, 
                                                      maxent01s, maxent01v, maxent01j, maxent01y, MoreArgs = list(na.rm = TRUE))
  maxprob_as_function_of_ancsize_and_decsize = matrix(data = maxprob_as_function_of_ancsize_and_decsize, 
                                                      nrow = nrow(maxent01s), ncol = ncol(maxent01s))
  maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize > 
                                               0] = 1
  maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize <= 
                                               0] = 0
  max_minsize_as_function_of_ancsize = apply(X = maxprob_as_function_of_ancsize_and_decsize, 
                                             MARGIN = 1, FUN = maxsize)
  tmpca_1 = rep(1, (numstates - 1))
  tmpcb_1 = rep(1, (numstates - 1))
  COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs = tmpca_1, 
                                                                 Rcpp_rightprobs = tmpcb_1, l = l, s = s, v = v, j = j, 
                                                                 y = y, dmat = dispersal_multipliers_matrix, maxent01s = maxent01s, 
                                                                 maxent01v = maxent01v, maxent01j = maxent01j, maxent01y = maxent01y, 
                                                                 max_minsize_as_function_of_ancsize = max_minsize_as_function_of_ancsize, 
                                                                 printmat = printmat)
  Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar = COO_weights_columnar)
  Qmat
  rowSums(Qmat)
  max(rowSums(Qmat))
  COO_weights_columnar
  Rsp_rowsums
  returned_mats = NULL
  returned_mats$states_list = states_list
  returned_mats$ranges_list = ranges_list
  returned_mats$spPmat_inputs = spPmat_inputs
  returned_mats$areas_list = areas_list
  returned_mats$areanames = areanames
  returned_mats$dmat = dmat
  returned_mats$Qmat = Qmat
  returned_mats$COO_weights_columnar = COO_weights_columnar
  returned_mats$Rsp_rowsums = Rsp_rowsums
  return(returned_mats)
}
