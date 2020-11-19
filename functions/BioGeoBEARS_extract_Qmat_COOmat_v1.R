

# Get the Q matrix and cladogenesis mode from the
# BioGeoBEARS_results_object "res"
get_Qmat_COOmat_from_res <- function(res, numstates=NULL, include_null_range=TRUE)
	{
	
	#######################################################
	# Load the model object
	#######################################################
	inputs = res$inputs
	output_model_object = res$output
	
	# Should the optim run be printed?
	print_optim = inputs$print_optim


	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(inputs$geogfn))
	
	
	# Should we do optimx speedup?
	speedup = inputs$speedup
	

	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes
	
	# Calculate the number of states, if needed
	if (is.null(numstates))
		{
		numstates = numstates_from_numareas(numareas=length(areas), maxareas=inputs$max_range_size, include_null_range=include_null_range)
		}
	
	# Change the names to tipranges@df:
	# this doesn't make sense if areas_list is 0-based indexes
	names(tipranges@df) = areas_list

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.na(inputs$max_range_size))
		{
		if (is.null(inputs$states_list))
			{
			# Maximum range size is all areas
			max_range_size = length(areas)
			} else {
			# If not NA
			# Get max rangesize from states list
			max_range_size = max(sapply(X=inputs$states_list, FUN=length), na.rm=TRUE)
			}
		} else {
		# Maximum range size hard-coded
		max_range_size = inputs$max_range_size
		}
	max_numareas = max_range_size
	
	#######################################################
	# Check that no tips have larger ranges than you allowed
	#######################################################
	TF = (rowSums(dfnums_to_numeric(tipranges@df))) > max_range_size
	if (sum(TF, na.rm=TRUE) > 0)
		{
		cat("\n\nERROR: Tips with ranges too big:\n", sep="")
		print(dfnums_to_numeric(tipranges@df)[TF, ])
		cat("\n\nCheck your input geography file!\n", sep="")
		txt = paste("ERROR: Some tips (listed above) have range sizes larger than ", max_range_size, sep="")
		stop(txt)
		}





	# Take the list of areas, and get list of possible states
	# (the user can manually input states if they like)
	if (is.null(BioGeoBEARS_run_object$states_list))
		{
		states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
		states_list
		#BioGeoBEARS_run_object$states_list = states_list
		#inputs$states_list = states_list
		} else {
		states_list = BioGeoBEARS_run_object$states_list
		#BioGeoBEARS_run_object$states_list = states_list
		#inputs$states_list = states_list
		}

	if (is.na(BioGeoBEARS_run_object$force_sparse))
		{
		if (length(states_list) > 128)
			{
			force_sparse = TRUE
			cat("\nNote: force_sparse being set to TRUE, as length(states_list) > 128\n", sep="")
			} else {
			force_sparse = FALSE
			}
		} else {
		force_sparse = BioGeoBEARS_run_object$force_sparse
		}

	if (force_sparse == TRUE)
		{
		cat("\nNote: force_sparse is set to TRUE; length(states_list)=", length(states_list), "\n", sep="")
		}






	#######################################################
	# Load the phylogenetic tree
	#######################################################
	trfn = np(inputs$trfn)
	phy = read.tree(file=trfn)



	#######################################################
	# Read the stratification/distances input files, if any
	#######################################################
	inputs = readfiles_BioGeoBEARS_run(inputs=inputs)



	#######################################################
	# FROM CALC_LOGLIKE_SP_FOR_OPTIM()
	#######################################################


	# Set the dispersal and extinction rate
	d = output_model_object@params_table["d","est"]
	e = output_model_object@params_table["e","est"]
	a = output_model_object@params_table["a","est"]


	#######################################################
	#######################################################
	# Do branch-length exponentiation if desired
	#######################################################
	#######################################################
	b = output_model_object@params_table["b","est"]
	# Modify the edge.lengths
	phy$edge.length = phy$edge.length ^ b
	# Make sure this doesn't duplicate a previous "^b", e.g.
	# the summarization step in bears_optim_run

	#######################################################
	#######################################################
	# Do distance-dependence and dispersal multipliers matrix
	#######################################################
	#######################################################
	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	areas = areas_list
	
	# If there is a distance matrix, use the first one 
	# (non-stratified analysis, here)
	if ( (is.null(BioGeoBEARS_run_object$list_of_distances_mats) == FALSE))
		{
		distances_mat = BioGeoBEARS_run_object$list_of_distances_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
		}

	# Get the exponent on distance, apply to distances matrix
	x = output_model_object@params_table["x","est"]
	dispersal_multipliers_matrix = distances_mat ^ x

	# Apply manual dispersal multipliers, if any
	# If there is a manual dispersal multipliers matrix, use the first one 
	# (non-stratified analysis, here)
	if ( (is.null(BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats) == FALSE))
		{
		manual_dispersal_multipliers_matrix = BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		manual_dispersal_multipliers_matrix = matrix(1, nrow=length(areas), ncol=length(areas))
		}
	
	# Apply element-wise
	dispersal_multipliers_matrix = dispersal_multipliers_matrix * manual_dispersal_multipliers_matrix

	#######################################################
	# multiply parameter d by dispersal_multipliers_matrix
	#######################################################
	dmat = dispersal_multipliers_matrix * matrix(d, nrow=length(areas), ncol=length(areas))
	amat = dispersal_multipliers_matrix * matrix(a, nrow=length(areas), ncol=length(areas))




	#######################################################
	#######################################################
	# Do area-dependence and extinction multipliers list
	#######################################################
	#######################################################
	if ( (is.null(BioGeoBEARS_run_object$list_of_area_of_areas) == FALSE))
		{
		area_of_areas = BioGeoBEARS_run_object$list_of_area_of_areas[[1]]
		} else {
		# Default is all areas effectively equidistant
		area_of_areas = rep(1, length(areas))
		}
		
	# Get the exponent on extinction, apply to extinction modifiers	
	u = output_model_object@params_table["u","est"]
	extinction_modifier_list = area_of_areas ^ (1 * u)
	
	# Apply to extinction rate
	elist = extinction_modifier_list * rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	# someday we'll have to put "a" (anagenic range-switching) in here...
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, amat, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	#######################################################
	# Cladogenic model
	#######################################################
	j = output_model_object@params_table["j","est"]
	ysv = output_model_object@params_table["ysv","est"]
	ys = output_model_object@params_table["ys","est"]
	v = output_model_object@params_table["v","est"]
	y = output_model_object@params_table["y","est"]
	s = output_model_object@params_table["s","est"]
	sum_SPweights = y + s + j + v


	maxent_constraint_01 = output_model_object@params_table["mx01","est"]
	
	# Text version of speciation matrix	
	maxent_constraint_01v = output_model_object@params_table["mx01v","est"]
	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = output_model_object@params_table["mx01s","est"]
	maxent01v_param = output_model_object@params_table["mx01v","est"]
	maxent01j_param = output_model_object@params_table["mx01j","est"]
	maxent01y_param = output_model_object@params_table["mx01y","est"]


	# Cladogenesis model inputs
	spPmat_inputs = NULL

	# Note that this gets the dispersal multipliers matrix, which is applied to 
	# e.g. the j events, NOT the dmat above which is d*dispersal_multipliers_matrix
	spPmat_inputs$dmat = dispersal_multipliers_matrix

	states_indices = states_list
	states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
	spPmat_inputs$l = states_indices
	spPmat_inputs$s = s
	spPmat_inputs$v = v
	spPmat_inputs$j = j
	spPmat_inputs$y = y
	spPmat_inputs$maxent01s_param = maxent01s_param
	spPmat_inputs$maxent01v_param = maxent01v_param
	spPmat_inputs$maxent01j_param = maxent01j_param
	spPmat_inputs$maxent01y_param = maxent01y_param





	#######################################################
	# From calc_loglike_sp()
	#######################################################

	# defaults
	calc_ancprobs = TRUE
	cppSpMethod = 3
	printmat = FALSE

	
	# Fix "l" (which is states_indices, i.e. a list of lists of state_indices)
	# so that there is no "null geographic range", i.e. no "_", no "-", no "NA"
	if ( is.null(spPmat_inputs)==FALSE )
		{
		spPmat_inputs$l[spPmat_inputs$l == c("_")] = NULL
		spPmat_inputs$l[spPmat_inputs$l == c("-")] = NULL
		spPmat_inputs$l[spPmat_inputs$l == c("-1")] = NULL
		#spPmat_inputs$l[spPmat_inputs$l == c(-1)] = NULL
		}



	# (above) Fix "l" (which is states_indices, i.e. a list of lists of state_indices)
	# so that there is no "null geographic range", i.e. no "_", no "-", no "NA"
	l = spPmat_inputs$l		# states_indices
	s = spPmat_inputs$s
	v = spPmat_inputs$v
	j = spPmat_inputs$j
	y = spPmat_inputs$y


	dmat = spPmat_inputs$dmat
	
	# Take the max of the indices of the possible areas, and add 1
	# numareas = max(unlist(spPmat_inputs$l), na.rm=TRUE) + 1 # old, bogus
	numareas = max(sapply(X=spPmat_inputs$l, FUN=length), na.rm=TRUE) + 0
	
	maxent01s_param = spPmat_inputs$maxent01s_param
	maxent01v_param = spPmat_inputs$maxent01v_param
	maxent01j_param = spPmat_inputs$maxent01j_param
	maxent01y_param = spPmat_inputs$maxent01y_param
	
	maxent01s = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01s_param, NA_val=0)
	maxent01v = relative_probabilities_of_vicariants(max_numareas=numareas, maxent_constraint_01v=maxent01v_param, NA_val=0)
	maxent01j = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01j_param, NA_val=0)
	maxent01y = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01y_param, NA_val=0)

	# You really need a list of sizes here:
	
	# Matrix of probs for each ancsize
	maxprob_as_function_of_ancsize_and_decsize = mapply(FUN=max, maxent01s, maxent01v, maxent01j, maxent01y, MoreArgs=list(na.rm=TRUE))
	maxprob_as_function_of_ancsize_and_decsize = matrix(data=maxprob_as_function_of_ancsize_and_decsize, nrow=nrow(maxent01s), ncol=ncol(maxent01s))
	maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize > 0] = 1
	maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize <= 0] = 0
	
	# Now, go through, and make a list of the max minsize for each decsize
	max_minsize_as_function_of_ancsize = apply(X=maxprob_as_function_of_ancsize_and_decsize, MARGIN=1, FUN=maxsize)


	tmpca_1 = rep(1, (numstates-1))
	tmpcb_1 = rep(1, (numstates-1))


	COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=printmat)
	
	# combine with C++ function
	# This causes an error with spPmat=NULL; spPmat_inputs=NULL; use_cpp=TRUE; sparse=FALSE
	# i.e. gives 16 states with a 0 on the end, rather than 15 states
	#Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar, numstates=numstates)
	
	# This gives 15 states
	Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar)


	Qmat
	rowSums(Qmat)	# yep, they sum to 0
	max(rowSums(Qmat))
	COO_weights_columnar
	Rsp_rowsums
	
	
	returned_mats = NULL
	returned_mats$Qmat = Qmat
	returned_mats$COO_weights_columnar = COO_weights_columnar
	returned_mats$Rsp_rowsums = Rsp_rowsums
	
	return(returned_mats)
	}








# Get the Q matrix and cladogenesis mode from the
# BioGeoBEARS_run_object
get_Qmat_COOmat_from_BioGeoBEARS_run_object <- function(BioGeoBEARS_run_object, numstates=NULL, include_null_range=TRUE)
	{
	#######################################################
	# Load the model object
	#######################################################
	inputs = BioGeoBEARS_run_object
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	
	# Should the optim run be printed?
	print_optim = inputs$print_optim


	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(inputs$geogfn))
	
	
	# Should we do optimx speedup?
	speedup = inputs$speedup
	

	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes

	# Change the names to tipranges@df:
	# this doesn't make sense if areas_list is 0-based indexes
	names(tipranges@df) = areas_list



	# Calculate the number of states, if needed
	if (is.null(numstates))
		{
		numstates = numstates_from_numareas(numareas=length(areas), maxareas=inputs$max_range_size, include_null_range=include_null_range)
		}




	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.na(inputs$max_range_size))
		{
		if (is.null(inputs$states_list))
			{
			# Maximum range size is all areas
			max_range_size = length(areas)
			} else {
			# If not NA
			# Get max rangesize from states list
			max_range_size = max(sapply(X=inputs$states_list, FUN=length), na.rm=TRUE)
			}
		} else {
		# Maximum range size hard-coded
		max_range_size = inputs$max_range_size
		}
	max_numareas = max_range_size
	
	#######################################################
	# Check that no tips have larger ranges than you allowed
	#######################################################
	TF = (rowSums(dfnums_to_numeric(tipranges@df))) > max_range_size
	if (sum(TF, na.rm=TRUE) > 0)
		{
		cat("\n\nERROR: Tips with ranges too big:\n", sep="")
		print(dfnums_to_numeric(tipranges@df)[TF, ])
		cat("\n\nCheck your input geography file!\n", sep="")
		txt = paste("ERROR: Some tips (listed above) have range sizes larger than ", max_range_size, sep="")
		stop(txt)
		}





	# Take the list of areas, and get list of possible states
	# (the user can manually input states if they like)
	if (is.null(BioGeoBEARS_run_object$states_list))
		{
		states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
		states_list
		#BioGeoBEARS_run_object$states_list = states_list
		#inputs$states_list = states_list
		} else {
		states_list = BioGeoBEARS_run_object$states_list
		#BioGeoBEARS_run_object$states_list = states_list
		#inputs$states_list = states_list
		}

	if (is.na(BioGeoBEARS_run_object$force_sparse))
		{
		if (length(states_list) > 128)
			{
			force_sparse = TRUE
			cat("\nNote: force_sparse being set to TRUE, as length(states_list) > 128\n", sep="")
			} else {
			force_sparse = FALSE
			}
		} else {
		force_sparse = BioGeoBEARS_run_object$force_sparse
		}

	if (force_sparse == TRUE)
		{
		cat("\nNote: force_sparse is set to TRUE; length(states_list)=", length(states_list), "\n", sep="")
		}






	#######################################################
	# Load the phylogenetic tree
	#######################################################
	trfn = np(inputs$trfn)
	phy = read.tree(file=trfn)



	#######################################################
	# Read the stratification/distances input files, if any
	#######################################################
	inputs = readfiles_BioGeoBEARS_run(inputs=inputs)



	#######################################################
	# FROM CALC_LOGLIKE_SP_FOR_OPTIM()
	#######################################################


	# Set the dispersal and extinction rate
	d = BioGeoBEARS_model_object@params_table["d","est"]
	e = BioGeoBEARS_model_object@params_table["e","est"]
	a = BioGeoBEARS_model_object@params_table["a","est"]


	#######################################################
	#######################################################
	# Do branch-length exponentiation if desired
	#######################################################
	#######################################################
	b = BioGeoBEARS_model_object@params_table["b","est"]
	# Modify the edge.lengths
	phy$edge.length = phy$edge.length ^ b
	# Make sure this doesn't duplicate a previous "^b", e.g.
	# the summarization step in bears_optim_run

	#######################################################
	#######################################################
	# Do distance-dependence and dispersal multipliers matrix
	#######################################################
	#######################################################
	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	areas = areas_list
	
	# If there is a distance matrix, use the first one 
	# (non-stratified analysis, here)
	if ( (is.null(BioGeoBEARS_run_object$list_of_distances_mats) == FALSE))
		{
		distances_mat = BioGeoBEARS_run_object$list_of_distances_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
		}

	# Get the exponent on distance, apply to distances matrix
	x = BioGeoBEARS_model_object@params_table["x","est"]
	dispersal_multipliers_matrix = distances_mat ^ x

	# Apply manual dispersal multipliers, if any
	# If there is a manual dispersal multipliers matrix, use the first one 
	# (non-stratified analysis, here)
	if ( (is.null(BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats) == FALSE))
		{
		manual_dispersal_multipliers_matrix = BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		manual_dispersal_multipliers_matrix = matrix(1, nrow=length(areas), ncol=length(areas))
		}
	
	# Apply element-wise
	dispersal_multipliers_matrix = dispersal_multipliers_matrix * manual_dispersal_multipliers_matrix

	#######################################################
	# multiply parameter d by dispersal_multipliers_matrix
	#######################################################
	dmat = dispersal_multipliers_matrix * matrix(d, nrow=length(areas), ncol=length(areas))
	amat = dispersal_multipliers_matrix * matrix(a, nrow=length(areas), ncol=length(areas))




	#######################################################
	#######################################################
	# Do area-dependence and extinction multipliers list
	#######################################################
	#######################################################
	if ( (is.null(BioGeoBEARS_run_object$list_of_area_of_areas) == FALSE))
		{
		area_of_areas = BioGeoBEARS_run_object$list_of_area_of_areas[[1]]
		} else {
		# Default is all areas effectively equidistant
		area_of_areas = rep(1, length(areas))
		}
		
	# Get the exponent on extinction, apply to extinction modifiers	
	u = BioGeoBEARS_model_object@params_table["u","est"]
	extinction_modifier_list = area_of_areas ^ (1 * u)
	
	# Apply to extinction rate
	elist = extinction_modifier_list * rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	# someday we'll have to put "a" (anagenic range-switching) in here...
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, amat, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	#######################################################
	# Cladogenic model
	#######################################################
	j = BioGeoBEARS_model_object@params_table["j","est"]
	ysv = BioGeoBEARS_model_object@params_table["ysv","est"]
	ys = BioGeoBEARS_model_object@params_table["ys","est"]
	v = BioGeoBEARS_model_object@params_table["v","est"]
	y = BioGeoBEARS_model_object@params_table["y","est"]
	s = BioGeoBEARS_model_object@params_table["s","est"]
	sum_SPweights = y + s + j + v


	maxent_constraint_01 = BioGeoBEARS_model_object@params_table["mx01","est"]
	
	# Text version of speciation matrix	
	maxent_constraint_01v = BioGeoBEARS_model_object@params_table["mx01v","est"]
	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = BioGeoBEARS_model_object@params_table["mx01s","est"]
	maxent01v_param = BioGeoBEARS_model_object@params_table["mx01v","est"]
	maxent01j_param = BioGeoBEARS_model_object@params_table["mx01j","est"]
	maxent01y_param = BioGeoBEARS_model_object@params_table["mx01y","est"]


	# Cladogenesis model inputs
	spPmat_inputs = NULL

	# Note that this gets the dispersal multipliers matrix, which is applied to 
	# e.g. the j events, NOT the dmat above which is d*dispersal_multipliers_matrix
	spPmat_inputs$dmat = dispersal_multipliers_matrix

	states_indices = states_list
	states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
	spPmat_inputs$l = states_indices
	spPmat_inputs$s = s
	spPmat_inputs$v = v
	spPmat_inputs$j = j
	spPmat_inputs$y = y
	spPmat_inputs$maxent01s_param = maxent01s_param
	spPmat_inputs$maxent01v_param = maxent01v_param
	spPmat_inputs$maxent01j_param = maxent01j_param
	spPmat_inputs$maxent01y_param = maxent01y_param





	#######################################################
	# From calc_loglike_sp()
	#######################################################

	# defaults
	calc_ancprobs = TRUE
	cppSpMethod = 3
	printmat = FALSE

	
	# Fix "l" (which is states_indices, i.e. a list of lists of state_indices)
	# so that there is no "null geographic range", i.e. no "_", no "-", no "NA"
	if ( is.null(spPmat_inputs)==FALSE )
		{
		spPmat_inputs$l[spPmat_inputs$l == c("_")] = NULL
		spPmat_inputs$l[spPmat_inputs$l == c("-")] = NULL
		spPmat_inputs$l[spPmat_inputs$l == c("-1")] = NULL
		#spPmat_inputs$l[spPmat_inputs$l == c(-1)] = NULL
		}



	# (above) Fix "l" (which is states_indices, i.e. a list of lists of state_indices)
	# so that there is no "null geographic range", i.e. no "_", no "-", no "NA"
	l = spPmat_inputs$l		# states_indices
	s = spPmat_inputs$s
	v = spPmat_inputs$v
	j = spPmat_inputs$j
	y = spPmat_inputs$y


	dmat = spPmat_inputs$dmat
	
	# Take the max of the indices of the possible areas, and add 1
	# numareas = max(unlist(spPmat_inputs$l), na.rm=TRUE) + 1 # old, bogus
	numareas = max(sapply(X=spPmat_inputs$l, FUN=length), na.rm=TRUE) + 0
	
	maxent01s_param = spPmat_inputs$maxent01s_param
	maxent01v_param = spPmat_inputs$maxent01v_param
	maxent01j_param = spPmat_inputs$maxent01j_param
	maxent01y_param = spPmat_inputs$maxent01y_param
	
	maxent01s = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01s_param, NA_val=0)
	maxent01v = relative_probabilities_of_vicariants(max_numareas=numareas, maxent_constraint_01v=maxent01v_param, NA_val=0)
	maxent01j = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01j_param, NA_val=0)
	maxent01y = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01y_param, NA_val=0)

	# You really need a list of sizes here:
	
	# Matrix of probs for each ancsize
	maxprob_as_function_of_ancsize_and_decsize = mapply(FUN=max, maxent01s, maxent01v, maxent01j, maxent01y, MoreArgs=list(na.rm=TRUE))
	maxprob_as_function_of_ancsize_and_decsize = matrix(data=maxprob_as_function_of_ancsize_and_decsize, nrow=nrow(maxent01s), ncol=ncol(maxent01s))
	maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize > 0] = 1
	maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize <= 0] = 0
	
	# Now, go through, and make a list of the max minsize for each decsize
	max_minsize_as_function_of_ancsize = apply(X=maxprob_as_function_of_ancsize_and_decsize, MARGIN=1, FUN=maxsize)


	tmpca_1 = rep(1, (numstates-1))
	tmpcb_1 = rep(1, (numstates-1))


	COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=printmat)
	
	# combine with C++ function
	# This causes an error with spPmat=NULL; spPmat_inputs=NULL; use_cpp=TRUE; sparse=FALSE
	# i.e. gives 16 states with a 0 on the end, rather than 15 states
	#Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar, numstates=numstates)
	
	# This gives 15 states
	Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar)


	Qmat
	rowSums(Qmat)	# yep, they sum to 0
	max(rowSums(Qmat))
	COO_weights_columnar
	Rsp_rowsums
	
	
	returned_mats = NULL
	returned_mats$Qmat = Qmat
	returned_mats$COO_weights_columnar = COO_weights_columnar
	returned_mats$Rsp_rowsums = Rsp_rowsums
	
	return(returned_mats)
	}