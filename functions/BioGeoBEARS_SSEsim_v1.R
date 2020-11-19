
#######################################################
# Set up the default inputs to SSEsim
#######################################################
SSEsim_setup_inputs <- function(SSEmodel=NULL, BioGeoBEARS_run_object=NULL)
	{
	# This object holds the SSE sim inputs
	SSEsim_inputs = NULL
	
	# Set up the Birth Rate, Death Rate, etc.
	# These are the parameters for the tree, and for 
	# dependence of speciation/extinction on states
	if (is.null(SSEmodel))
		{
		SSEmodel = NULL
		
		# Birth rate (lambda)
		# This is the ML estimate under Yule on the Psychotria tree
		SSEmodel$brate = 0.3289132
		
		# Death rate (omega)
		SSEmodel$drate = SSEmodel$brate/3
		
		# Exponent on rangesize multiplier
		# Positive exponent means larger ranges have
		# increased speciation rates
		SSEmodel$rangesize_b_exponent = 1	

		# Exponent on rangesize multiplier
		# More negative exponents mean larger ranges have
		# decreases extinction rates
		SSEmodel$rangesize_d_exponent = -1

		# To get the values out:		
# 		brate = SSEmodel$brate
# 		drate = SSEmodel$drate
# 		rangesize_b_exponent = SSEmodel$rangesize_b_exponent
# 		rangesize_d_exponent = SSEmodel$rangesize_d_exponent

		}
	
	
	# Set up the BioGeoBEARS model (the model of range evolution)
	if (is.null(BioGeoBEARS_run_object))
		{
		# Use a default BioGeoBEARS_run_object, but specify some different parameters
		BioGeoBEARS_run_object = define_BioGeoBEARS_run()
		BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
		include_null_range = BioGeoBEARS_run_object$include_null_range

		#######################################################
		# Define the areas and states
		#######################################################
		# Get geographic ranges at tips
		tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(BioGeoBEARS_run_object$geogfn))

		# Get the list of geographic areas
		areas = getareas_from_tipranges_object(tipranges)
		areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes
		areanames = areas
		areas
		areas_list

		max_range_size = length(areas)	# if Psychotria, this is 4
		BioGeoBEARS_run_object$max_range_size = max_range_size
		
		states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
		states_list

		state_indices_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=BioGeoBEARS_run_object$max_range_size, include_null_range=TRUE)
		state_indices_0based

		# Get the ranges
		ranges_list = areas_list_to_states_list_new(areas=areas, maxareas=length(areas),
		include_null_range=include_null_range, split_ABC=FALSE)
		ranges_list
		ranges = unlist(ranges_list)
		ranges
		

		###############################################
		# Default params for SSEsim
		###############################################
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = 0.03
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = 0.03

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = 0.03
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = 0.03

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.1

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","init"] = 0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","est"] = 0

		# Update linked parameters
		BioGeoBEARS_run_object$BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_run_object$BioGeoBEARS_model_object)
		BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object

		} else {
		# Read the input files, if any
		BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
		include_null_range = BioGeoBEARS_run_object$include_null_range

		#######################################################
		# Define the areas and states
		#######################################################
		# Get geographic ranges at tips
		tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(BioGeoBEARS_run_object$geogfn))

		# Get the list of geographic areas
		areas = getareas_from_tipranges_object(tipranges)
		areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes
		areanames = areas
		areas
		areas_list

		max_range_size = BioGeoBEARS_run_object$max_range_size
		
		states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
		states_list

		state_indices_0based = states_list

		# Get the ranges
		ranges_list = areas_list_to_states_list_new(areas=areas, maxareas=length(areas),
		include_null_range=include_null_range, split_ABC=FALSE)
		ranges_list
		ranges = unlist(ranges_list)
		ranges
		


		# Update linked parameters
		BioGeoBEARS_run_object$BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_run_object$BioGeoBEARS_model_object)
		BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
		} # end is.null(BioGeoBEARS_run_object)
		

	
	#######################################################
	# Read in the SSE params
	#######################################################
	brate = SSEmodel$brate
	drate = SSEmodel$drate
	rangesize_b_exponent = SSEmodel$rangesize_b_exponent
	rangesize_d_exponent = SSEmodel$rangesize_d_exponent


	#######################################################
	# Use the anagenetic and cladogenetic parameters to make a Qmat and COOmat
	#######################################################	
	# Set the dispersal and extinction rate
	d = BioGeoBEARS_model_object@params_table["d","est"]
	e = BioGeoBEARS_model_object@params_table["e","est"]
	a = BioGeoBEARS_model_object@params_table["a","est"]
	
	# More inputs
	force_sparse = BioGeoBEARS_run_object$force_sparse
	areas = areas_list
	
	# Calculate the dispersal_multipliers_matrix
	dispersal_multipliers_matrix = dispersal_multipliers_matrix_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object)
	
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
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, amat, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)


	#######################################################
	# Get the cladogenesis matrix
	#######################################################
	spPmat_inputs = spPmat_inputs_from_BioGeoBEARS_model_object(BioGeoBEARS_run_object, states_list, dispersal_multipliers_matrix)
	
	COOmat_Rsp_rowsums = spPmat_inputs_to_COO_weights_columnar(spPmat_inputs, cppSpMethod=3, numstates_in_cladogenesis_matrix=length(states_list), printmat=FALSE)
	COO_weights_columnar = COOmat_Rsp_rowsums$COO_weights_columnar
	Rsp_rowsums = COOmat_Rsp_rowsums$Rsp_rowsums
	
	
	# Store the SSEsim inputs
	SSEsim_inputs$BioGeoBEARS_run_object = BioGeoBEARS_run_object
	SSEsim_inputs$Qmat = Qmat
	SSEsim_inputs$COO_weights_columnar = COO_weights_columnar
	SSEsim_inputs$Rsp_rowsums = Rsp_rowsums
	SSEsim_inputs$state_indices_0based = state_indices_0based
	SSEsim_inputs$ranges = ranges
	SSEsim_inputs$areanames = areanames
	SSEsim_inputs$SSEmodel = SSEmodel

	# To get them back out
# 	SSEmodel = SSEsim_inputs$SSEmodel
# 	Qmat = SSEsim_inputs$Qmat
# 	COO_weights_columnar = SSEsim_inputs$COO_weights_columnar
# 	Rsp_rowsums = SSEsim_inputs$Rsp_rowsums
#	state_indices_0based = SSEsim_inputs$state_indices_0based
#	ranges = SSEsim_inputs$ranges

	return(SSEsim_inputs)
	}



SSEsim_run <- function(SSEsim_inputs, rootstate=2, time_stop=1000, taxa_stop=1000, seed=54321, printlevel=2, testwd="~")
	{
	defaults='
	rootstate = 2
	time_stop=100
	taxa_stop=50
	seed=1
	printlevel=4
	testwd=testwd
	'
	
	# Get the model parameters
	Qmat = SSEsim_inputs$Qmat
	COO_weights_columnar = SSEsim_inputs$COO_weights_columnar
	Rsp_rowsums = SSEsim_inputs$Rsp_rowsums
	state_indices_0based = SSEsim_inputs$state_indices_0based
	ranges = SSEsim_inputs$ranges
	SSEmodel = SSEsim_inputs$SSEmodel
	
	numstates = length(state_indices_0based)

	#######################################################
	# Set up the parameters of a full forward simulation
	#######################################################
	brate = SSEmodel$brate
	drate = SSEmodel$drate
	rangesize_b_exponent = SSEmodel$rangesize_b_exponent
	rangesize_d_exponent = SSEmodel$rangesize_d_exponent

	# SSE rates for states
	# Set the birthrate to be a function of the number of areas in each range
	range_sizes = unlist(lapply(X=state_indices_0based, FUN=length))
	brates = brate * (range_sizes ^ rangesize_b_exponent)
	brates

	# Set the deathrate to be a function of the number of areas in each range
	drates = drate * (range_sizes ^ rangesize_d_exponent)
	drates

	# Get the rates from the Qmat
	# Source: https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/R/simulation.R?root=proteinevoutk20
	rates_de = -diag(Qmat) #exponential rate of waiting, diagonal of Qmat
	rates_de


	#######################################################
	# Run the forward simulation!
	#######################################################
	# Set the seed
	seed_input = set.seed(seed)
	
	
	if (time_stop == 0 & taxa_stop == 0)
		{
		stop("Must have stopping criterion\n")
		}

	# Number of tries
	trynum = 0
	
	# Run this while loop until you hit a successful simulation
	while (1)
		{
		trynum = trynum + 1
		
		if (printlevel >= 1)
			{
		cat("\n\n=================================\nTree simulation attempt #", trynum, "\n\n=================================\n\n", sep="")
			}

		###########################################################
		# Set up edge.ranges and begin tree		
		###########################################################
		# "edge" is a 2x2 matrix to start
		# rows are lineages, col#1 is parent node, col#2 is daughter node
		edge <- rbind(c(1, 2), c(1, 3))
		edge_length <- rep(NA, 2)
	
		# have another matrix giving zone 1 (tropics) or zone 2 (temperate)
		# start in the tropics for now
		# edge "state" is the state of the parent node
		edge.zone <- c(1, 1)
	
		# Simulate the first split
		# edge "range" contains the geographical range
		initial_state_1based = rootstate		# Start simulation in Kauai
	
		COO_probs_columnar = COO_weights_columnar
		daughter_states = given_a_starting_state_simulate_split(index_Qmat_0based_of_starting_state=(initial_state_1based-1), COO_probs_columnar=COO_probs_columnar, numstates=numstates)
	
		# Convert to 1-based
		daughter_states = daughter_states + 1
	
		edge_ranges <- rbind(daughter_states[1], daughter_states[2])
	
		# Create array to hold branch lengths
		stem_depth = numeric(2)
		alive_TF = rep(TRUE, 2)
		t = 0.0
		next_node = 4
		###########################################################


		tmp_edgelength = 0
		repeatnum = 0
		de_event_num = 0
		cl_event_num = 0

		# Record anagenetic and cladogenetic events
		table_of_range_change_events = NULL
		table_of_cladogenetic_events = NULL
		
		# A flag to stop if you hit a taxa_stop
		stopped_due_to_taxon_count = FALSE
	
		# Repeat until one of the break statements is reached
		repeat
			{
			repeatnum = repeatnum+1
			popsize = sum(alive_TF)

			if (printlevel >= 2)
				{
				cat("\n")
				cat("Treebuilding step #", repeatnum, ", # species living: ", popsize, "\n", sep="")
				#print(" ")
				cat(paste("Time: ", format(t, digits=3), " my, #alive = ", sum(alive_TF), "\n", sep=""))
				}
		
			# Density-dependence...
			#d = original_d + ((popsize/k) * (original_b-original_d))
			#print(paste("birth rate=b=", b, ", death rate=d=", d))
			#d = (events$top[14]-events$top[7])
			#b = 1-(d)
			#cat("birth events total: ", b, ", death rate=d=", d, "\n", sep="")

			# change things with 0 range to dead
	# 		for (i in 1:nrow(edge.ranges))
	# 			{
	# 			tmp_ranges = edge.ranges[i, ]
	# 			if (sum(tmp_ranges) == 0)
	# 				{
	# 				alive[i] = FALSE
	# 				}
	# 			}
		
		
			# Stop the simulation if everything is dead
			if (sum(alive_TF) == 0)
				{
				break
				}


			# Get waiting time to the next event (scaling factor * number alive); 
			# depends on number of lineages, perhaps times area
		
			# The rates depend on the states at each tip
			# Go through each tip
			current_ntips = length(edge_ranges[alive_TF])
			rates_of_events_per_tip = matrix(data=NA, nrow=current_ntips, ncol=3)
			for (j in 1:current_ntips)
				{
				rates_of_events_per_tip[j,1] = rates_de[edge_ranges[alive_TF][j]]
				rates_of_events_per_tip[j,2] = brates[edge_ranges[alive_TF][j]]
				rates_of_events_per_tip[j,3] = drates[edge_ranges[alive_TF][j]]
				}
		
		
			#rate = sf * get_rate(alive, edge_ranges, rate_calc="per_lineage")
			rate = sum(rates_of_events_per_tip)
			rate
		
			# dt is the amount of time until the next event on the tree
			# (the branches with non-events will be extended)
			# Higher rates result in shorter waiting times
			dt <- rexp(n=1, rate = rate)
		
			# This is the total time since start
			t <- t + dt
		
			# Stop if out of time before the next event
			time_stop_hit = FALSE
			if (time_stop)
				{
					#print(dt)
					#print(edge.zone[alive_TF]==1)
					#print(edge.zone[alive_TF]==2)
					if (t >= time_stop)
						{
						t <- time_stop
						time_stop_hit = TRUE
						break
						}
				}

			# If you didn't hit the time barrier after the last speciation event,
			# stop the simulation if you've reached the # taxa_stop
			# You don't want to stop instantly, because then you will have
			# 2 zero-length branches
			# Instead, stop after (dt/1), i.e. when the next event happens
			# (without actually doing that next event; just extend the branches;
			# do this at the end of the while loop, for both time_stop and taxa_stop)
			if (taxa_stop)
				{
				#print(paste("Stop? #alive=", sum(alive_TF), "taxa_stop=", taxa_stop, sep=""))
				if (sum(alive_TF) >= taxa_stop)
					{
					break					
					}
				}


		
			# Choose which event happened
			event_probs = c(rates_of_events_per_tip) / rate
			event_probs

			# Which event happened?
			eventnum = sample(x=1:length(event_probs), size=1, replace=FALSE, prob=event_probs)
			eventnum
		
			eventnums = matrix(data=1:length(event_probs), nrow=current_ntips, ncol=3, byrow=FALSE)
			event_TF = eventnums == eventnum
			eventnums
			event_TF
		
			event_rownums = matrix(data=1:current_ntips, nrow=current_ntips, ncol=3, byrow=FALSE)
			event_colnums = matrix(data=1:3, nrow=current_ntips, ncol=3, byrow=TRUE)
			event_rownums
			event_colnums
		
			rownum = event_rownums[event_TF]
			colnum = event_colnums[event_TF]
		
			# Get the tip to change
			tip_to_change = rownum
		
			# Get event type
			event_type = colnum
		
			# List the edges to modify
			e <- matrix(edge[alive_TF, ], ncol = 2)
		
			edgenums_alive = (1:nrow(edge))[alive_TF]
			edgenums_alive
		
			# Which edge had the event?
			edgenum_w_event = edgenums_alive[tip_to_change]
		
			# Get the parent
			parent <- edge[edgenum_w_event, 2]
		
			#######################################################
			# Make the changes, based on the event
			#######################################################
		
			# Anagenetic range expansion/contraction event
			# (can result in extinction if ranges of size 1 drop to 0)
			if (event_type == 1)
				{
				if (printlevel >= 2)
					{
					cat("Event type: #1 (anagenetic range-change)\n", sep="")
					}
			
				starting_range = edge_ranges[edgenum_w_event,]
			
				# Get the probabilities of new ranges
				# (zeroing out the diagonal, since we know the 
				# range doesn't stay the same)
				probs_of_new_ranges = Qmat[starting_range, ]
				probs_of_new_ranges[probs_of_new_ranges < 0] = 0
				probs_of_new_ranges = probs_of_new_ranges / sum(probs_of_new_ranges)
			
				new_range = sample(x=1:numstates, size=1, replace=FALSE, prob=probs_of_new_ranges)
				edge_ranges[edgenum_w_event, ] = new_range
			
				# If the range contraction led to extinction:
				if (new_range == 1)
					{
					alive_TF[edgenum_w_event] <- FALSE
					}
			
				# Also keep track of this event
				de_event_num = de_event_num + 1
				tmprow = c(edgenum_w_event, t, dt, starting_range, new_range)
				cmdstr = paste("de", de_event_num, " = tmprow", sep="")
				eval(parse(text=cmdstr))
			
				if (printlevel >= 2)
					{
					print(tmprow)
					}
			
				# Add row to the table
				cmdstr = paste("table_of_range_change_events = rbind(table_of_range_change_events, de", de_event_num, ")", sep="")
				eval(parse(text=cmdstr))
				} # end anagenetic range-changing event
		
		
			# Speciation event
			# (requires sampling an cladogenetic range-changing event)
			if (event_type == 2)
				{
				if (printlevel >= 2)
					{
					cat("Event type: #2 (cladogenesis, including range-scenario)\n", sep="")
					}

				# Speciation event
				# take the lineage and temporarily "kill" it
				#parental_zone <- edge.zone[alive_TF][random_lineage]
				alive_TF[edgenum_w_event] <- FALSE

				# assign two new daughter nodes
				edge <- rbind(edge, c(parent, next_node), c(parent, next_node + 1))

				# move the "next node"
				next_node <- next_node + 2
			
				# Add two new living lineages
				alive_TF <- c(alive_TF, TRUE, TRUE)

				# Add the range of the parents to the daughters
				# (direct inheritance)
				#parent.range <- edge.ranges[edge_to_change,]
				# Copy the parent range to each daughter
				#edge.ranges <- rbind(edge.ranges, parent.range)
				#edge.ranges <- rbind(edge.ranges, parent.range)
			
				# Get the parent range
				parent_range = edge_ranges[edgenum_w_event,]
			
				# Sample a cladogenetic range-changing event
				daughter_states = given_a_starting_state_simulate_split(index_Qmat_0based_of_starting_state=(parent_range-1), COO_probs_columnar=COO_probs_columnar, numstates)
			
				# Convert to 1-based
				daughter_states = daughter_states + 1
			
				# Add the two daughter ranges to edge_ranges:
				edge_ranges <- rbind(edge_ranges, daughter_states[1], daughter_states[2])
	
				# add to stem depth array
				stem_depth <- c(stem_depth, t, t)
			
				# add to edge length array
				#x <- which(edge[, 2] == parent)
				#edge_length[x] = t - stem_depth[x]
				edge_length = c(edge_length, NA, NA)


				# Write the event
				if (printlevel >= 2)
					{
					event_txt = paste(parent_range, " -> ", daughter_states[1], ", ", daughter_states[2], "\n", sep="")
					cat(event_txt)
					}
			
				# Add to the cladogenesis event table
				cl_event_num = cl_event_num + 1
				tmprow = c(edgenum_w_event, parent_range, daughter_states[1], daughter_states[2])
				cmdstr = paste("cl", cl_event_num, " = tmprow", sep="")
				eval(parse(text=cmdstr))
			
				cmdstr = paste("table_of_cladogenetic_events = rbind(table_of_cladogenetic_events, cl", cl_event_num, ")", sep="")
				eval(parse(text=cmdstr))
			
				} # end cladogenesis event

			# Extinction event (lineage-wide event; extinctions can also
			# occur through range-contraction)
			if (event_type == 3)
				{
				if (printlevel >= 2)
					{
					cat("Event type: #3 (lineage extinction)\n", sep="")
					}
				
				# Get the parent edge
				starting_range = edge_ranges[edgenum_w_event,]

				# Lineage extinction event
				# take the lineage and *permanently* "kill" it
				alive_TF[edgenum_w_event] <- FALSE
			
				# The range nevertheless gets copied up
				new_range = starting_range
				edge_ranges[edgenum_w_event, ] = new_range
			
				} # end extinction event
		
			# Update the edge lengths on everything that was alive during this step
			# update to edge length array

			# List the edges to modify
			#e <- matrix(edge[alive_TF, ], ncol = 2)
			x <- which(edge[, 2] == parent)
			edge_length[x] <- t - stem_depth[x]		
	
			} #end repeat

		# break unless you need to delete extinct lineages (?)
		#if (return.all.extinct == T | sum(alive_TF) > 1) 

		# Don't break if everything died!!
		if (sum(alive_TF) == 0)
			{
			pass_txt = "\nThis simulation failed i.e. died-out.\n"
			cat(pass_txt)
			} else {

			# Don't break if it was a time-stop hit!
			if (time_stop_hit == TRUE)
				{
				pass_txt = paste("\n\nSimulation trynum#", trynum, " failed as it didn't hit taxa_stop=", taxa_stop, " within ", time_stop, " my; trying again.\n", sep="")
				cat(pass_txt)
				} else {
				# Only break if you got a success!!
				# After you get a successful simulation, 
				# add the remaining amount of time until the next event
				edge_length[alive_TF] <- t - stem_depth[alive_TF]
			
				# You got a successful simulation, so exit
				success = TRUE
				break
				}
			}

		# If failure, print the output
		num_alive = sum(alive_TF)
		d = SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"]
		e = SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"]
		j = SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"]
		
		txtnums = adf2(matrix(data=c(trynum, num_alive, t, d, e, j, SSEmodel), nrow=1))
		names(txtnums) = c("try", "num_alive", "t", "d", "e", "j", "brate", "drate", "b_exp", "d_exp")
		print(txtnums)
		
		# Brian_crazy
 		#if (trynum > 100)
 		if (trynum > 5000)
 			{
 			simtr = NA
 			# (This would cause problems later on, though, let's skip for now)
 			success = FALSE
 			break
 			}
		} # end while loop

	# raw edge matrix from simulation
	edges_from_sim = edge


	
	# Error check, i.e. if tree failed
	if (success == FALSE)
		{
		cat("\n\nSimulation of tree failed after many tries\n\n", sep="")
		# Convert the simulation into a real tree
		# Get the branch lengths
		edge_length[alive_TF] <- t - stem_depth[alive_TF]

		# Store the original edges
		original_edges = edge


		# change internal node numbers to (-1, -2, ... , -n)
		n = -1
		for (i in 1:max(edge))
			{
			if (any(edge[, 1] == i))
				{
				edge[which(edge[, 1] == i), 1] = n
				edge[which(edge[, 2] == i), 2] = n
				n = n - 1
				}
			}

		# re-number tips (everything that's left)
		edge[edge > 0] <- 1:sum(edge > 0)

		# Make an edge conversion table
		edge_conversion_table = cbind(original_edges, edge)
		edge_conversion_table = as.data.frame(edge_conversion_table)
		names(edge_conversion_table) = c("orig_ancnode", "orig_decnode", "ancnode", "decnode")
		edge_conversion_table

		# Create numbers for tip labels
		tip.label <- 1:sum(edge > 0)

		# change arrays "edge" and "tip.label" from numeric to character
		mode(edge) <- "character"
		mode(tip.label) <- "character"

		# Create the "phylo" object
		obj <- list(edge = edge, edge.length = edge_length, tip.label = tip.label)

		# call it a class phylo object
		class(obj) <- "phylo"

		# replace the old object with a valid APE object
		# This is the "raw" simulation, and the node numbers
		# are still screwed up compared to default APE
		obj2 <- old2new.phylo(obj)

		
		SSEsim_results = NULL
		SSEsim_results$states_list = SSEsim_inputs$state_indices_0based
		SSEsim_results$ranges_list = SSEsim_inputs$ranges
		SSEsim_results$areanames = SSEsim_inputs$areanames
		SSEsim_results$rootstate = rootstate
		SSEsim_results$rootnode = length(obj2$tip.label) + 1
		SSEsim_results$trynum = trynum
		SSEsim_results$simtr = obj2
		SSEsim_results$edge_conversion_table = NA
		SSEsim_results$edge_ranges = NA
		SSEsim_results$table_of_cladogenetic_events_translated = NA
		SSEsim_results$table_of_range_change_events_translated = NA
		SSEsim_results$success = FALSE
	
		return(SSEsim_results)
		} # end error check




	# Convert the table of range-changing (anagenetic) events to 
	# a data.frame
	# Error check, i.e. if no events
	if ((length(table_of_range_change_events) == 0) || (length(table_of_range_change_events) == 1))
		{
		# Fill with NAs
		table_of_range_change_events = adf2(matrix(data=NA, ncol=5))
		names(table_of_range_change_events) = c("edgenum_w_event", "t", "dt", "starting_range", "new_range")
		table_of_range_change_events_PROBLEM = TRUE
		} else {
		table_of_range_change_events = adf2(matrix(data=table_of_range_change_events, ncol=5))
		names(table_of_range_change_events) = c("edgenum_w_event", "t", "dt", "starting_range", "new_range")
		table_of_range_change_events		
		table_of_range_change_events_PROBLEM = FALSE
		}

	if ((length(table_of_cladogenetic_events) == 0) || (length(table_of_cladogenetic_events) == 1))
		{
		# Fill with NAs
		table_of_cladogenetic_events = NA
		table_of_cladogenetic_events_PROBLEM = TRUE
		} else {
		table_of_cladogenetic_events = data.frame(table_of_cladogenetic_events)
		names(table_of_cladogenetic_events) = c("edgenum_w_event", "parent_range", "Left_state", "Right_state")
		table_of_cladogenetic_events_PROBLEM = FALSE
		}









	# Convert the simulation into a real tree
	edge_ranges

	# Get the branch lengths
	edge_length[alive_TF] <- t - stem_depth[alive_TF]

	# Store the original edges
	original_edges = edge


	# change internal node numbers to (-1, -2, ... , -n)
	n = -1
	for (i in 1:max(edge))
		{
		if (any(edge[, 1] == i))
			{
			edge[which(edge[, 1] == i), 1] = n
			edge[which(edge[, 2] == i), 2] = n
			n = n - 1
			}
		}

	# re-number tips (everything that's left)
	edge[edge > 0] <- 1:sum(edge > 0)

	# Make an edge conversion table
	edge_conversion_table = cbind(original_edges, edge)
	edge_conversion_table = as.data.frame(edge_conversion_table)
	names(edge_conversion_table) = c("orig_ancnode", "orig_decnode", "ancnode", "decnode")
	edge_conversion_table

	# Create numbers for tip labels
	tip.label <- 1:sum(edge > 0)

	# change arrays "edge" and "tip.label" from numeric to character
	mode(edge) <- "character"
	mode(tip.label) <- "character"

	# Create the "phylo" object
	obj <- list(edge = edge, edge.length = edge_length, tip.label = tip.label)

	# call it a class phylo object
	class(obj) <- "phylo"

	# replace the old object with a valid APE object
	# This is the "raw" simulation, and the node numbers
	# are still screwed up compared to default APE
	obj2 <- old2new.phylo(obj)


	# Set a temporary directory for the newick file
	orig_wd = getwd()
	
	#testwd = "/drives/SkyDrive/_________thesis/_doc2/ch2_submission/2014-02-11_reviews/testsim/"
	setwd(testwd)

	# Convert the simulated tree node labels to default APE nodelabels
	# Write the tree out and read it back in
	trfn = "obj2.newick"
	write.tree(obj2, file=trfn)
	obj4 = read.tree(file=trfn)

	sim_tip_nodenums = as.numeric(obj4$tip.label)
	ape_tip_nodenums = 1:length(obj4$tip.label)
	tip_nodenums = cbind(sim_tip_nodenums, ape_tip_nodenums)
	tip_nodenums

	sim_int_nodenums = get_lagrange_nodenums(obj2)
	ape_int_nodenums = get_lagrange_nodenums(obj4)
	int_nodenums = cbind(sim_int_nodenums[,1], ape_int_nodenums[,1])
	int_nodenums
	int_nodenums = int_nodenums[order(int_nodenums[,2]), ]
	int_nodenums

	translation_table = rbind(tip_nodenums, int_nodenums)
	translation_table = adf2(translation_table)
	names(translation_table) = c("sim_nodenums", "ape_nodenums")
	translation_table



	# The edge conversion table has to convert negative tipnums
	# to obj2 nodenums, this is numtips-(negative nodenum)
	# (see old2new.phylo)
	ntips = length(obj2$tip.label)
	TF = edge_conversion_table < 0
	edge_conversion_table2 = edge_conversion_table
	edge_conversion_table2[TF] = ntips - as.numeric(edge_conversion_table[TF])

	edge_conversion_table = cbind(edge_conversion_table, edge_conversion_table2[,3:4], edge_conversion_table2[,3:4])
	names(edge_conversion_table) = c("orig_ancnode", "orig_decnode", "ancnode", "decnode", "sim_ancnode", "sim_decnode", "ape_ancnode", "ape_decnode")

	edge_conversion_table


	# Translate obj2 (sim nodenums) into obj3 (ape nodenums)
	obj3 = obj2
	obj4_tiplabels = as.numeric(obj4$tip.label)
	
	# Change the node numbers to match default APE node numbers
	for (i in 1:nrow(translation_table))
		{
		# Translate nodes in edges
		TF = obj2$edge == translation_table$sim_nodenums[i]
		obj3$edge[TF] = translation_table$ape_nodenums[i]

		# Translate tipnames
		TF = obj4_tiplabels == translation_table$sim_nodenums[i]
		TF
		obj3$tip.label[TF] = translation_table$ape_nodenums[i]
		
		# Translate edge_conversion_table
		TF = edge_conversion_table$sim_ancnode == translation_table$sim_nodenums[i]
		edge_conversion_table$ape_ancnode[TF] = translation_table$ape_nodenums[i]

		TF = edge_conversion_table$sim_decnode == translation_table$sim_nodenums[i]
		edge_conversion_table$ape_decnode[TF] = translation_table$ape_nodenums[i]
		
		# Translate raw sim node numbers into 
		#TF = edge_conversion_table$sim_ancnode == translation_table$sim_nodenums[i]
		#raw_simnode = unique(edge_conversion_table$orig_ancnode[TF])
		#simnode = unique(edge_conversion_table$sim_ancnode[TF])
		}



	# Plot the raw simtree, and the APE simtree
	if (printlevel >= 4)
		{
		plot(obj2, label.offset=0.2)
		axisPhylo()
		title("Simulation tree: raw node labels")
		tiplabels()
		nodelabels()

		plot(obj3, label.offset=0.2)
		axisPhylo()
		title("Simulation tree: APE node labels")
		tiplabels()
		nodelabels()

		cbind(obj2$edge, obj3$edge)
		}

	# Relabel the tips in obj3 so it says "sp1", "sp2", etc.
	obj3$tip.label = paste("sp", obj3$tip.label, sep="")
	

	
	
	# Translate the simulation edgenums_w_events into APE anc and decscendent nodes
	ape_ancnode = obj3$edge[table_of_cladogenetic_events$edgenum_w_event,1]
	ape_decnode = obj3$edge[table_of_cladogenetic_events$edgenum_w_event,2]
	table_of_cladogenetic_events_translated = cbind(ape_ancnode, ape_decnode, table_of_cladogenetic_events)
	names(table_of_cladogenetic_events_translated)[3] = "sim_edgenum_w_event"
	table_of_cladogenetic_events_translated

	if (table_of_range_change_events_PROBLEM == FALSE)
		{
		table_of_range_change_events_translated = table_of_range_change_events
		ape_ancnode = obj3$edge[table_of_range_change_events$edgenum_w_event,1]
		ape_decnode = obj3$edge[table_of_range_change_events$edgenum_w_event,2]
		table_of_range_change_events_translated = cbind(ape_ancnode, ape_decnode, table_of_range_change_events)
		names(table_of_range_change_events_translated)[3] = "sim_edgenum_w_event"
		table_of_range_change_events_translated
		} else {
		table_of_range_change_events_translated = table_of_range_change_events
		}

	
	#######################################################
	# Store the simulation results, tables, etc.
	#######################################################
	SSEsim_results = NULL
	
	# We will need the states_list and ranges_list later
	SSEsim_results$states_list = SSEsim_inputs$state_indices_0based
	SSEsim_results$ranges_list = SSEsim_inputs$ranges
	SSEsim_results$areanames = SSEsim_inputs$areanames

	SSEsim_results$rootstate = rootstate
	SSEsim_results$rootnode = length(obj3$tip.label) + 1
	
	# Number of tries to get a successful simulation
	SSEsim_results$trynum = trynum

	# Simulated tree with APE node numbers
	# (should save and read the same, thankfully)
	SSEsim_results$simtr = obj3

	# The ranges at the ends of the branches/edges, from the 
	# original simulation.
	# The edges remain in the original order
	# See the edge conversion table; the nodes holding these edge
	# ranges are edge_conversion_table$sim_decnode
	SSEsim_results$edge_conversion_table = edge_conversion_table
	SSEsim_results$edge_ranges = edge_ranges
	
	SSEsim_results$table_of_cladogenetic_events_translated = table_of_cladogenetic_events_translated
	SSEsim_results$table_of_range_change_events_translated = table_of_range_change_events_translated
	
	
	SSEsim_results$success = success
	
	# Double-check the output
	# cbind(obj3$edge, obj2$edge, obj3$edge.length, obj2$edge.length, edge_conversion_table, edge_ranges)
	
	# Return to the original working directory
	setwd(orig_wd)
	
	return(SSEsim_results)
	}






# Take the SSEsim results and convert them to:
# 1. Newick and geography files WITH extinct tips
# 2. Newick and geography files WITHOUT extinct tips
# 3. True state probabilities for each node in complete tree 
# 4. True state probabilities for each node in observed tree
# 5. DE event counts
# 6. Cladogenetic event counts
# 7. Cladogenetic event counts preserved in observed tree
SSEsim_to_files <- function(SSEsim_results, simdir, fossils_older_than=0.001, printlevel=2, write_files=TRUE)
	{
	defaults='
	simdir = "/drives/SkyDrive/_________thesis/_doc2/ch2_submission/2014-02-11_reviews/testsim/"
	fossils_older_than=0.001
	printlevel=4
	write_files=TRUE
	'
	
	# Working directories
	orig_wd = getwd()
	setwd(simdir)
	
	simtr = SSEsim_results$simtr
	states_list = SSEsim_results$states_list
	ranges_list = SSEsim_results$ranges_list
	areanames = SSEsim_results$areanames
	numstates = length(states_list)

	rootstate = SSEsim_results$rootstate
	rootnode = SSEsim_results$rootnode
	
	# Label the original nodes of the full tree, so we can link
	# them to the observed tree
	ntips = length(simtr$tip.label)
	tipnums = 1:ntips
	nodenums = (ntips+1):(ntips+simtr$Nnode)
	nodenums
	all_nodenums = c(tipnums, nodenums)
	numnodes = length(all_nodenums)
	
	simtr$node.label = paste("fulltr_node", nodenums, sep="")
	
	# Get the tips to drop
	simtr_table = prt(simtr, fossils_older_than=fossils_older_than, printflag=FALSE)
	tips_to_drop_TF = simtr_table$fossils[1:ntips]
	num_fossil_tips = sum(tips_to_drop_TF)

	if (num_fossil_tips > 0)
		{
		tips_to_drop = simtr$tip.label[tips_to_drop_TF]
		simtr_observed = ape::drop.tip(phy=simtr, tip=tips_to_drop, trim.internal=TRUE)
		} else {
		# No tips to drop
		simtr_observed = simtr
		}
	
	# Plot the tree
	if (printlevel >= 4)
		{
		plot(simtr_observed, label.offset=0.2)
		axisPhylo()
		title("Simulated tree (observed)")
		nodelabels()
		tiplabels()
		}

	# Get the nodenumbers of the original full tree, which correspond to the
	# ordered node numbers of the subset tree
	list_of_full_nodelabels_in_subset_tree = c(simtr_observed$tip.label, simtr_observed$node.label)
	list_of_full_nodenums_in_subset_tree = gsub(pattern="sp", replacement="", x=list_of_full_nodelabels_in_subset_tree)
	list_of_full_nodenums_in_subset_tree = gsub(pattern="fulltr_node", replacement="", x=list_of_full_nodenums_in_subset_tree)
	list_of_full_nodenums_in_subset_tree = as.numeric(list_of_full_nodenums_in_subset_tree)
	list_of_full_nodenums_in_subset_tree
	
	list_of_full_tipnums_in_subset_tree = list_of_full_nodenums_in_subset_tree[1:length(simtr_observed$tip.label)]
	list_of_full_tipnums_in_subset_tree
	
	# Label the cladogenetic events
	table_of_cladogenetic_events_translated = SSEsim_results$table_of_cladogenetic_events_translated
	cladogenetic_event_labels = label_table_of_cladogenetic_events(table_of_cladogenetic_events_translated, states_list, ranges_list)
	table_of_cladogenetic_events_translated = cbind(table_of_cladogenetic_events_translated, cladogenetic_event_labels)
	table_of_cladogenetic_events_translated
	
	# Label the anagenetic events
	table_of_range_change_events_translated = SSEsim_results$table_of_range_change_events_translated
	
	if (!is.na(table_of_range_change_events_translated[1]))
		{
		anagenetic_event_labels = label_table_of_anagenetic_events(table_of_range_change_events_translated, states_list, ranges_list)
	
		table_of_range_change_events_translated = cbind(table_of_range_change_events_translated, anagenetic_event_labels)
		}

	#######################################################
	# Make tipranges object for complete, and observed trees
	#######################################################
	edge_ranges = SSEsim_results$edge_ranges
	edges_table = cbind(SSEsim_results$edge_conversion_table, edge_ranges)
	edges_table2 = edges_table[order(edges_table$ape_decnode), ]
	edges_table2
	
	tiplabels_for_tipranges = simtr$tip.label
	tipranges_simulated = edges_table_to_tipranges_object(edges_table2, states_list, ntips, tiplabels_for_tipranges, areanames, addval=0)	

	# Get the (TRUE!) probabilities of the ancestral states
	# on the complete tree
	ancstate_probs_simfull = matrix(data=0, nrow=numnodes, ncol=numstates)
	for (i in 1:numnodes)
		{
		tmpstate_index = edges_table2$edge_ranges[i]
		tmpnodenum_index = edges_table2$ape_decnode[i]
		ancstate_probs_simfull[tmpnodenum_index, tmpstate_index] = 1
		}
	# Put in the root state
	ancstate_probs_simfull[rootnode, rootstate] = 1
	ancstate_probs_simfull
	rowSums(ancstate_probs_simfull)
	

	
	# Reduce the event lists based on what went extinct
	internal_nodes_that_are_observed_TF = simtr$node.label %in% simtr_observed$node.label
	internal_nodes_that_are_observed_TF
	internal_nodes_that_went_extinct_TF = internal_nodes_that_are_observed_TF == FALSE
	if (sum(internal_nodes_that_went_extinct_TF) > 0)
		{
		internal_nodes_that_went_extinct_labels = simtr$node.label[internal_nodes_that_went_extinct_TF]
		internal_nodes_that_went_extinct_labels
		internal_nodenums_that_went_extinct = as.numeric(gsub(pattern="fulltr_node", replacement="", x=internal_nodes_that_went_extinct_labels))
		
		tipnums_that_went_extinct = as.numeric(gsub(pattern="sp", replacement="", x=tips_to_drop))
		
		internal_nodenums_that_went_extinct
		tipnums_that_went_extinct
		
		nodenums_extinct = sort(c(tipnums_that_went_extinct, internal_nodenums_that_went_extinct))
		
		# Now, subset the de_events table and the cladogenetic events table
		
		# cladogenetic
		TF = table_of_cladogenetic_events_translated$ape_decnode %in% nodenums_extinct == FALSE
		table_of_cladogenetic_events_observed = table_of_cladogenetic_events_translated[TF, ]
		table_of_cladogenetic_events_observed
		
		# anagenetic
		TF = table_of_range_change_events_translated$ape_decnode %in% nodenums_extinct == FALSE
		table_of_range_change_events_observed = table_of_range_change_events_translated[TF, ]
		table_of_range_change_events_observed

		# subset ancestral probabilities
		rownums = 1:nrow(ancstate_probs_simfull)
		keepTF = (rownums %in% nodenums_extinct) == FALSE
		ancstate_probs_simobs = ancstate_probs_simfull[list_of_full_nodenums_in_subset_tree, ]
		
		# subset tipranges
		tipranges_observed = tipranges_simulated
		tipranges_observed@df = tipranges_simulated@df[list_of_full_tipnums_in_subset_tree, ]
		tipranges_observed
	
		
		
		} else {
		# Nothing went extinct, just use these
		table_of_cladogenetic_events_observed = table_of_cladogenetic_events_translated
		table_of_range_change_events_observed = table_of_range_change_events_translated
		ancstate_probs_simobs = ancstate_probs_simfull
		tipranges_observed = tipranges_simulated
		}


	
	#######################################################
	# Get some summary statistics
	#######################################################
	simstats = NULL
	
	num_fossil_tips = sum(tips_to_drop_TF)
	num_extinct_nodes = simtr$Nnode - simtr_observed$Nnode

	y_actual = sum(table_of_cladogenetic_events_translated$event_type == "sympatry (y)")
	s_actual = sum(table_of_cladogenetic_events_translated$event_type == "subset (s)")
	v_actual = sum(table_of_cladogenetic_events_translated$event_type == "vicariance (v)")
	j_actual = sum(table_of_cladogenetic_events_translated$event_type == "founder (j)")

	y_observed = sum(table_of_cladogenetic_events_observed$event_type == "sympatry (y)")
	s_observed = sum(table_of_cladogenetic_events_observed$event_type == "subset (s)")
	v_observed = sum(table_of_cladogenetic_events_observed$event_type == "vicariance (v)")
	j_observed = sum(table_of_cladogenetic_events_observed$event_type == "founder (j)")

	d_actual = sum(table_of_range_change_events_translated$event_type == "expansion (d)")
	e_actual = sum(table_of_range_change_events_translated$event_type == "contraction (e)")
	a_actual = sum(table_of_range_change_events_translated$event_type == "range-switching (a)")

	d_observed = sum(table_of_range_change_events_observed$event_type == "expansion (d)")
	e_observed = sum(table_of_range_change_events_observed$event_type == "contraction (e)")
	a_observed = sum(table_of_range_change_events_observed$event_type == "range-switching (a)")
	
	simstats = c(num_fossil_tips, num_extinct_nodes, y_actual, s_actual, v_actual, j_actual, y_observed, s_observed, v_observed, j_observed, d_actual, e_actual, a_actual, d_observed, e_observed, a_observed)

	simstats = adf2(matrix(data=simstats, nrow=1))
	names(simstats) = c("num_fossil_tips", "num_extinct_nodes", "y_actual", "s_actual", "v_actual", "j_actual", "y_observed", "s_observed", "v_observed", "j_observed", "d_actual", "e_actual", "a_actual", "d_observed", "e_observed", "a_observed")
	
	
	
	#######################################################
	# Output object
	#######################################################
	SSEsim_results_processed = NULL
	
	SSEsim_results_processed$simtr = simtr
	SSEsim_results_processed$simtr_observed = simtr_observed

	SSEsim_results_processed$ancstate_probs_simfull = ancstate_probs_simfull
	SSEsim_results_processed$ancstate_probs_simobs = ancstate_probs_simobs

	SSEsim_results_processed$tipranges_simulated = tipranges_simulated
	SSEsim_results_processed$tipranges_observed = tipranges_observed

	SSEsim_results_processed$table_of_cladogenetic_events_translated = table_of_cladogenetic_events_translated
	SSEsim_results_processed$table_of_cladogenetic_events_observed = table_of_cladogenetic_events_observed

	SSEsim_results_processed$table_of_range_change_events_translated = table_of_range_change_events_translated
	SSEsim_results_processed$table_of_range_change_events_observed = table_of_range_change_events_observed

	SSEsim_results_processed$simstats = simstats
	SSEsim_results_processed$simdir = simdir

	if (write_files == TRUE)
		{
		write.tree(phy=simtr, file="simtr_complete.newick")
		write.tree(phy=simtr_observed, file="simtr_observed.newick")
		write.tree(phy=simtr_observed, file="tree.newick")
		
		save_tipranges_to_LagrangePHYLIP(tipranges=tipranges_simulated, lgdata_fn="geog_sim_complete.txt")
		save_tipranges_to_LagrangePHYLIP(tipranges=tipranges_observed, lgdata_fn="geog_sim_observed.txt")
		save_tipranges_to_LagrangePHYLIP(tipranges=tipranges_observed, lgdata_fn="geog.data")
		}
	
	# Return to original working directory
	setwd(orig_wd)

	return(SSEsim_results_processed)
	} # end SSEsim_to_files()



# Label the recorded anagenetic events as d, e...
label_table_of_anagenetic_events <- function(table_of_range_change_events_translated, states_list, ranges_list)
	{
	starting_ranges = table_of_range_change_events_translated$starting_range
	new_ranges = table_of_range_change_events_translated$new_range

	numevents = nrow(table_of_range_change_events_translated)
	event_type = rep(NA, numevents)
	event_txt = rep(NA, numevents)

	for (i in 1:numevents)
		{
		# Text describing the event
		event_txt[i] = paste(ranges_list[[starting_ranges[i]]], "->", ranges_list[[new_ranges[i]]], sep="")
		
		# "dispersal" (range expansion)
		if (length(states_list[[starting_ranges[i]]]) < length(states_list[[new_ranges[i]]]))
			{
			event_type[i] = "expansion (d)"
			next()
			}
		
		# "extinction" (range contraction / local extirpation)
		if (length(states_list[[starting_ranges[i]]]) > length(states_list[[new_ranges[i]]]))
			{
			event_type[i] = "contraction (e)"
			next()
			}

		if (length(states_list[[starting_ranges[i]]]) == length(states_list[[new_ranges[i]]]))
			{
			# Contraction will result in the same length, check for this
			if (ranges_list[[new_ranges[i]]] == "_")
				{
				event_type[i] = "contraction (e)"
				next()				
				} else {
				event_type[i] = "range-switching (a)"
				next()
				}
			
			}
		} # end for-loop
	
	anagenetic_event_labels = as.data.frame(cbind(event_type, event_txt))
	anagenetic_event_labels
	
	return(anagenetic_event_labels)
	} # end label_table_of_anagenetic_events()


# Label the recorded cladogenetic events as j, v, s, y...
label_table_of_cladogenetic_events <- function(table_of_cladogenetic_events_translated, states_list, ranges_list)
	{
	parent_ranges = table_of_cladogenetic_events_translated$parent_range
	Left_states = table_of_cladogenetic_events_translated$Left_state
	Right_states = table_of_cladogenetic_events_translated$Right_state

	numevents = nrow(table_of_cladogenetic_events_translated)
	event_type = rep(NA, numevents)
	event_txt = rep(NA, numevents)
	
	for (i in 1:numevents)
		{
		# Text describing the event
		event_txt[i] = paste(ranges_list[[parent_ranges[i]]], "->", ranges_list[[Left_states[i]]], ",",  ranges_list[[Right_states[i]]], sep="")
		
		# Sympatry
		if (Left_states[i] == Right_states[i])
			{
			event_type[i] = "sympatry (y)"
			next()
			}
		
		parent_indices = sort(states_list[[parent_ranges[i]]])
		Left_indices = states_list[[Left_states[i]]]
		Right_indices = states_list[[Right_states[i]]]
		desc_merged = sort(c(Left_indices, Right_indices))
		desc_indices = sort(unique(c(Left_indices, Right_indices)))
		
		# Founder-events
		if (length(desc_indices) > length(parent_indices))
			{
			event_type[i] = "founder (j)"
			next()			
			}
		
		# Vicariance
		if ((length(parent_indices) == length(desc_merged)) && (all(parent_indices == desc_merged) == TRUE))
			{
			event_type[i] = "vicariance (v)"
			next()						
			} else {
			event_type[i] = "subset (s)"
			next()			
			}
		
		} # end for-loop
		
		
	cladogenetic_event_labels = as.data.frame(cbind(event_type, event_txt))
	cladogenetic_event_labels
	
	return(cladogenetic_event_labels)
	} # end label_table_of_cladogenetic_events()



edges_table_to_tipranges_object <- function(edges_table2, states_list, ntips, tiplabels_for_tipranges, areanames, addval=0)
	{
	defaults = '
	tiplabels = paste("sp", 1:ntips, sep="")
	areanames = c("K", "O", "M", "H")
	addval=0
	'
	binary_table = matrix(data=0, nrow=ntips, ncol=length(areanames))

	for (i in 1:ntips)
		{
		# Get the index of the range, from the simulation (APE nodenums)
		tmp_range_index = edges_table2$edge_ranges[i]
		
		# Get the 1-based indexes of the areas
		area_indices_1based = states_list[[tmp_range_index]] + 1 + addval
		
		# Convert to binary array
		binary_table[i, area_indices_1based] = 1
		}

	binary_table

	# Convert to tipranges object
	tmpdf2 = adf2(data.matrix(binary_table))
	names(tmpdf2) = areanames
	rownames(tmpdf2) = tiplabels_for_tipranges
	tmpdf2
	
	tipranges_simulated = define_tipranges_object(tmpdf=tmpdf2)
	tipranges_simulated
	
	return(tipranges_simulated)
	}





plot_simDEC_DECJ_inferences <- function(simdir=NULL, sim_params_Rdata_table="/simdata/BGB/dej_params.Rdata", modelnull="DEC", modelalt="DEC+J", makePDFs=TRUE, openPDFs=TRUE)
	{
	defaults='
	simdir = "/simdata/BGB/ps0011_sim003/"
	sim_params_Rdata_table="/simdata/BGB/dej_params.Rdata"
	modelnull="DEC"
	modelalt="DEC+J"
	makePDFs=TRUE
	openPDFs=TRUE
	'
	
	# If you just want the stats, not the plots
	if (makePDFs == FALSE)
		{
		juststats = TRUE
		openPDFs = FALSE
		} else {
		juststats = FALSE
		}
	
	
	#######################################################
	# Get whatever directory your simulations are in
	#######################################################
	if (is.null(simdir))
		{
		cat("\n\nWARNING: no simulation directory specified, using default on Nick Matzke's computer.\n\n")
		simdir = "/simdata/BGB/ps0011_sim003/"
		}
	

	# Set the working directory to the simulation directory
	setwd(simdir)
	
	# Process the directory to get the row of the simulation parameters
	words = strsplit(simdir, "/")[[1]]
	dir_txt = words[length(words)]

	# Get the row of parameter values
	words = strsplit(dir_txt, "_")[[1]]
	word = words[1]
	param_iter = as.numeric(gsub(pattern="ps", replacement="", x=word))
	param_iter


	# Get the parameters for this run
	load(file="/simdata/BGB/dej_params.Rdata")
	nums = 1:nrow(dej_params)
	dej_params = cbind(nums, dej_params)
	dim(dej_params)
	param_vals = dej_params[param_iter,]
	param_vals




	#######################################################
	# Read saved results
	#######################################################
	# Read the results
	resDEC_fn = "DEC_inf.Rdata"
	load(resDEC_fn)
	resDEC = res

	resDECJ_fn = "DECJ_inf.Rdata"
	load(resDECJ_fn)
	resDECj = res


	SSEsim_results_fn = "SSEsim_results_processed.Rdata"
	load(SSEsim_results_fn)
	SSEsim_results_processed




	# Inputs for graphics and stats
	tr = read.tree(resDEC$inputs$trfn)
	tipranges = getranges_from_LagrangePHYLIP(resDEC$inputs$geogfn)
	BioGeoBEARS_run_object = resDEC$inputs



	# Results table
	restable = NULL
	teststable = NULL

	
	if (makePDFs==TRUE)
		{
		pdffn = paste(dir_txt, "_DEC_vs_DECj_SSE.pdf", sep="")
		pdf(pdffn, width=8.5, height=11)
		} # end makePDFs
	
	#######################################################
	# Plot ancestral states - DEC
	#######################################################
	analysis_titletxt = paste("BioGeoBEARS ", modelnull, " on ", dir_txt, " d=", param_vals$d, " e=", param_vals$e, " j=", param_vals$j, " brate=", param_vals$brate, " drate=", param_vals$drate, " bexp=", param_vals$b_exp, " dexp=", param_vals$d_exp, sep="")

	# Setup
	results_object = resDEC
	scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

	# States
	res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges, juststats=juststats)
	row.names(res1) = modelnull
	
	# Pie chart
	plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges, juststats=juststats)

	#######################################################
	# Plot ancestral states - DECJ
	#######################################################
	analysis_titletxt = paste("BioGeoBEARS ", modelalt, " on ", dir_txt, " d=", param_vals$d, " e=", param_vals$e, " j=", param_vals$j, " brate=", param_vals$brate, " drate=", param_vals$drate, " bexp=", param_vals$b_exp, " dexp=", param_vals$d_exp, sep="")

	# Setup
	results_object = resDECj
	scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

	# States
	res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges, juststats=juststats)
	row.names(res2) = modelalt
	
	# Pie chart
	plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges, juststats=juststats)

	if (makePDFs==TRUE)
		{
		dev.off()
		if (openPDFs==TRUE)
			{
			cmdstr = paste("open ", pdffn, sep="")
			system(cmdstr)
			}
		}


	#######################################################
	# Stats
	#######################################################
	# We have to extract the log-likelihood differently, depending on the 
	# version of optim/optimx
	if (BioGeoBEARS_run_object$use_optimx == TRUE)
		{
		# Using optimx() results
		if (packageVersion("optimx") < 2013)
			{
			# optimx 2012
			LnL_2 = as.numeric(resDEC$optim_result$fvalues)
			LnL_1 = as.numeric(resDECj$optim_result$fvalues)
			} else {
			# optimx 2013
			LnL_2 = as.numeric(resDEC$optim_result$value)
			LnL_1 = as.numeric(resDECj$optim_result$value)
			} # end optimx 2012 vs. 2013
		} else {
		# Using optim() results
		LnL_2 = as.numeric(resDEC$optim_result$value)
		LnL_1 = as.numeric(resDECj$optim_result$value)
		} # end optim vs. optimx

	numparams1 = 3
	numparams2 = 2
	stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
	row.names(stats) = paste(modelnull, "_v_", modelalt, sep="")
	stats

	res1
	res2

	rbind(res1, res2)
	tmp_tests = conditional_format_table(stats)

	restable = rbind(restable, res1, res2)
	teststable = rbind(teststable, tmp_tests)


	restable
	teststable

	cat("\nParameter values for this simulation:\n")
	print(param_vals)
	
	
	# Assemble output list
	simDEC_DECJ_inf_results = NULL	
	simDEC_DECJ_inf_results$SSEsim_results_processed = SSEsim_results_processed
	simDEC_DECJ_inf_results$res1 = res1
	simDEC_DECJ_inf_results$res2 = res2
	simDEC_DECJ_inf_results$tests = tmp_tests
	
	return(simDEC_DECJ_inf_results)
	}




hist_sim_v_inf_vals <- function(inferred_vals, trueval)
	{
	h1 = hist(inferred_vals, plot=FALSE)
	xticks = c(0, pretty(c(0, h1$mids), n=1))
	yticks = c(0, pretty(c(0,h1$counts), n=1))
	xlims = c(0, max(xticks))
	ylims = c(0, max(yticks))
	plot(h1$mids, h1$counts, pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n", xlim=xlims, ylim=ylims)
	plot(h1, add=TRUE)
	
	abline(v=trueval, lwd=2, col="blue", lty="dashed")
	box()
	# x-axis tick marks
	axis(side=1, at=xticks, labels=TRUE, tick=TRUE)
	# y-axis tick marks (not really needed for histograms)
	#axis(side=1, at=xticks, labels=TRUE, tick=TRUE)
	
	return(h1)
	}
















#######################################################
# Make plots of DEC vs. DEC+J parameter inference and model choice
#######################################################



plot_inf_accuracy_vs_SSEsims_v2 <- function(plot_inputs, xwidth_manual=NULL)
	{
	

defaults='
titletxt1 = bquote(paste("100 SSE simulations (", lambda, "=0.3, ", mu, "=0.3, ", alpha, "=1, ", omega, "=-1)", sep=""))
titletxt2 = paste("DEC (white) vs. DEC+J (grey) inference", sep="")
pdffn = "DEC_DECJ_inf_boxplots_SSE_b03d03bx1dx-1.pdf"

simtype = "SSE"
segwidth_truth = 0.05
plot_inputs$segwidth_truth = segwidth_truth

simtype = "SSE"
brate_plot = 0.3
drate_plot = 0.3
b_exp_plot = 1
d_exp_plot = -1

plot_inputs = NULL
plot_inputs$titletxt1 = NULL
plot_inputs$pdffn = pdffn
plot_inputs$doPDF = TRUE
plot_inputs$siminf_stats_good = siminf_stats_good
plot_inputs$simtype = simtype
plot_inputs$brate_plot = brate_plot
plot_inputs$drate_plot = drate_plot
plot_inputs$b_exp_plot = b_exp_plot
plot_inputs$d_exp_plot = d_exp_plot
plot_inputs$uniq_params_txt = uniq_params_txt
plot_inputs$segwidth_truth = segwidth_truth

'# end defaults

titletxt1 = plot_inputs$titletxt1
pdffn = plot_inputs$pdffn
doPDF = plot_inputs$doPDF
siminf_stats_good = plot_inputs$siminf_stats_good
brate_plot = plot_inputs$brate_plot
drate_plot = plot_inputs$drate_plot
b_exp_plot = plot_inputs$b_exp_plot
d_exp_plot = plot_inputs$d_exp_plot
uniq_params_txt = plot_inputs$uniq_params_txt
segwidth_truth = plot_inputs$segwidth_truth
simtype = plot_inputs$simtype

if (is.null(titletxt1))
	{
	titletxt1 = bquote(paste(.(simtype), " simulations (", lambda, "=", .(brate_plot), ", ",  mu, "=", .(drate_plot), ", ", alpha, "=", .(b_exp_plot), ", ", omega, "=", .(d_exp_plot), ")", sep=""))
	}



#######################################################
# Spacer between models
#######################################################
xspace_value = 0.5

if (is.null(xwidth_manual))
	{
	xwidth = uniq_params_num + (uniq_params_num-1) * xspace_value
	} else {
	xwidth = xwidth_manual
	}

#######################################################
# Make plots of DEC vs. DEC+J parameter inference and model choice
#######################################################
if(plot_inputs$doPDF == TRUE)
	{
	pdf(pdffn, height=11, width=8.5)
	}

par(oma = c(5,4,4,1))
par(mar=c(0,2,1,1))
nums = 1:5
mat = matrix(data=nums, nrow=5, ncol=5, byrow=FALSE)
mat

layout(mat)

i=1
ymax = 0.25
yaxis_ticks = c(0, 0.05, 0.1, 0.15, 0.2, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote(italic("d")), side=2, las=1, line=4)

# Plot the page title
mtext(text=titletxt1, side=3, line=0, cex=1.2, outer=TRUE)
#mtext(text=titletxt2, side=3, line=-0.5, cex=0.8, outer=TRUE)

# Plot a legend
xpos = 1
points(x=xpos, y=0.925*ymax, pch=18, cex=2, col="black")
segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=0.925*ymax, y1=0.925*ymax, lwd=2, col="black")
text(x=xpos+1.5*segwidth_truth, y=0.925*ymax, labels="true parameter value", pos=4, offset=1)
points(x=xpos, y=0.825*ymax, pch=22, cex=3, col="black", bg="white")
text(x=xpos+1.5*segwidth_truth, y=0.825*ymax, labels="DEC inference", pos=4, offset=1)
points(x=xpos, y=0.725*ymax, pch=22, cex=3, col="black", bg="lightgray")
text(x=xpos+1.5*segwidth_truth, y=0.725*ymax, labels="DEC+J inference", pos=4, offset=1)



# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	siminf_stats_sim = siminf_stats_subset[TF,]
	siminf_stats_sim
	
	
	dtrue = unique(siminf_stats_sim$d)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=siminf_stats_sim$d_DECinf, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=siminf_stats_sim$d_DECJinf, at=xpos, add=TRUE, col="lightgray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	points(x=xpos, y=dtrue, pch=18, cex=2, col="black")
	segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=dtrue, y1=dtrue, lwd=2, col="black")
	
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through d truth/inference plots



i=1
ymax = 0.25
yaxis_ticks = c(0, 0.05, 0.1, 0.15, 0.2, ymax)

plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote(italic("e")), side=2, las=1, line=4)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	siminf_stats_sim = siminf_stats_subset[TF,]
	siminf_stats_sim
	
	
	etrue = unique(siminf_stats_sim$e)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=siminf_stats_sim$e_DECinf, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=siminf_stats_sim$e_DECJinf, at=xpos, add=TRUE, col="lightgray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	points(x=xpos, y=etrue, pch=18, cex=2, col="black")
	segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=etrue, y1=etrue, lwd=2, col="black")
	
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through e truth/inference plots



i=1
ymax = 0.5
yaxis_ticks = pretty(c(0, ymax))
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote(italic("j")), side=2, las=1, line=4)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	siminf_stats_sim = siminf_stats_subset[TF,]
	siminf_stats_sim
	
	
	jtrue = unique(siminf_stats_sim$j)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=siminf_stats_sim$j_DECinf, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=siminf_stats_sim$j_DECJinf, at=xpos, add=TRUE, col="lightgray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	points(x=xpos, y=jtrue, pch=18, cex=2, col="black")
	segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=jtrue, y1=jtrue, lwd=2, col="black")

	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through j truth/inference plots





#######################################################
# Accuracy
#######################################################
i=1
ymax = 1
yaxis_ticks =c(0, 1/15, 0.25, 0.5, 0.75, ymax)
yaxis_ticks_txt =c("0", "1/15", "0.25", "0.5", "0.75", "1")

plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks_txt, las=1)
mtext(text="Proportion correct", side=2, las=3, line=3.5, cex=1)

# Line representing random guess among states
abline(h=1/15, lty="dotted", col="black", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	
	# Calculate accuracy
	MLstate_TF_ancstate_accuracies_DEC_subset = ancstate_accuracies_DEC_subset[TF,99:147]
	MLstate_accuracy_DEC = apply(X=MLstate_TF_ancstate_accuracies_DEC_subset, MARGIN=1, FUN=mean, na.rm=TRUE)

	MLstate_TF_ancstate_accuracies_DECJ_subset = ancstate_accuracies_DECJ_subset[TF,99:147]
	MLstate_accuracy_DECJ = apply(X=MLstate_TF_ancstate_accuracies_DECJ_subset, MARGIN=1, FUN=mean, na.rm=TRUE)
	
	
	jtrue = unique(siminf_stats_sim$j)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=MLstate_accuracy_DEC, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=MLstate_accuracy_DECJ, at=xpos, add=TRUE, col="lightgray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	#points(x=i, y=jtrue, pch="*", cex=4, col="darkblue")

	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through accuracy boxplots








#######################################################
# Model choice -- with AICc
#######################################################
i=1
ymax = 130
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 50, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote(paste(Delta, "AICc", sep=""))
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3.5, cex=1)

# Dashed red line at the significance level
abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	siminf_stats_sim = siminf_stats_subset[TF,]
	siminf_stats_sim
	
	LnL_DEC = siminf_stats_sim$LnL_DEC
	LnL_DECJ = siminf_stats_sim$LnL_DECJ
	
	# Boxplot of the log-likelihood advantage of DEC+J
	LnL_advantage = LnL_DECJ-LnL_DEC
	# Min on this plot is 0.1
	# Also a few optimization issues at high d, e values
	# (could be datasets with no signal left; functions will misfire if LnL DEC>LnL DEC+J)
	TF = LnL_advantage < 0
	LnL_DECJ[TF] = LnL_DEC[TF]
	
	if (sum(TF) > 0)
		{
		cat("Note: In ", sum(TF), " cases, LnL DEC > LnL DEC+J.\n", sep="")
		# This may be very minor, or may indicate a failure of optimx ML routine
		# to find the true optimum. For now, in these cases we are setting LnL 
		# DEC+J to equal LnL DEC.
		}
	
	# Calculate DeltaAICc 
	AICc_DEC = getAICc(LnL=LnL_DEC, numparams=2, samplesize=100)
	AICc_DECJ = getAICc(LnL=LnL_DECJ, numparams=3, samplesize=100)
	delta_AICc_vals = AICc_DEC - AICc_DECJ

	# Plot AICc weights
	# DEC or DEC+J
	xpos = extra_xspace + i + 0.0
	
	# Plot one box (grey or white)
	# AICc weights look stupid
	jtrue = unique(siminf_stats_sim$j)
	if (jtrue > 0)
		{
		# width makes fatter boxes
		b1 = boxplot(x=delta_AICc_vals, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n", col="lightgray", border="black", width=3)
		} else {
		# width makes fatter boxes
		b1 = boxplot(x=delta_AICc_vals, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n", col="lightgray", border="black", width=3)
		}
	
	# Add "true model:" text
	if (i==1)
		{
		txtplot = "True model:"
		mtext(side=1, at=-1, text=txtplot, cex=0.55, las=1, line=0)
		}
	
	if (jtrue > 0)
		{
		txtplot = "+J"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		} else {
		txtplot = "DEC"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		}
	
		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot


if(plot_inputs$doPDF == TRUE)
	{
	dev.off()
	cmdstr = paste("open ", pdffn, sep="")
	system(cmdstr)
	}



	
	} # end plot_inf_accuracy_vs_SSEsims_v2()
























#######################################################
# Do a subset plot for publication
#######################################################

plot_inf_accuracy_vs_SSEsims_v3<- function(plot_inputs, xwidth_manual=NULL)
	{
	

defaults='
titletxt1 = bquote(paste("100 SSE simulations (", lambda, "=0.3, ", mu, "=0.3, ", alpha, "=1, ", omega, "=-1)", sep=""))
titletxt2 = paste("DEC (white) vs. DEC+J (grey) inference", sep="")
pdffn = "DEC_DECJ_inf_boxplots_SSE_b03d03bx1dx-1.pdf"

simtype = "SSE"
segwidth_truth = 0.05
plot_inputs$segwidth_truth = segwidth_truth

simtype = "SSE"
brate_plot = 0.3
drate_plot = 0.3
b_exp_plot = 1
d_exp_plot = -1

plot_inputs = NULL
plot_inputs$titletxt1 = NULL
plot_inputs$pdffn = pdffn
plot_inputs$doPDF = TRUE
plot_inputs$siminf_stats_good = siminf_stats_good
plot_inputs$simtype = simtype
plot_inputs$brate_plot = brate_plot
plot_inputs$drate_plot = drate_plot
plot_inputs$b_exp_plot = b_exp_plot
plot_inputs$d_exp_plot = d_exp_plot
plot_inputs$uniq_params_txt = uniq_params_txt
plot_inputs$segwidth_truth = segwidth_truth

'# end defaults

titletxt1 = plot_inputs$titletxt1
pdffn = plot_inputs$pdffn
doPDF = plot_inputs$doPDF
siminf_stats_good = plot_inputs$siminf_stats_good
brate_plot = plot_inputs$brate_plot
drate_plot = plot_inputs$drate_plot
b_exp_plot = plot_inputs$b_exp_plot
d_exp_plot = plot_inputs$d_exp_plot
uniq_params_txt = plot_inputs$uniq_params_txt
segwidth_truth = plot_inputs$segwidth_truth
simtype = plot_inputs$simtype

if (is.null(titletxt1))
	{
	titletxt1 = bquote(paste(.(simtype), " simulations (", lambda, "=", .(brate_plot), ", ",  mu, "=", .(drate_plot), ", ", alpha, "=", .(b_exp_plot), ", ", omega, "=", .(d_exp_plot), ")", sep=""))
	}



#######################################################
# Spacer between models
#######################################################
xspace_value = 0.5

if (is.null(xwidth_manual))
	{
	xwidth = uniq_params_num + (uniq_params_num-1) * xspace_value
	} else {
	xwidth = xwidth_manual
	}

#######################################################
# Make plots of DEC vs. DEC+J parameter inference and model choice
#######################################################
if(plot_inputs$doPDF == TRUE)
	{
	pdf(pdffn, height=6, width=6)
	}

par(oma = c(5,4,4,1))
par(mar=c(0,2,1,1))
nums = 1:4
mat = matrix(data=nums, nrow=4, ncol=4, byrow=FALSE)
mat

layout(mat)

i=1
ymax = 0.20
yaxis_ticks = c(0, 0.05, 0.1, 0.15, 0.2, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote(italic("d")), side=2, las=1, line=4)

# Plot the page title
mtext(text=titletxt1, side=3, line=0, cex=1.2, outer=TRUE)
#mtext(text=titletxt2, side=3, line=-0.5, cex=0.8, outer=TRUE)

# Plot a legend
xpos = 1
points(x=xpos, y=0.925*ymax, pch=18, cex=2, col="black")
segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=0.925*ymax, y1=0.925*ymax, lwd=2, col="black")
text(x=xpos+1.5*segwidth_truth, y=0.925*ymax, labels="true parameter value", pos=4, offset=1)
points(x=xpos, y=0.725*ymax, pch=22, cex=3, col="black", bg="white")
text(x=xpos+1.5*segwidth_truth, y=0.725*ymax, labels="DEC inference", pos=4, offset=1)
points(x=xpos, y=0.525*ymax, pch=22, cex=3, col="black", bg="lightgray")
text(x=xpos+1.5*segwidth_truth, y=0.525*ymax, labels="DEC+J inference", pos=4, offset=1)



# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	siminf_stats_sim = siminf_stats_subset[TF,]
	siminf_stats_sim
	
	
	dtrue = unique(siminf_stats_sim$d)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=siminf_stats_sim$d_DECinf, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=siminf_stats_sim$d_DECJinf, at=xpos, add=TRUE, col="lightgray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	points(x=xpos, y=dtrue, pch=18, cex=2, col="black")
	segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=dtrue, y1=dtrue, lwd=2, col="black")
	
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through d truth/inference plots




# Cut "e"




i=1
ymax = 0.5
yaxis_ticks = pretty(c(0, ymax))
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote(italic("j")), side=2, las=1, line=4)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	siminf_stats_sim = siminf_stats_subset[TF,]
	siminf_stats_sim
	
	
	jtrue = unique(siminf_stats_sim$j)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=siminf_stats_sim$j_DECinf, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=siminf_stats_sim$j_DECJinf, at=xpos, add=TRUE, col="lightgray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	points(x=xpos, y=jtrue, pch=18, cex=2, col="black")
	segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=jtrue, y1=jtrue, lwd=2, col="black")

	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through j truth/inference plots





#######################################################
# Accuracy
#######################################################
i=1
ymax = 1
yaxis_ticks =c(1/15, 0.25, 0.5, 0.75, ymax)
yaxis_ticks_txt =c("1/15", "0.25", "0.5", "0.75", "1")

plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks_txt, las=1)
mtext(text="Proportion correct", side=2, las=3, line=3.5, cex=1)

# Line representing random guess among states
abline(h=1/15, lty="dotted", col="black", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	
	# Calculate accuracy
	MLstate_TF_ancstate_accuracies_DEC_subset = ancstate_accuracies_DEC_subset[TF,99:147]
	MLstate_accuracy_DEC = apply(X=MLstate_TF_ancstate_accuracies_DEC_subset, MARGIN=1, FUN=mean, na.rm=TRUE)

	MLstate_TF_ancstate_accuracies_DECJ_subset = ancstate_accuracies_DECJ_subset[TF,99:147]
	MLstate_accuracy_DECJ = apply(X=MLstate_TF_ancstate_accuracies_DECJ_subset, MARGIN=1, FUN=mean, na.rm=TRUE)
	
	
	jtrue = unique(siminf_stats_sim$j)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=MLstate_accuracy_DEC, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=MLstate_accuracy_DECJ, at=xpos, add=TRUE, col="lightgray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	#points(x=i, y=jtrue, pch="*", cex=4, col="darkblue")

	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through accuracy boxplots








#######################################################
# Model choice -- with AICc
#######################################################
i=1
ymax = 130
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 50, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote(paste(Delta, "AICc", sep=""))
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3.5, cex=1)

# Dashed red line at the significance level
abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	siminf_stats_sim = siminf_stats_subset[TF,]
	siminf_stats_sim
	
	LnL_DEC = siminf_stats_sim$LnL_DEC
	LnL_DECJ = siminf_stats_sim$LnL_DECJ
	
	# Boxplot of the log-likelihood advantage of DEC+J
	LnL_advantage = LnL_DECJ-LnL_DEC
	# Min on this plot is 0.1
	# Also a few optimization issues at high d, e values
	# (could be datasets with no signal left; functions will misfire if LnL DEC>LnL DEC+J)
	TF = LnL_advantage < 0
	LnL_DECJ[TF] = LnL_DEC[TF]
	
	if (sum(TF) > 0)
		{
		cat("Note: In ", sum(TF), " cases, LnL DEC > LnL DEC+J.\n", sep="")
		# This may be very minor, or may indicate a failure of optimx ML routine
		# to find the true optimum. For now, in these cases we are setting LnL 
		# DEC+J to equal LnL DEC.
		}
	
	# Calculate DeltaAICc 
	AICc_DEC = getAICc(LnL=LnL_DEC, numparams=2, samplesize=100)
	AICc_DECJ = getAICc(LnL=LnL_DECJ, numparams=3, samplesize=100)
	delta_AICc_vals = AICc_DEC - AICc_DECJ

	# Plot AICc weights
	# DEC or DEC+J
	xpos = extra_xspace + i + 0.0
	
	# Plot one box (grey or white)
	# AICc weights look stupid
	jtrue = unique(siminf_stats_sim$j)
	if (jtrue > 0)
		{
		# width makes fatter boxes
		b1 = boxplot(x=delta_AICc_vals, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n", col="lightgray", border="black", width=3)
		} else {
		# width makes fatter boxes
		b1 = boxplot(x=delta_AICc_vals, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n", col="lightgray", border="black", width=3)
		}
	
	# Add "true model:" text
	if (i==1)
		{
		txtplot = "True model:"
		mtext(side=1, at=-0.25, text=txtplot, cex=0.55, las=1, line=0)
		}
	
	if (jtrue > 0)
		{
		txtplot = "+J"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		} else {
		txtplot = "DEC"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		}
	
		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot


if(plot_inputs$doPDF == TRUE)
	{
	dev.off()
	cmdstr = paste("open ", pdffn, sep="")
	system(cmdstr)
	}



	
	} # end plot_inf_accuracy_vs_SSEsims_v2()






















#######################################################
# Make plots of Tree Statistics for these simulations,
# to show variability in the simulated trees
#######################################################



#######################################################
# Make plots of keys statistics for simulated tree/event histories
#######################################################

plot_SSEsims_treestats_pt1 <- function(plot_inputs, xwidth_manual=NULL, sizes=c(0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4))
	{
	

defaults='
titletxt1 = bquote(paste("100 SSE simulations (", lambda, "=0.3, ", mu, "=0.3, ", alpha, "=1, ", omega, "=-1)", sep=""))
titletxt2 = paste("fossil (white) vs. observed (grey) tree stats", sep="")
pdffn = "simtree_stats_SSE_b03d03bx1dx-1.pdf"

simtype = "SSE"
segwidth_truth = 0.05
plot_inputs$segwidth_truth = segwidth_truth

simtype = "SSE"
brate_plot = 0.3
drate_plot = 0.3
b_exp_plot = 1
d_exp_plot = -1

plot_inputs = NULL
plot_inputs$titletxt1 = NULL
plot_inputs$pdffn = pdffn
plot_inputs$doPDF = TRUE
plot_inputs$siminf_stats_good = siminf_stats_good
plot_inputs$sim_tipstates_good = sim_tipstates_good
plot_inputs$ancstate_accuracies_DEC_good = ancstate_accuracies_DEC_good
plot_inputs$ancstate_accuracies_DECJ_good = ancstate_accuracies_DECJ_good
plot_inputs$simevent_counts_good = simevent_counts_good
plot_inputs$treeheights_incl_fossil_good = treeheights_incl_fossil_good

plot_inputs$simtype = simtype
plot_inputs$brate_plot = brate_plot
plot_inputs$drate_plot = drate_plot
plot_inputs$b_exp_plot = b_exp_plot
plot_inputs$d_exp_plot = d_exp_plot
plot_inputs$uniq_params_txt = uniq_params_txt
plot_inputs$segwidth_truth = segwidth_truth

sizes = c(0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4)
'# end defaults


# Load the SSE parameters, subset the data

titletxt1 = plot_inputs$titletxt1
pdffn = plot_inputs$pdffn
doPDF = plot_inputs$doPDF
siminf_stats_good = plot_inputs$siminf_stats_good
sim_tipstates_good = plot_inputs$sim_tipstates_good
ancstate_accuracies_DEC_good = plot_inputs$ancstate_accuracies_DEC_good
ancstate_accuracies_DECJ_good = plot_inputs$ancstate_accuracies_DECJ_good
simevent_counts_good = plot_inputs$simevent_counts_good
treeheights_incl_fossil_good = plot_inputs$treeheights_incl_fossil_good

brate_plot = plot_inputs$brate_plot
drate_plot = plot_inputs$drate_plot
b_exp_plot = plot_inputs$b_exp_plot
d_exp_plot = plot_inputs$d_exp_plot
uniq_params_txt = plot_inputs$uniq_params_txt
segwidth_truth = plot_inputs$segwidth_truth
simtype = plot_inputs$simtype

TF1 = siminf_stats_good$brate == brate_plot
TF2 = siminf_stats_good$drate == drate_plot
TF3 = siminf_stats_good$b_exp == b_exp_plot
TF4 = siminf_stats_good$d_exp == d_exp_plot
TF = (TF1 + TF2 + TF3 + TF4) == 4


# Do the subsetting
siminf_stats_subset = dfnums_to_numeric(siminf_stats_good[TF, ])
sim_tipstates_subset = sim_tipstates_good[TF, ]
ancstate_accuracies_DEC_subset = ancstate_accuracies_DEC_good[TF,]
ancstate_accuracies_DECJ_subset = ancstate_accuracies_DECJ_good[TF,]
simevent_counts_subset = simevent_counts_good[TF,]
treeheights_incl_fossil_subset = treeheights_incl_fossil_good[TF]


if (is.null(titletxt1))
	{
	titletxt1 = bquote(paste("Tree stats for ", .(simtype), " simulations (", lambda, "=", .(brate_plot), ", ",  mu, "=", .(drate_plot), ", ", alpha, "=", .(b_exp_plot), ", ", omega, "=", .(d_exp_plot), ")", sep=""))
	}
titletxt1


#######################################################
# Spacer between models
#######################################################
xspace_value = 0.5

if (is.null(xwidth_manual))
	{
	xwidth = uniq_params_num + (uniq_params_num-1) * xspace_value
	} else {
	xwidth = xwidth_manual
	}

#######################################################
# Make plots of Tree Statistics
#######################################################
if(plot_inputs$doPDF == TRUE)
	{
	pdf(pdffn, height=11, width=8.5)
	}

par(oma = c(5,4,4,1))
par(mar=c(0,2,1,1))
nums = 1:5
mat = matrix(data=nums, nrow=5, ncol=5, byrow=FALSE)
mat

layout(mat)


#######################################################
# Tree heights, all and observed
#######################################################

i=1
ymax = 150
yaxis_ticks = c(0, 50, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote("Tree height"), side=2, las=3, line=3)

# Plot the page title
mtext(text=titletxt1, side=3, line=0, cex=1.2, outer=TRUE)
#mtext(text=titletxt2, side=3, line=-0.5, cex=0.8, outer=TRUE)

# Plot a legend
xpos = 1
#points(x=xpos, y=0.925*ymax, pch=18, cex=2, col="black")
#segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=0.925*ymax, y1=0.925*ymax, lwd=2, col="black")
#text(x=xpos+1.5*segwidth_truth, y=0.925*ymax, labels="true parameter value", pos=4, offset=1)
points(x=xpos, y=0.925*ymax, pch=22, cex=3, col="black", bg="white")
text(x=xpos+1.5*segwidth_truth, y=0.925*ymax, labels="with fossils", pos=4, offset=1)
points(x=xpos, y=0.825*ymax, pch=22, cex=3, col="black", bg="gray")
text(x=xpos+1.5*segwidth_truth, y=0.825*ymax, labels="observed tree", pos=4, offset=1)



# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	observed_tree_ages = ancstate_accuracies_DEC_subset[TF,1]
	observed_tree_ages[c(1:5, 96:100)]
	actual_tree_ages = treeheights_incl_fossil_subset[TF]
	actual_tree_ages[c(1:5, 96:100)]
	
	
	#dtrue = unique(siminf_stats_sim$d)
	
	# Actual tree age (with fossils)
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=actual_tree_ages, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Observed tree ages (excluding fossils, extinct nodes, etc.)
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=observed_tree_ages, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	#points(x=xpos, y=dtrue, pch=18, cex=2, col="black")
	#segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=dtrue, y1=dtrue, lwd=2, col="black")
	
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through d truth/inference plots



#######################################################
# Number of tips, fossil and observed
#######################################################

i=1
ymax = 400
yaxis_ticks = c(0, 200, ymax)

plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text="# species", side=2, las=3, line=3)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	simevent_counts_sim
	

	# Line representing number of species observed
	#abline(h=50, lty="dotted", col="black", lwd=0.5)
	axis(side=2, at=50, labels="50 species\nobserved", las=1, cex=0.5)
	#mtext(side=2, at=50, text="observed: 50", las=1, line=1, cex=0.55, lty=1)


	#etrue = unique(siminf_stats_sim$e)
	
	# Number of tips in true tree, including fossils
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=50+simevent_counts_sim$num_fossil_tips, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Number of tips in observed tree (always 50)
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=50+simevent_counts_sim$num_fossil_tips-simevent_counts_sim$num_fossil_tips, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	#points(x=xpos, y=etrue, pch=18, cex=2, col="black")
	#segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=etrue, y1=etrue, lwd=2, col="black")
	
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through e truth/inference plots



#######################################################
# Tip range sizes
#######################################################
# Range sizes
#sizes = c(0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4)
#sizes


foo <- function(val, sizes)
	{
	sizes[val]
	}



i=1
ymax = 4
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 1, 2, 3, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote("range size"), side=2, las=3, line=3)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	sim_tipstates_sim = sim_tipstates_subset[TF,]
	sim_tipstates_sim2 = sim_tipstates_sim - 1
	
	sim_tipsizes = sapply(X=sim_tipstates_sim2, FUN=foo, sizes=sizes)
	sim_tipsizes2 = matrix(data=sim_tipsizes, nrow=nrow(sim_tipstates_sim2), ncol=ncol(sim_tipstates_sim2))

	meanSizes = rowMeans(sim_tipsizes2, na.rm=FALSE)


	
	#jtrue = unique(siminf_stats_sim$j)
	
	xpos = extra_xspace + i-0.0
	b1 = boxplot(x=meanSizes, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	#points(x=xpos, y=jtrue, pch=18, cex=2, col="black")
	#segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=jtrue, y1=jtrue, lwd=2, col="black")

	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through j truth/inference plots





#######################################################
# Number of range-change events
#######################################################

i=1
ymax = 400
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 200, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text=bquote("# events (total)"), side=2, las=3, line=3)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$y_actual + simevent_counts_sim$s_actual + simevent_counts_sim$v_actual + simevent_counts_sim$j_actual + simevent_counts_sim$d_actual + simevent_counts_sim$e_actual + simevent_counts_sim$a_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$y_observed + simevent_counts_sim$s_observed + simevent_counts_sim$v_observed + simevent_counts_sim$j_observed + simevent_counts_sim$d_observed + simevent_counts_sim$e_observed + simevent_counts_sim$a_observed
	
	
	#jtrue = unique(siminf_stats_sim$j)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	xpos = extra_xspace + i
	#points(x=xpos, y=jtrue, pch=18, cex=2, col="black")
	#segments(x0=xpos-segwidth_truth, x1=xpos+segwidth_truth, y0=jtrue, y1=jtrue, lwd=2, col="black")

	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through j truth/inference plots





#######################################################
# Num events - cladogenetic, nonsympatric
#######################################################
i=1
ymax = 200
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 100, ymax)

plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(0, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
mtext(text="# events (clad.)", side=2, las=3, line=3, cex=1)

# Line representing random guess among states
#abline(h=1/15, lty="dotted", col="black", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$s_actual + simevent_counts_sim$v_actual + simevent_counts_sim$j_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$s_observed + simevent_counts_sim$v_observed + simevent_counts_sim$j_observed
		
	#jtrue = unique(siminf_stats_sim$j)
	
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	# Plot the true value
	#points(x=i, y=jtrue, pch="*", cex=4, col="darkblue")




	# Plot one box (grey or white)
	# AICc weights look stupid
	siminf_stats_sim = siminf_stats_subset[TF,]
	jtrue = unique(siminf_stats_sim$j)
	xpos = extra_xspace + i+0.0
	
	# Add "true model:" text
	if (i==1)
		{
		txtplot = "True model:"
		mtext(side=1, at=-1, text=txtplot, cex=0.55, las=1, line=0)
		}
	
	if (jtrue > 0)
		{
		txtplot = "+J"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		} else {
		txtplot = "DEC"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		}
	



	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	} # End loop through accuracy boxplots







par(oma = c(5,4,4,1))
par(mar=c(0,2,1,1))
nums = 1:6
mat = matrix(data=nums, nrow=6, ncol=6, byrow=FALSE)
mat

layout(mat)



#######################################################
# Number of "d" events
#######################################################
i=1
ymax = 200
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote("# events (d)")
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3, cex=1)

# Dashed red line at the significance level
#abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]


	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$d_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$d_observed
		
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")





		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot





#######################################################
# 2nd page (counts of e, y, s, v, j)
#######################################################


titletxt1_cont = bquote(paste("Tree stats for ", .(simtype), " simulations, cont. (", lambda, "=", .(brate_plot), ", ",  mu, "=", .(drate_plot), ", ", alpha, "=", .(b_exp_plot), ", ", omega, "=", .(d_exp_plot), ")", sep=""))

titletxt1_cont



#######################################################
# Number of "e" events
#######################################################
i=1
ymax = 200
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote("# events (e)")
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3, cex=1)


# Plot the page title
mtext(text=titletxt1_cont, side=3, line=0, cex=1.2, outer=TRUE)


# Dashed red line at the significance level
#abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$e_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$e_observed
		
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")

		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot






#######################################################
# Number of "y" events
#######################################################
i=1
ymax = 200
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote("# events (y)")
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3, cex=1)

# Dashed red line at the significance level
#abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$y_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$y_observed
		
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")

		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot









#######################################################
# Number of "s" events
#######################################################
i=1
ymax = 200
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote("# events (s)")
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3, cex=1)

# Dashed red line at the significance level
#abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$s_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$s_observed
		
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")

		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot










#######################################################
# Number of "v" events
#######################################################
i=1
ymax = 200
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote("# events (v)")
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3, cex=1)

# Dashed red line at the significance level
#abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$v_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$v_observed
		
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")

		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot






#######################################################
# Number of "j" events
#######################################################
i=1
ymax = 200
#yaxis_ticks = pretty(c(0, ymax))
yaxis_ticks = c(0, 100, ymax)
plot(1:xwidth, 1:xwidth, pch=".", col="white", xlim=c(1,xwidth), ylim=c(-5, ymax), xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=2, at=yaxis_ticks, labels=yaxis_ticks, las=1)
xlabel_txt = bquote("# events (j)")
xlabel_txt
mtext(text=xlabel_txt, side=2, las=3, line=3, cex=1)

# Dashed red line at the significance level
#abline(h=0, lty="dashed", col="darkgrey", lwd=0.5)

# Reset extra_xspace
extra_xspace = 0
# Loop through param combinations
for (i in 1:uniq_params_num)
	{
	if (is.na(uniq_params_txt[i]))
		{
		next()
		}
	TF = siminf_stats_subset$params_txt == uniq_params_txt[i]
	simevent_counts_sim = simevent_counts_subset[TF,]
	
	# Total # of actual events:
	num_events_actual = simevent_counts_sim$j_actual
	# Total # of observed events:
	num_events_observed = simevent_counts_sim$j_observed
		
	xpos = extra_xspace + i-0.25
	b1 = boxplot(x=num_events_actual, at=xpos, add=TRUE, outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")
	
	xpos = extra_xspace + i+0.25
	b1 = boxplot(x=num_events_observed, at=xpos, add=TRUE, col="gray", border="black", outline=FALSE, medlwd=1, medlty="dotted", xaxt="n", yaxt="n")



	# Plot one box (grey or white)
	# AICc weights look stupid
	siminf_stats_sim = siminf_stats_subset[TF,]
	jtrue = unique(siminf_stats_sim$j)
	xpos = extra_xspace + i+0.0
	
	# Add "true model:" text
	if (i==1)
		{
		txtplot = "True model:"
		mtext(side=1, at=-1, text=txtplot, cex=0.55, las=1, line=0)
		}
	
	if (jtrue > 0)
		{
		txtplot = "+J"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		} else {
		txtplot = "DEC"
		mtext(side=1, at=xpos, text=txtplot, cex=0.55, las=1, line=0)
		}
	

		
	# Add to the extra_xspace
	extra_xspace = extra_xspace + xspace_value
	
	} # End loop through model-choice plot





if(plot_inputs$doPDF == TRUE)
	{
	dev.off()
	cmdstr = paste("open ", pdffn, sep="")
	system(cmdstr)
	}



	
	} # end plot_SSEsims_treestats_pt1()













