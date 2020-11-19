
#######################################################
# Simulate and infer DNA under standard and cladogenetic models
#######################################################
library(BioGeoBEARS)
source('/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_DNA_cladogenesis_sim_v1.R', chdir = TRUE)


# Kimura 1980 (Kimura 2-parameter
# rate of:

# transitions	(a1, A<->G and C<->T)
# transversions	(a2, A<->C, A<->T, C<->G, G<->T)

A = c("-", "a2", "a1", "a2")
C = c("a2", "-", "a2", "a1")
G = c("a1", "a2", "-", "a2")
T = c("a2", "a1", "a2", "-")

Qmat_txt = adf2(rbind(A, C, G, T))
names(Qmat_txt) = c("A", "C", "G", "T")
Qmat_txt

a_params = adf2(matrix(data=c(0.25, 0.25), ncol=1))
names(a_params) = "a"
rownames(a_params) = c("a1", "a2")
a_params



# The "jump" equivalent would be a cladogenetic model

ancstate1 = c("A", "A", "A", "A", "C", "C", "C", "C", "G", "G", "G", "G", "T", "T", "T", "T")
Lstate1 = c("A", "A", "A", "A", "C", "C", "C", "C", "G", "G", "G", "G", "T", "T", "T", "T")
Rstate1 = c("A", "C", "G", "T", "A", "C", "G", "T", "A", "C", "G", "T", "A", "C", "G", "T")

ancstate2 = c("A", "A", "A", "C", "C", "C", "G", "G", "G", "T", "T", "T")
Lstate2 = c("C", "G", "T", "A", "G", "T", "A", "C", "T", "A", "C", "G")
Rstate2 = c("A", "A", "A", "C", "C", "C", "G", "G", "G", "T", "T", "T")

ancstate = c(ancstate1, ancstate2)
Lstate = c(Lstate1, Lstate2)
Rstate = c(Rstate1, Rstate2)

cladogenesis_txt = adf2(cbind(ancstate, Lstate, Rstate))
tmp_rownames = paste(ancstate, "->", Lstate, ",", Rstate, sep="")
rownames(cladogenesis_txt) = tmp_rownames

# For cladogenesis, we will also have 2 parameters, controlling
# rate of transitions (c1) and transversions (c2)
#
# transitions	(c1, A<->G and C<->T)
# transversions	(c2, A<->C, A<->T, C<->G, G<->T)
#
# These are rates relative to the "stay the same" rate (e.g. A->A,A),
# which will always get weight 1. As with BioGeoBEARS, the weights will be 
# added and then divided by the sum of weights to produce the 
# 
# Transitions of type A->G,T are disallowed (probability 0)
# 
rates = c(	"-", "c2", "c1", "c2", 
			"c2", "-", "c2", "c1", 
			"c1", "c2", "-", "c2", 
			"c2", "c1", "c2", "-", 
			      "c2", "c1", "c2", 
			"c2",       "c2", "c1", 
			"c1", "c2",       "c2", 
			"c2", "c1", "c2"       )

type = rates
type[rates=="c1"] = "transition"
type[rates=="c2"] = "transversion"
type[rates=="-"] = "none"


cladogenesis_txt = cbind(cladogenesis_txt, rates, type)
cladogenesis_txt


# Calculate the anagenetic and cladogenetic transition probability matrices
a1 = 0.1
a2= 0.02
Qmat = DNA_anagenesis_to_Qmat(a1, a2, Qmat_txt)
Qmat

c1 = 0
c2 = 0
cladogenesis_probs = DNA_cladogenesis_to_inheritance_condprobs(c1, c2, cladogenesis_txt)
cladogenesis_probs


# Generate a Yule tree
# Simulate trees under the Yule (pure-birth) process
library(geiger) # for sim.bdtree()
birthrate = 0.3
deathrate = 0
numtips = 50
seedset = set.seed(54321)

# extinct=TRUE doesn't matter for a pure-birth process
yuletree = sim.bdtree(b=birthrate, d=deathrate, stop="taxa", n=numtips, extinct=TRUE)
ntips = length(yuletree$tip.label)
num_internal_nodes = yuletree$Nnode
numnodes = ntips + num_internal_nodes
nodenums = c(1:numnodes)
tip_nodenums = 1:ntips
internal_nodenums = (ntips+1):numnodes


simulated_states_by_node = sim_DNA_clado(phy=yuletree, Qmat, cladogenesis_probs, starting_state=1)
simulated_states_by_node

simulated_states_by_node_txt = as.character(simulated_states_by_node)
simulated_states_by_node_txt[simulated_states_by_node_txt=="1"] = "A"
simulated_states_by_node_txt[simulated_states_by_node_txt=="2"] = "C"
simulated_states_by_node_txt[simulated_states_by_node_txt=="3"] = "G"
simulated_states_by_node_txt[simulated_states_by_node_txt=="4"] = "T"

node_simstates = simulated_states_by_node_txt[internal_nodenums]
tip_simstates = simulated_states_by_node_txt[tip_nodenums]

plot(yuletree)
title("Tree simulated under the pure-birth process, DNA model: K2p")
axisPhylo()
nodelabels(node_simstates)
tiplabels(tip_simstates)



tip_condlikes_of_data_on_each_state = tip_simstate_nums_to_tip_condlikes_of_data_on_each_state(simulated_states_by_node, numtips, numstates=NULL)
rowSums(tip_condlikes_of_data_on_each_state)


LnL = calc_loglike_for_optim(params, phy, Qmat_txt, cladogenesis_txt, tip_condlikes_of_data_on_each_state)
LnL



params = c(a1, a2, c1, c2)
params_lower = c(0, 0, 0, 0)
params_upper = c(2, 2, 2, 2)

optim_result = optim(par=params, fn=calc_loglike_for_optim, phy=yuletree, Qmat_txt=Qmat_txt, cladogenesis_txt=cladogenesis_txt, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, method="L-BFGS-B", lower=params_lower, upper=params_upper, control=c(fnscale=-1))

optim_result$value 
LnL


total_loglikelihood = calc_DNA_loglike(tip_condlikes_of_data_on_each_state, phy, Qmat, cladogenesis_probs, returnwhat="loglike")
total_loglikelihood

