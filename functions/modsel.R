# modsel function to perform model selection after a BioGeoBEARS run on several clades.

# Input:
# results_list is a list of results objects from BioGeoBEARS, 
# named with the model: DEC, DECJ, DIVA, DIVAJ, BAYAREA, BAYAREAJ.
# Each element of the list is a list with the results for the different clades.

# treedata_list is a list of phylogenetic trees in treedata format (treeio package), 
# for which the BioGeoBEARS models in results list have been applied.

# Output: 
# - restable_list is a list of tables with the parameters of each model for each clade.
# - table_bestmodels is a table indicating the best-fit model for each clade. 
# - BESTMODEL_consensus is a list of matrices with the probabilities of each node for each clade. 
# - Qmat_consensus, COOmat_consensus and rootstate_consensus are lists of elements for each clade 
# that are needed to perform the simulations.

modsel <- function(results_list, treedata_list){
  
  restable_list <<- vector("list", length(treedata_list))
  names(restable_list) <<- names(treedata_list)
  
  BESTMODELconsensus <<- vector("list", length(treedata_list))
  names(BESTMODELconsensus) <<- names(treedata_list)
  
  table_bestmodels <<- data.frame(genera=names(treedata_list), bestmodel=NA)
  rownames(table_bestmodels) <<- names(treedata_list)
  
  Qmat_consensus <<- vector("list", length(treedata_list))
  names(Qmat_consensus) <<- names(treedata_list)
  COOmat_consensus <<- vector("list", length(treedata_list))
  names(COOmat_consensus) <<- names(treedata_list)
  rootstate_consensus <<- vector("list", length(treedata_list))
  names(rootstate_consensus) <<- names(treedata_list)
  
  for (i in 1:length(results_list[['DEC']])){
    
    restable <- data.frame()
    
    if ("DEC" %in% names(results_list)){
      resDEC <- results_list[['DEC']][[i]]
      LnL_DEC = get_LnL_from_BioGeoBEARS_results_object(resDEC)
      numparams_DEC = 2
      cols_DEC = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, 
                                                                returnwhat="table", 
                                                                addl_params=c("j"), 
                                                                paramsstr_digits=4)
      restable <- rbind(restable, cols_DEC)
    }
    
    if ("DECJ" %in% names(results_list)){
      resDECj <- results_list[['DECJ']][[i]]
      LnL_DECJ = get_LnL_from_BioGeoBEARS_results_object(resDECj)
      numparams_DECJ = 3
      cols_DECJ = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, 
                                                                 returnwhat="table", 
                                                                 addl_params=c("j"), 
                                                                 paramsstr_digits=4)
      restable <- rbind(restable, cols_DECJ)
    }
    
    if ("DIVA" %in% names(results_list)){
      resDIVALIKE <- results_list[['DIVA']][[i]]
      LnL_DIVA = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
      numparams_DIVA = 2
      cols_DIVA = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
      restable <- rbind(restable, cols_DIVA)
    }
    
    if ("DIVAJ" %in% names(results_list)){
      resDIVALIKEj <- results_list[['DIVAJ']][[i]]
      LnL_DIVAJ = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)
      numparams_DIVAJ = 3
      cols_DIVAJ = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
      restable <- rbind(restable, cols_DIVAJ)
    }
    
    if ("BAYAREA" %in% names(results_list)){
      resBAYAREALIKE <- results_list[['BAYAREA']][[i]]
      LnL_BAYAREA = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
      numparams_BAYAREA = 2
      cols_BAYAREA = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
      restable <- rbind(restable, cols_BAYAREA)
    }
    
    if ("BAYAREAJ" %in% names(results_list)){
      resBAYAREALIKEj <- results_list[['BAYAREAJ']][[i]]
      LnL_BAYAREAJ = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)
      numparams_BAYAREAJ = 3
      cols_BAYAREAJ = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
      restable <- rbind(restable, cols_BAYAREAJ)
    }
    
    # Table with LnL
    row.names(restable) = names(results_list)
    restable
    # Add AIC and AIC weight
    AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
    restable = cbind(restable, AICtable)
    #restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
    #restable_AIC_rellike
    
    
    if (length(treedata_list[[i]]@phylo$tip.label) >= 4){
      
      # Add AICc and AICc weight for groups with 4 or more species
      samplesize = length(treedata_list[[i]]@phylo$tip.label)
      AICctable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
      restable = cbind(restable, AICctable)
      
      restable_list[[i]] <<- restable
      
      # Model selection based on min(AICc)
      BESTMODELconsensus[[i]] <<- 
        results_list[[rownames(restable)[restable$AICc == min(restable$AICc)]]][[i]]$ML_marginal_prob_each_state_at_branch_top_AT_node
      
      table_bestmodels[names(treedata_list)[i], ]$bestmodel <<- 
        rownames(restable)[restable$AICc == min(restable$AICc)]
      
      Qmat_consensus[[i]] <<- 
        get_Qmat_COOmat_from_res_nexus(results_list[[rownames(restable)[restable$AICc == min(restable$AICc)]]][[i]])$Qmat
      
      COOmat_consensus[[i]] <<- 
        get_Qmat_COOmat_from_res_nexus(results_list[[rownames(restable)[restable$AICc == min(restable$AICc)]]][[i]])$COO_weights_columnar
      
      rootstate_consensus[[i]] <<- 
        which.max(results_list[[rownames(restable)[restable$AICc == min(restable$AICc)]]][[i]]$relative_probs_of_each_state_at_bottom_of_root_branch)
      
      
      # For groups with less than 4 species, AICc calculation returns an error. 
      # So we use AIC instead.
    } else {
      restable_list[[i]] <<- restable
      
      BESTMODELconsensus[[i]] <<- 
        results_list[[rownames(restable)[restable$AIC == min(restable$AIC)]]][[i]]$ML_marginal_prob_each_state_at_branch_top_AT_node
      
      table_bestmodels[names(treedata_list)[i], ]$bestmodel <<- 
        rownames(restable)[restable$AIC == min(restable$AIC)]
      
      Qmat_consensus[[i]] <<- 
        get_Qmat_COOmat_from_res_nexus(results_list[[rownames(restable)[restable$AIC == min(restable$AIC)]]][[i]])$Qmat
      
      COOmat_consensus[[i]] <<- 
        get_Qmat_COOmat_from_res_nexus(results_list[[rownames(restable)[restable$AIC == min(restable$AIC)]]][[i]])$COO_weights_columnar
      
      rootstate_consensus[[i]] <<- 
        which.max(results_list[[rownames(restable)[restable$AIC == min(restable$AIC)]]][[i]]$relative_probs_of_each_state_at_bottom_of_root_branch)
      
    }
    
    
    
  }
  
  r <- list(table_bestmodels = table_bestmodels, 
            BESTMODELconsensus = BESTMODELconsensus, 
            restable_list = restable_list,
            Qmat_consensus = Qmat_consensus, 
            COOmat_consensus = COOmat_consensus, 
            rootstate_consensus = rootstate_consensus)
  
  return(r)
  
}
