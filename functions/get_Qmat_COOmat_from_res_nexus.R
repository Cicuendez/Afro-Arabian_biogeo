get_Qmat_COOmat_from_res_nexus <- function (res, numstates = NULL, include_null_range = TRUE, max_range_size = NULL, 
          timeperiod_i = 1) 
{
  setup = "\n\t# Set up a fake results object\n\tres = NULL\n\tres$inputs = define_BioGeoBEARS_run()\n\tres$outputs = res$inputs$BioGeoBEARS_model_object\n\tres\n\t\n\t# Produces stop error\n\tget_Qmat_COOmat_from_res(res, numstates=NULL, include_null_range=TRUE, max_range_size=NULL)\n\t\n\t# Works:\n\t# Get the matrices (trivial case, non-stratified Psychotria example)\n\treturned_mats = get_Qmat_COOmat_from_res(res, numstates=NULL, max_range_size=4)\n\treturned_mats\n\tnames(returned_mats)\n\t\n\t# Also works:\n\tres$inputs$max_range_size = 4\n\treturned_mats = get_Qmat_COOmat_from_res(res, numstates=NULL, max_range_size=NULL)\n\treturned_mats\n\tnames(returned_mats)\n\t"
  if (isblank_TF(res$inputs$max_range_size) == FALSE) {
    max_range_size = res$inputs$max_range_size
  }
  else {
    if (isblank_TF(max_range_size) == FALSE) {
      res$inputs$max_range_size = max_range_size
    }
  }
  if (is.null(max_range_size)) {
    txt = "STOP ERROR in get_Qmat_COOmat_from_res(): either max_range_size or res$inputs$max_range_size must be specified."
    cat("\n\n")
    cat(txt)
    stop(txt)
    cat("\n\n")
  }
  BioGeoBEARS_run_object = res$inputs
  BioGeoBEARS_model_object = res$output
  returned_mats = get_Qmat_COOmat_from_BioGeoBEARS_run_object_nexus(BioGeoBEARS_run_object = BioGeoBEARS_run_object, 
                                                              BioGeoBEARS_model_object = BioGeoBEARS_model_object, 
                                                              numstates = numstates, include_null_range = include_null_range, 
                                                              timeperiod_i = timeperiod_i)
  return(returned_mats)
}
