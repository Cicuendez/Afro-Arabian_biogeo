readfiles_BioGeoBEARS_run_nexus <- function (inputs) 
{
  if (is.character(inputs$timesfn)) {
    inputs$timeperiods = read_times_fn(inputs)
  }
  if (is.character(inputs$distsfn)) {
    inputs$list_of_distances_mats = read_distances_fn(inputs)
  }
  if (is.character(inputs$envdistsfn)) {
    inputs$list_of_envdistances_mats = read_envdistances_fn(inputs = NULL, 
                                                            envdistsfn = inputs$envdistsfn)
  }
  if (is.character(inputs$dispersal_multipliers_fn)) {
    inputs$list_of_dispersal_multipliers_mats = read_dispersal_multipliers_fn(inputs)
  }
  if (is.character(inputs$area_of_areas_fn)) {
    inputs$list_of_area_of_areas = read_area_of_areas_fn(inputs)
  }
  if (is.character(inputs$areas_allowed_fn)) {
    inputs$list_of_areas_allowed_mats = read_areas_allowed_fn(inputs)
  }
  if (is.character(inputs$areas_adjacency_fn)) {
    inputs$list_of_areas_adjacency_mats = read_areas_adjacency_fn(inputs)
  }
  if (inputs$use_detection_model == TRUE) {
    phy = check_trfn_nexus(trfn = inputs$trfn)
    if (is.character(inputs$detects_fn)) {
      inputs$detects_df = read_detections(inputs$detects_fn, 
                                          OTUnames = NULL, areanames = NULL, tmpskip = 0, 
                                          phy = phy)
    }
    if (is.character(inputs$detects_fn)) {
      inputs$controls_df = read_controls(inputs$controls_fn, 
                                         OTUnames = NULL, areanames = NULL, tmpskip = 0, 
                                         phy = phy)
    }
  }
  return(inputs)
}
