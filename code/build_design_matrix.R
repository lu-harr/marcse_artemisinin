#' Build up the design matrix from covariate rasters and an optiona set of
#' coordinates
#'
#' @title Build Covariate Design Matrix
#' @param covariates A multilayer SpatRaster with the values of the covariates
#' @param coords An optional matrix of coordinates at which to extract the
#'   covariate values and construct the design matrix
#' @param scale Should the covariates be scaled (to mean zero, standard
#'   deviation 1) before extraction?
#' @return
#' @author Nick Golding; edited by Lucy Harrison to account for temporally varying covariates
#' @export
build_design_matrix <- function(covariates, coords = NULL, scale = TRUE, 
                                temporal_var = TRUE, temporal_range) {
  
  # optionally scale covariates
  if (scale) {
    covariates <- scale(covariates)
  }
  
  if (temporal_var){
    # don't want to cellFromXY loads of times but the functionality is missing to XYZ
    # pre-treat years ... this is super clunky and I do not like it ... come back
    message(names(coords))
    coords$year_truncated <- case_when(coords$year < min(temporal_range) ~ min(temporal_range),
                                       coords$year > max(temporal_range) ~ max(temporal_range),
                                        .default = coords$year)
    
    years <- unique(coords$year_truncated)
    covs <- rep(NA, nrow(coords)) # should probably be a df for more than one covariate
    
    for (year in years){
      # we're assuming that all of the covariates are labelled with a year 
      # (and none are labelled with no year)
      # and we're assuming that all years have the same number of covariates
      
      coord_idx <- which(coords$year_truncated == year)
      lyr_idx <- grep(as.character(year), names(covariates))
      cell_ids <- terra::cellFromXY(covariates[[lyr_idx]], 
                                    as.data.frame(coords[coord_idx, c("x","y")]))
      covs[coord_idx] <- terra::extract(covariates[[lyr_idx]], cell_ids) %>%
        unlist()
    }
  }
  
  else {
    # if coords is provided use those cells, otherwise use all cells
    if (is.null(coords)) {
      cell_ids <- terra::cells(covariates)
      
    } else {
      
      cell_ids <- terra::cellFromXY(covariates, coords[,c("x", "y")])
    }
    
    # extract, pad with intercept dummy, and return
    covs <- terra::extract(covariates, cell_ids)
  }
  
  cbind(Intercept = 1, 
        year = if ("scaled_year" %in% names(coords)) {coords$scaled_year} else {coords$year}, # bad 
        pfpr = as.matrix(covs))
  # add year as a covt?
  
}
