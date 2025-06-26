#' Build up the design matrix from covariate rasters and an optiona set of
#' coordinates
#'
#' @title Build Covariate Design Matrix
#' @param covariates A multilayer SpatRaster with the values of the covariates
#' @param coords An optional matrix of coordinates at which to extract the
#'   covariate values and construct the design matrix
#' @param scale Should the covariates be scaled (to mean zero, standard
#'   deviation 1) before extraction?
#' @param temporal_var Is there temporal variance in the `covariates`?
#' @param temporal_range Vector of numerics with length equal or greater than two
#' @param degs_to_rads Boolean: are coords in degrees? Would user like to convert 
#' them to radians? Requires columns `x` and `y` in `coords`.
#' @param obs_cols Optional argument: names of columns to include in design matrix, 
#' a subset of full set returned by default
#' @return
#' @author Nick Golding; edited by Lucy Harrison to account for temporally 
#' varying covariates et al
#' @export

scale_years <- function(temporal_range){
  years <- do.call(seq, as.list(temporal_range))
  scaled_years <- scale(years)
  # fiddle with output of scale:
  scaled_years <- lapply(years, function(x){scaled_years[x - min(temporal_range) + 1]})
  names(scaled_years) <- years
  
  scaled_years
} 


build_design_matrix <- function(covariates, 
                                coords = NULL, 
                                scale = TRUE, 
                                temporal_var = FALSE, 
                                temporal_range = NULL,
                                degs_to_rads = FALSE,
                                obs_cols = NULL) {
  
  # optionally scale covariates
  if (scale) {
    covariates <- scale(covariates)
  }
  
  if (length(intersect(c("x", "y"), names(coords))) != 2){
    message("`x` and `y` not in `names(coords)`")
    return(FALSE)
  }
  
  if (temporal_var){
    # this is little silly
    temporal_range <- range(temporal_range)
    
    # don't want to cellFromXY loads of times but the functionality is missing to XYZ
    # pre-treat years ... this is super clunky and I do not like it ... come back
    coords$year_truncated <- case_when(coords$year < min(temporal_range) ~ min(temporal_range),
                                       coords$year > max(temporal_range) ~ max(temporal_range),
                                        .default = coords$year)
    
    # scale years after they are truncated ........ 
    scaled_years <- scale_years(temporal_range)
    #coords$year_scaled <- scaled_years[coords$year_truncated - min(temporal_range) + 1]
    coords$year_scaled <- unlist(scaled_years[as.character(coords$year_truncated)])
    
    years <- unique(coords$year_truncated)
    covs <- rep(NA, nrow(coords)) 
    # should be a df for more than one covariate ..
    
    for (year in years){
      
      coord_idx <- which(coords$year_truncated == year)
      lyr_idx <- grep(as.character(year), names(covariates))
      
      if (length(lyr_idx) == 0){
        
        # we're assuming that all of the covariates are labelled with a year 
        # (and none are labelled with no year)
        # and we're assuming that all years have the same number of covariates
        message(paste0("Warning: zero matches to year ", year))
      }
      cell_ids <- terra::cellFromXY(covariates[[lyr_idx]], 
                                    as.data.frame(coords[coord_idx, c("x","y")]))
      covs[coord_idx] <- terra::extract(covariates[[lyr_idx]], cell_ids) %>%
        unlist()
    }
  }
  
  else {
    scaled_years = NULL
    
    # if coords is provided use those cells, otherwise use all cells
    if (is.null(coords)) {
      cell_ids <- terra::cells(covariates)
      
    } else {
      
      cell_ids <- terra::cellFromXY(covariates, coords[,c("x", "y")])
    }
    
    # extract, pad with intercept dummy, and return
    covs <- terra::extract(covariates, cell_ids) %>%
      unlist()
  }
  
  if (degs_to_rads){
    coords$x_rd <- degrees_to_radians(coords$x)
    coords$y_rd <- degrees_to_radians(coords$y)
  }
  
  df <- cbind(intercept = 1, 
                coords,
                pfpr = as.matrix(covs))
  
  if (!is.null(obs_cols)){
    df <- df[[obs_cols]]
  }
  
  return(list(df = df,
              scaled_years = scaled_years))
}
