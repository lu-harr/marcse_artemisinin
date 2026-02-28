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
#' @param temporal_range Vector of numerics with length equal or greater than 
#' two - the years represented in the covariate set. User can specify so they can 
#' be more or less extreme than what is represented in data
#' @param degs_to_rads Boolean: are coords in degrees? Would user like to convert 
#' them to radians? Requires columns `x` and `y` in `coords`.
#' @param obs_cols Optional argument: names of columns to include in design matrix, 
#' a subset of full set returned by default
#' @return
#' @author Nick Golding; edited by Lucy Harrison to account for temporally 
#' varying covariates et al
#' @export
build_design_matrix <- function(covariates, 
                                coords = NULL, 
                                scale = TRUE, 
                                temporal_var = FALSE, 
                                temporal_covt_range = NULL,
                                degs_to_rads = FALSE,
                                obs_cols = NULL,
                                buffer = 0) {
  
  # optionally scale covariates
  if (scale) {
    covariates <- scale(covariates)
  }
  
  if (length(intersect(c("x", "y"), names(coords))) != 2){
    message("`x` and `y` not in `names(coords)`")
    return(FALSE)
  }
  
  if (temporal_var){
    # don't want to cellFromXY loads of times but the functionality is missing to XYZ
    
    temporal_covt_range <- range(temporal_covt_range)
    
    # scale years before they're truncated !
    scaled_years <- scale_years(range(coords$year))
    coords$year_scaled <- unlist(scaled_years[as.character(coords$year)])
    
    coords$year_truncated <- case_when(coords$year < min(temporal_covt_range) ~ min(temporal_covt_range),
                                       coords$year > max(temporal_covt_range) ~ max(temporal_covt_range),
                                       .default = coords$year)
    
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
      
      if(length(coord_idx) == 1){
        covs[coord_idx] <- terra::extract(covariates[[lyr_idx]], 
                                          as.data.frame(coords[coord_idx, c("x","y")]), 
                                          search_radius = buffer, 
                                          ID = FALSE)[1,1] # don't ask :/
      } else {
        covs[coord_idx] <- terra::extract(covariates[[lyr_idx]], 
                                          as.data.frame(coords[coord_idx, c("x","y")]), 
                                          search_radius = buffer, 
                                          ID = FALSE) %>%
          dplyr::select(names(covariates)[[lyr_idx]]) %>%
          unlist()
      }
      
      
    }
  }
  
  else {
    scaled_years = NULL
    
    # if coords is provided use those cells, otherwise use all cells
    if (is.null(coords)) {
      coords <- crds(covariates)
      
    }
    
    # extract, pad with intercept dummy, and return
    covs <- terra::extract(covariates, 
                           coords[,c("x", "y")],
                           search_radius = buffer,
                           ID = FALSE) %>%
      dplyr::select(names(covariates)) %>%
      unlist()
    message("using buffer: haven't tested this before ...")
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



scale_years <- function(temporal_range){
  years <- do.call(seq, as.list(temporal_range))
  scaled_years <- scale(years)
  # fiddle with output of scale:
  scaled_years <- lapply(years, function(x){scaled_years[x - min(temporal_range) + 1]})
  names(scaled_years) <- years
  
  scaled_years
} 


# why did I set it up like that originally ??
extended_scaled_years <- function(scaled_years, extra_years){
  # this is the dumbest code I'm going to write for a while
  existing_years <- as.numeric(names(scaled_years))
  roc <- scaled_years[[2]] - scaled_years[[1]]
  intercept <- scaled_years[[1]] - roc * existing_years[1]
  to_add <- as.list(intercept + roc * extra_years)
  names(to_add) <- extra_years  
  
  scaled_years <- c(scaled_years, to_add)
}

# as in:
# extended_scaled_years(scaled_years, c(2000, 2027))




