
############ RUN JULIA ############

#' Run Julia version of CIRCUITSCAPE from R
#'
#' Execute julia CS from R
#'
#' @param JULIA_HOME Path to the folder containing the Julia binary (See Details). 
#' @param rast Accepts a RasterLayer object or SpatRaster. 
#' @param input_locs Provide a \code{\link[sp]{SpatialPoints}} object, \code{sf} point object, or \code{sf} polygon object.
#' @param regions If providing a point file for `input_locs`, specify short-circuit regions that overlap the points. Provide as \code{sf} polygon object. If `regions` are provided, your location points will not be used directly, but a random point from within each region will be used. (Default = NULL)
#' @param field (Default = NULL). Can be used to designate the identity of points or polygons when converting them to raster objects.If providing a polygon (SpatVector or sf POLYGON) this must be specified as the field name to get patch ID values from.
#' @param CurrentMap Logical. If TRUE, the cumulative current resistance map will be generated during the CS run (Default = TRUE)
#' @param cumulative_map_only Only write individual pairwise current maps, or also write each pairwise current map? (Default = TRUE)
#' @param VoltageMaps Logical. If TRUE, voltage maps will be produced.
#' @param scenario Specify 'pairwise' (Default), 'one-to-all', 'all-to-one', or 'advanced'. See CIRCUITSCAPE guide for details on usage.
#' @param Neighbor_Connect Specify 4 or 8 neighbor connection scheme (Default = 8)
#' @param output_dir Directory where CS results will be written. Only specify if Circuitscape current map outputs are requested. It is critical that there are NO SPACES in the directory, as this may cause the function to fail.
#' @param output_name Name to be applied to output files. Defaults to using the raster layer names.
#' @param output Specify either "matrix" or "raster". "matrix" will return the pairwise resistance matrix, while "raster" (Default) will return a \code{spatRaster} or \code{RasterLayer} object of the cumulative current map. The raster map can only be returned if \code{CurrentMap=TRUE}
#' @param pairs_to_include Default is NULL. If you wish to use the advanced CIRCUITSCAPE setting mode to include or exclude certain pairs of sample locations, provide a vector consisting of 1 (keep) or 0 (drop) for each pairwise observation (see example). This is an option if you do not want to assess all pairwise observations. 
#' @param variable_source When using 'one-to-all', 'all-to-one', or 'advanced' modes. Provide a two-column data frame or matrix with the first column being point/patch ID and the second column being the input source strength. Must be specified in conjunction with `variable_ground` in advanced mode. Can be specified by itself in 'one-to-all' and 'all-to-one'.
#' @param variable_ground Required when using 'advanced' mode. Provide a two-column data frame or matrix with the first column being point/patch ID and the second column being in the input source strength. Must be specified in conjunction with `variable_source` in advanced mode.
#' @param parallel (Logical; Default = FALSE) Do you want to run CIRCUITSCAPE in parallel? Only consider when dealing with rasters >2 million cells and if using Mac or Linux OS.
#' @param cores If `parallel = TRUE`, how many cores should be used for parallel processing? 
#' @param cholmod (Logical; Default = TRUE). Should the cholmod solver be used?  
#' @param precision (Logical; Default = FALSE). Should experimental single precision method be used? See details. 
#' @param is_resistance Default = TRUE. Is the landscape represented as a resistance (TRUE) or conductance (FALSE) surface?
#' @param focal_node_current_zero Set focal node current to 0 (Default = FALSE). May help with visualization of current flows in mapping.
#' @param max_map Should maximum current map be created and saved (Default = FALSE)?
#' @param remove_files Remove temporary files. Default is TRUE. 
#' @param silent Printing of output from Circuitscape will be suppressed unless an error occurs (Default = TRUE)
#' 
#' @return Function will save all output files created by CIRCUITSCAPE. Can choose to have a Full square distance matrix or Raster objects of the current maps returned to R. Raster object will be SpatRaster or RasterLayer, depending upon input format.
#' 
#' @details There is extensive documentation for Circuitscape here: https://docs.circuitscape.org/Circuitscape.jl/latest/ . Not all features and functionality have been built into this R function.
#' @examples 
#' ## To do
#' 
#' @usage run_cs(JULIA_HOME = NULL,
#' rast,
#' input_locs = NULL,
#' regions = NULL,
#' field = NULL,
#' CurrentMap = TRUE,
#' VoltageMaps = FALSE,
#' Neighbor_Connect = 8,
#' full_mat = FALSE,
#' output_dir = NULL,
#' output_name = NULL,
#' output = "raster",
#' parallel = FALSE, 
#' cores = NULL,
#' cholmod = TRUE,
#' precision = FALSE, 
#' is_resistance = TRUE,
#' focal_node_current_zero = FALSE,
#' max_map = FALSE,
#' remove_files = TRUE,
#' silent = TRUE)

#' @export
#' @author Bill Peterman <peterman.73@@osu.edu>
#' 
#' @examples  
#' ## Not run:
#' ## *** TO BE COMPLETED *** ##
#' 
#' ## End (Not run)
run_cs <- function(JULIA_HOME = NULL,
                   rast,
                   input_locs = NULL,
                   regions = NULL,
                   field = NULL,
                   CurrentMap = TRUE,
                   cumulative_map_only = TRUE,
                   VoltageMaps = FALSE,
                   scenario = 'pairwise',
                   Neighbor_Connect = 8,
                   output_dir = NULL,
                   output_name = NULL,
                   output = "raster",
                   pairs_to_include = NULL,
                   variable_source = NULL,
                   variable_ground = NULL,
                   parallel = FALSE,
                   cores = NULL,
                   cholmod = TRUE,
                   precision = FALSE,
                   is_resistance = TRUE,
                   focal_node_current_zero = FALSE,
                   max_map = FALSE,
                   remove_files = TRUE,
                   silent = TRUE) {
  
  wd <- getwd()
  
  rast_class <- class(rast)
  rast_crs <- try(crs(rast), silent = T)
  
  if(!exists('rast')){
    stop("\nMissing variable `rast`: Provide a raster surface!")
  }
  
  if(is.null(JULIA_HOME)){
    stop('\n`JULIA_HOME` must be specified')
  }
  
  if(is.null(input_locs)) {
    stop("\n`input_locs` must be specified!")
  }
  
  if(class(rast) == 'SpatRaster'){
    rast <- raster::raster(rast)
  } 
  asc.dir <- NULL
  
  if(isTRUE(cholmod)){
    solver <- "cholmod"
  } else {
    solver <- NULL
  }
  
  
  if(any(class(input_locs) == 'SpatialPointsDataFrame')){
    use_poly <- FALSE
    if(!is.null(field)){
      field_id <- input_locs[[field]]
      input_locs_rast <- try(raster::rasterize(x = input_locs, 
                                               raster = rast,
                                               field = field_id),
                             silent = T)
      if(class(input_locs_rast) == 'try-error'){
        # cat('\nSpecified `field` does not exist!\nManually setting patch values...check results carefully!!!\n')
        
        input_locs_rast <- raster::rasterize(coordinates(input_locs), 
                                             rast)
        
      }
    } else {
      input_sp <- as(input_locs, 'Spatial')
      input_locs_rast <- raster::rasterize(coordinates(input_locs), 
                                           rast)
    }
  } 
  
  
  if(any(class(input_locs) == 'SpatialPoints')){
    use_poly <- FALSE
    input_locs_rast <- raster::rasterize(input_locs, 
                                         rast)
  } 
  
  if(any(class(input_locs) == 'sf')){
    if(st_geometry_type(input_locs)[[1]] == 'POLYGON' |
       st_geometry_type(input_locs)[[1]] == 'MULTIPOLYGON'){
      use_poly <- TRUE
      if(!is.null(field)){
        region_rast <- input_locs_rast <- try(fasterize::fasterize(input_locs, 
                                                                   raster = rast,
                                                                   field = field), 
                                              silent = T)
        if(class(input_locs_rast) == 'try-error'){
          cat('\nSpecified `field` does not exist!\nManually setting patch values...check results carefully!!!\n')
          
          input_locs$id <- 1:length(input_locs$geometry)
          field_id <- 'id'
          
          region_rast <- input_locs_rast <- try(fasterize::fasterize(input_locs, 
                                                                     raster = rast,
                                                                     field = field_id),
                                                silent = T)
        }
      } else {
        cat('\nManually setting patch values...check results carefully!!!\n')
        
        input_locs$id <- 1:length(input_locs$geometry)
        field_id <- 'id'
        region_rast <- input_locs_rast <- try(fasterize::fasterize(input_locs, 
                                                                   raster = rast,
                                                                   field = field_id))
      }
    } else { ## sf Points
      use_poly <- FALSE
      if(!is.null(field)){
        field_id <- input_locs[[field]]
        input_sp <- as_Spatial(input_locs)
        input_sp <- as(input_sp, 'SpatialPoints')
        
        input_locs_rast <- try(raster::rasterize(input_sp, 
                                                 rast,
                                                 field = field_id))
        if(class(input_locs_rast) == 'try-error'){
          cat('\nIssues with specified `field`!\nManually setting point values...check results carefully!!!\n')
          
          input_locs_rast <- raster::rasterize(input_sp, 
                                               rast)
        }
      } else {
        input_sp <- as_Spatial(input_locs)
        input_sp <- as(input_sp, 'SpatialPoints')
        
        input_locs_rast <- raster::rasterize(input_sp, 
                                             rast)
      } 
    } 
  } ## End sf statement
  
  
  if(any(class(input_locs) == 'SpatVector')){
    if(terra::geomtype(input_locs) == 'polygons'){
      use_poly <- TRUE
      ## Convert to raster
      input_sf <- sf::st_as_sf(input_locs)
      if(!is.null(field)){
        
        region_rast <- input_locs_rast <- try(fasterize::fasterize(input_sf, 
                                                                   raster = rast,
                                                                   field = field),
                                              silent = T)
        
        if(class(input_locs_rast) == 'try-error'){
          cat('\nSpecified `field` does not exist!\nManually setting patch values...check results carefully!!!')
          
          input_sf$id <- 1:length(input_sf$geometry)
          field_id <- 'id'
          
          region_rast <- input_locs_rast <- try(fasterize::fasterize(input_sf, 
                                                                     raster = rast,
                                                                     field = field_id))
        } 
      } else {
        input_sf$id <- 1:length(input_sf$geometry)
        field_id <- 'id'
        
        region_rast <- input_locs_rast <- try(fasterize::fasterize(input_sf, 
                                                                   raster = rast,
                                                                   field = field_id))
      }
    } else {
      use_poly <- FALSE
      input_sp <- as(input_locs, 'Spatial')
      input_sp <- as(input_sp, 'SpatialPoints')
      input_locs_rast <- raster::rasterize(input_sp, 
                                           rast)
    }
  }## End SpatVector  
  
  ## Polygon region
  if(!is.null(regions)){
    use_poly <- TRUE
    # sf_class <- grep('sf', class(regions))
    if(isFALSE('sf' %in% class(regions))){
      stop('`regions` variable is defined but is not a simple features polygon of class `sf`. Try again!\n')
    }
    region_pts <- st_sample(regions, rep(1, length(regions$geometry)))
    
    region_sp <- as_Spatial(region_pts)
    region_sp <- as(region_sp, 'SpatialPoints')
    input_locs_rast <- raster::rasterize(region_sp, 
                                         rast)
    
    if(!is.null(field)){
      region_rast <- try(fasterize::fasterize(regions, 
                                              raster = rast,
                                              field = field), 
                         silent = T)
      if(class(region_rast) == 'try-error'){
        cat('\nSpecified `field` does not exist!\nManually setting patch values...check results carefully!!!\n')
        
        regions$id <- 1:length(regions$geometry)
        field_id <- 'id'
        
        region_rast <- try(fasterize::fasterize(regions, 
                                                raster = rast,
                                                field = field_id),
                           silent = T)
      }
    } else {
      regions$id <- 1:length(regions$geometry)
      field_id <- 'id'
      
      region_rast <- try(fasterize::fasterize(regions, 
                                              raster = rast,
                                              field = field_id),
                         silent = T)
    }
  }
  
  loc_rast_made <- try(exists('input_locs_rast'), silent = T)
  if(class(loc_rast_made) == 'try-error'){
    stop('\nFailed to create input location raster!\n')
  }
  
  if (class(rast)[1] == 'SpatRaster') {
    rast <- raster(rast)
    asc.dir <- NULL
  } else {
    rast <- rast
    asc.dir <- NULL
  }
  
  
  if(is.null(output_dir)) {   
    if(Sys.info()[['sysname']] == "Windows") {
      output_dir <- paste0(tempdir(),"/")
      output_dir <- normalizePath(output_dir, winslash = "/")
      output_dir <- paste0(output_dir,'/')
    } else {
      output_dir <- paste0(normalizePath(tempdir()),"/")
    }
  } else {
    output_dir <- paste0(normalizePath(output_dir),'/')
    
  }
  
  temp_rast <- rm.rast <- tempfile(pattern = "raster_", 
                                   tmpdir = tempdir(),
                                   fileext = ".asc") 
  
  tmp.name <- basename(temp_rast) %>% strsplit(., '.asc') %>% unlist()
  
  temp_loc <- rm.loc <- tempfile(pattern = "loc_", 
                                 tmpdir = tempdir(),
                                 fileext = ".asc") 
  
  loc.name <- basename(temp_loc) %>% strsplit(., '.asc') %>% unlist()
  
  temp_poly <- rm.poly <- tempfile(pattern = "poly_", 
                                   tmpdir = tempdir(),
                                   fileext = ".asc") 
  
  poly.name <- basename(temp_poly) %>% strsplit(., '.asc') %>% unlist()
  
  temp_ini <- rm.ini <- tempfile(pattern = "cs_run_", 
                                 tmpdir = tempdir(),
                                 fileext = ".ini") 
  
  
  # Write Rasters for CS ----------------------------------------------------
  ## Conductance / Resistance
  writeRaster(
    x = rast,
    filename = temp_rast,
    overwrite = TRUE
  )
  
  ## Point raster
  writeRaster(
    x = input_locs_rast,
    filename = temp_loc,
    overwrite = TRUE
  )  
  
  ## Region raster
  if(isTRUE(use_poly)){
    writeRaster(
      x = region_rast,
      filename = temp_poly,
      overwrite = TRUE
    )
  }
  
  
  # Run Write ini -------------------------------------------------------------
  if(output == "raster"){
    if(isFALSE(CurrentMap)){
      cat("\nRaster output requested; Setting `CurrentMap = TRUE`")
      CurrentMap <- TRUE
    }
  }
  
  
  if(!is.null(variable_source)){
    if(isTRUE(use_poly)){
      var_source_rast <- reclassify(input_locs_rast, as.matrix(variable_source))
      var_source_df <- rasterToPoints(var_source_rast)
      variable_source <- var_source_df[,c(3,1,2)]
    }
    if(is.null(variable_ground)){
      VARSOURCE_T <- TRUE
    }
    if(is.null(variable_ground) & scenario == 'advanced'){
      stop('Must also specify `variable_ground` file with advanced mode.')
    }
  } else {
    VARSOURCE_T <- FALSE
  }
  
  
  if(isTRUE(focal_node_current_zero)){
    FOCALNODE <- "set_focal_node_currents_to_zero = True"
  } else {
    FOCALNODE <- "set_focal_node_currents_to_zero = False"
  }
  
  
  if(isTRUE(max_map)){
    MAXMAP <- "write_max_cur_maps = True"
  } else {
    MAXMAP <- "write_max_cur_maps = False"
  }
  
  
  if(is.null(variable_source)){
    VARSOURCE <- 'use_variable_source_strengths = FALSE'
    SOURCEFILE <- 'variable_source_file = (Browse for file)'
    GRNDFILE <- 'ground_file = (Browse for file)'
    VARSOURCEFILE <- 'variable_source_file = (Browse for file)'
  }
  
  if(isTRUE(VARSOURCE_T)){
    temp_src <- rm.src <- tempfile(pattern = "source_", 
                                   tmpdir = tempdir(),
                                   fileext = ".txt") 
   
    write.table(variable_source,
                temp_src,
                sep = '\t',
                row.names = F, col.names = F)
    
    VARSOURCE <- 'use_variable_source_strengths = True'
    VARSOURCEFILE <- paste0('variable_source_file = ', temp_src)
    GRNDFILE <- 'ground_file = (Browse for file)'
    SOURCEFILE <- 'variable_source_file = (Browse for file)'
    
  } 
  
  if(!is.null(variable_ground)){
    if(ncol(variable_source != 2)){
      stop('Must specify a two-column matrix or data frame for variable sources.')
    }
    temp_src <- rm.src <- tempfile(pattern = "source_", 
                                   tmpdir = tempdir(),
                                   fileext = ".txt") 
    
    write.table(variable_source,
                temp_src,
                sep = '\t',
                row.names = F, col.names = F)
    
    SOURCEFILE <- paste0('variable_source_file = ', temp_src)
    
    if(ncol(variable_ground != 2)){
      stop('Must specify a two-column matrix or data frame for variable grounds.')
    }
    temp_grnd <- rm.grnd <- tempfile(pattern = "ground_", 
                                     tmpdir = tempdir(),
                                     fileext = ".txt") 
    
    write.table(variable_ground,
                temp_grnd,
                sep = '\t',
                row.names = F, col.names = F)
    
    GRNDFILE <- paste0('variable_ground_file = ', temp_grnd)
  }
  
  if(scenario == 'pairwise'){
    SCENARIO <- "scenario = pairwise"
    
  } else if(scenario == 'one-to-all'){
    SCENARIO <- "scenario = one-to-all"
    
  } else if(scenario == 'all-to-one'){
    SCENARIO <- "scenario = all-to-one"
    
  } else if (scenario == 'advanced'){
    SCENARIO <- "scenario = advanced"
    
  } else {
    SCENARIO <- "scenario = pairwise"
    cat('\nDefaulting to "pairwise" mode!\n')
  }
  
  if (CurrentMap == FALSE) {
    File.name <- names(rast)
    MAP = "write_cum_cur_map_only = False"
    CURRENT.MAP = "write_cur_maps = 0"
    
  } else {
    if(is.null(output_dir)){
      File.name <- tmp.name
      # tmp.name <- File.name
    } else {
      if(!is.null(output_name)){
        File.name <- output_name
        tmp.name <- File.name 
      } else {
        File.name <- names(rast)
        tmp.name <- File.name 
      }
    }
    if(isTRUE(cumulative_map_only)){
      MAP = "write_cum_cur_map_only = True"
    } else {
      MAP = "write_cum_cur_map_only = False"
    }
    CURRENT.MAP = "write_cur_maps = True"
  }
  
  if(VoltageMaps == FALSE){
    VOLTAGE <- "write_volt_maps = False"
  } else {
    VOLTAGE <- "write_volt_maps = True"
  }
  
  # Modify and write Circuitscape.ini file
  if(Neighbor_Connect == 4) {
    connect <- "True"
  } else {
    connect <- "False"
  }
  
  if(isTRUE(use_poly)){
    POLYGON_FILE <- paste0("polygon_file = ", temp_poly)
  } else {
    POLYGON_FILE <- paste0("polygon_file = (Browse for polygon file)") # temp_loc
  }
  
  ## Turned off because slower to calculate specific pairs in circuitscape
  # if(is.null(jl.inputs$pairs_to_include)) {
  PAIRS_TO_INCLUDE <-
    paste0("included_pairs_file = (Browse for a file with pairs to include or exclude)")
  PAIRS <- paste0("use_included_pairs = False")
  
  suppressWarnings(OUT <- paste0("output_file = ", normalizePath(paste0(output_dir, tmp.name, ""))))
  suppressWarnings(BATCH <- normalizePath(paste0(output_dir, tmp.name, ".ini")))
  
  write_ini(
    BATCH = temp_ini,
    OUT = OUT,
    HABITAT = paste0("habitat_file = ", temp_rast),
    LOCATION.FILE = paste0("point_file = ", temp_loc),
    CONNECTION = paste0("connect_four_neighbors_only = ", connect),
    GRNDFILE = GRNDFILE,
    VARSOURCE = VARSOURCE,
    SOURCEFILE = SOURCEFILE,
    VARSOURCEFILE = VARSOURCEFILE,
    MAP = MAP,
    VOLTAGE = VOLTAGE,
    CURRENT.MAP = CURRENT.MAP,
    PAIRS_TO_INCLUDE = PAIRS_TO_INCLUDE,
    PAIRS = PAIRS,
    PARALLELIZE = parallel,
    CORES = cores,
    solver = solver,
    precision = precision,
    POLYGON_FILE = POLYGON_FILE,
    USE_POLYGONS = use_poly,
    SCENARIO = SCENARIO,
    FOCALNODE = FOCALNODE,
    MAXMAP = MAXMAP,
    silent = silent,
    is_resistance = is_resistance
  )
  
  
  # rm(list = c('input_locs', 'input_locs_rast', 'rast'))
  
  # Run CIRCUITSCAPE.jl -----------------------------------------------------
  
  # out <- try(JuliaCall::julia_call('compute', normalizePath(
  #   paste0(output_dir, tmp.name, ".ini")
  # ))[-1,-1], silent = T)
  
  JuliaCall::julia_setup(JULIA_HOME = JULIA_HOME)
  try(JuliaCall::julia_library("Circuitscape"), silent = T)
  
  out <- try(JuliaCall::julia_call('compute', temp_ini)[-1,-1], silent = T)
  
  if(any(class(out) == 'try-error')){
    cat(readChar(temp_ini, file.info(temp_ini)$size))
    cat("\n\n")
    if(isTRUE(remove_files)){
      try(unlink(rm.loc, force = TRUE), silent = T)
      try(unlink(rm.rast, force = TRUE), silent = T)
      try(unlink(rm.ini, force = TRUE), silent = T)
      try(unlink(rm.src, force = TRUE), silent = T)
      try(unlink(rm.grnd, force = TRUE), silent = T)
      try(unlink(rm.poly, force = TRUE), silent = T)
      try(unlink(OUT, force = TRUE), silent = T)
    }
    
    stop('\nCIRCUITSCAPE encountered an error!\nCheck inputs printed above carefully\n\n')
  }
  
  if(wd != getwd()) {
    setwd(wd)
  }
  
  if (any(out[lower.tri(out)] == 0)) {
    message(
      "\n Zero values were generated by CIRCUITSCAPE \n Check point file to see if multiple points share the same raster cell!"
    )
  }
  
  
  # Delete files ------------------------------------------------------------
  if(isTRUE(remove_files)){
    try(unlink(rm.loc, force = TRUE), silent = T)
    try(unlink(rm.rast, force = TRUE), silent = T)
    try(unlink(rm.ini, force = TRUE), silent = T)
    try(unlink(rm.src, force = TRUE), silent = T)
    try(unlink(rm.grnd, force = TRUE), silent = T)
    try(unlink(rm.poly, force = TRUE), silent = T)
    try(unlink(OUT, force = TRUE), silent = T)
  }
  
  
  # Return outputs -----------------------------------------------------------
  
  if (output == "raster") {
    rast_cs <- try(raster(normalizePath(paste0(output_dir, '/', File.name, "_cum_curmap.asc"))),
                   silent = T)
    
    names(rast_cs) <- File.name
    
    if(class(rast_crs) != 'character'){
      crs(rast_cs) <- rast_crs
    }
    
    if(isTRUE(VoltageMaps)){
      volt_list <- list.files(path = normalizePath(paste0(output_dir, '/')),
                              pattern = paste0(File.name, "_voltmap"),
                              all.files = T,
                              full.names = T)
      volt_stack <- stack(sapply(volt_list, raster))
      volt_names <- basename(unlist(volt_list))
      names(volt_stack) <- sub('.asc', "", volt_names)
      
      # if(rast_class == 'SpatRaster'){
      #   volt_stack <- suppressWarnings(terra::rast(volt_stack))
      # }
      rast_cs <- stack(rast_cs, volt_stack)
    }
    
    if(isFALSE(cumulative_map_only)){
      cur_list <- list.files(path = normalizePath(paste0(output_dir, '/')),
                             pattern = paste0(File.name, "_curmap_"),
                             all.files = T,
                             full.names = T)
      cur_stack <- stack(sapply(cur_list, raster))
      cur_names <- basename(unlist(cur_list))
      names(cur_stack) <- sub('.asc', "", cur_names)
      
      rast_cs <- stack(rast_cs, cur_stack)
      
      # if(rast_class == 'SpatRaster'){
      #   cur_stack <- suppressWarnings(terra::rast(cur_stack))
      # }
    }
    
    if(rast_class == 'SpatRaster'){
      rast_cs <- suppressWarnings(terra::rast(rast_cs))
    }
    
    return(rast_cs)
  }
  
  
  return(out)
} # End Julia function

write_ini <- function(BATCH,
                      OUT,
                      HABITAT,
                      LOCATION.FILE,
                      CONNECTION = 'connect_four_neighbors_only = FALSE',
                      CURRENT.MAP = "write_cur_maps = 0",
                      CUMONLY = 'write_cum_cur_map_only = True',
                      VOLTAGE = "write_volt_maps = False",
                      MAP = "write_cum_cur_map_only = False",
                      GRNDFILE,
                      VARSOURCE,
                      SOURCEFILE,
                      VARSOURCEFILE,
                      PARALLELIZE = FALSE,
                      CORES = NULL,
                      solver = 'cholmod',
                      precision = FALSE,
                      is_resistance = TRUE,
                      PAIRS_TO_INCLUDE = "included_pairs_file = (Browse for a file with pairs to include or exclude)",
                      POLYGON_FILE,
                      USE_POLYGONS = FALSE,
                      PAIRS = "use_included_pairs = False",
                      SCENARIO = "scenario = pairwise",
                      FOCALNODE = "set_focal_node_currents_to_zero = False",
                      MAXMAP = "write_max_cur_maps = False",
                      silent = NULL){
  if(PARALLELIZE == TRUE && !is.null(CORES)) {
    PARALLELIZE <- "parallelize = True"
    CORES <- paste0("max_parallel = ", CORES)
  } else {
    PARALLELIZE = "parallelize = False"
    CORES <- "max_parallel = 0"
  }
  
  if(!is.null(solver)) {
    solver <- "solver = cholmod"
  } else {
    solver <- "solver = cg+amg"
  }
  
  if(isTRUE(silent)) {
    LOG <- "log_level = critical"
  } else {
    LOG <- "log_level = INFO"
  }
  
  if(isTRUE(precision == 'single')) {
    precision <- "precision = single"
  } else {
    precision <- "precision = None"
  } 
  
  if(isTRUE(is_resistance)) {
    RESISTANCE <- "connect_using_avg_resistances = False"
    HABITAT_RES <- "habitat_map_is_resistances = True"
  } else {
    RESISTANCE <- "connect_using_avg_resistances = False"
    HABITAT_RES <- "habitat_map_is_resistances = False"
  }
  
  sink(BATCH)
  cat("[Options for advanced mode]")
  cat("\n")
  cat("ground_file_is_resistances = True")
  cat("\n")
  cat("remove_src_or_gnd = rmvsrc")
  cat("\n")
  cat(GRNDFILE)
  cat("\n")
  cat("use_unit_currents = False")
  cat("\n")
  cat(SOURCEFILE)
  cat("\n")
  cat("use_direct_grounds = False")
  cat("\n")
  cat("\n")
  cat("[Mask file]")
  cat("\n")
  cat("mask_file = (Browse for a raster mask file)")
  cat("\n")
  cat("use_mask = False")
  cat("\n")
  cat("\n")
  cat("[Calculation options]")
  cat("\n")
  cat("low_memory_mode = False")
  cat("\n")
  cat(PARALLELIZE)
  cat("\n")
  cat(CORES)
  cat("\n")
  cat(solver)
  cat("\n")
  cat("print_timings = False")
  cat("\n")
  cat("preemptive_memory_release = False")
  cat("\n")
  cat("print_rusages = False")
  cat("\n")
  cat("\n")
  cat("[Short circuit regions (aka polygons)]")
  cat("\n")
  cat(POLYGON_FILE)
  cat("\n")
  cat(paste0('use_polygons = ', USE_POLYGONS))
  # cat("polygon_file = (Browse for a short-circuit region file)")
  # cat("\n")
  # cat("use_polygons = False")
  cat("\n")
  cat("\n")
  cat("[Options for one-to-all and all-to-one modes]")
  cat("\n")
  cat(VARSOURCE)
  cat("\n")
  cat(VARSOURCEFILE)
  cat("\n")
  cat("\n")
  cat("[Output options]")
  cat("\n")
  cat("set_null_currents_to_nodata = True")
  cat("\n")
  cat(FOCALNODE)
  cat("\n")
  cat("set_null_voltages_to_nodata = True")
  cat("\n")
  cat("compress_grids = False")
  cat("\n")
  cat(VOLTAGE)
  cat("\n")
  cat(CURRENT.MAP)
  cat("\n")
  cat(OUT)
  cat("\n")
  cat(MAP)
  cat("\n")
  cat("log_transform_maps = False")
  cat("\n")
  cat(MAXMAP)
  cat("\n")
  cat("\n")
  cat("[Options for reclassification of habitat data]")
  cat("\n")
  cat("reclass_file = (Browse for file with reclassification data)")
  cat("\n")
  cat("use_reclass_table = False")
  cat("\n")
  cat("\n")
  cat("[Logging Options]")
  cat("\n")
  cat(LOG)
  cat("\n")
  cat("log_file = None")
  cat("\n")
  cat("profiler_log_file = None")
  cat("\n")
  cat("screenprint_log = False")
  cat("\n")
  cat("\n")
  cat("[Options for pairwise and one-to-all and all-to-one modes]")
  cat("\n")
  cat(PAIRS_TO_INCLUDE)
  cat("\n")
  cat(PAIRS)
  cat("\n")
  cat(LOCATION.FILE)
  cat("\n")
  cat("\n")
  cat("[Connection scheme for raster habitat data]")
  cat("\n")
  cat(RESISTANCE)
  cat("\n")
  cat(CONNECTION)
  cat("\n")
  cat("\n")
  cat("[Habitat raster or graph]")
  cat("\n")
  cat(HABITAT_RES)
  cat("\n")
  cat(HABITAT)
  cat("\n")
  cat("\n")
  cat("[Circuitscape mode]")
  cat("\n")
  cat("data_type = raster")
  cat("\n")
  cat(SCENARIO)
  cat("\n")
  cat(precision)
  sink()
}
