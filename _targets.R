#Small changes to implement before next run
#L22 functions: capitalize S in spain
#L56 and L158 functions: catch metacols regardless of capitalization
#Standardize country names across all input datasets
#Standardize metacols across all input datasets
#site names in "data/wp1/Results_present_period_final/data/Genal/Genal_sampling_sites_ReachIDs.csv" are incorrect
#Need to know how sites were snapped to network "C:\DRYvER_wp2\WP2 - Predicting biodiversity changes in DRNs\Coordinates\Shapefiles with sites moved to the river network\Croatia_near_coords.shp"
#3s flow acc for Europe: https://data.hydrosheds.org/file/hydrosheds-v1-acc/eu_acc_3s.zip

#Make sure that biological data are standardized by area to get densities
#Get ancillary catchment data

library(rprojroot)
rootdir <- rprojroot::find_root(has_dir('R'))
setwd(rootdir)

source('src/packages.R')
source("src/functions.R")

hydromod_present_dir <- 'data'
resdir <- 'results'

tar_option_set(format = "qs")

#--------------------------  Define targets plan -------------------------------
list(
  #------------------------------- Define paths ----------------------------------
  #Path to modeled hydrological data
  tar_target(
    hydromod_paths_dt,
    define_hydromod_paths(in_hydromod_dir = hydromod_present_dir)
  ),
  
  #------------------------------- Read in data ----------------------------------
  #Read reach data
  tar_target(
    reaches_dt,
    lapply(hydromod_paths_dt$country, function(in_country) {
      hydromod_paths_dt[country==in_country, 
                        fread(reaches_path, sep="\t", skip=1, 
                              header=T, drop=c(17,18)) %>%
                          .[-c(1,2,3),] %>%
                          .[ID != '# end of reach.par',] %>%
                          .[, lapply(.SD, as.numeric)] %>%
                          setnames(gsub('[-]', '_', names(.)))
      ] %>%
        .[, country := in_country]
    }) %>% rbindlist(use.names=T, fill=T) 
    #Reach that flows to 9999 is outlet in Czech, Finland, France, and Hungary
    #Reach that flows to 0 is outlet in Spain and Croatia
  ),

  #Subset river network shapefiles to only keep sections within which there
  #are sampling site-reaches (based on sub-catchment file given by countries)
  tar_target(
    network_sub_gpkg_list,
    subset_network(in_hydromod_paths_dt = hydromod_paths_dt,
                   out_dir = file.path('results', 'gis'),
                   overwrite = T)
  )
  ,
  
  #Clean networks
  tar_target(
    network_clean_gpkg_list,
    lapply(names(network_sub_gpkg_list), function(in_country) {
      #Set size of simplifying radius to remove loops. See function
      clustering_dist_list <- list(Croatia=50, Czech=60, Finland=50, 
                                   France=50, Hungary=50, Spain=40)
      
      out_net_path <- clean_network(
        rivnet_path = network_sub_gpkg_list[[in_country]],
        idcol = 'cat',
        node_clustering_dist = clustering_dist_list[[in_country]],
        min_segment_length = 20,
        outdir = file.path(resdir, 'gis'),
        save_gpkg = TRUE, 
        return_path = TRUE)
      
      #Manual corrections
      if (in_country == 'Czech') {
        clean_net <- st_read(out_net_path) %>%
          .[!(.[['cat']] %in% c(40112, 40106, 40478)),]
        st_write(clean_net, out_net_path, append=F)
      } else if (in_country == 'Croatia') {
        clean_net <- st_read(out_net_path) %>%
          .[!(.[['UID']] %in% c(1102, 2988, 457, 458, 462, 464)),] 
        #Change startpoint of disconnected line where there is a site 
        #(startpoint, because the line dir is reversed)
        line_to_edit <- clean_net[clean_net$UID==651,]$geom
        clean_net[clean_net$UID==651,]$geom <- st_sfc(st_linestring(
          rbind(
            c(X=596968.793, Y=4899716.907),
            st_coordinates(line_to_edit)[-1, c('X', 'Y')]
          )))
        
        st_write(clean_net, out_net_path, append=F)
      }
      return(out_net_path)
    }) %>% setNames(names(network_sub_gpkg_list)) 
  )
  ,
  
  tar_target(
    network_directed_gpkg_list,
    lapply(names(network_clean_gpkg_list), function(in_country) {
      #Set size of simplifying radius to remove loops. See function
      outlet_uid_list <- list(Croatia = 463, Czech = 4, Finland = 682,
                              France = 1, Hungary = 5, Spain = 86)
      
      out_net_path <- direct_network(
        rivnet_path = network_clean_gpkg_list[[in_country]],
        idcol = 'UID',
        outletid = outlet_uid_list[[in_country]],
        outdir = file.path(resdir, 'gis'), 
        save_gpkg = TRUE) 
      
      return(out_net_path)
    }) %>% setNames(names(network_clean_gpkg_list))
  )
  ,
  
  
  tar_target(
    network_nocomplexconf_gpkg_list,
    lapply(names(network_directed_gpkg_list), function(in_country) {
      fix_complex_confluences(
        rivnet_path = network_directed_gpkg_list[[in_country]], 
        max_node_shift = 5,
        out_path= file.path(resdir, 'gis', 
                            paste0(tolower(in_country)
                                   , '_river_network_nocomplexconf',
                                   format(Sys.time(), "%Y%m%d"), '.gpkg'
                            ))
      )
    }) %>% setNames(names(network_directed_gpkg_list))
  )
  ,
  
  #Compute strahler stream order
  tar_target(
    network_strahler,
    lapply(names(network_nocomplexconf_gpkg_list), function(in_country) {
      assign_strahler_order(
        in_rivnet = network_nocomplexconf_gpkg_list[[in_country]], 
        idcol = 'UID')
    }) %>% setNames(names(network_nocomplexconf_gpkg_list))
  )
  ,
  
  #Re-assign correct IDs to match with hydrological data
  tar_target(
    network_ssnready_gpkg_list,
    lapply(names(network_nocomplexconf_gpkg_list), function(in_country) {
      out_net_path <- reassign_netids(
        rivnet_path = network_nocomplexconf_gpkg_list[[in_country]], 
        strahler_dt = network_strahler[[in_country]], 
        in_reaches_hydromod_dt = reaches_dt[country==in_country,], 
        outdir = file.path(resdir, 'gis'),
        country = in_country
      )
      return(out_net_path)
    }) %>% setNames(names(network_nocomplexconf_gpkg_list))
  )
  ,
  
  #Copy gpkg to shapefiles for sharing
  tar_target(
    network_ssnready_shp_list,
    lapply(network_ssnready_gpkg_list, function(path) {
      old_cols <- c('UID', 'strahler', 'length_uid', 'cat_cor', 'from', 'to',
                    'to_reach_shpcor', 'to_reach_hydromod', 
                    'hydromod_shpcor_match', 'geom')
      new_cols <- c('UID', 'strahler', 'length_m', 'cat', 'from', 'to',
                    'to_cat_shp', 'to_cat_mod', 'mod_match', 'geom')
      lyr <- st_read(path)[,old_cols]
      names(lyr) <- new_cols
      lyr$UID <- seq_along(lyr$UID)
      lyr$strahler <- as.integer(lyr$strahler)
      lyr$cat <- as.integer(lyr$cat)
      lyr$to_cat_shp <- as.integer(lyr$to_cat_shp)
      lyr$to_cat_mod <- as.integer(lyr$to_cat_mod)
      write_sf(lyr, paste0(tools::file_path_sans_ext(path), '.shp'))
    })
  )



