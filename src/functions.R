#-------------- utility functions ----------------------------------------------
#------ dist_proj  -------------------------------------------------
# Define standard two-point equidistance projection for a given bounding box
#https://gis.stackexchange.com/questions/313721/automatically-get-an-adequate-projection-system-based-on-a-bounding-box
## distance projection (tpeqd - two-point equidistant) with projection parameters 
## derived from feature extent
dist_proj <- function(x) {
  bb <- sf::st_bbox(x)
  paste0("+proj=tpeqd +lat_1=",
         bb[2],
         " +lon_1=", 
         bb[1],
         " +lat_2=",
         bb[4], 
         " +lon_2=", 
         bb[3],
         " +x_0=0",
         " +y_0=0",
         " +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
}

#------ snap sites to nearest segment ---------------------------------------------
snap_points_inner <- function(in_pts,
                              in_target,
                              sites_idcol,
                              attri_to_join=NULL
) {
  #Snap points (fastest custom way in R, it seems):
  #first computing a line between site and snapping place on nearest segment
  sitesnap_l <- terra::nearest(in_pts, in_target, centroids = F, lines = T)
  values(sitesnap_l) <- values(in_pts)
  sitesnap_l$snap_dist_m <- perim(sitesnap_l)
  
  #convert the line to a point (the line's end point)
  sitesnap_p <- terra::as.points(sitesnap_l) %>%
    .[duplicated(values(.)[, sites_idcol]),]
  
  #Join attributes of nearest line to that point
  if (!is.null(attri_to_join)) {
    if (attri_to_join == 'all') { 
      sitesnap_p[, names(in_target)] <- terra::nearby(
        sitesnap_p, in_target, k=1, centroids=FALSE)[,'k1'] %>% #Could grab the nth nearest or place a distance limit
        as.data.frame(in_target)[.,] 
    } else {
      sitesnap_p[, attri_to_join] <- terra::nearby(
        sitesnap_p, in_target, k=1, centroids=FALSE)[,'k1'] %>%
        as.data.frame(in_target)[., attri_to_join] 
    }
  }
  
  return(sitesnap_p)
}

#------ split_sp_line----------------------------------------------------------
# Original author: Miguel Porto
# #From https://github.com/miguel-porto/fix-streams
# splits a Lines object (or coordinate matrix) in n segments of length length
#(starting in the begining) plus the remaining segment (what is left)

split_sp_line <- function(line, n, length, debug = F) {
  # splits a Lines object (or coordinate matrix) in n segments of length length (starting in the begining) plus the remaining segment (what is left)
  if (debug) plot(line)
  coo <- sp::coordinates(line)
  
  if (inherits(coo, "list")) {
    if (length(coo) > 1) {
      stop("Multiple lines not allowed")
    } else {
      coo <- coo[[1]]
      if (!inherits(coo, "matrix")) {
        if (!inherits(coo, "list")) stop("Invalid line object")
        if (length(coo) > 1) stop("Multiple lines not allowed")
        coo <- coo[[1]]
      }
    }
  } else {
    if (!inherits(coo, "matrix")) stop("Invalid line object")
  }
  
  pieces <- list()
  accum <- 0
  i <- 1
  remainder <- 0
  newcoords <- matrix(nc = 2, nr = 0)
  
  repeat {
    newcoords <- rbind(newcoords, coo[i, ])
    v <- c(coo[i + 1, ] - coo[i, ])
    hyp <- sqrt(sum(v^2))
    accum <- accum + hyp
    if (accum > length && length(pieces) < n) { # cut in this segment
      oript <- coo[i, ]
      # 			accum=hyp
      repeat {
        newpt <- oript + v / hyp * (length - remainder)
        newcoords <- rbind(newcoords, newpt)
        pieces <- c(pieces, list(newcoords))
        newcoords <- matrix(newpt, nr = 1)
        remainder <- 0
        if (debug) points(newpt[1], newpt[2], pch = 19)
        accum <- accum - length
        if (accum < length || length(pieces) >= n) break
        oript <- newpt
      }
      remainder <- accum
    } else {
      remainder <- accum
    }
    
    i <- i + 1
    if (i >= dim(coo)[1]) {
      newcoords <- rbind(newcoords, coo[dim(coo)[1], ])
      pieces <- c(pieces, list(newcoords))
      break
    }
  }
  
  if (debug) {
    for (i in 1:length(pieces)) {
      pieces[[i]][, 1] <- pieces[[i]][, 1] - 20 * i
      lines(pieces[[i]], col = "red")
      points(pieces[[i]][1, 1], pieces[[i]][1, 2])
    }
  }
  
  return(pieces)
}

#------ fix_confluences_inner ---------------------------------------------------------
# Fixes complex confluences in stream networks
# Original author: Miguel Porto (only cosmetic changes were performed)
# From https://github.com/miguel-porto/fix-streams/blob/master/fix-streams.r
# All nodes which have >2 streams flowing to it are corrected. The outermost streams' end vertices
# are adjusted by "step" meters along the downgoing stream. New nodes are created at suitable places,
# and existing lines suitably split.
#*** Requires a SpatialLinesDataFrame with the proper FROM_NODE and TO_NODE fields. 
# The network is assumed to be correct in all other aspects, there is no error checking.
###### USAGE EXAMPLE
# rios=readOGR("streams_Pt.shp","streams_Pt")
# correctedshp=fix.streams(rios,step=10)
# writeOGR(correctedshp,"streams_corrected.shp","streams_corrected","ESRI Shapefile")

fix_confluences_inner <- function(shp, from = "FROM_NODE", to = "TO_NODE", 
                                  step = 10, fields_to_keep=NULL) {
  # step is the desired length (in map units) by which the river sinks are adjusted (separated) downstream.
  pieces <- list()
  probrivers <- list()
  removeindexes <- integer(0)
  CRS <- shp@proj4string
  
  # find multiple confluence nodes
  nv <- table(shp@data[, to])
  mc <- as.numeric(names(nv[nv > 2])) # these are the nodes with >2 rivers flowing to	them
  maxnode <- max(c(shp@data[, to], 
                   shp@data[, from])) # max node ID, for creating new nodes
  
  if (length(mc) > 0) {
    cat(length(mc), "complex confluences found.\n")
    cat("Cutting lines and tweaking vertices...\n")
    flush.console()
    
    for (i in mc) { # for each problematic node
      # msrc=shp[shp@data[,to] %in% i,]	# get source rivers flowing to it
      privers <- which(shp@data[, to] == i) # get source rivers flowing to it (problematic rivers)
      sinkriver <- which(shp@data[, from] == i) # get sink river
      msrc <- shp[privers, ]
      mto <- shp[sinkriver, ]
      delta <- sum(shp@data[, to] %in% i) - 2 # how many problematic rivers flow to there? leave only two, the others correct
      newnodes <- c(mto@data[, from], (maxnode + 1):(maxnode + delta), mto@data[, to]) # the IDs of the nodes that will be created (the first remains the same for the two "good" rivers)
      
      coo <- sp::coordinates(mto)[[1]][[1]] # coordinates of the sink river
      
      # order the rivers by their angle, so that the outermost rivers are first adjusted (alternating the side)
      v1 <- coo[1, ] - coo[2, ]
      v2 <- matrix(nc = 2, nr = length(msrc))
      for (j in 1:length(msrc)) { #For each source river
        tmp <- sp::coordinates(msrc[j, ])[[1]][[1]] #Get coordinates of all of its nodes
        v2[j, ] <- tmp[dim(tmp)[1] - 1, ] - tmp[dim(tmp)[1], ] #Compute diff in lon and lat between the second to last and last nodes
      }
      ang <- atan2(v2[, 2], v2[, 1]) - atan2(v1[2], v1[1]) #Compute angle between source rivers and sink river
      ang <- (ang + pi) %% (2 * pi) - pi
      # ang[ang > 0] <- ang[ang > 0] %% pi
      # ang[ang < 0] <- -((-ang[ang < 0]) %% pi)
      names(ang) <- privers
      ang <- ang[order(abs(ang), decreasing = T)]
      
      if (sum(ang > 0) > sum(ang < 0)) {
        angneg <- c(as.numeric(names(ang[ang < 0])), rep(NA, sum(ang > 0) - sum(ang < 0)))
        angpos <- as.numeric(names(ang[ang >= 0]))
      } else {
        angpos <- c(as.numeric(names(ang[ang >= 0])), rep(NA, sum(ang < 0) - sum(ang > 0)))
        angneg <- as.numeric(names(ang[ang < 0]))
      }
      angs <- matrix(c(angpos, angneg), nc = 2)
      privers <- na.omit(as.vector(t(angs)))[1:delta]
      
      # split sink river in delta pieces plus the remainder
      tmplines <- split_sp_line(line = coo, 
                                n = delta, 
                                length = step,
                                debug = F)
      
      if (length(tmplines) <= delta) stop("You must decrease step: river ", sinkriver)
      
      for (j in 1:length(tmplines)) {
        # cut sink river into pieces, as many as necessary
        rid <- runif(1, 10^6, 10^7) # random ID for the piece
        # create a new piece with the j'th (step+1) vertices of the sink river
        piece <- sp::SpatialLines(list(sp::Lines(list(sp::Line(tmplines[[j]])), rid)), 
                                  proj4string = CRS)
        newdata <- mto@data
        newdata[1, ] <- NA
        newdata[1, from] <- newnodes[j]
        newdata[1, to] <- newnodes[j + 1]
        rownames(newdata) <- rid
        
        #Re-assign kept data 
        if (!is.null(fields_to_keep)) {
          newdata[[fields_to_keep]] <- mto[[fields_to_keep]]
        }
        
        pieces <- c(pieces, 
                    list(sp::SpatialLinesDataFrame(piece, newdata, match = F))) # save pieces for later use
        
        if (j > 1) {
          # now change coords of problematic rivers
          pri <- privers[j - 1] # pick the j'th problematic river
          tmp <- sp::coordinates(shp[pri, ])[[1]][[1]]
          tmp[dim(tmp)[1], ] <- tmplines[[j]][1, ] # change the coordinate of the last vertex of problematic river
          tmp1 <- sp::SpatialLines(list(sp::Lines(list(sp::Line(tmp)), 
                                                  shp[pri, ]@lines[[1]]@ID)), 
                                   proj4string = CRS) # keep same ID (original will be removed)
          tmp1 <- sp::SpatialLinesDataFrame(tmp1, shp[pri, ]@data)
          tmp1@data[1, to] <- newnodes[j]
          probrivers <- c(probrivers, list(tmp1)) # collect new rivers to replace old
        }
      }
      
      removeindexes <- c(removeindexes, c(privers[1:delta], sinkriver))
      maxnode <- maxnode + delta
    }
    
    cat("Now reassembling shape...\n")
    flush.console()
    newlines <- pieces[[1]]
    
    for (j in 2:length(pieces)) {
      newlines <- rbind(newlines, pieces[[j]])
    }
    
    for (j in 1:length(probrivers)) {
      newlines <- rbind(newlines, probrivers[[j]])
    }
    
    # remove all problematic + sink rivers
    newshp <- shp[-removeindexes, ]
    newshp <- rbind(newshp, newlines)
  } else {
    print("No complex confluences found. Returning input network unchanged.")
    newshp <- shp
  }
  
  return(newshp)
}

#-------------- workflow functions ---------------------------------------------
# path_list = tar_read(bio_data_paths)
# in_metadata_edna <- tar_read(metadata_edna)

#------ define_hydromod_paths --------------------------------------------------
#in_hydromod_dir <- hydromod_present_dir

#List data paths for hydrological data
define_hydromod_paths <- function(in_hydromod_dir) {
  hydro_drn_paths_dt <- data.table(
    country = c("Croatia", "Czech", "Finland", "France",  "Hungary", "Spain"),
    catchment = c("Butiznica", "Velicka", "Lepsamaanjoki", "Albarine", "Bukkosdi", "Genal"), 
    all_sims_filename = c(
      "Butiznica_2022-12-15_option0.nc", #_run8_final
      "Velicka_2023-02-01_option0.nc",
      "Lepsamaanjoki_2022-12-16_option0.nc",
      'Albarine_2022-12-16_option0.nc',
      "Bukkosdi_2022-12-16_option0.nc", #_run8_final
      "Genal_2023-01-18_option0.nc"
    ),
    sel_sim_filename = c(
      "Butiznica_2022-12-15_option0_run8_final.nc",
      "Velicka_2023-02-01_option0_run9_final.nc",
      "Lepsamaanjoki_2022-12-16_option0_run20_final.nc",
      "Albarine_2022-12-16_option0_run3_final.nc",
      "Bukkosdi_2022-12-16_option0_run8_final.nc", 
      "Genal_2023-01-18_option0_run15_final.nc"
    )
  ) %>%
    .[, `:=`(all_sims_path = file.path(in_hydromod_dir, catchment, 
                                       "Results_present_period", all_sims_filename),
             sel_sim_path = file.path(in_hydromod_dir, catchment, 
                                      "Results_present_period", sel_sim_filename),
             catchment_path = file.path(in_hydromod_dir, catchment,
                                        "watershed_small_catchment.shp"),
             network_path = file.path(in_hydromod_dir, catchment,
                                      "river_network.shp"),
             reaches_path = file.path(in_hydromod_dir, catchment, "reach.par"),
             sites_reachids = file.path(in_hydromod_dir, catchment,
                                        paste0(catchment, 
                                               '_sampling_sites_ReachIDs.csv')))]
}

#------ subset_network -----------------------------------------------------------
# in_hydromod_paths_dt <- tar_read(hydromod_paths_dt)
# out_dir <- file.path('results', 'gis')
# overwrite = FALSE

subset_network <- function(in_hydromod_paths_dt, out_dir, overwrite=FALSE) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  in_hydromod_paths_dt[
    , network_sub_path := file.path(
      out_dir, 
      paste0(tolower(country), '_river_network_sub_',
             format(Sys.time(), "%Y%m%d"),'.gpkg'))]
  
  in_hydromod_paths_dt[, {
    if (!file.exists(network_sub_path) | overwrite) {
      terra::vect(network_path) %>%
        .[relate(terra::vect(catchment_path), .,"intersects")[1,],] %>% #Contains removes some segments
        terra::writeVector(filename = network_sub_path,
                           overwrite = T)
    }
  }, by=country]
  
  return(in_hydromod_paths_dt[, stats::setNames(network_sub_path, country)])
}

#------ clean_network ------------------------------
# in_country <- 'Croatia'
# rivnet_path <- tar_read(network_sub_gpkg_list)[[in_country]]
# idcol <- 'cat'
# node_clustering_dist = 50
# min_segment_length = 20
# outdir = file.path(resdir, 'gis')
# save_gpkg = TRUE

clean_network <- function(rivnet_path, idcol, 
                          node_clustering_dist,
                          min_segment_length = 20,
                          outdir=NULL, save_gpkg=FALSE, 
                          return_path=FALSE) {
  #Read input network
  rivnet <- st_read(rivnet_path) %>%
    st_cast("LINESTRING") %>%
    #Make sure that the geometry column is equally named regardless 
    #of file format (see https://github.com/r-spatial/sf/issues/719)
    st_set_geometry('geometry') 
  
  #Preformat basic network
  sfnet_ini <- rivnet %>%
    as_sfnetwork %>%
    activate("edges") %>%
    filter(!edge_is_multiple()) %>% #Keep shortest of edges that connect the same pair of nodes
    filter(!edge_is_loop()) #Remove obvious loops: edges that start and end at the same node
  
  #------------------ Split lines at intersections -----------------------------
  #Get confluence nodes (nodes of third degree: with at least 3 intersecting edges)
  #2nd degree nodes are pseudonodes and 1st degree nodes are dangling
  splitting_nodes <- activate(sfnet_ini, nodes) %>%
    mutate(degree = igraph::degree(.)) %>%
    filter(degree >= 3) %>%
    st_as_sf("nodes")
  
  #Visualize interactively
  # ggplotly(
  #   ggplot(rivnet) +
  #     geom_sf() +
  #     geom_sf(data=splitting_nodes, color='red')
  # )
  #Write split nodes to double check
  #st_write(splitting_nodes, 'split_nodes.gpkg')
  
  #Simplify network by first fully dissolving and then splitting at confluences
  rivnet_agg <- sf::st_union(rivnet) %>%
    sf::st_line_merge(.) %>%
    lwgeom::st_split(., splitting_nodes)  %>%
    sf::st_collection_extract("LINESTRING")
  
  #------------------ Remove loops by clustering nearby nodes--------------------
  sfnet <- as_sfnetwork(rivnet_agg)
  
  node_coords <- sfnet %>%
    activate("nodes") %>%
    st_coordinates()
  
  # Cluster the nodes with the DBSCAN spatial clustering algorithm.
  # We set eps = 40 such that:
  # Nodes within a distance of 40 m from each other will be in the same cluster.
  # We set minPts = 1 such that:
  # A node is assigned a cluster even if it is the only member of that cluster.
  clusters = dbscan::dbscan(node_coords, 
                            eps = node_clustering_dist, minPts = 1)$cluster
  
  # Add the cluster information to the nodes of the network.
  clustered = sfnet %>%
    activate("nodes") %>%
    mutate(cls = clusters)
  
  #contracts groups of nodes based on cluster number as a grouping variable. 
  #The geometry of each contracted node is the centroid of the original 
  #group members’ geometries. Moreover, the geometries of the edges that start
  #or end at a contracted node are updated such that their boundaries match
  #the new node geometries
  sfnet_clustered <- tidygraph::convert( #
    clustered,
    sfnetworks::to_spatial_contracted,
    cls,
    simplify = TRUE
  ) %>%
    tidygraph::convert(sfnetworks::to_spatial_smooth) #Remove pseudo nodes
  
  #------------------ Trim spits (segments with dangle point under minimum length ) ------------
  #Convert back to sf
  rivnet_clustered <- st_as_sf(sfnet_clustered, "edges")
  
  #Re-aggregate and split the network
  rivnet_clustered_agg <- sf::st_union(rivnet_clustered) %>%
    sf::st_line_merge(.) %>%
    lwgeom::st_split(., splitting_nodes)  %>%
    sf::st_collection_extract("LINESTRING") %>%
    st_sf(geometry=.)
  
  #Compute segment length
  rivnet_clustered_agg$length <- as.numeric(st_length(rivnet_clustered_agg))
  
  #Identify dangle points
  agg_endpts <- c(lwgeom::st_startpoint(rivnet_clustered_agg), 
                  lwgeom::st_endpoint(rivnet_clustered_agg)) 
  
  danglepts <- agg_endpts[!(duplicated(agg_endpts) | 
                              duplicated(agg_endpts, fromLast=TRUE)),] %>%
    st_sf %>%
    mutate(dangle = 'dangle')
  
  #Remove spits under min length
  rivnet_clustered_aggsub <- st_join(x = rivnet_clustered_agg,
                                     y = danglepts, 
                                     join = st_intersects,
                                     left = TRUE,
                                     all.x = TRUE) %>%
    mutate(dangle = replace_na(dangle, "connected")) %>%
    dplyr::filter(!((dangle == 'dangle') &
                      (length < min_segment_length))) %>%
    select(-c(length, dangle))
  
  #------------------ Re-assign original IDs to all segments ---------------
  #Split back into component linestrings
  rivnet_endpts <- c(lwgeom::st_startpoint(rivnet), 
                     lwgeom::st_endpoint(rivnet))
  rivnet_clustered_resplit <- lwgeom::st_split(rivnet_clustered_aggsub, 
                                               rivnet_endpts)  %>%
    sf::st_collection_extract("LINESTRING") %>%
    dplyr::distinct(.) #Remove duplicate lines
  #.[-unlist(st_equals(., retain_unique=TRUE)),] 
  
  rivnet_clustered_resplit$UID <- seq_along(rivnet_clustered_resplit$geometry)
  
  #Join those that have not moved
  rivnet_clustered_joinini <- st_join(
    rivnet_clustered_resplit,
    rivnet,
    join = st_equals,
    suffix = c(".ini", ".clustered"),
    left = TRUE,
    largest = FALSE
  )
  
  #Convert those without match to points every 10 meters
  rivnet_nomatch_pts <- rivnet_clustered_joinini %>%
    .[is.na(rivnet_clustered_joinini[[idcol]]),] %>%
    st_line_sample(density = 1/10) %>%
    st_sf(geometry = ., crs = st_crs(rivnet_clustered_joinini))
  rivnet_nomatch_pts$UID <- rivnet_clustered_joinini[
    is.na(rivnet_clustered_joinini[[idcol]]),][['UID']]
  rivnet_nomatch_pts <- st_cast(rivnet_nomatch_pts, 'POINT')
  
  #For each point, get id of nearest line in initial river network
  nearest_segix <- rivnet_nomatch_pts %>%
    st_nearest_feature(., rivnet)
  rivnet_nomatch_pts$nearest_id <- rivnet[nearest_segix,][[idcol]]
  
  #Compute distance to that nearest line
  rivnet_nomatch_pts$nearest_dist <- st_distance(
    rivnet_nomatch_pts, rivnet[nearest_segix,], 
    by_element=TRUE)
  
  #Keep line for which the sum of the inverse distance to the points is greatest              
  rivnet_nomatch_selid <- as.data.table(rivnet_nomatch_pts) %>%
    .[, list(inverse_dist_sum = sum(1/nearest_dist)), 
      by=.(UID, nearest_id)] %>%
    .[, list(nearest_id=.SD[which.max(inverse_dist_sum), 
                            nearest_id]), 
      by=UID] 
  
  #Fill NAs with those
  rivnet_clustered_joinall <- merge(rivnet_clustered_joinini,
                                    rivnet_nomatch_selid, 
                                    by='UID', all.x=T)
  rivnet_clustered_joinall[[idcol]] <- dplyr::coalesce(
    rivnet_clustered_joinall[[idcol]],
    rivnet_clustered_joinall$nearest_id)
  
  #------------------ Write out results ------------------------------------------
  out_net <- rivnet_clustered_joinall[
    , c('UID', names(rivnet)[names(rivnet) != 'geometry'])]
  
  out_path <- file.path(outdir,
                        paste0(tools::file_path_sans_ext(basename(rivnet_path)),
                               '_clean',
                               format(Sys.time(), "%Y%m%d"),
                               '.gpkg')
  )
  
  if (save_gpkg) {
    st_write(out_net, out_path, append=F)
  }
  
  return(ifelse(return_path, out_path, out_net))
}




#------ direct_network -----------------------------------------
# Define helper functions
get_endpoint <- function(line) st_coordinates(line)[nrow(st_coordinates(line)), ]
get_startpoint <- function(line) st_coordinates(line)[1, ]

direct_network_inner <- function(segment, in_network, idcol, visited = NULL) {
  # #Reverse upstream segments recursively
  # visited <- NULL
  # segment <- rivnet[rivnet[[idcol]] == 22,] 
  # in_network = rivnet
  
  #print(segment[[idcol]])
  visited <- c(visited, segment[[idcol]])
  
  # Find connected segments 
  inseg_startpoint <- get_startpoint(segment$geometry)
  
  connected <- in_network[!(in_network[[idcol]] %in% visited),] %>%
    .[as.vector(st_intersects(., lwgeom::st_startpoint(segment$geometry), 
                              sparse=F)),]
  
  # ggplotly(ggplot(in_network[!(in_network[[idcol]] %in% visited),]) +
  #            geom_sf() +
  #            geom_sf(data=segment, color='red') +
  #            geom_sf(data=lwgeom::st_startpoint(segment$geometry)))
  
  # Reverse and recurse
  for (i in seq_len(nrow(connected))) {
    seg_id <- connected[[idcol]][i]
    seg_geom <- connected$geometry[i]
    # Reverse if not aligned
    if (!(all.equal(get_endpoint(seg_geom), inseg_startpoint) == TRUE)) {
      in_network[in_network[[idcol]]==seg_id,]$geometry <- st_reverse(seg_geom)
    }
    # Recursively process upstream
    in_network <- direct_network_inner(
      segment = in_network[in_network[[idcol]]==seg_id,], 
      in_network, idcol, visited)
  }
  
  return(in_network)
}


# outlet_uid_list <- list(Croatia = 458,
#                         Czech = 4,
#                         Finland = 682,
#                         France = 1,
#                         Hungary = 5,
#                         Spain = 86
# )
# 
# in_country <- 'Croatia'
# rivnet_path <- tar_read(network_clean_gpkg_list)[[in_country]]
# idcol <- 'UID'
# outletid <- outlet_uid_list[[in_country]]
# outdir = file.path(resdir, 'gis')
# save_gpkg = TRUE

direct_network <- function(rivnet_path, idcol,
                           outletid, outdir=NULL, 
                           save_gpkg=FALSE) {
  #Read input network
  rivnet <- st_read(rivnet_path) %>%
    st_cast("LINESTRING") %>%
    #Make sure that the geometry column is equally named regardless 
    #of file format (see https://github.com/r-spatial/sf/issues/719)
    st_set_geometry('geometry') 
  
  #Interactively plot network
  # (ggplot(rivnet) +
  #   geom_sf(size = 0.1,
  #           arrow = arrow(angle = 30,
  #                         length = unit(0.075, 'inches'),
  #                         ends = "last",
  #                         type = "closed")))
  
  #Identify dangle points
  agg_endpts <- c(lwgeom::st_startpoint(rivnet), 
                  lwgeom::st_endpoint(rivnet)) 
  
  danglepts <- agg_endpts[!(duplicated(agg_endpts) | 
                              duplicated(agg_endpts, fromLast=TRUE)),] %>%
    st_sf %>%
    mutate(dangle = 'dangle')
  
  #Make sure outlet segment is in the right direction
  outlet_endpt <- lwgeom::st_endpoint(
    rivnet[rivnet[[idcol]] == outletid,])
  
  if (!(outlet_endpt %in% danglepts$geometry)) {
    rivnet[rivnet[[idcol]] == outletid, 'geometry'] <- st_reverse(
      rivnet[rivnet[[idcol]] == outletid, 'geometry'])
  }
  
  outlet_seg <- rivnet[rivnet[[idcol]] == outletid,]
  
  # Apply to the entire network
  out_net <- direct_network_inner(segment = outlet_seg, 
                                  in_network = rivnet, 
                                  idcol = idcol,
                                  visited = NULL)
  
  #------------------ Write out results ------------------------------------------
  out_path <- file.path(outdir,
                        paste0(tools::file_path_sans_ext(basename(rivnet_path)),
                               '_directed',
                               format(Sys.time(), "%Y%m%d"),
                               '.gpkg')
  )
  
  if (save_gpkg) {
    st_write(out_net, out_path, append=F)
  }
  
  return(out_path)
}

#------ fix_complex_confluences ------------------------------------------------
# in_country <- 'Croatia'
# rivnet_path = tar_read(network_directed_gpkg_list)[[in_country]]
# outdir = file.path(resdir, 'gis')
# max_node_shift = 5

#Network must be directed & topologically correct aside from the complex confluences
fix_complex_confluences <- function(rivnet_path, max_node_shift = 5,
                                    outdir=NULL, out_path=NULL) {
  # Read input network
  rivnet <- st_read(rivnet_path) %>%
    st_cast("LINESTRING") %>%
    st_set_geometry("geometry")
  
  #Compute from-to fields
  net_fromto <- as_sfnetwork(rivnet) %>%
    activate(edges) %>%
    as.data.table %>%
    .[, .(from, to, UID)]
  
  rivnet_fromto <- merge(rivnet, net_fromto, by='UID', all.x=T)
  
  #Run fix streamlines from Miguel Porto
  rivnet_fixed <- fix_confluences_inner(shp = as_Spatial(rivnet_fromto), 
                                        from = "from",
                                        to = "to", 
                                        step = max_node_shift,
                                        fields_to_keep = 'cat') %>%
    .[, !(names(.) %in% c('from', 'to'))] %>%
    st_as_sf
  
  rivnet_fixed[is.na(rivnet_fixed$UID),]$UID <- max(rivnet_fixed$UID, na.rm=T) + 
    seq_len(sum(is.na(rivnet_fixed$UID)))
  
  #------------------ Write out results ------------------------------------------
  if (is.null(out_path)) {
    out_path <- file.path(outdir,
                          paste0(tools::file_path_sans_ext(basename(rivnet_path)),
                                 '_conflufixed',
                                 format(Sys.time(), "%Y%m%d"),
                                 '.gpkg')
    )
  }
  
  st_write(rivnet_fixed, out_path, append=F)
  
  return(out_path)
}

#------ assign_strahler_order --------------------------------------------------
# in_country <- 'Croatia'
# in_rivnet = network_ssnready_gpkg_list[[in_country]]
# idcol = 'UID'

assign_strahler_order <- function(in_rivnet, idcol, verbose = F) {
  if (is.character(in_rivnet)) {
    rivnet <- st_read(in_rivnet)
    
    #Compute from-to fields
    net_fromto <- as_sfnetwork(rivnet) %>%
      activate(edges) %>%
      as.data.table %>%
      .[, c('from', 'to', idcol), with=F] %>%
      merge(.[, list(nsource = .N), by=to], 
            by.x='from', by.y='to', all.x=T) %>%
      .[is.na(nsource), nsource := 0]
    
    #Join to spatial network
    rivnet_fromto <- merge(rivnet, net_fromto, by=idcol)
    
    #Convert to data.table for speed and syntax
    rivnet_fromto_dt <- as.data.table(rivnet_fromto)
    
  } else if (is.data.table(in_rivnet)) {
    rivnet_fromto_dt <- merge(in_rivnet,
                              in_rivnet[, list(nsource = .N), by=to], 
                              by.x='from', by.y='to', all.x=T) %>%
      .[is.na(nsource), nsource := 0]
  }
  
  #Assign strahler 1 to lines with no source line
  rivnet_fromto_dt[nsource == 0, strahler := 1]
  
  # Compute Strahler order iteratively
  #while (sum(is.na(rivnet_fromto_dt$strahler)) > 417) {
  while (any(is.na(rivnet_fromto_dt$strahler))) {
    if (verbose) { print(sum(is.na(rivnet_fromto_dt$strahler)))}
    # Identify lines whose sources' Strahler orders are all assigned
    rivnet_fromto_dt <-  merge(rivnet_fromto_dt, rivnet_fromto_dt[, .(to, strahler)], 
                               by.x='from', by.y='to', suffixes = c('_down', '_up'),
                               all.x=T)
    
    #Extend strahler order downstream for consecutive sections on the same segment 
    #(only one source section)
    rivnet_fromto_dt[is.na(strahler_down) & !is.na(strahler_up) & nsource == 1,
                     strahler_down := strahler_up]
    
    #Flag as eligible those segments for which all upstream segments have
    #a strahler order
    rivnet_fromto_dt[is.na(strahler_down) & nsource >= 2 & !is.na(strahler_up),
                     eligible := ((.N>=2) & .N==nsource), by=idcol]
    
    #Identify max upstream strahler and number of upstream sections with that strahler order
    rivnet_fromto_dt[eligible & !is.na(eligible), 
                     max_strahler_u := max(strahler_up),
                     by=idcol]
    rivnet_fromto_dt[eligible & !is.na(eligible),
                     n_max_strahler_u := .SD[strahler_up==max_strahler_u, .N],
                     by=idcol]
    
    rivnet_fromto_dt[eligible & !is.na(eligible),
                     strahler_down := fifelse(
                       n_max_strahler_u >= 2,
                       max_strahler_u + 1,
                       max_strahler_u),
                     by = idcol]
    
    setnames(rivnet_fromto_dt, 'strahler_down', 'strahler')
    
    rivnet_fromto_dt <- rivnet_fromto_dt[
      !duplicated(get(idcol)), 
      -c('strahler_up', 'eligible', 'max_strahler_u', 'n_max_strahler_u'), 
      with=F]
  }
  
  return(rivnet_fromto_dt[
    , c(idcol, 'from', 'to', 'strahler', 'nsource'), with=F])
}

#------ reassign_netids --------------------------------------------------
#Problem: 'cat' or 'reach_id' in the shapefile, associated with lines across 
#confluences. Makes it impossible to correctly use network hydrology data.

#There are three spatial units in this code:
#Segment sections: uniquely identified by UID, these are the individual lines
#                 contained in the shapefile generated by WP1 after hydrological
#                 modeling.
#Segments: uniquely identifed by UID_fullseg, representing lines extending between 
#           two confluences. They often contain multiple sections.
#Cat: uniquely identified by cat (or reach_id), these are unique reaches as modeled
#     by the hydrological model, but not well represented/assigned to segment
#     sections  in the shapefile. Some may be omitted by the shapefile, 
#     many span multiple segment sections and even multiple segments across 
#     confluences and Strahler stream order.

#The goal of this function is to re-assign a topologically/hydrologically
#logical cat to each segment section to be able to join the hydrological model 
#outputs to the shapefile.

#Helper functions
#For those segments where only one cat remains, assign cat to the entire segment
#at least a given proportion of the length of the full segment
assign_singlecat_to_seg <- function(net_dt, length_threshold = 1/3) {
  # Identify segments with only one remaining section
  #with a cat (length_uid_fullsegper == 1) and where the actual length of that
  #section is at least a given proportion of the length of the full segment
  #including all sections.
  cat_toassign <- net_dt[
    is.na(cat_cor) & (length_uid_fullsegper == 1) &
      (length_uid / length_fullseg > length_threshold),
    .(UID_fullseg, cat)
  ] %>%
    setnames('cat', 'cat_cor')
  
  #Assign those cat_cor to all sections within those segments
  net_dt[is.na(cat_cor), 
         cat_cor := cat_toassign[.SD, on = 'UID_fullseg', cat_cor]]
}

#For segment sections with cat==NAs and cat_cor==NA, 
#get cat_cor from upstream section in the same segment
get_nearby_cat <- function(in_dt, source_direction, 
                           strahler_list = seq(1,15)) {
  if (source_direction == 'upstream') {
    source_col <- 'to' #Column of section to get cat_cor from
    sink_col <- 'from' #Column of section to assign cat_cor to
  } else if (source_direction == 'downstream') {
    source_col <- 'from' #Column of section to get cat_cor from
    sink_col <- 'to' #Column of section to assign cat_cor to
  }
  
  # Identify sections to fill
  sections_tofill <- in_dt[
    , which(is.na(cat) & is.na(cat_cor) & (strahler %in% strahler_list))]
  
  # Prepare cat_cor values to assign
  cat_cor_toassign <- in_dt[
    !is.na(cat_cor) & 
      (strahler %in% strahler_list) & 
      (get(source_col) %in% in_dt[sections_tofill, unique(get(sink_col))]), 
    c(source_col, 'UID_fullseg', 'cat_cor'), with=F] %>%
    setnames(source_col, sink_col)
  
  # Assign cat_cor to sections_tofill using a join
  in_dt[sections_tofill,
        cat_cor := cat_cor_toassign[
          .SD, on=c(sink_col, 'UID_fullseg'), cat_cor]]
}

#Remove pseudo-nodes among segments sections provided equal attributes
#Re-compute from-to and length_uid
remove_pseudonodes <- function(in_net, equal_cols = FALSE, 
                               summarise_attributes ='first') {
  out_net <- as_sfnetwork(in_net) %>%
    activate(edges) %>%
    convert(to_spatial_smooth, 
            require_equal = equal_cols,
            summarise_attributes = summarise_attributes) %>%
    st_as_sf() %>%
    mutate(length_uid = as.numeric(st_length(.)))
  out_net
}

#Parameters
# in_country <- 'Spain'
# rivnet_path <- tar_read(network_nocomplexconf_gpkg_list)[[in_country]]
# strahler_dt <- tar_read(network_strahler)[[in_country]]
# in_reaches_hydromod_dt <- tar_read(reaches_dt)[country==in_country,]

reassign_netids <- function(rivnet_path, strahler_dt, 
                            in_reaches_hydromod_dt, outdir,
                            country = NULL, in_ext='gpkg') {
  #---------- Prepare data -------------------------------------------------------
  #Read network and join with hydromod data
  rivnet <- st_read(rivnet_path) %>%
    #Make sure that the geometry column is equally named regardless 
    #of file format (see https://github.com/r-spatial/sf/issues/719)
    st_set_geometry('geometry') %>%
    merge(strahler_dt, by='UID') 
  
  #Remove pseudonodes to define full segments between confluences
  rivnet_fullseg <- as_sfnetwork(rivnet) %>%
    activate("edges") %>%
    convert(to_spatial_smooth) %>%
    st_as_sf() 
  
  #Give these full segments between confluences a separate ID: UID_fullseg
  rivnet_fullseg$UID_fullseg <- seq_len(nrow(rivnet_fullseg))
  
  #Intersect with original network to check length of overlap for each 'cat'
  rivnet_fullseg_inters <- st_intersection(
    rivnet, rivnet_fullseg[, c('geometry', 'UID_fullseg')]) %>%
    .[st_geometry_type(.) != 'POINT',] %>%
    st_cast('MULTILINESTRING') %>%
    st_cast('LINESTRING')
  
  #Remove pseudo-nodes among segments sections of the same cat
  #and compute length of segment sections (UID)
  rivnet_fullseg_inters_nopseudo <- remove_pseudonodes(
    in_net = rivnet_fullseg_inters, 
    equal_cols = c("cat", "UID_fullseg"), 
    summarise_attributes='first')
  
  #Create data.table to work with
  rivnet_inters_dt <- as.data.table(rivnet_fullseg_inters_nopseudo) 
  
  #Back up cat before editing it
  rivnet_inters_dt[, cat_copy := cat]
  
  #Compute other lengths
  rivnet_inters_dt[, length_cat := sum(length_uid), by=cat]
  rivnet_inters_dt[, length_fullseg := sum(length_uid), by=UID_fullseg]
  
  #Compute the percentage length of that cat ("hydromod" reach) that is 
  #represented by that segment section (identified by UID)
  rivnet_inters_dt[, length_uid_catper := length_uid/length_cat]
  #Check whether this is the largest segment section for that cat 
  rivnet_inters_dt[, length_uid_catmax := (length_uid == max(length_uid)), by=cat]
  #Compute the percentage length of that full segment that is represented by
  #that segment section (identified by UID)
  rivnet_inters_dt[, length_uid_fullsegper := length_uid/length_fullseg]
  
  #Computer number of full segments that a cat overlaps
  rivnet_inters_dt[, n_seg_overlap := length(unique(UID_fullseg)), by=cat]
  
  #Format network data from hydrological model
  reaches_hydromod_format <- in_reaches_hydromod_dt[
    , .(ID, to_reach, length, strahler)] %>%
    setnames(c('ID_hydromod', 'to_reach_hydromod', 
               'length_hydromod', 'strahler_hydromod')) %>%
    .[, `:=`(from = ID_hydromod, to = to_reach_hydromod)] %>%
    merge( #Compute actual strahler order
      assign_strahler_order(in_rivnet = ., idcol = 'ID_hydromod', verbose = F),
      by = 'ID_hydromod'
    ) %>%
    setnames(c('strahler', 'nsource'),
             c('strahler_hydromod_recalc', 'nsource_hydromod'))
  
  #---------- Initial assignment of correct cat -------------------------------------------------
  #For first order segments where the cat is only represented by that segment,
  #keep that cat for this section of the segment (see cat==2938 for Croatia)
  rivnet_inters_dt[n_seg_overlap == 1 & strahler == 1, 
                   cat_cor := cat]
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  
  #For full segments that only have one cat associated with them, confirm cat
  rivnet_inters_dt[length_uid == length_fullseg, 
                   cat_cor := cat]
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  
  #For cats that are already assigned, remove their "cat" from other segments
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #Re-compute length_uid_fullsegper, excluding segment sections whose initial
  #cat was assigned elsewhere
  rivnet_inters_dt[!is.na(cat), 
                   length_uid_fullsegper := length_uid/sum(length_uid),
                   by=UID_fullseg]
  
  #For those segments where only one cat remains, assign cat to the entire segment
  assign_singlecat_to_seg(rivnet_inters_dt)
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  
  #For cats that are already assigned, remove their "cat" from other segments
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #For segment sections with cat==NAs and cat_cor==NA in first order streams, 
  #get cat_cor from upstream segment in same segment
  get_nearby_cat(rivnet_inters_dt, 
                 source_direction = 'upstream', 
                 strahler_list = 1) 
  
  #Re-compute length_uid_fullsegper, excluding segment sections whose initial
  #cat was assigned elsewhere
  rivnet_inters_dt[!is.na(cat), 
                   length_uid_fullsegper := length_uid/sum(length_uid),
                   by=UID_fullseg]
  
  #For those segments where only one cat remains, assign cat to the entire segment
  assign_singlecat_to_seg(rivnet_inters_dt)
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  
  #For cats that are already assigned, remove their "cat" from other segments
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #Computer number of full segments that a cat overlaps
  rivnet_inters_dt[, n_seg_overlap := length(unique(UID_fullseg)), by=cat]
  
  #For first and second order segments where the cat is now only represented 
  #by that segment, keep that cat for this section of the segment
  rivnet_inters_dt[n_seg_overlap == 1 & strahler <= 2 & is.na(cat_cor), 
                   cat_cor := cat]
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  
  #For cats that are already assigned, remove their "cat" from other segments
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #For segment section with cat==NAs in first order streams, assign cat_cor of
  #upstream segment
  get_nearby_cat(rivnet_inters_dt, 
                 source_direction = 'upstream', 
                 strahler_list = 1) 
  
  #Re-compute length_uid_fullsegper, excluding segment sections whose initial
  #cat was assigned elsewhere
  rivnet_inters_dt[!is.na(cat), 
                   length_uid_fullsegper := length_uid/sum(length_uid),
                   by=UID_fullseg]
  
  #For those segments where only one cat remains, assign cat to the entire segment
  assign_singlecat_to_seg(rivnet_inters_dt)
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  
  #For cats that are already assigned, remove their "cat" from other segments
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #For sections that represent over 75% of the total length of that cat, assign cat_cor
  rivnet_inters_dt[length_uid_catper>0.7 & is.na(cat_cor), cat_cor := cat]
  
  #For cats that are already assigned, remove their "cat" from other segments
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #Re-compute length_uid_fullsegper, excluding segment sections whose initial
  #cat was assigned elsewhere
  rivnet_inters_dt[!is.na(cat), 
                   length_uid_fullsegper := length_uid/sum(length_uid),
                   by=UID_fullseg]
  
  #For those segments where only one cat remains, assign cat to the entire segment
  assign_singlecat_to_seg(rivnet_inters_dt)
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  
  #For cats that are already assigned, remove their "cat" from other segments
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #For segment sections in strahler==1 & n_seg_overlap > 2, cat := NA 
  rivnet_inters_dt[strahler == 1 & n_seg_overlap > 2 & is.na(cat_cor), 
                   cat := NA] 
  
  #if main representative of cat, assign cat_cor
  rivnet_inters_dt[!is.na(cat_cor), cat := cat_cor]
  rivnet_inters_dt[, length_cat := sum(length_uid), by=cat]
  rivnet_inters_dt[, length_uid_catper := length_uid/length_cat, by=cat]
  rivnet_inters_dt[length_uid_catper > 0.7 & is.na(cat_cor), cat_cor := cat]
  
  #For cats that are already assigned, remove their "cat" from other segments
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #Re-compute length_uid_fullsegper, excluding segment sections whose initial
  #cat was assigned elsewhere
  rivnet_inters_dt[!is.na(cat_cor), 
                   length_uid_fullsegper := length_uid/sum(length_uid),
                   by=UID_fullseg]
  
  ##For those segments where only one cat remains, assign cat to the entire segment
  #regardless of the proportion of the segment it represents
  cat_cor_toassign <-  rivnet_inters_dt[
    is.na(cat_cor) & length_uid_fullsegper == 1,
    .(UID_fullseg, cat)] %>%
    setnames('cat', 'cat_cor')
  rivnet_inters_dt[is.na(cat_cor), 
                   cat_cor := cat_cor_toassign[.SD, on='UID_fullseg', 
                                               cat_cor]]
  
  #For the remaining segments where all cats are NAs (have been assigned elsewhere), 
  #assign the cat of the section representing the highest percent of the segment
  #These are usually near loops that have been removed 
  rivnet_inters_dt[is.na(cat_cor), NAlength := sum(length_uid),
                   by=UID_fullseg]
  rivnet_inters_dt[(NAlength==length_fullseg), 
                   cat_cor := fifelse(length_uid == max(length_uid), cat_copy, NA), 
                   by=UID_fullseg]
  
  #For all remaining sections, remove their "cat" 
  rivnet_inters_dt[is.na(cat_cor), cat := NA]
  
  for (i in 1:3) {
    #For segment section with cat==NAs and cat_cor==NA, 
    #assign cat_cor of upstream segment of same strahler order
    get_nearby_cat(rivnet_inters_dt, source_direction = 'upstream')
    
    #For is.na(cat) & is.na(cat_cor), 
    #assign downstream cat_cor of same UID_fullseg
    get_nearby_cat(rivnet_inters_dt, source_direction = 'downstream') 
  }
  
  #---------- Merge back with spatial data and remove pseudo nodes ---------------
  #Merge processed attributes with pre-formatted spatial network
  rivnet_catcor <- merge(
    rivnet_fullseg_inters_nopseudo, 
    rivnet_inters_dt[, .(UID, cat_cor, cat_copy)],
    by='UID') %>%
    .[, c('UID', 'cat_cor', 'cat_copy', 'strahler', 
          'nsource', 'UID_fullseg', 'geometry')]
  
  #Remove pseudo-nodes among segments sections of the same cat_cor
  #re-compute length uid
  rivnet_catcor_smooth <- remove_pseudonodes(
    in_net = rivnet_catcor, 
    equal_cols = c("cat_cor", "UID_fullseg"))
  
  
  #---------- Adjust results based on topology from hydrological model network ---
  #Merge with hydrological model network topology data
  rivnet_catcor_hydromod <- merge(
    rivnet_catcor_smooth,
    reaches_hydromod_format[
      , .(ID_hydromod, to_reach_hydromod, length_hydromod, 
          strahler_hydromod_recalc, nsource_hydromod)],
    by.x = 'cat_cor', by.y = 'ID_hydromod',
    all.x = T
  ) %>%
    as.data.table
  
  if (rivnet_catcor_hydromod[is.na(cat_cor), .N] == 0) {
    #Remove cats that are not supposed to be downstream of any other cat
    #when there are sections of the correct stream order on that segment
    rivnet_catcor_hydromod[, n_cats_fullseg := length(unique(cat_cor)),
                           by=UID_fullseg]
    
    rivnet_catcor_hydromod[n_cats_fullseg > 1 & 
                             nsource_hydromod == 0 & nsource > 0,
                           `:=`(cat = NA, cat_cor  = NA)]
    
    #assign cat_cor of upstream segment of same strahler order
    for (i in 1:3) {
      get_nearby_cat(in_dt=rivnet_catcor_hydromod, source_direction = 'upstream')
    }
    get_nearby_cat(in_dt=rivnet_catcor_hydromod, source_direction = 'downstream')
    
    rivnet_catcor_hydromod <- st_as_sf(rivnet_catcor_hydromod) %>%
      remove_pseudonodes(equal_cols = c("cat_cor", "UID_fullseg")) %>%
      as.data.table
  }
  
  #---------- Check results and correct based on hydrological model topology --------
  #Compare downstream segment cat
  to_reach_shpcor <-  rivnet_catcor_hydromod[
    , .(from, cat_cor)] %>%
    setnames(c('to', 'to_reach_shpcor'))
  
  rivnet_catcor_hydromod <- merge(rivnet_catcor_hydromod, 
                                  to_reach_shpcor, by='to', all.x=T) %>%
    .[, hydromod_shpcor_match := (to_reach_hydromod == to_reach_shpcor)]
  
  #If two upstream sections both should point to the same reach but do not
  #correct if downstream cat is in network topology but missing in shapefile
  #or if the assigned cat_cor is duplicated
  dupli_catcor <- rivnet_catcor_hydromod[, .N, by=cat_cor][N>=2,]$cat_cor
  missing_catcor <- reaches_hydromod_format[
    !(ID_hydromod %in% unique(rivnet_catcor_hydromod$cat_cor)),]$ID_hydromod
  
  #Remove first-order sections with no upstream sections whose downstream cat_cor
  #and strahler order does not match network topology from hydrological model
  rivnet_catcor_hydromod <- rivnet_catcor_hydromod[
    !(hydromod_shpcor_match==FALSE & nsource == 0 & strahler_hydromod_recalc>1),]
  
  #Identify pairs of reaches that are both upstream of the wrong reach
  double_downstream_mismatch_correct <- rivnet_catcor_hydromod[
    (to_reach_hydromod != to_reach_shpcor),
    list(n1=.N, to, to_reach_shpcor), by=to_reach_hydromod] %>%
    .[n1>=2, list(n2=.N, to_reach_hydromod, to_reach_shpcor), by=to] %>%
    .[n2>=2 & ((to_reach_hydromod %in% missing_catcor) |
                 (to_reach_shpcor %in% dupli_catcor)),] %>%
    unique %>%
    setnames(c('to_reach_shpcor', 'to'), c('cat_cor', 'from'))
  
  rivnet_catcor_hydromod[from %in% double_downstream_mismatch_correct$from &
                           cat_cor %in% double_downstream_mismatch_correct$cat_cor
                         , cat_cor := double_downstream_mismatch_correct[.SD, on=c('cat_cor', 'from'), 
                                                                         to_reach_hydromod]]
  
  #---------- Implement a few manual corrections  ---------------------------------
  if (country == 'Croatia') {
    rivnet_catcor_manual <- rivnet_catcor_hydromod %>%
      filter(!(UID %in% c(489, 90))) %>% #cat_cor 2082, and 2480, respectively) %>%
      mutate(cat_cor = case_match(
        UID,
        103 ~ 1176,  #UID 103 (catcor 1832) -> catcor 1776
        116 ~ 1832, #UID 116 (catcor 1836) -> catcor 1832
        31 ~ 1788, #UID 31 (catcor 1792) -> catcor 1788
        1045 ~ 2898, #UID 1045 (cat_cor 2908) -> catcor 289
        .default = cat_cor
      )
      )
  } else if (country == 'Czech') {
    rivnet_catcor_manual <- rivnet_catcor_hydromod %>%
      mutate(cat_cor = case_match(
        UID,
        298 ~ 40232, #UID 298 (catcor 40126) -> catcor 1776
        158 ~ 40258, #UID 158 (40260) -> catcor 40258
        .default = cat_cor
      )
      )
  } else if (country == 'Finland') {
    rivnet_catcor_manual <- rivnet_catcor_hydromod %>%
      mutate(cat_cor = case_match(
        UID,
        538 ~ 1049400, #UID 538 (catcor 1051800) -> catcor 1049400
        392 ~ 1069001, #UID 392 (catcor 1086800) -> catcor 1069001
        397 ~ 1069001, #UID 397 (catcor 1086600) -> catcor 1069001 
        .default = cat_cor
      )
      )
  } else if (country == 'France') {
    rivnet_catcor_manual <- rivnet_catcor_hydromod %>%
      filter(!(UID %in% c(529, 630, 642, 991))) %>% #2489000, 2475800, 2476200, 2475800
      mutate(cat_cor = case_match(
        UID,
        205 ~ 2422600, #UID 205  (catcor 2497800) -> catcor 2422600
        693 ~ 2444000, #UID 693 (catcor 2465400) -> catcor 2444000  
        805 ~ 2467600, #UID 805 (catcor 2467800) -> catcor 2467600
        500 ~ 2485400, #UID 500 (catcor 2485600) -> catcor 2485400
        .default = cat_cor
      )
      )
  } else if (country == 'Hungary') {
    rivnet_catcor_manual <- rivnet_catcor_hydromod %>%
      mutate(cat_cor = case_match(
        UID,
        51 ~ 650800,  #UID 51 (catcor 651000) -> catcor 650800
        331 ~ 652001,  #UID 331 (catcor 652200) -> catcor 652001
        379 ~ 673000,  #UID 379 (catcor 673400) -> catcor 673000
        .default = cat_cor
      )
      )
  } else if (country == 'Spain') {
    rivnet_catcor_manual <- rivnet_catcor_hydromod 
    rivnet_catcor_manual[rivnet_catcor_manual$UID == 25, 'cat_cor'] <- 5426 #UID 25 (catcor 5428) -> catcor 5426
  }
  
  
  #------ Final processing -------------------------------------------------------
  #Remove pseudo-nodes among segments sections of the same cat_cor
  #re-compute length uid, and re-assign from-to
  out_rivnet<- remove_pseudonodes(
    in_net = st_as_sf(rivnet_catcor_manual), 
    equal_cols = c("cat_cor", "UID_fullseg"))
  
  #re-assign downstream segment cat
  to_reach_shpcor_dt <-  as.data.table(out_rivnet)[
    , .(from, cat_cor)] %>%
    setnames(c('to', 'to_reach_shpcor'))
  
  out_rivnet <- merge(
    as.data.table(out_rivnet)[, -c('to_reach_shpcor', 'to_reach_hydromod'), with=F], 
    to_reach_shpcor_dt, by='to', all.x=T) %>%
    merge(reaches_hydromod_format[, .(ID_hydromod, to_reach_hydromod)],
          by.x='cat_cor', by.y='ID_hydromod', all.x=T) %>%
    .[, hydromod_shpcor_match := (to_reach_hydromod == to_reach_shpcor)]
  
  #---------- Write out results ------------------------------------------------------
  out_path <- file.path(outdir,
                        paste0(tools::file_path_sans_ext(basename(rivnet_path)),
                               '_reided',
                               format(Sys.time(), "%Y%m%d"),
                               '.', in_ext)
  )
  
  #Export results to gpkg
  write_sf(st_as_sf(out_rivnet)[
    , c('UID', 'strahler','length_uid', 'UID_fullseg', 'cat_cor', 'from', 'to',
        'to_reach_shpcor', 'to_reach_hydromod', 'hydromod_shpcor_match')],
    out_path)
  
  return(out_path)
}
