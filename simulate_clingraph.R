#main worker function for simulating two-group (or more) graph structure
simulate_clingraph <- function(gobj, npergroup=100, groupmanipulation=NULL, edgenoise=0, proportional=FALSE) {
  # groupmanipulation syntax/structure: list of lists
  # [[1]]
  #   $edgeshift_weight_mean_bw: mean shift in weights across edges for this group
  #   $edgeshift_weight_sd_bw: variation in magnitude of weight shift across replications/subjects
  #   $edgeshift_weight_sd_wi: variation in magnitude of weight shift across nodes within subject
  #   $group_name: name of the group
  #   $nodemanipulation: list of nodes to be manipulated
  #
  # nodemanipulation structure: list of nodes (embedded within groupmanipulation)
  # [[1]]
  #   $name: name of node: must match vertex names in graph
  #   $edgeshift_weight_mean_bw: nodal mean shift in edge weights
  #   $edgeshift_weight_sd_bw: nodal variation in magnitude of weight shift across replications
  #   $edgeshift_weight_sd_wi: nodal variation in magnitude of weight shift across edges connected to this node
  
  require(foreach)
  require(doSNOW)
  
  setDefaultClusterOptions(master="localhost")
  clusterobj <- makeSOCKcluster(4) #4 cpus
  registerDoSNOW(clusterobj)
  
  on.exit(try(stopCluster(clusterobj)))
  
  weighted <- FALSE
  if ("weight" %in% edge_attr_names(gobj)) { weighted = TRUE }
  
  stopifnot(weighted==TRUE) #currently working on weighted variant
  
  #handle empty group structure (single group sim): basically just useful to look at edge noise across replications
  #expected group structure
  if (is.null(groupmanipulation)) {
    groupmanipulation <- list(list(group_name="all", edgeshift_weight_mean_bw=0, edgeshift_weight_sd_wi=0, edgeshift_weight_sd_bw=0))    
  }
  
  #node manipulations are treated as property of group specification
  which_nodemanip <- sapply(groupmanipulation, function(group) { !is.null(group$nodemanipulation) })
  #has_nodemanip <- any(which_nodemanip)
  
  allsims <- list()
  #build sims
  for (g_iter in 1:length(groupmanipulation)) {
    #Gaussian variation in strength of shift across replications for a group
    repshift <- groupmanipulation[[g_iter]]$edgeshift_weight_mean_bw + rnorm(npergroup, 0, groupmanipulation[[g_iter]]$edgeshift_weight_sd_bw) #mean shift + between-rep variation
    
    groupsims <- foreach(j = iter(1:npergroup), .packages="igraph") %dopar% {
      g_rep <- gobj
      
      #simulate Gaussian variation on edges 
      if (weighted) {
        w <- E(gobj)$weight #vector of master weights
        
        if (edgenoise > 0) { w <- w + rnorm(length(w), 0, edgenoise) } #add Gaussian variation in edge strength over replications
        
        w <- w + rnorm(length(w), repshift[j], groupmanipulation[[g_iter]]$edgeshift_weight_sd_wi) #add Gaussian variation with mean = between and variation around that mean = within
        E(g_rep)$weight <- w #weights with shift and noise  
      }
      
      g_rep #return
    }

    #only apply node manipulations if specified for this group
    if (which_nodemanip[g_iter]) {      
      # setup an edgelist of shifts
    
      cat("NB: Shifts in edge strength are applied in order of the node list.\n")
      cat("Furthermore, once a shift has been applied to an edge, that edge is immune to further shifts (no compounding)\n")
      
      edgesmanipulated <- c()
      for (node in groupmanipulation[[g_iter]]$nodemanipulation) {
        #node specification includes mean shift, variation in shift between replications (bw), and variation in shift within incident edges (wi)
        #obtain edge sequence for connections incident to this node
        node_edges <- E(gobj)[inc(node$name)]
        node_edges <- difference(node_edges, edgesmanipulated) #pull out any edges already tweaked
        edgesmanipulated <- c(edgesmanipulated, node_edges) #add edges to edge sequence for this node
        
        noderepshift <- node$edgeshift_weight_mean_bw + rnorm(npergroup, 0, node$edgeshift_weight_sd_bw) #mean shift + between-rep variation
        
        #tweak this node for each replication (loop over replications)
        groupsims <- lapply(1:npergroup, function(j) {
              #add Gaussian variation with mean = replication between and variation around that mean = within
              E(groupsims[[j]])[node_edges]$weight <- E(groupsims[[j]])[node_edges]$weight + rnorm(length(node_edges), noderepshift[j], node$edgeshift_weight_sd_wi)
              return(groupsims[[j]])
            })
      }
    }
    
    allsims[[ groupmanipulation[[g_iter]]$group_name ]] <- groupsims
  }
  
  return(allsims)
  
}


sanitycheck_sims <- function(g_positive) {
  #check that variation in graphs is as expected
  
  edgeshift_weight_sd_wi=.06
  edgeshift_weight_mean_bw=0.5
  edgeshift_weight_sd_bw=.25
  
  simlist <- simulate_clingraph(g_positive, npergroup=100, edgenoise=0.0, #need 0 edgenoise for sanity check
      groupmanipulation=list(
          list(group_name="patients", edgeshift_weight_mean_bw=edgeshift_weight_mean_bw,
              edgeshift_weight_sd_wi=edgeshift_weight_sd_wi, edgeshift_weight_sd_bw=edgeshift_weight_sd_bw)))
  
  cat("Running sanity checks on weight shifts\n\n")
  #sanity checks: are the means and SDs as expected for patient group?
  #these checks will only hold up under 0 edge noise (or for huge groups in the asymptote)
  patients <- simlist[["patients"]]
  
  #per-replication average deviation between ground_truth and average weight
  repdevmeans <- sapply(patients, function(replication) {
        mean(E(replication)$weight - E(g_positive)$weight)
      })
  
  #between-replication mean shift: should be close to edgeshift_weight_mean_bw in the asymptote (VERIFIED 30Sep2016)
  cat("Expected mean shift in group:", edgeshift_weight_mean_bw, "\n")
  cat("Observed mean shift:", mean(repdevmeans), "\n\n")
  
  #across-replication variation in mean shift (sd of mean shifts): should be close to edgeshift_weight_sd_bw (VERIFIED 30Sep2016)
  cat("Expected variation in shift across replications (bw):", edgeshift_weight_sd_bw, "\n")
  cat("Observed variation:", sd(repdevmeans), "\n\n")
  
  #per-replication (within) edge variation from ground truth (across nodes)
  repdevsds <- sapply(patients, function(replication) {
        sd(E(replication)$weight - E(g_positive)$weight)
      })
  
  #overall sd across nodes (within): this should be close to edgeshift_weight_sd_wi in the asymptote (VERIFIED 30Sep2016)
  cat("Expected variation in shift across nodes (wi):", edgeshift_weight_sd_wi, "\n")
  cat("Observed variation:",   mean(repdevsds), "\n\n")
 
  simlist <- simulate_clingraph(g_positive, npergroup=100, edgenoise=0.4) #no group structure
  
  #convert to nconnections x nreps matrix
  weightmat <- sapply(simlist[["all"]], function(rep) { E(rep)$weight })
  
  edgesds <- apply(weightmat, 1, sd)
  
  cat("Expected edge noise (additive for each node): 0.4\n")
  cat("Observed edge noise:", mean(edgesds), "\n")
}