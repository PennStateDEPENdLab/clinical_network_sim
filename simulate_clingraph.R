#setup for project (just get stuff loaded)
suppressMessages(library(dplyr))
suppressMessages(library(igraph))

gsetupComplete <- FALSE
if (!gsetupComplete) {
  adjmat <- as.matrix(read.table("10895_bb264_corrmat.corr.1D"))
  adjmat_ztrans <- as.matrix(read.table("10895_bb264_corrmat.zval.1D")) #from AFNI
  
  # zval was computed from @ROI_Corr_Mat
  # 3dcalc   -a ${init_path}/${outfile}.corr+orig \
  #         -expr "min(0.5*log((1+a)/max(1-a,0.00001))*sqrt(${N}-3),15)" \
  #         -prefix ${init_path}/${outfile}.zval  >& /dev/null
  #
  # So the calc is: 0.5*log((1 + r)/(1 - r)) [usual fisher r->z] multiplied by sqrt(N - 3), which appears to be the SE of the statistic
  # why they multiply by SE instead of divide is beyond me
  
  #fisher r->z (aka atanh)
  
  #range(adjmat_noz[lower.tri(adjmat_noz)])
  #range(atanh(adjmat_noz[lower.tri(adjmat_noz)]))
  #range(adjmat[lower.tri(adjmat)])
  
  #cor(atanh(adjmat_noz[lower.tri(adjmat_noz)]),adjmat[lower.tri(adjmat)])
  #range(atanh(adjmat_noz[lower.tri(adjmat_noz)]))
  #head(cbind(atanh(adjmat_noz[lower.tri(adjmat_noz)]),adjmat[lower.tri(adjmat)]), n=50)
  
  #plot(cbind(atanh(adjmat_noz[lower.tri(adjmat_noz)]),adjmat[lower.tri(adjmat)]))
  adjmat <- atanh(adjmat)
  
  g_orig <- graph.adjacency(adjmat, mode="undirected", weighted=TRUE, diag=FALSE)
  g_positive <- delete.edges(g_orig, which(E(g_orig)$weight < 0)) #remove edges with negative weights 
  
  #mean z-val
  mean(E(g_orig)$weight)
  mean(E(g_positive)$weight)
  #hist(E(g_positive)$weight)
  
  #read in and support modularity
  mod <- gdata::read.xls(file.path(getMainDir(), "clinical_brainnetwork_sim", "Consensus264.xls"), header=TRUE, skip=1)
  mod <- mod %>% select(ROI, X.7, X.12) %>% rename(modnum=X.7, modlabel=X.12) 
  V(g_positive)$module <- mod$modnum
  V(g_positive)$module_name <- mod$modlabel
  
  V(g_orig)$module <- mod$modnum
  V(g_orig)$module_name <- mod$modlabel
  
  densities <- seq(.01, .20, by=.01)
  
  #threshold top 1-20% of connections
  #these will be positive, but it is not the top 1% of *positive* connections, for example
  #rather it's 264*263/2*.01
  allg_density <- lapply(densities, function(d, g) {
    thresh <- quantile(E(g)$weight, 1-d)
    
    gthresh <- delete.edges(g, which(E(g)$weight <= thresh)) #return binary weighted density-thresholded graph
    gthresh <- remove.edge.attribute(gthresh, "weight")
    return(gthresh)
  }, g_orig)
  
  names(allg_density) <- paste0(densities*100, "%")
  
  #lapply(allg_density, function(g_orig) { length(E(g_orig)) })
  gsetupComplete <- TRUE
}


#nodal graph measure worker function
gmeasures <- function(g) {
  deg <- degree(g)
  evc <- evcent(g)$vector
  s <- strength(g)
  bet <- betweenness(g, normalize=TRUE)*100
  locclust <- transitivity(g, type="local")
  pr <- page.rank(g, algo="prpack")$vector
  d <- graph.density(g) #add observed density (graph-level) to nodal stats for further analysis
  if (!is.null(E(g)$weight)) {
    meanFC <- mean(E(g)$weight, na.rm=TRUE)
  } else {
    meanFC <- NA
  }
  
  #return a 264 x nmetrics matrix for analysis/reduction            
  data.frame(node=V(g)$name, density=d, degree=deg, evcent=evc, betweenness=bet,
             locclust=locclust, pagerank=pr, strength=s, meanFC=meanFC)
}

#run the nodal measures function on a three-level graph list: densities, groups, graphs
nodal_measures_glist <- function(threshlist, ncores=4) {
  lapply(threshlist, function(grouplist) {
    lapply(grouplist, function(group) {
      mclapply(group, function(gthresh) {
        gmeasures(gthresh)
      }, mc.cores=ncores)
    })      
  })
} 

summarize_nodal_measures <- function(nlist) {
  suppressMessages(require(dplyr))
  thresh_summaries <- lapply(nlist, function(grouplist) {
    groupsummaries <- lapply(grouplist, function(group) {
      #here, group is a list of subjects where each is a 264 x measures data.frame
      recode <- lapply(1:length(group), function(g) {
        group[[g]]$replication <- g #add replication
        group[[g]]$nodeNum <- 1:nrow(group[[g]]) #add node number
        return(group[[g]])
      })            
      df <- do.call(rbind, recode)
      df %>% group_by(node) %>% summarize_at(vars(degree, evcent, betweenness, locclust, pagerank), funs(mean=mean, sd=sd))
    })
    groupsummaries <- lapply(1:length(grouplist), function(G) { groupsummaries[[G]]$group <- names(grouplist)[G]; return(groupsummaries[[G]])} )
    groupsummaries <- bind_rows(groupsummaries)
  })
  
  thresh_summaries <- lapply(1:length(thresh_summaries), function(D) { thresh_summaries[[D]]$thresh <- names(nlist)[D]; return(thresh_summaries[[D]])} )
  bind_rows(thresh_summaries) #return a single data frame
}

bind_nodal_measures <- function(nlist) { 
  by_threshold <- lapply(nlist, function(rlist) {
    ll <- lapply(rlist, function(group) {
      #here, group is a list of subjects where each is a 264 x measures data.frame
      recode <- lapply(1:length(group), function(g) {
        group[[g]]$replication <- g #add replication
        group[[g]]$nodeNum <- 1:nrow(group[[g]]) #add replication
        return(group[[g]])
      })            
      df <- do.call(rbind, recode)
      return(df)
    })
    
    ll$controls$group <- "controls"
    ll$patients$group <- "patients"
    return(do.call(rbind, ll))
  })
  
  by_threshold <- lapply(1:length(by_threshold), function(D) { by_threshold[[D]]$thresh <- names(nlist)[D]; return(by_threshold[[D]])} )
  bind_rows(by_threshold)
  
}

#operates on a data.frame created by bind_nodal_measures
group_lm <- function(nodal_combined, f=formula(degree ~ group), analyze_nodes=NULL) {
  if (!is.null(analyze_nodes)) {
    #filter to just a subset of nodes
    nodal_combined <- filter(nodal_combined, node %in% analyze_nodes)
  }
  
  effects <- nodal_combined %>% group_by(thresh, nodeNum) %>%
    do({
      mm <- lm(f, .)
      partab <- tidy(mm)
      pfilt <- filter(partab, term=="grouppatients")
      pfilt
    }) %>% ungroup()
  return(effects)
}

#add a nodeType to the data.frame (e.g., output of group_lm) to differentiate positive and negative targets
flag_nodes <- function(nodal_combined, pos, neg) {
  if (is.character(pos)) { pos <- as.numeric(sub("V", "", pos, fixed=TRUE))} #tolerate V<[0-9]> syntax
  if (is.character(neg)) { neg <- as.numeric(sub("V", "", neg, fixed=TRUE))}
  nodal_combined %>% mutate(nodeType = case_when(
    nodeNum %in% pos ~ "Positive Target",
    nodeNum %in% neg ~ "Negative Target",
    TRUE ~ "Not Targeted"
  ))
}

#simple function to apply thresholding to a weighted graph, either based on density or strength threshold
threshold_graph <- function(g, w, method="density", rmweights=TRUE) {
  stopifnot(is_igraph(g))
  stopifnot(is.numeric(w))
  
  #Obtains desired density given graph diameter
  if (method=="density") {
    weights <- sort(E(g)$weight, decreasing=TRUE)
    threshold <- weights[length(V(g))*(length(V(g))-1)/2 * w]
  } else {
    threshold <- w #just use strength-based threshold
  }
  
  gthresh <- delete.edges(g, which(E(g)$weight < threshold))
  if (rmweights) {
    E(gthresh)$w_hidden <- E(gthresh)$weight #retain weight attribute for additional analyses (but remove "weight" to avoid igraph thinking it's a weighted graph)
    gthresh <- remove.edge.attribute(gthresh, "weight")
  }
  gthresh$wthresh <- threshold #copy threshold weight into object for tracking  
  if (method=="density") { gthresh$target_density <- w } #copy density into object for tracking  
  return(gthresh)  
}

threshold_glist <- function(glist, thresholds, method="density", ncores=4, ...) {
  allg_thresh <- lapply(thresholds, function(t, glist) {
    lapply(glist, function(group) {
      mclapply(group, function(g) {
        return(threshold_graph(g, t, method=method, ...))
      }, mc.cores=ncores)
    })
  }, glist)
  
  names(allg_thresh) <- as.character(paste(method, thresholds, sep="_"))
  return(allg_thresh)
}

global.measures <- function(g) {
  suppressMessages(require(igraph))
  #assumes binary
  Cglob <- 
    Lglob <-  #average path
    Lglob <- mean_distance(g) #average path
  
  data.frame(
    Cglob=transitivity(g, type="global"), #global clustering/transitivity
    Lglob=mean_distance(g), #average shortest path
    Dcent=centr_degree(g)$centralization, #global centralization according to degree
    Dens=edge_density(g)
  )
}

#adaptation of brainGraph within_module_deg_z_score function
#mine is clunkier, but easier to read
#it also returns the within and between degree, as well as between z stats.
#upon further testing, this is painfully slow
#i've adapted the igraph code below to retain the matrix multiplication approach to get degree
#wibw_module_degree_slow <- function(gobj) {
#  nnodes <- length(V(gobj))
#  modules <- sort(unique(V(gobj)$module))
#  #Ki = number of edges from node i to nodes in the same module
#  #Ksi = average number of edges from one node to another in a module
#  #sigmaKsi = standard deviation of within-module degree distribution
#  Ki <- Ksi <- sigmaKsi <- rep(0, nnodes)
#  Bi <- Bsi <- sigmaBsi <- rep(0, nnodes) #same for between
#  module <- rep(0, nnodes) #module tracking
#  for (m in modules) {
#    ncomm <- V(gobj)[V(gobj)$module == m]
#    module[V(gobj)$module == m] <- m #super redundant, but safe! :-)
#    othernodes <- V(gobj)[difference(V(gobj), ncomm)] #all others not in ncomm
#    #all_edges <- E(gobj)[inc(ncomm)]
#    #within_edges <- E(gobj)[ncomm %--% ncomm] ## within nodes of this module
#    #between_edges <- difference(all_edges, within_edges) ## the rest
#    
#    for (n in 1:length(ncomm)) {
#      Ki[ncomm[n]] <- length(E(gobj)[ncomm[n] %--% ncomm])
#      Bi[ncomm[n]] <- length(E(gobj)[ncomm[n] %--% othernodes])
#    }
#    
#    Ksi[ncomm] <- mean(Ki[ncomm])
#    sigmaKsi[ncomm] <- sd(Ki[ncomm])
#    Bsi[ncomm] <- mean(Bi[ncomm])
#    sigmaBsi[ncomm] <- sd(Bi[ncomm])
#  }
#  z_within <- (Ki - Ksi)/sigmaKsi
#  z_between <- (Bi - Bsi)/sigmaBsi
#  return(data.frame(module=module, Ki=Ki, Ksi=Ksi, sigmaKsi=sigmaKsi, z_within=z_within, Bi=Bi, Bsi=Bsi, sigmaBsi=sigmaBsi, z_between=z_between))
#}

#adapted function from igraph
wibw_module_degree <- function (g) {
  stopifnot(is_igraph(g))
  memb <- V(g)$module
  modules <- unique(memb)
  N <- length(modules)
  A <- as_adj(g, sparse = FALSE, names = FALSE)
  z_within <- Ki <- rep(0, nrow(A))
  z_between <- Bi <- rep(0, nrow(A))
  Ksi <- sigKsi <- rep(0, nrow(A))
  Bsi <- sigBsi <- rep(0, nrow(A))
  for (S in modules) {
    x_wi <- A[memb == S, ] %*% (memb == S) #from S to S
    x_bw <- A[memb == S, ] %*% (memb != S) #from S to other nodes
    Ki[memb == S] <- x_wi; Ksi[memb==S] <- mean(x_wi); sigKsi[memb==S] <- sd(x_wi)
    Bi[memb == S] <- x_bw; Bsi[memb==S] <- mean(x_bw); sigBsi[memb==S] <- sd(x_bw)
  }
  z_within <- (Ki - Ksi)/sigKsi
  z_within <- ifelse(!is.finite(z_within), 0, z_within)
  z_between <- (Bi - Bsi)/sigBsi
  z_between <- ifelse(!is.finite(z_between), 0, z_between)
  return(data.frame(node=names(V(g)), module=memb, Ki=Ki, Ksi=Ksi, sigKsi=sigKsi, z_within=z_within, Bi=Bi, Bsi=Bsi, sigBsi=sigBsi, z_between=z_between))
}


#ex-Gaussian simulation function that accepts a total SD parameter, sigma,
#which is divvied up between the Gaussian and exponential.
#by default, the data are also recentered at zero.
rexgauss <- function(n, mu, sigma, beta, recenter=0) {
  #beta is the 'survival parameter' that is 1/tau (rate)
  #beta defines both the mean and standard deviation for an exponential distribution
  #exgauss: rnorm(n, mu, sigma) + rexp(n, 1/tau)
  
  #variance of exponential is beta^2 (aka 1/tau^2)
  #mean of exponential is beta (aka 1/tau)
  #here we've specified tau in terms of the target mean, not the "rate"
  #hence, the sd of the exponential piece is technically sqrt(beta) (aka 1/((1/tau)^2)) ) since we invert tau to get back to a rate
  if (sigma^2 - beta^2 < 0) { stop("Cannot simulate rexgauss data with these settings. Negative Gaussian variance") }
  sigma_gauss <- sqrt(sigma^2 - beta^2) #here the sd of the Gaussian part is the sqrt of the overall variance minus the variance of the exponential part
  vec <- rnorm(n, mu, sigma_gauss) + rexp(n, 1/beta)
  if (!is.null(recenter)) { 
    vec <- vec - mean(vec) + recenter #recenter at target
  }  
  return(vec)
}

#rexgauss parameter testing
#n <- 10000
#sigma <- .2
#mu <- 5
#beta <- .1 #sd of exp is equal to beta, and sd of gaussian is sigma^2 - beta^2
#recenter <- 0
#
##the ratio of sigma / beta determines how exp versus gaussian it is
##to the extent that sigma is much larger than beta, the distribution will be mostly gaussian
##mu is not a parameter of interest in the shape because it just shifts around the location of the gaussian part
#
#x <- rexgauss(n, mu, sigma, beta, recenter=NULL)
#
#hist(x)
#moments::skewness(x)
#
#sigma_vec <- seq(0.1, 1, by=0.1)
#sbratio_vec <- seq(1, 5, by=0.1)
#
##because ex Gaussian adds an exponential variate onto a Gaussian variate, it always adds positive values onto the Gaussian, regardless of mean
##NB: an exponential distribution has skewness of ~2.
#
#m <- array(NA, dim=c(length(sigma_vec), length(sbratio_vec), n))
#skew_arr <- array(NA, dim=c(length(sigma_vec), length(sbratio_vec)), dimnames=list(sigma=sigma_vec, sbratio=sbratio_vec))
#for (i in 1:length(sigma_vec)) {
#  for (j in 1:length(sbratio_vec)) {
#    m[i,j,] <- rexgauss(n, mu=0, sigma=sigma_vec[i], beta=sigma_vec[i]/sbratio_vec[j], recenter=0)
#    skew_arr[i,j] <- moments::skewness(m[i,j,])
#  }
#}
#
#melted <- reshape2::melt(m)
#setwd("~/")
#library(ggplot2)
#library(dplyr)
#pdf("rexpgauss_test.pdf", width=20, height=15)
##ggplot(melted, aes(x=value)) + geom_histogram() + facet_grid(Var1 ~ Var2, scales="free")
#ggplot(filter(melted, Var1==5), aes(x=value)) + geom_histogram() + facet_wrap(~ Var2, scales="free")
#dev.off()
#
##lattice::histogram(m[1,1,])
#hist(m[1,3,], breaks=20)
#moments::skewness(m[1,3,])

#conclusion: regardless of overall sigma, sbratio of 1 is pure exponential, 1.1 is skewness 1.5, 1.2 is skewness ~1.2.
#for now, go with 1.2 ratio since this gives a reasonably long positive tail without shifting the center of the distribution far from
#zero.

#main worker function for simulating two-group (or more) graph structure
simulate_clingraph <- function(gobj, npergroup=100, compute_brainGraph=FALSE, groupmanipulation=NULL, 
                               edgenoise=list(sd=0, dist="gaussian"), proportional=FALSE, ncpus=4, partial=FALSE, seed=NULL) {
  # groupmanipulation syntax/structure: list of lists
  # [[1]]
  #   $edgeshift_weight_bw_mean: mean shift in weights across edges for this group
  #   $edgeshift_weight_bw_sd: variation in magnitude of weight shift across replications/subjects
  #   $edgeshift_weight_wi_sd: variation in magnitude of weight shift across nodes within subject
  #   $group_name: name of the group
  #   $nodemanipulation: list of nodes to be manipulated
  #
  # nodemanipulation structure: list of nodes (embedded within groupmanipulation)
  # [[1]]
  #   $name: name of node: must match vertex names in graph
  #   $edgeshift_weight_bw_mean: nodal mean shift in edge weights
  #   $edgeshift_weight_bw_sd: nodal variation in magnitude of weight shift across replications
  #   $edgeshift_weight_wi_sd: nodal variation in magnitude of weight shift across edges connected to this node
  
  suppressMessages(require(foreach))
  #require(doSNOW)
  suppressMessages(require(doParallel))
  #suppressMessages(require(doRNG))
  suppressMessages(require(corpcor))
  
  #setDefaultClusterOptions(master="localhost")
  if (ncpus > 1) {
    clusterobj <- makePSOCKcluster(ncpus) #4 cpus
    #registerDoSNOW(clusterobj)
    registerDoParallel(clusterobj)
    #  if (!is.null(seed)) {
    #    registerDoRNG(seed=seed)    
    #  }
    cat("parallel")
    on.exit(try(stopCluster(clusterobj)))    
  } else {
    cat("sequential")
    registerDoSEQ()
  }
  
  weighted <- FALSE
  if ("weight" %in% edge_attr_names(gobj)) { weighted = TRUE }
  
  stopifnot(weighted==TRUE) #currently working on weighted variant
  
  #handle empty group structure (single group sim): basically just useful to look at edge noise across replications
  #expected group structure
  if (is.null(groupmanipulation)) {
    groupmanipulation <- list(list(group_name="all", edgeshift_weight_bw_mean=0, edgeshift_weight_wi_sd=0, edgeshift_weight_bw_sd=0))    
  }
  
  #node manipulations are treated as property of group specification
  which_nodemanip <- sapply(groupmanipulation, function(group) { !is.null(group$nodemanipulation) })
  #has_nodemanip <- any(which_nodemanip)
  
  allsims <- list()
  #build sims
  for (g_iter in 1:length(groupmanipulation)) {
    #check noise model for weight shifts between subjects: default Gaussian
    if (is.null(groupmanipulation[[g_iter]]$edgeshift_weight_dist_bw)) {
      groupmanipulation[[g_iter]]$edgeshift_weight_dist_bw <- "gaussian"
    }
    
    if (groupmanipulation[[g_iter]]$edgeshift_weight_dist_bw == "gaussian") {
      #mean shift + between-rep variation
      #Gaussian variation in strength of shift across replications for a group
      repshift <- rnorm(npergroup, groupmanipulation[[g_iter]]$edgeshift_weight_bw_mean, groupmanipulation[[g_iter]]$edgeshift_weight_bw_sd) #mean shift + between-rep variation
    } else if (groupmanipulation[[g_iter]]$edgeshift_weight_dist_bw == "ex-gaussian") {
      repshift <- rexgauss(n=npergroup, mu=0, sigma=groupmanipulation[[g_iter]]$edgeshift_weight_bw_sd, 
                           beta=groupmanipulation[[g_iter]]$edgeshift_weight_bw_sd/groupmanipulation[[g_iter]]$edgeshift_weight_bw_sbratio, #sigma/beta ratio
                           recenter=groupmanipulation[[g_iter]]$edgeshift_weight_bw_mean) #reset mean to target
    }
    
    groupsims <- foreach(j = iter(1:npergroup), .packages="igraph", .export=c("rexgauss")) %dopar% {
      g_rep <- gobj
      
      #simulate (ex-)Gaussian variation on edges 
      if (weighted) {
        w <- E(gobj)$weight #vector of master weights
        
        if (edgenoise$sd > 0) {
          #add Gaussian variation in edge strength over replications
          if (edgenoise$dist == "gaussian") { 
            enoise <- rnorm(length(w), 0, edgenoise$sd)
          } else if (edgenoise$dist == "ex-gaussian") {
            enoise <- rexgauss(length(w), mu=0, sigma=edgenoise$sd, beta=edgenoise$sd/edgenoise$sbratio, recenter=0)
          }
          w <- w + enoise
          E(g_rep)$edgenoise <- enoise #store noise as edge attribute
        } else {
          E(g_rep)$edgenoise <- 0 #no edge noise
        }
        
        #add within-replication variation with mean = between and variation around that mean = within
        if (is.null(groupmanipulation[[g_iter]]$edgeshift_weight_dist_wi)) {
          groupmanipulation[[g_iter]]$edgeshift_weight_dist_wi <- "gaussian" #default to Gaussian variation 
        }
        
        if (groupmanipulation[[g_iter]]$edgeshift_weight_dist_wi == "gaussian") {
          w <- w + rnorm(length(w), repshift[j], groupmanipulation[[g_iter]]$edgeshift_weight_wi_sd)
        } else if (groupmanipulation[[g_iter]]$edgeshift_weight_dist_wi == "ex-gaussian") {
          w <- w + rexgauss(length(w), mu=0, sigma=groupmanipulation[[g_iter]]$edgeshift_weight_wi_sd, 
                            beta=groupmanipulation[[g_iter]]$edgeshift_weight_wi_sd/groupmanipulation[[g_iter]]$edgeshift_weight_bw_sbratio, recenter=repshift[j])
        }
        E(g_rep)$weight <- w #weights with shift and noise  
      }
      
      g_rep #return
    }
    
    edgesmanipulated <- c() #running tally of edges that have already been shifted
    
    #only apply node manipulations if specified for this group
    if (which_nodemanip[g_iter]) {      
      # setup an edgelist of shifts
      
      cat("  NB: Shifts in edge strength are applied in order of the node list.\n")
      cat("  Furthermore, once a shift has been applied to an edge, that edge is immune to further shifts (no compounding)\n")
      cat("  This also applies to network/community shifts, which fall after this in the pipeline\n\n")
      
      for (node in groupmanipulation[[g_iter]]$nodemanipulation) {
        #node specification includes mean shift, variation in shift between replications (bw), and variation in shift within incident edges (wi)
        #obtain edge sequence for connections incident to this node
        node_edges <- E(gobj)[inc(node$name)]
        node_edges <- difference(node_edges, edgesmanipulated) #pull out any edges already tweaked
        edgesmanipulated <- c(edgesmanipulated, node_edges) #add edges to edge sequence for this node
        
        if (is.null(node$edgeshift_weight_bw_dist)) { node$edgeshift_weight_bw_dist <- "gaussian" }
        if (is.null(node$edgeshift_weight_wi_dist)) { node$edgeshift_weight_wi_dist <- "gaussian" }
        if (node$edgeshift_weight_bw_dist == "gaussian") {
          noderepshift <- rnorm(npergroup, node$edgeshift_weight_bw_mean, node$edgeshift_weight_bw_sd) #mean shift + between-rep variation  
        } else if (node$edgeshift_weight_bw_dist == "ex-gaussian") {
          noderepshift <- rexgauss(npergroup, mu=0, recenter=node$edgeshift_weight_bw_mean, 
                                   sigma=node$edgeshift_weight_bw_sd, beta=node$edgeshift_weight_bw_sd/node$edgeshift_weight_bw_sbratio) #mean shift + between-rep variation
        }
        
        #tweak this node for each replication (loop over replications)
        groupsims <- lapply(1:npergroup, function(j) {
          #add Gaussian variation with mean = replication between and variation around that mean = within
          if (node$edgeshift_weight_wi_dist == "gaussian") {
            E(groupsims[[j]])[node_edges]$weight <- E(groupsims[[j]])[node_edges]$weight + rnorm(length(node_edges), noderepshift[j], node$edgeshift_weight_wi_sd)
          } else if (node$edgeshift_weight_wi_dist == "ex-gaussian") {
            E(groupsims[[j]])[node_edges]$weight <- E(groupsims[[j]])[node_edges]$weight + 
              rexgauss(length(node_edges), mu=0, recenter=noderepshift[j], sigma=node$edgeshift_weight_wi_sd, beta=node$edgeshift_weight_wi_sd/node$edgeshift_weight_wi_sbratio)
          }
          return(groupsims[[j]])
        })
      }
    }
    #for (node in groupmanipulation[[g_iter]]$nodemanipulation) {
    for (network in groupmanipulation[[g_iter]]$networkmanipulation) {
      #network specification includes:
      # a) between-network variability
      #   1) mean shift between replications *between* networks (edgeshift_weight_bw_bwn_mean),
      #   2) variation in mean shift between replications *between* networks (edgeshift_weight_bw_bwn_sd),
      #   3) variation in shift within a given replication for *between*-network connections (edgeshift_weight_wi_bwn_sd)
      # b) within-network variability
      #   4) mean shift between replications *within* networks (edgeshift_weight_bw_win_mean),
      #   5) variation in mean shift between replications *between* networks (edgeshift_weight_bw_win_sd),
      #   6) variation in shift within a given replication for *between*-network connections (edgeshift_weight_wi_win_sd)
      
      if (is.null(network$edgeshift_weight_bw_bwn_dist)) { network$edgeshift_weight_bw_bwn_dist <- "gaussian" }
      if (is.null(network$edgeshift_weight_wi_bwn_dist)) { network$edgeshift_weight_wi_bwn_dist <- "gaussian" }
      if (is.null(network$edgeshift_weight_bw_win_dist)) { network$edgeshift_weight_bw_win_dist <- "gaussian" }
      if (is.null(network$edgeshift_weight_wi_win_dist)) { network$edgeshift_weight_wi_win_dist <- "gaussian" }
      
      #nodes in the community
      ncomm <- V(gobj)[V(gobj)$module %in% network$module] #allow $module to be vector
      
      all_edges <- E(gobj)[inc(ncomm)]
      all_edges <- difference(all_edges, edgesmanipulated) #pull out nodal manipulations previously applied
      within_edges <- E(gobj)[ncomm %--% ncomm] # %--% syntax
      between_edges <- difference(all_edges, within_edges)
      
      #handle between-network connectivity modulation
      #compute distribution of between network variation across replications
      if (network$edgeshift_weight_bw_bwn_dist == "gaussian") {
        bwn_repshift <- rnorm(npergroup, network$edgeshift_weight_bw_bwn_mean, network$edgeshift_weight_bw_bwn_sd) #mean shift + between-rep variation in *between network* strength 
      } else if (network$edgeshift_weight_bw_bwn_dist == "ex-gaussian") {
        bwn_repshift <- rexgauss(npergroup, mu=0, recenter=network$edgeshift_weight_bw_bwn_mean, 
                                 sigma=network$edgeshift_weight_bw_bwn_sd, beta=network$edgeshift_weight_bw_bwn_sd/network$edgeshift_weight_bw_bwn_sbratio) #mean shift + between-rep variation
      }
      
      #add variability in between-network connectivity
      groupsims <- lapply(1:npergroup, function(j) {
        #add Gaussian variation with mean = replication between and variation around that mean = within
        if (network$edgeshift_weight_wi_bwn_dist == "gaussian") {
          E(groupsims[[j]])[between_edges]$weight <- E(groupsims[[j]])[between_edges]$weight + rnorm(length(between_edges), bwn_repshift[j], network$edgeshift_weight_wi_bwn_sd)
        } else if (network$edgeshift_weight_wi_bwn_dist == "ex-gaussian") {
          E(groupsims[[j]])[between_edges]$weight <- E(groupsims[[j]])[between_edges]$weight + 
            rexgauss(length(between_edges), mu=0, recenter=bwn_repshift[j], sigma=network$edgeshift_weight_wi_bwn_sd, beta=network$edgeshift_weight_wi_bwn_sd/network$edgeshift_weight_wi_bwn_sbratio)
        }
        return(groupsims[[j]])
      })
      
      #handle within-network connectivity modulation
      #compute distribution of between network variation across replications
      if (network$edgeshift_weight_bw_win_dist == "gaussian") {
        win_repshift <- rnorm(npergroup, network$edgeshift_weight_bw_win_mean, network$edgeshift_weight_bw_win_sd) #mean shift + between-rep variation in *between network* strength 
      } else if (network$edgeshift_weight_bw_win_dist == "ex-gaussian") {
        win_repshift <- rexgauss(npergroup, mu=0, recenter=network$edgeshift_weight_bw_win_mean, 
                                 sigma=network$edgeshift_weight_bw_win_sd, beta=network$edgeshift_weight_bw_win_sd/network$edgeshift_weight_bw_win_sbratio) #mean shift + between-rep variation
      }
      
      #add variability in within-network connectivity
      groupsims <- lapply(1:npergroup, function(j) {
        #add Gaussian variation with mean = replication between and variation around that mean = within
        if (network$edgeshift_weight_wi_win_dist == "gaussian") {
          E(groupsims[[j]])[within_edges]$weight <- E(groupsims[[j]])[within_edges]$weight + rnorm(length(within_edges), win_repshift[j], network$edgeshift_weight_wi_win_sd)
        } else if (network$edgeshift_weight_wi_win_dist == "ex-gaussian") {
          E(groupsims[[j]])[within_edges]$weight <- E(groupsims[[j]])[within_edges]$weight + 
            rexgauss(length(within_edges), mu=0, recenter=win_repshift[j], sigma=network$edgeshift_weight_wi_win_sd, 
                     beta=network$edgeshift_weight_wi_win_sd/network$edgeshift_weight_wi_win_sbratio)
        }
        return(groupsims[[j]])
      })
      
    }
    
    #handle conversion to partial correlation matrix
    #because RS data are rank deficient, use 
    if (partial) {
      groupsims <- foreach(g=iter(groupsims), .packages=c("corpcor", "Matrix")) %dopar% {
        #groupsims <- lapply(groupsims, function(g) {
        #transform back into correlation metric for cor -> pcor
        #as_adj is much slower than the matrix subset approach below
        #amat <- tanh(as_adj(groupsims[[1]], sparse=FALSE, names=FALSE, attr="weight"))
        #diag(vv) <- 1
        
        vv2 <- tanh(g[]) #use subsetting syntax to force conversion to matrix
        diag(vv2) <- 1 #put 1 back on diagonal for inversion
        
        #problems with non-positive definite matrix
        #use the nearPD algorithm to compute the nearest positive definite matrix
        #eigen(vv2)$values #so many negative!!
        vpos <- as.matrix(Matrix::nearPD(vv2, corr=TRUE)$mat) #$mat has the raw matrix
        
        #m <- -pseudoinverse(vpos$mat)
        #diag(m) <- -diag(m)
        #xx2 <- cov2cor(m)
        
        vpartial <- corpcor::cor2pcor(vpos) #compute this on the positive definite matrix (otherwise will blow up)
        
        #painfully slow
        #adjmat <- EBICglasso(S=vv2, n=550, gamma=0.5, penalize.diagonal=FALSE, nlambda=100,
        #    lambda.min.ratio = 0.01, returnAllResults = FALSE, checkPD = FALSE, 
        #    countDiagonal = FALSE, refit = FALSE)
        #cor(as.vector(vpos), as.vector(vpartial)) #correlation of partial matrix with original
        #cor(as.vector(vpos), as.vector(vv2)) #correlation of pos def matrix with original
        
        #repopulate edge weights with partial correlation
        #g2 <- g
        #ee <- E(g2)$weight #test that edges can be updated by lower.tri (i.e., that edge list is in that order)
        #ee2 <- g2[]
        #ee2 <- ee2[lower.tri(ee2)]
        #identical(ee, ee2) #yes, we can populate the graph using lower triangle of the partial correlation matrix
        #E(g2)$weight <- vpartial[lower.tri(vpartial)]
        #g3 <- g2[]
        #diag(g3) <- 1
        #all.equal(vpartial, as.matrix(g3))
        #range(vpartial - as.matrix(g3)) #identical within numerical precision
        
        E(g)$weight <- vpartial[lower.tri(vpartial)]
        return(g)
      }
    }
    
    groupsims <- lapply(groupsims, function(g) { g$group <- groupmanipulation[[g_iter]]$group_name; return(g) }) #add global group attribute
    allsims[[ groupmanipulation[[g_iter]]$group_name ]] <- groupsims
  }
  
  return(allsims)
  
}

sanitycheck_sims <- function(g_positive) {
  #check that variation in graphs is as expected
  
  edgeshift_weight_wi_sd=.06
  edgeshift_weight_bw_mean=0.5
  edgeshift_weight_bw_sd=.25
  npergroup = 200
  
  simlist <- simulate_clingraph(g_positive, npergroup=npergroup, edgenoise=list(sd=0, dist="gaussian"), #need 0 edgenoise for sanity check ###edgenoise=list(sd=0, dist="ex-gaussian", sbratio=1.1), #need 0 edgenoise for sanity check 
                                groupmanipulation=list(
                                  list(group_name="patients", edgeshift_weight_bw_mean=edgeshift_weight_bw_mean,
                                       edgeshift_weight_wi_sd=edgeshift_weight_wi_sd, edgeshift_weight_bw_sd=edgeshift_weight_bw_sd)))
  
  cat("Running sanity checks on weight shifts\n\n")
  #sanity checks: are the means and SDs as expected for patient group?
  #these checks will only hold up under 0 edge noise (or for huge groups in the asymptote)
  patients <- simlist[["patients"]]
  
  #per-replication average deviation between ground_truth and average weight
  repdevmeans <- sapply(patients, function(replication) { mean(E(replication)$weight - E(g_positive)$weight) })
  
  #between-replication mean shift: should be close to edgeshift_weight_bw_mean in the asymptote (VERIFIED 30Sep2016)
  cat("  Expected mean shift in group: ", edgeshift_weight_bw_mean, "\n", sep="")
  cat("  Observed mean shift: ", round(mean(repdevmeans), 3), "\n\n", sep="")
  
  #across-replication variation in mean shift (sd of mean shifts): should be close to edgeshift_weight_bw_sd (VERIFIED 30Sep2016)
  cat("  Expected variation in mean shift between replications (bw): ", edgeshift_weight_bw_sd, "\n")
  cat("  Observed variation: ", round(sd(repdevmeans), 3), "\n\n")
  
  #per-replication (within) edge variation from ground truth (across nodes)
  repdevsds <- sapply(patients, function(replication) { sd(E(replication)$weight - E(g_positive)$weight) })
  
  #overall sd across nodes (within): this should be close to edgeshift_weight_wi_sd in the asymptote (VERIFIED 30Sep2016)
  cat("  Expected variation in shift across nodes within a replication (wi): ", edgeshift_weight_wi_sd, "\n", sep="")
  cat("  Observed variation: ", round(mean(repdevsds), 3), "\n\n", sep="")
  
  edgesd <- 0.4
  simlist <- simulate_clingraph(g_positive, npergroup=npergroup, edgenoise=list(sd=edgesd, dist="ex-gaussian", sbratio=1.1)) #edgenoise=list(sd=0.4, dist="gaussian")) #no group structure
  
  #verify distribution of added noise compared to ground truth
  #plot(hist(E(simlist[[1]][[1]])$weight - E(g_positive)$weight))
  #moments::skewness(E(simlist[[1]][[1]])$weight - E(g_positive)$weight)
  #hist(E(simlist[[1]][[1]])$weight)
  #moments::skewness(E(simlist[[1]][[1]])$edgenoise)
  #hist(E(g_positive)$weight)
  
  #convert to nconnections x nreps matrix
  weightmat <- sapply(simlist[["all"]], function(rep) { E(rep)$weight })
  
  edgesds <- apply(weightmat, 1, sd)
  
  cat("  Expected edge noise (additive for each node): ", edgesd, "\n", sep="")
  cat("  Observed edge noise: ", round(mean(edgesds), 3), "\n", sep="")
  
  cat("\n==================\n\n")
  cat("Testing node-level manipulations\n\n")
  
  #sanity check of nodal manipulations
  v1edgeMean <- 20; v1edgeSDWi <- 7; v1edgeSDbw <- 5
  V189edgeMean <- 100; V189edgeSDWi <- 25; V189edgeSDbw <- 5
  
  allsims <- simulate_clingraph(g_positive, npergroup=200, edgenoise=list(sd=0.0, dist="gaussian"), #just for sanity check
                                groupmanipulation=list(
                                  controls=list(group_name="control", 
                                                edgeshift_weight_bw_mean=0, edgeshift_weight_wi_sd=0, edgeshift_weight_bw_sd=0, #global shifts in weights
                                                nodemanipulation=list( #nodal shifts in (weighted) edges incident to a node
                                                  list(name="V1", edgeshift_weight_bw_mean=v1edgeMean, edgeshift_weight_wi_sd=v1edgeSDWi, edgeshift_weight_bw_sd=v1edgeSDbw),
                                                  list(name="V189", edgeshift_weight_bw_mean=V189edgeMean, edgeshift_weight_wi_sd=V189edgeSDWi, edgeshift_weight_bw_sd=V189edgeSDbw)
                                                )
                                  ),
                                  patients=list(group_name="patients", edgeshift_weight_bw_mean=0.5, edgeshift_weight_wi_sd=.06, edgeshift_weight_bw_sd=.25)))
  
  cat("  Shifting V1 by M =", v1edgeMean, ", SDwi = ", v1edgeSDWi, ", SDbw = ", v1edgeSDbw, "\n\n", sep="")
  
  controls <- allsims[["control"]]
  
  #per-replication average deviation between ground_truth and average weight
  v1repdevmeans <- sapply(controls, function(replication) { mean(E(replication)[inc("V1")]$weight - E(g_positive)[inc("V1")]$weight) })
  
  cat("  Expected mean shift in V1 edges: ", v1edgeMean, "\n", sep="")
  cat("  Observed mean shift in V1 edges across replications: ", round(mean(v1repdevmeans), 3), "\n\n", sep="")
  
  #across-replication variation in mean shift (sd of mean shifts): should be close to edgeshift_weight_bw_sd (VERIFIED 30Sep2016)
  cat("  Expected variation in shift across replications (bw):", v1edgeSDbw, "\n", sep="")
  cat("  Observed variation:", round(sd(v1repdevmeans), 3), "\n\n", sep="")
  
  #per-replication (within) edge variation from ground truth (across nodes)
  repdevsds <- sapply(controls, function(replication) { sd(E(replication)[inc("V1")]$weight - E(g_positive)[inc("V1")]$weight) })
  
  #overall sd across nodes (within): this should be close to edgeshift_weight_wi_sd in the asymptote (VERIFIED 30Sep2016)
  cat("  Expected variation in shift across nodes (wi):", v1edgeSDWi, "\n")
  cat("  Observed variation:", round(mean(repdevsds), 3), "\n\n")
  
  cat("  After V1, shifting V189 by M = ", V189edgeMean, ", SDwi = ", V189edgeSDWi, ", SDbw = ", V189edgeSDbw, "\n", sep="")
  
  #per-replication average deviation between ground_truth and average weight
  V189repdevmeans <- sapply(controls, function(replication) { mean(E(replication)[inc("V189")]$weight - E(g_positive)[inc("V189")]$weight) })
  
  cat("  Expected mean shift in V189 edges: ", V189edgeMean, "\n", sep="")
  cat("  Observed mean shift in V189 edges across replications: ", mean(V189repdevmeans), "\n\n", sep="")
  
  #across-replication variation in mean shift (sd of mean shifts): should be close to edgeshift_weight_bw_sd (VERIFIED 30Sep2016)
  cat("  Expected variation in shift across replications (bw):", V189edgeSDbw, "\n")
  cat("  Observed variation:", sd(V189repdevmeans), "\n\n")
  
  #per-replication (within) edge variation from ground truth (across nodes)
  V189repdevsds <- sapply(controls, function(replication) { sd(E(replication)[inc("V189")]$weight - E(g_positive)[inc("V189")]$weight) })
  
  #overall sd across nodes (within): this should be close to edgeshift_weight_wi_sd in the asymptote (VERIFIED 30Sep2016)
  cat("  Expected variation in shift across nodes (wi): ", V189edgeSDWi, "\n", sep="")
  cat("  Observed variation: ", round(mean(V189repdevsds), 3), "\n\n", sep="")
  
  cat("  Checking that the V1 -- V189 connection is not dialed up to the V189 level since it was tweaked first in the sequence\n")
  cat("  Change should be something on the order of 20-25, not 100-150\n")
  
  v1_V189 <- mean(sapply(controls, function(replication) { replication["V1","V189"] })) #even more compact subset notation (igraph internally treats this as a weighted adjacency matrix)
  cat("  Observed mean connection between V1 and V189: ", v1_V189, "\n\n")
  
  ###TEST NETWORK/MODULE LEVEL MANIPULATIONS
  cat("\n==================\n\n")
  cat("Testing network-level manipulations\n\n")
  
  edgeshift_weight_bw_bwn_mean=10
  edgeshift_weight_bw_bwn_sd=.5
  edgeshift_weight_wi_bwn_sd=.7
  
  edgeshift_weight_bw_win_mean=50
  edgeshift_weight_bw_win_sd=5
  edgeshift_weight_wi_win_sd=12
  
  testcommunity <- simulate_clingraph(g_positive, npergroup=100, edgenoise=list(sd=0.0, dist="gaussian"),
                                      groupmanipulation=list(
                                        controls=list(group_name="controls", 
                                                      edgeshift_weight_bw_mean=0, edgeshift_weight_wi_sd=0, edgeshift_weight_bw_sd=0
                                        ),         #global shifts in weights
                                        patients=list(group_name="patients", edgeshift_weight_bw_mean=0.0, edgeshift_weight_wi_sd=0, edgeshift_weight_bw_sd=0, #global shifts in weights
                                                      networkmanipulation=list( #network-level shifts in (weighted) edges between and within a network
                                                        list(module=1, edgeshift_weight_bw_bwn_mean=edgeshift_weight_bw_bwn_mean, edgeshift_weight_bw_bwn_sd=edgeshift_weight_bw_bwn_sd, edgeshift_weight_wi_bwn_sd=edgeshift_weight_wi_bwn_sd,
                                                             edgeshift_weight_bw_win_mean=edgeshift_weight_bw_win_mean, edgeshift_weight_bw_win_sd=edgeshift_weight_bw_win_sd, edgeshift_weight_wi_win_sd=edgeshift_weight_wi_win_sd)
                                                      )
                                        )
                                      )
  )
  
  patients <- testcommunity[["patients"]]
  
  #per-replication average deviation between ground_truth and average weight
  m1nodes <- V(g_positive)[V(g_positive)$module == 1]
  
  all_edges <- E(g_positive)[inc(m1nodes)]
  #all_edges <- difference(all_edges, edgesmanipulated) #pull out nodal manipulations previously applied
  within_edges <- E(g_positive)[m1nodes %--% m1nodes] # %--% syntax denotes edges between one vertex set and another
  between_edges <- difference(all_edges, within_edges)
  
  #mean shift in between-network connectivity for module 1 across replications (VERIFIED 15Apr2017)
  m1bwn_meanshift <- sapply(patients, function(replication) { mean(E(replication)[between_edges]$weight - E(g_positive)[between_edges]$weight) })
  
  cat("  Expected mean shift in between-network connectivity for module 1 (edgeshift_weight_bw_bwn_mean): ", edgeshift_weight_bw_bwn_mean, "\n", sep="")
  cat("  Observed variation: ", round(mean(m1bwn_meanshift), 3), "\n\n", sep="")
  
  cat("  Expected variation between replications in between-network connectivity for module 1 (edgeshift_weight_bw_bwn_sd): ", edgeshift_weight_bw_bwn_sd, "\n", sep="")
  cat("  Observed variation: ", round(sd(m1bwn_meanshift), 3), "\n\n")
  
  #per-replication between-network edge variation (SD) from ground truth
  m1bwn_repsds <- sapply(patients, function(replication) { sd(E(replication)[between_edges]$weight - E(g_positive)[between_edges]$weight) })
  
  cat("  Expected variation within replications in between-network connectivity for module 1 (edgeshift_weight_wi_bwn_sd): ", edgeshift_weight_wi_bwn_sd, "\n", sep="")
  cat("  Observed variation: ", round(mean(m1bwn_repsds), 3), "\n\n", sep="")
  
  
  #mean shift in within-network connectivity for module 1 across replications (VERIFIED 15Apr2017)
  m1win_meanshift <- sapply(patients, function(replication) { mean(E(replication)[within_edges]$weight - E(g_positive)[within_edges]$weight) })
  
  cat("  Expected mean shift in within-network connectivity for module 1 (edgeshift_weight_bw_win_mean): ", edgeshift_weight_bw_win_mean, "\n", sep="")
  cat("  Observed variation: ", round(mean(m1win_meanshift), 3), "\n\n", sep="")
  
  cat("  Expected variation between replications in within-network connectivity for module 1 (edgeshift_weight_bw_win_sd): ", edgeshift_weight_bw_win_sd, "\n", sep="")
  cat("  Observed variation: ", round(sd(m1win_meanshift), 3), "\n\n")
  
  #per-replication between-network edge variation (SD) from ground truth
  m1win_repsds <- sapply(patients, function(replication) { sd(E(replication)[within_edges]$weight - E(g_positive)[within_edges]$weight) })
  
  cat("  Expected variation within replications in within-network connectivity for module 1 (edgeshift_weight_wi_win_sd): ", edgeshift_weight_wi_win_sd, "\n", sep="")
  cat("  Observed variation: ", round(mean(m1win_repsds), 3), "\n\n", sep="")
  
  #some leftover checks
  #  v1edgeReps <- do.call(rbind, lapply(controls, function(g) { E(g)[inc("V1")]$weight}))
  #  v1edgeMeans <- colMeans(v1edgeReps)
  #  colSDs()
  #  
  #  
  ##check on nodal tweaks to controls
  #  firstcontrol = allsims[["control"]][[1]]
  #  E(firstcontrol)[inc("V1")] #edges incident to vertex V1
  #  E(firstcontrol)[inc("V1")]$weight #edges weights incident to vertex V1 -- yup those are huge (given +20 mean)
  #  sd(E(firstcontrol)[inc("V1")]$weight) #wi node variation in strength for one replication: Should match edgeshift_weight_wi_sd
  #  
  ##verify that V1 -- V189 is not double tweaked (check the sequential processing of node tweaks): VERIFIED 30Sep2016
  #  E(firstcontrol)[inc("V189")] #which nodes are connected to 189 
  #  sort(E(firstcontrol)[inc("V189")]$weight) #edges incident to vertex V189 -- yup those are huge (given +100 mean)
  #  E(firstcontrol)[inc("V189")]$weight #edges incident to vertex V189 -- yup those are huge (given +100 mean)
  #  sd(E(firstcontrol)[inc("V1")]$weight) #SD of edges incident to V1. should correspond to edgeshift_weight_wi_sd for V1
  #  sd(E(firstcontrol)[inc("V189")]$weight) #wi node variation in strength for one replication. should correspond to edgeshift_weight_wi_sd for V2
  #  
  #  E(firstcontrol)["V1" %--% "V189"]$weight #weight of the specific connection that should be weaker (~20-30) since V1 -- V189 was tweaked first
  #  firstcontrol["V1","V189"] #even more compact subset notation (igraph internally treats this as a weighted adjacency matrix) 
  #  g_positive["V1","V189"] #even more compact subset notation (igraph internally treats this as a weighted adjacency matrix) 
  #  
  
}


## grabbed from braingraph, adapted not to require their rewiring function or the package itself (given RGtk2 dependency)

#' Calculate graph small-worldness
#'
#' This function will calculate the characteristic path length and clustering
#' coefficient, which are used to calculate small-worldness.
#'
#' @param g The graph (or list of graphs) of interest
#' @param rand List of (lists of) equivalent random graphs (output from
#' \code{\link{sim.rand.graph.par}})
#' @export
#'
#' @return A data frame with the following components:
#' \item{density}{The range of density thresholds used.}
#' \item{N}{The number of random graphs that were generated.}
#' \item{Lp}{The characteristic path length.}
#' \item{Cp}{The clustering coefficient.}
#' \item{Lp.rand}{The mean characteristic path length of the random graphs with
#' the same degree distribution as g.}
#' \item{Cp.rand}{The mean clustering coefficient of the random graphs with
#' the same degree distribution as g.}
#' \item{Lp.norm}{The normalized characteristic path length.}
#' \item{Cp.norm}{The normalized clustering coefficient.}
#' \item{sigma}{The small-world measure of the graph.}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Watts D.J., Strogatz S.H. (1998) \emph{Collective dynamics of
#' 'small-world' networks}. Nature, 393:440-442.

small.world <- function(g, rand) {
  require(data.table)
  if (is_igraph(g)) g <- list(g)  # Single graph at a single density
  Lp <- vapply(g, function(x) graph_attr(x, 'Lp'), numeric(1))
  Cp <- vapply(g, function(x) graph_attr(x, 'Cp'), numeric(1))
  densities <- vapply(g, function(x) graph_attr(x, 'density'), numeric(1))
  
  if (is_igraph(rand[[1]])) {
    Lp.rand <- mean(vapply(rand, function(x) graph_attr(x, 'Lp'), numeric(1)))
    Cp.rand <- mean(vapply(rand, function(x) graph_attr(x, 'Cp'), numeric(1)))
    N <- length(rand)
  } else {
    if (length(rand[[1]]) == 1) {  # If there's 1 rand graph for each density
      Lp.rand <- vapply(rand,
                        function(x) vapply(x, function(y) graph_attr(y, 'Lp'), numeric(1)),
                        numeric(1))
      Cp.rand <- vapply(rand,
                        function(x) vapply(x, function(y) graph_attr(y, 'Cp'), numeric(1)),
                        numeric(1))
      N <- 1
    } else {
      Lp.rand <- colMeans(vapply(rand,
                                 function(x) vapply(x, function(y) graph_attr(y, 'Lp'), numeric(1)),
                                 numeric(length(rand[[1]]))))
      Cp.rand <- colMeans(vapply(rand,
                                 function(x) vapply(x, function(y) graph_attr(y, 'Cp'), numeric(1)),
                                 numeric(length(rand[[1]]))))
      N <- lengths(rand)
    }
  }
  Cp.norm <- Cp / Cp.rand
  Lp.norm <- Lp / Lp.rand
  sigma <- Cp.norm / Lp.norm
  return(data.table(density=densities, N, Lp, Cp, Lp.rand, Cp.rand, Lp.norm,
                    Cp.norm, sigma))
}