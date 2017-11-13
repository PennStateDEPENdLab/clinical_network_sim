#use 10895 data processed as part of possum motion
setwd(file.path(getMainDir(), "clinical_brainnetwork_sim"))
source("simulate_clingraph.R")

#create graphs from adjacency matrices
library(igraph)

mean(adjmat[lower.tri(adjmat)])
par(mfrow=c(1,2)); hist(adjmat[lower.tri(adjmat)], main="from zstat"); hist(adjmat_noz[lower.tri(adjmat_noz)], main="no z")
range(adjmat[lower.tri(adjmat)])
which(adjmat[lower.tri(adjmat)]==15)

#original graphs with correlation weights 

#for basic checks, compute degree, eigenvector centrality, local clustering, betweenness, and pagerank at each density threshold
all_centrality <- lapply(allg_density, function(gthresh) {
      gmeasures(gthresh)
    })

#number of edges
sapply(allg_density, function(gthresh) { length(E(gthresh)) })

#verify that functions work as expected
sink(paste0("sanity_checks_", format(Sys.Date(), "%d%b%Y"), ".txt"))
sanitycheck_sims(g_positive)
sink()

testsim1 <- simulate_clingraph(g_positive, npergroup=100, edgenoise=list(sd=0.05, dist="gaussian"),
    groupmanipulation=list(
        controls=list(group_name="control", 
            edgeshift_weight_bw_mean=0, edgeshift_weight_wi_sd=0, edgeshift_weight_bw_sd=0
        ),         #global shifts in weights
        patients=list(group_name="patient", edgeshift_weight_bw_mean=0.0, edgeshift_weight_wi_sd=0, edgeshift_weight_bw_sd=0, #global shifts in weights
            nodemanipulation=list( #nodal shifts in (weighted) edges incident to a node
                list(name="V215", edgeshift_weight_bw_mean=.3, edgeshift_weight_wi_sd=.03, edgeshift_weight_bw_sd=.03), #D ACC
                list(name="V208", edgeshift_weight_bw_mean=.4, edgeshift_weight_wi_sd=.04, edgeshift_weight_bw_sd=.04),  #L IFG
                list(name="V209", edgeshift_weight_bw_mean=.4, edgeshift_weight_wi_sd=.04, edgeshift_weight_bw_sd=.04)   #R IFG
            )
        )
    )
)

# testsim2 <- simulate_clingraph(g_positive, npergroup=100, edgenoise=list(sd=0.05, dist="gaussian"),
#     groupmanipulation=list(
#         controls=list(group_name="control", 
#             edgeshift_weight_bw_mean=0, edgeshift_weight_wi_sd=0, edgeshift_weight_bw_sd=0
#         ),         #global shifts in weights
#         patients=list(group_name="patient", edgeshift_weight_bw_mean=0.0, edgeshift_weight_wi_sd=0, edgeshift_weight_bw_sd=0, #global shifts in weights
#             networkmanipulation=list( #network-level shifts in (weighted) edges between and within a network
#                 list(module=1, edgeshift_weight_bw_bwn_mean=.3, edgeshift_weight_bw_bwn_sd=.03, edgeshift_weight_wi_bwn_sd=.03,
#                     edgeshift_weight_bw_win_mean=.3, edgeshift_weight_bw_win_sd=.03, edgeshift_weight_wi_win_sd=.03)
#             ),
#             nodemanipulation=list( 
#                 list(name="V215", edgeshift_weight_bw_mean=.3, edgeshift_weight_wi_sd=.03, edgeshift_weight_bw_sd=.03), #D ACC
#                 list(name="V208", edgeshift_weight_bw_mean=.4, edgeshift_weight_wi_sd=.04, edgeshift_weight_bw_sd=.04),  #L IFG
#                 list(name="V209", edgeshift_weight_bw_mean=.4, edgeshift_weight_wi_sd=.04, edgeshift_weight_bw_sd=.04)   #R IFG
#             )
#         )
#     )
# )



str(all_centrality[[1]][["control"]][[1]])




#what edges are present at a given density for controls (unaltered) and patients (nodeal increases). Whack a node
#plot these on brain

#tr


#test <- lapply(1:100,function(x) rnorm(100000))
#system.time(x <- lapply(test,function(x) loess.smooth(x,x)))
#system.time(x <- mclapply(test,function(x) loess.smooth(x,x), mc.cores=4))




#global increase in edge strength in patients, then compare partial versus full correlation




