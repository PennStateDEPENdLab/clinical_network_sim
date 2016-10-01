#use 10895 data processed as part of possum motion
setwd(file.path(getMainDir(), "clinical_network_sim"))
source("simulate_clingraph.R")
adjmat <- as.matrix(read.table("10895_bb264_corrmat.zval.1D"))
mean(adjmat[lower.tri(adjmat)])
hist(adjmat[lower.tri(adjmat)])

#create graphs from adjacency matrices
library(igraph)

#original graphs with correlation weights 
g <- graph.adjacency(adjmat, mode="undirected", weighted=TRUE, diag=FALSE)
g_positive <- delete.edges(g, which(E(g)$weight < 0)) #remove edges with negative weights 

#mean z-val
mean(E(g)$weight)
mean(E(g_positive)$weight)

hist(E(g_positive)$weight)

densities <- seq(.01, .20, by=.01)

allg_density <- lapply(densities, function(d, g) {
      thresh <- quantile(E(g)$weight, 1-d)
      
      gthresh <- delete.edges(g, which(E(g)$weight <= thresh)) #return binary weighted density-thresholded graph
      gthresh <- remove.edge.attribute(gthresh, "weight")
      return(gthresh)
    }, g)

names(allg_density) <- paste0(densities*100, "%")

gmeasures <- function(g) {
  deg <- degree(g)
  evc <- evcent(g)$vector
  bet <- betweenness(g, normalize=TRUE)*100
  locclust <- transitivity(g, type="local")
  pr <- page.rank(g, algo="prpack")$vector
  data.frame(degree=deg, evcent=evc, betweenness=bet, locclust=locclust, pagerank=pr) #return a 264 x nmetrics matrix for analysis/reduction            
}

#for basic checks, compute degree, eigenvector centrality, local clustering, betweenness, and pagerank at each density threshold
all_centrality <- lapply(allg_density, function(gthresh) {
      gmeasures(gthresh)
    })

#number of edges
sapply(allg_density, function(gthresh) { length(E(gthresh)) })

sanitycheck_sims(g_positive)

allsims <- simulate_clingraph(g_positive, npergroup=100, edgenoise=0.0, #just for sanity check
    groupmanipulation=list(
        controls=list(group_name="control", 
            edgeshift_weight_mean_bw=0, edgeshift_weight_sd_wi=0, edgeshift_weight_sd_bw=0, #global shifts in weights
            nodemanipulation=list( #nodal shifts in (weighted) edges incident to a node
                list(name="V1", edgeshift_weight_mean_bw=20, edgeshift_weight_sd_wi=10, edgeshift_weight_sd_bw=5),
                list(name="V189", edgeshift_weight_mean_bw=100, edgeshift_weight_sd_wi=10, edgeshift_weight_sd_bw=5)
            )
        ),
        patients=list(group_name="patients", edgeshift_weight_mean_bw=0.5, edgeshift_weight_sd_wi=.06, edgeshift_weight_sd_bw=.25)))

#check on nodal tweaks to controls
firstcontrol = allsims[["control"]][[1]]
E(firstcontrol)[inc("V1")] #edges incident to vertex V1
E(firstcontrol)[inc("V1")]$weight #edges weights incident to vertex V1 -- yup those are huge (given +20 mean)
sd(E(firstcontrol)[inc("V1")]$weight) #wi node variation in strength for one replication

#verify that V1 -- V189 is not double tweaked (check the sequential processing of node tweaks): VERIFIED 30Sep2016
E(firstcontrol)[inc("V189")] #which nodes are connected to 189 
sort(E(firstcontrol)[inc("V189")]$weight) #edges incident to vertex V189 -- yup those are huge (given +100 mean)
sd(E(firstcontrol)[inc("V1")]$weight) #wi node variation in strength for one replication

E(firstcontrol)["V1" %--% "V189"]$weight #weight of the specific connection that should be weaker (~20-30) since V1 -- V189 was tweaked first
firstcontrol["V1","V189"] #even more compact subset notation (igraph internally treats this as a weighted adjacency matrix) 