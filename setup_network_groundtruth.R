#use 10895 data processed as part of possum motion
library(dplyr)
setwd(file.path(getMainDir(), "clinical_brainnetwork_sim"))
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
  data.frame(node=V(g)$name, degree=deg, evcent=evc, betweenness=bet, locclust=locclust, pagerank=pr) #return a 264 x nmetrics matrix for analysis/reduction            
}

#for basic checks, compute degree, eigenvector centrality, local clustering, betweenness, and pagerank at each density threshold
all_centrality <- lapply(allg_density, function(gthresh) {
      gmeasures(gthresh)
    })

#number of edges
sapply(allg_density, function(gthresh) { length(E(gthresh)) })

#verify that functions work as expected
sanitycheck_sims(g_positive)


testsim1 <- simulate_clingraph(g_positive, npergroup=100, edgenoise=0.05,
    groupmanipulation=list(
        controls=list(group_name="control", 
            edgeshift_weight_mean_bw=0, edgeshift_weight_sd_wi=0, edgeshift_weight_sd_bw=0
        ),         #global shifts in weights
        patients=list(group_name="patient", edgeshift_weight_mean_bw=0.0, edgeshift_weight_sd_wi=0, edgeshift_weight_sd_bw=0, #global shifts in weights
            nodemanipulation=list( #nodal shifts in (weighted) edges incident to a node
                list(name="V215", edgeshift_weight_mean_bw=.3, edgeshift_weight_sd_wi=.03, edgeshift_weight_sd_bw=.03), #D ACC
                list(name="V208", edgeshift_weight_mean_bw=.4, edgeshift_weight_sd_wi=.04, edgeshift_weight_sd_bw=.04),  #L IFG
                list(name="V209", edgeshift_weight_mean_bw=.4, edgeshift_weight_sd_wi=.04, edgeshift_weight_sd_bw=.04)   #R IFG
            )
        )
    )
)

#now look at between group differences
#for a moment, just look at 10% density
densities <- c(.05, .06, .07, .08, .09, .1)
library(parallel)
allg_density <- lapply(densities, function(d, simlist) {
      lapply(simlist, function(group) {
            mclapply(group, function(g) {
                  #g is a replication igraph object
                  #thresh <- quantile(E(g)$weight, 1-d) #this takes quantile of non-negative edges, not overall density based on graph diameter
                  #gthresh <- delete.edges(g, which(E(g)$weight <= thresh)) #return binary weighted density-thresholded graph
                  #gthresh <- remove.edge.attribute(gthresh, "weight")
                  
                  #related approach, but will obtain desired density given graph diameter
                  weights <- sort(E(g)$weight, decreasing=TRUE)
                  threshold <- weights[length(V(g))*(length(V(g))-1)/2 * d]
                  gthresh <- delete.edges(g, which(E(g)$weight < threshold))
                  gthresh <- remove.edge.attribute(gthresh, "weight")
                  
                  return(gthresh)
                }, mc.cores=4)
          })
    }, testsim1)

names(allg_density) <- as.character(densities)

x <- allg_density[["0.1"]][["control"]][[2]]
graph.density(x)
is.directed(x)

#for basic checks, compute degree, eigenvector centrality, local clustering, betweenness, and pagerank at each density threshold
all_centrality <- lapply(allg_density, function(dlist) {
      lapply(dlist, function(group) {
            mclapply(group, function(gthresh) {
                  gmeasures(gthresh)            
                })
          })      
    })

str(all_centrality[[1]][["control"]][[1]])
str(all_centrality[[1]][["patient"]][[1]])

#compare target nodes
subset(all_centrality[["0.05"]][["control"]][[1]], node %in% c("V207", "V208", "V209", "V215"))
subset(all_centrality[["0.05"]][["patient"]][[1]], node %in% c("V207", "V208", "V209", "V215"))


#get mean and sd centrality
all_agg <- lapply(all_centrality, function(dlist) {
      lapply(dlist, function(group) {
            #here, group is a list of subjects where each is a 264 x measures data.frame
            recode <- lapply(1:length(group), function(g) {
                  group[[g]]$replication <- g #add replication
                  group[[g]]$nodeNum <- 1:nrow(group[[g]]) #add replication
                  return(group[[g]])
                })            
            df <- do.call(rbind, recode)
            df %>% group_by(node) %>% summarize_at(vars(degree, evcent, betweenness, locclust, pagerank), funs(mean=mean, sd=sd))
          })
    })

#take 0.1 for a sec
xx <- all_centrality[["0.1"]]

patients <- do.call(rbind, xx[["patient"]])
patients$group <- "patient"
controls <- do.call(rbind, xx[["control"]])
controls$group <- "control"

both <- rbind(patients, controls)

library(broom)
res <- lapply(split(both, both$node), function(df) {
      if (sd(df$degree) < .01) { return(NULL)}
      mm <- lm(degree ~ group, df)
      partab <- tidy(mm)
      pfilt <- filter(partab, term=="grouppatient") %>% mutate(node = df$node[1])
      return(pfilt)
    })

resdf <- do.call(rbind, res)
filter(resdf, estimate < 0)
filter(resdf, estimate > 0)

as_tibble(lapply(resdf, function(col) { 
        if (is.numeric(col)) { return(round(col, 3))} else { return(col) }

library(gridExtra)
pdf("degreediff.pdf", height=11, width=8.5)
grid.table(resdf[1:50,])
dev.off()

sum(filter(resdf, estimate < 0)$estimate) 
sum(filter(resdf, estimate > 0)$estimate)


library(broom)
res <- lapply(split(both, both$node), function(df) {
      if (sd(df$degree) < .01) { return(NULL)}
      mm <- lm(evcent ~ group, df)
      partab <- tidy(mm)
      pfilt <- filter(partab, term=="grouppatient") %>% mutate(node = df$node[1])
      return(pfilt)
    })

resdf <- do.call(rbind, res)
filter(resdf, estimate < 0)
filter(resdf, estimate > 0)

sum(filter(resdf, estimate < 0)$estimate) 
sum(filter(resdf, estimate > 0)$estimate)


#now should be in a position to compare nodal statistics by group after thresholding
d10 <- all_agg[["0.09"]]
d10[["control"]]$group <- "control"
d10[["patient"]]$group <- "patient"
d10bind <- do.call(rbind, d10)



library(ggplot2)
pdf("deg diff test.pdf", width=12, height=40)
ggplot(d10bind, aes(x=factor(node), y=degree_mean, ymin=degree_mean-degree_sd, ymax=degree_mean+degree_sd, fill=group)) + 
    geom_bar(stat="identity", position="dodge", width=0.5) + 
    geom_errorbar(position="dodge", width=0.5) +theme_bw(base_size=12) + coord_flip()
    #facet_wrap(~node, ncol=1) + coord_flip() + theme_bw(base_size=12)
dev.off()




str(all_centrality[[1]][["control"]][[1]])



#what edges are present at a given density for controls (unaltered) and patients (nodeal increases). Whack a node
#plot these on brain

#tr


#test <- lapply(1:100,function(x) rnorm(100000))
#system.time(x <- lapply(test,function(x) loess.smooth(x,x)))
#system.time(x <- mclapply(test,function(x) loess.smooth(x,x), mc.cores=4))




#global increase in edge strength in patients, then compare partial versus full correlation



