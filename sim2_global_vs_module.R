setwd(file.path(getMainDir(), "clinical_brainnetwork_sim"))
source("simulate_clingraph.R") #loads adjmat into workspace
library(brainGraph)
library(parallel)
library(tidyr)
library(dplyr)
library(broom)
##Sim 2: insensitivity of global metrics to modular changes

#group 1: increase C-O task control and salience nodes
#in Power 264, these are 47-60 [C-O] and 174-181,186-202 [Salience]

#group 2: increase DMN patient nodes
#in Power 264, these are 74-83, 86-131

#if (!file.exists("sim2_global_with_braingraph.RData")) {  
#  comsim1 <- simulate_clingraph(g_positive, npergroup=50, edgenoise=list(sd=0.2, dist="gaussian"),
#    groupmanipulation=list(
#      controls=list(
#        group_name="controls", edgeshift_weight_bw_mean=0, edgeshift_weight_wi_sd=0.2, edgeshift_weight_bw_sd=0.2, #global shifts in weights
#        networkmanipulation=list( #increase C-O [3] and Salience networks [9] 
#          list(module=c(3,9), 
#            edgeshift_weight_bw_bwn_mean=.25, 
#            edgeshift_weight_bw_bwn_sd=.1, 
#            edgeshift_weight_wi_bwn_sd=.2,
#            edgeshift_weight_wi_bwn_dist="ex-gaussian",
#            edgeshift_weight_wi_bwn_sbratio=1.2, #pos skew of between-network weight shifts within replication 
#            edgeshift_weight_bw_win_mean=.25, #start with equal shifts in win and win
#            edgeshift_weight_bw_win_sd=.1, 
#            edgeshift_weight_wi_win_sd=.2,
#            edgeshift_weight_wi_win_dist="ex-gaussian",
#            edgeshift_weight_wi_win_sbratio=1.2 #pos skew of within-network weight shifts within replication
#          )
#        )
#      ),
#      patients=list(
#        group_name="patients", edgeshift_weight_bw_mean=0.0, edgeshift_weight_wi_sd=0.2, edgeshift_weight_bw_sd=0.2, #global shifts in weights
#        networkmanipulation=list( #increase DMN [5] 
#          list(module=5, 
#            edgeshift_weight_bw_bwn_mean=.25, 
#            edgeshift_weight_bw_bwn_sd=.1, 
#            edgeshift_weight_wi_bwn_sd=.2,
#            edgeshift_weight_wi_bwn_dist="ex-gaussian",
#            edgeshift_weight_wi_bwn_sbratio=1.2, #pos skew of between-network weight shifts within replication 
#            edgeshift_weight_bw_win_mean=.25, #start with equal shifts in win and win
#            edgeshift_weight_bw_win_sd=.1, 
#            edgeshift_weight_wi_win_sd=.2,
#            edgeshift_weight_wi_win_dist="ex-gaussian",
#            edgeshift_weight_wi_win_sbratio=1.2 #pos skew of within-network weight shifts within replication
#          )
#        )
#      )
#    )
#  )
  
  comsim1 <- simulate_clingraph(g_positive, npergroup=50, edgenoise=list(sd=0.2, dist="gaussian"),
    groupmanipulation=list(
      controls=list(
        group_name="controls", edgeshift_weight_bw_mean=0, edgeshift_weight_wi_sd=0.2, edgeshift_weight_bw_sd=0.2, #global shifts in weights
        networkmanipulation=list( #increase C-O [3] and Salience networks [9] 
          list(module=c(8,12), #dial up F-P and DAN # c(3,9), 
            edgeshift_weight_bw_bwn_mean=.1, 
            edgeshift_weight_bw_bwn_sd=0.05, 
            edgeshift_weight_wi_bwn_sd=0.05,
            edgeshift_weight_wi_bwn_dist="gaussian",
            edgeshift_weight_bw_win_mean=.2, #start with equal shifts in bwn and win
            edgeshift_weight_bw_win_sd=.1, 
            edgeshift_weight_wi_win_sd=.1,
            edgeshift_weight_wi_win_dist="gaussian"
          )
        )
      ),
      patients=list(
        group_name="patients", edgeshift_weight_bw_mean=0.0, edgeshift_weight_wi_sd=0.2, edgeshift_weight_bw_sd=0.2, #global shifts in weights
        networkmanipulation=list( #increase DMN [5] 
          list(module=c(5), #DMN is 5, 12 is dorsal attention, 8 is fronto-parietal 
            edgeshift_weight_bw_bwn_mean=.1, 
            edgeshift_weight_bw_bwn_sd=0.05, 
            edgeshift_weight_wi_bwn_sd=0.05,
            edgeshift_weight_wi_bwn_dist="gaussian",
            edgeshift_weight_bw_win_mean=0.2, #start with equal shifts in bwn and win
            edgeshift_weight_bw_win_sd=0.1, 
            edgeshift_weight_wi_win_sd=0.1,
            edgeshift_weight_wi_win_dist="gaussian"
          )
        )
      )
    )
  )  
  
  #length(comsim1[[1]])
  
  #this is very slow if the set_brainGraph_attr is in place 
  #densities <- c(.06, .08, .1, .12)
  densities <- seq(.05, .25, .025)
  
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
              #gthresh <- set_brainGraph_attr(gthresh, modality="fmri") #add various statistics expected by brainGraph
              gthresh$Cp <- transitivity(gthresh, type='localaverage')
              gthresh$Lp <- mean_distance(gthresh)
              gthresh$density <- edge_density(gthresh)
              return(gthresh)
            }, mc.cores=4)
        })
    }, comsim1)
  
  names(allg_density) <- as.character(densities)
  
  #str(allg_density[[1]][[1]][[1]])
  #allg_density[[1]][[1]][[1]]$group
  
#  save(file="sim2_global_with_braingraph.RData", allg_density)
#} else {
#  load(file="sim2_global_with_braingraph.RData")
#}

#for basic checks, compute degree, eigenvector centrality, local clustering, betweenness, and pagerank at each density threshold
all_centrality <- lapply(allg_density, function(dlist) {
    lapply(dlist, function(group) {
        mclapply(group, function(gthresh) {
            gmeasures(gthresh)            
          })
      })      
  })

all_global <- lapply(allg_density, function(dlist) {
    gden <- lapply(dlist, function(group) {
        glist <- do.call(rbind,mclapply(group, function(gthresh) {
              global.measures(gthresh)            
            }))
        
        glist$replication <- 1:nrow(glist)
        return(glist)
      })
    
    gden <- lapply(1:length(gden), function(d) { gden[[d]]$group <- names(gden)[d]; return(gden[[d]]) })
    gden <- do.call(rbind, gden)
  })


#first participant in lowest density from control group
# library(foreach)
# library(doSNOW)
# clus <- makeSOCKcluster(4)
# system.time(rand1 <- sim.rand.graph.par(allg_density[[1]][[1]][[1]], N=1e1, clustering=F))
# stopCluster(clus)



#g <- allg_density[[1]][[1]][[1]]
#
##rewire 
#ff <- lapply(1:100, function(it) { 
#    r <- rewire(g, with = keeping_degseq(niter = ecount(g) * 10))
#    r$Cp <- transitivity(r, type='localaverage')
#    r$Lp <- mean_distance(r)
#    return(r)
#  })


rewire_glist <- function(glist) {
  lapply(glist, function(dlist) {
      glist <- lapply(dlist, function(group) {
          dd <- mclapply(group, function(g) {
              rlist <- lapply(1:100, function(it) { #should be 100
                  r <- rewire(g, with = keeping_degseq(niter = ecount(g) * 10)) #should be 10
                  r$Cp <- transitivity(r, type='localaverage')
                  r$Lp <- mean_distance(r)
                  return(r)
                })
              return(small.world(g, rlist))
            }, mc.cores=4)
          dd <- bind_rows(dd)
        })
      
      #add group as column
      for (g in 1:length(glist)) { glist[[g]]$group <- names(glist)[g] }
      
      #combine into single data.frame    
      bind_rows(glist)
    })
}

test2 <- rewire_glist(allg_density)
test2 <- lapply(test2, function(d) { d$id <- 1:nrow(d); return(d) }) #add ID for mlm
test3 <- bind_rows(test2)

#vv <- small.world(g, ff)
#
#x <- rewire(allg_density[[1]][[1]][[1]], keeping_degseq(niter = 20))
#
#small.world(allg_density[[1]][[1]][[1]], rand1)
#
#system.time(rand1.cl <- sim.rand.graph.par(allg_density[[1]][[1]][[1]], N=1e1, max.iters=1e3))


#all_global <- do.call(rbind, lapply(1:length(all_global), function(d) {
#      all_global[[d]]$density <- as.numeric(names(all_global)[d])
#      return(all_global[[d]])
#    }))


#test uber script

#analysis_random_graphs(allg_density[[1]], N=1e1, #covars.dti,
#  savedir='/Users/michael/Data_Analysis/clinical_brainnetwork_sim', clustering=F)

library(ggplot2)
#ggplot(all_global, aes(x=factor(density), y=Cglob, color=group)) + geom_boxplot()

#ggplot(all_global, aes(x=factor(density), y=Lglob, color=group)) + geom_boxplot()

#pdf("clustering_pathlength_ratio.pdf", width=8, height=5)
#ggplot(all_global, aes(x=group, y=Cglob/Lglob)) + geom_boxplot() + facet_wrap(~density, scales="free") + theme_bw(base_size=12)
#dev.off()
#
#str(all_global[[1]])

library(tidyr)
library(purrr)

pdf("sim2_small_worldness.pdf", width=7, height=4)
test3 %>% filter(density > .05) %>% mutate(density=density*100) %>% group_by(density, group) %>% summarize(m=median(sigma), se=plotrix::std.error(sigma), sd=sd(sigma), Lower=quantile(sigma, .10), Upper=quantile(sigma, .90)) %>%
  ggplot(., aes(x=density, y=m, ymin=Lower, ymax=Upper, color=group)) + geom_crossbar(position=position_dodge(width=1.3), width=1) +
  scale_color_brewer("Group", palette="Dark2") + theme_bw(base_size = 15) + xlab("Density (%)") + ylab(expression(paste("Small worldness (", italic(sigma), ")")))
dev.off()
  
test3 %>% filter(density > .06) %>% group_by(density) %>% nest() %>% mutate(model=map(data, ~ t.test(sigma ~ group, data=.))) %>% 
  unnest(model %>% purrr::map(broom::glance)) %>% summarize(mean(statistic), sd(statistic), min(p.value))

mm <- lmer(sigma ~ density * group + (1|id), filter(test3, density > .05))
car::Anova(mm)

#all_global %>% mutate(SmW=Cglob/Lglob) %>% group_by(density) %>% nest() %>% mutate(model=map(data, ~ t.test(SmW ~ group, data=.))) %>% 
#  unnest(model %>% purrr::map(broom::glance))
#
#all_global %>% mutate(SmW=Cglob/Lglob) %>% group_by(density) %>% nest() %>% mutate(model=map(data, ~ t.test(Cglob ~ group, data=.))) %>% 
#  unnest(model %>% purrr::map(broom::glance))

#look at within and between modules connectivity

wibw_z <- lapply(allg_density, function(dlist) {
    gden <- lapply(dlist, function(group) {
        glist <- lapply(group, function(gthresh) {
            V(gthresh)$module[V(gthresh)$module==-1] <- 14 #recode -1 as 14 to get this function to work
            zstats <- wibw_module_degree(gthresh)
            #zstats_doublecheck <- wibw_module_degree_slow(gthresh) #Apr2017: slow, but easy, code is identical to faster matrix-based code
            #zz <- within_module_deg_z_score(gthresh, V(gthresh)$module) #wibw_module_degree also checks with brainGraph function
            #zz[V(gthresh)$module == 14]
            #degree(gthresh)
          })
    
        glist <- do.call(rbind, lapply(1:length(glist), function(rep) { glist[[rep]]$replication <- rep; return(glist[[rep]]) }))
        return(glist)
      })
    
    gden <- lapply(1:length(gden), function(d) { gden[[d]]$group <- names(gden)[d]; return(gden[[d]]) })
    gden <- do.call(rbind, gden)
  })

wibw_z <- do.call(rbind, lapply(1:length(wibw_z), function(d) {
      wibw_z[[d]]$density <- names(wibw_z)[d]
      return(wibw_z[[d]])
    }))

wibw_z$id <- with(wibw_z, paste0(group, replication))

wibw_z$module <- factor(wibw_z$module, levels=1:14, labels=paste0("M", 1:14))

#z-score per module for clarity
wibw_z <- wibw_z %>% group_by(density, module) %>%
  mutate(Ki_z=(Ki - mean(Ki))/sd(Ki), Bi_z=(Bi - mean(Bi))/sd(Bi)) %>% ungroup()


wibw_z %>% group_by(density, module) %>% do(psych::describe(.$Ki_z))

m2 <- lmer(Ki_z ~ density*module*group + (1|id), filter(wibw_z, density > .06))
car::Anova(m2)
cm <- lmerCellMeans(m2)
cm$Ki_z <- as.vector(cm$Ki_z)

cmwi <- filter(cm, density %in% c("0.1", "0.2")) %>% rename(stat=Ki_z) %>% mutate(struct="Within-Module")

ggplot(filter(cm, density %in% c("0.1", "0.2")), aes(x=module, y=Ki_z, ymin=Ki_z-se, ymax=Ki_z+se, color=group)) + geom_pointrange() + facet_wrap(~density) + coord_flip()

m3 <- lmer(Bi_z ~ density*module*group + (1|id), filter(wibw_z, density > .06))
car::Anova(m3)
cm3 <- lmerCellMeans(m3)
cm3$Bi_z <- as.vector(cm3$Bi_z)

cmbw <- filter(cm3, density %in% c("0.1", "0.2")) %>% rename(stat=Bi_z) %>% mutate(struct="Between-Module")

cmboth <- bind_rows(cmwi, cmbw) %>% filter(module %in% c("M5", "M8", "M12", "M7")) %>% mutate(module=ordered(module, levels=c("M5", "M8", "M12", "M7"), labels=c("DMN", "FPN", "DAN", "Visual")),
  density=factor(density, levels=c("0.1", "0.2"), labels=c("10% Density", "20% Density")))

pdf("wibw_deg_cellmeans.pdf", width=6, height=4)
ggplot(cmboth, aes(x=module, y=stat, ymin=tlo, ymax=thi, color=group)) + geom_pointrange(position=position_dodge(width=0.4)) + facet_grid(struct~density) + coord_flip(ylim=c(-0.9, 0.9)) +
  theme_bw(base_size=14) + xlab("Module") + ylab(expression(paste("Degree (", italic(z), ")"))) + scale_color_brewer("Group", palette="Dark2") + 
  theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.x=element_blank(),
    axis.title.x=element_text(margin=margin(t=10)), axis.title.y=element_text(margin=margin(r=12))) +
  scale_y_continuous(breaks = pretty_breaks(n = 3))
dev.off()


#means and SDs
options(dplyr.print_max = 1e9)
wibw_z %>% group_by(module, density, group) %>% select(Ki, Bi) %>% 
  summarize_all(funs(mean, sd))

pdf("within_module_shift.pdf", width=15, height=7)
for (d in as.character(densities)) {
  df <- filter(wibw_z, density==d)
  g <- ggplot(df, aes(x = group, y=Ki)) + facet_wrap(~module, scales="free") + stat_summary(fun.data = mean_cl_boot) + theme_bw(base_size=15) + ggtitle(paste0("density: ", d))
  plot(g)
} 
dev.off()

#z score Ki statistics to put onto same scale across modules




pdf("between_module_shift.pdf", width=15, height=7)
for (d in as.character(densities)) {
  df <- filter(wibw_z, density==d)
  g <- ggplot(df, aes(x = group, y=Bi)) + facet_wrap(~module, scales="free") + stat_summary(fun.data = mean_cl_boot) + theme_bw(base_size=15) + ggtitle(paste0("density: ", d))
  plot(g)
} 
dev.off()


#mutate(SmW=Cglob/Lglob) %>% group_by(density) %>% nest() %>% mutate(model=map(data, ~ t.test(SmW ~ group, data=.))) %>% 
#  unnest(model %>% purrr::map(broom::glance))
#
#all_global %>% mutate(SmW=Cglob/Lglob) %>% group_by(density) %>% nest() %>% mutate(model=map(data, ~ t.test(Cglob ~ group, data=.))) %>% 
#  unnest(model %>% purrr::map(broom::glance))
