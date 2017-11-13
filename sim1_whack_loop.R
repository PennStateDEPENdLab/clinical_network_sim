setwd(file.path(getMainDir(), "clinical_brainnetwork_sim")) #should be where this script lives already
source("simulate_clingraph.R")
suppressMessages(library(tidyverse))
suppressMessages(library(igraph))
suppressMessages(library(knitr))
suppressMessages(library(broom))
suppressMessages(library(parallel))
suppressMessages(library(doParallel))

library(ggplot2)
library(dplyr)

whack_increase <- .14
whack_decrease <- -.04
neg_edgeshift_weight_wi_sd <- .02 #within replication variation in edge shift for a node
neg_edgeshift_weight_bw_sd <- .02 #between replication variation in the magnitude of the mean weight shift for a node
pos_edgeshift_weight_wi_sd <- .07 #within replication variation in edge shift for a node
pos_edgeshift_weight_bw_sd <- .07 #between replication variation in the magnitude of the mean weight shift for a node

npos <- 3
nneg <- 3
ncomparator <- 3 #untargeted nodes
nnnodes <- 264

#thresholds for testing
densities_desired <- seq(.05, .25, .025)
rs_desired <- seq(.2, .4, .025)

replist <- list()
#for (s in 1:100) {
clusterobj <- makePSOCKcluster(4) #4 cpus
registerDoParallel(clusterobj)

# replist <- foreach(s=1:100, .packages=c("igraph", "tidyverse", "broom"),
#                    .noexport=c("adjmat", "adjmat_ztrans", "g_orig", "allg_density")) %dopar% {
for (s in 1:1) {                     
  #select positive and negative targets at random
  towhack <- paste0("V", sample(264, npos + nneg + ncomparator))
  poswhack <- towhack[1:npos]
  negwhack <- towhack[(npos+1):(npos+nneg)]
  comparators <- towhack[(npos+nneg+1):(npos+nneg+ncomparator)]

  #poswhack <- c("V215", "V208", "V209")
  # poswhack <- c("V200", "V218", "V230")
  # negwhack <- c("V19", "V63", "V77")
  # comparators <- c("V22", "V25", "V39")
  # towhack <- c(poswhack, negwhack, comparators)
  
  #build list of node manipulations
  nodemanip <- list()
  for (n in poswhack) {
    nodemanip[[length(nodemanip) + 1]] <- list(name=n, edgeshift_weight_bw_mean=whack_increase, 
                                               edgeshift_weight_wi_sd=pos_edgeshift_weight_wi_sd, 
                                               edgeshift_weight_bw_sd=pos_edgeshift_weight_bw_sd)
  }
  
  for (n in negwhack) {
    nodemanip[[length(nodemanip) + 1]] <- list(name=n, edgeshift_weight_bw_mean=whack_decrease, 
                                               edgeshift_weight_wi_sd=neg_edgeshift_weight_wi_sd, 
                                               edgeshift_weight_bw_sd=neg_edgeshift_weight_bw_sd)
  }
  
  #simulate data for this replication
  whack <- simulate_clingraph(
    g_positive, npergroup=50, ncpus=1, edgenoise=list(sd=0.2, dist="gaussian"),
    groupmanipulation=list(
      controls=list(
        group_name="controls",
        edgeshift_weight_bw_mean=0, 
        edgeshift_weight_wi_sd=0.2, 
        edgeshift_weight_bw_sd=0.2 #global shifts in weights (needed to get reasonable b/w subs variability)
      ),
      patients=list(
        group_name="patients",
        edgeshift_weight_bw_mean=0,
        edgeshift_weight_wi_sd=0.2,
        edgeshift_weight_bw_sd=0.2, #global shifts in weights
        nodemanipulation=nodemanip
      )
    )
  )
  
  #create a density-thresholded version of this replication
  allg_selective_bydensity <- threshold_glist(whack, densities_desired, method="density", ncores = 1)
  
  # dlist <- lapply(allg_selective_bydensity, function(d) { 
  #   gout <- lapply(d, function(g) {
  #     sapply(g, function(s) {
  #       s$wthresh
  #     })
  #   })
  #   bind_rows(gout)
  # })
  # 
  # plot(sapply(dlist, function(d) {
  #   mean(colMeans(d))
  # }))

  #compute nodal measures and form into a data.frame
  nodal_bydensity <- bind_nodal_measures(nodal_measures_glist(allg_selective_bydensity, ncores = 1))
  
  #compute regression of degree on group for positive and negative targets, as well as comparator nodes
  degree_effects_bydensity <- group_lm(nodal_bydensity, analyze_nodes = towhack) %>% 
    flag_nodes(pos=poswhack, neg=negwhack)
  
  summdens <- degree_effects_bydensity %>% group_by(nodeType, thresh) %>% summarize(pval = mean(p.value), tval=mean(statistic)) %>% mutate(method="density")
  
  #check on range of densities
  #nodal_bydensity  %>% group_by(thresh, group) %>% summarize(mean(density), min(density), max(density))
  
  #same approach for strength
  allg_bystrength <- threshold_glist(whack, rs_desired, method="strength", ncores=1)
  
  nodal_bystrength <- bind_nodal_measures(nodal_measures_glist(allg_bystrength, ncores=1))
  degree_effects_bystrength <- group_lm(nodal_bystrength, analyze_nodes = towhack) %>% 
    flag_nodes(pos=poswhack, neg=negwhack)
  
  summstrength <- degree_effects_bystrength %>% group_by(nodeType, thresh) %>% summarize(pval=mean(p.value), tval=mean(statistic)) %>% mutate(method="strength")
  
  #range of densities under strength threshold
  #nodal_bystrength %>% group_by(thresh, group) %>% summarize(mean(density), min(density), max(density))
  
  #what about introducing density as a covariate
  degree_effects_bystrength_dcov <- group_lm(nodal_bystrength, analyze_nodes = towhack, f=formula(degree ~ group + density)) %>% 
    flag_nodes(pos=poswhack, neg=negwhack)
  
  summstrength_dcov <- degree_effects_bystrength_dcov %>% group_by(nodeType, thresh) %>% summarize(pval=mean(p.value), tval=mean(statistic)) %>% mutate(method="strength_dcov")
  
  #weighted analysis
  weighted_whack <- threshold_glist(whack, thresholds = 0, method="weight", rmweights=FALSE, ncores=1)
  
  nodal_weighted <- bind_nodal_measures(nodal_measures_glist(weighted_whack, ncores=1))
  weighted_analysis <- group_lm(nodal_weighted, analyze_nodes = towhack) %>%
    flag_nodes(pos=poswhack, neg=negwhack)
  
  summweighted <- weighted_analysis %>% group_by(nodeType) %>% summarize(pval=mean(p.value), tval=mean(statistic)) %>% mutate(method="weighted")
  
  repsummaries <- bind_rows(summdens, summstrength, summstrength_dcov, summweighted)
  repsummaries$replication <- s

  #replist[[s]] <- repsummaries
  repsummaries
}

stopCluster(clusterobj)

save(allreps, allreps_agg, replist,
     whack_increase, whack_decrease, neg_edgeshift_weight_wi_sd, neg_edgeshift_weight_bw_sd,
     pos_edgeshift_weight_wi_sd, pos_edgeshift_weight_bw_sd, file="whack_100reps_p14increase_p04decrease_p2fcshift.RData")

#load("whack_50reps_p1increase_p01decrease.RData")
#load("whack_50reps_p1increase_p01decrease_largevar.RData")
#load("whack_50reps_04increase_01decrease.RData")
#load("whack_50reps_10increase.RData")
#load(file="whack_100reps_p14increase_p02decrease.RData")
load(file="whack_100reps_p14increase_p04decrease.RData") #this is the one to use for paper
#load("whack_100reps_p14increase_p04decrease_p2fcshift.RData") #larger global FC spread

allreps <- bind_rows(replist) %>% ungroup()

allreps_agg <- allreps %>% group_by(thresh, nodeType, method) %>%
  summarize(#m_p = mean(pval), m_t=mean(tval), 
    m_p=median(pval), m_t=median(tval),
    mean_p=mean(pval), mean_t=mean(tval),
    se_p=plotrix::std.error(pval), se_t=plotrix::std.error(tval),
    sd_p=sd(pval), sd_t=sd(tval),
    # Lower=Hmisc::smean.cl.boot(tval, conf.int = .99)[2],
    # Upper=Hmisc::smean.cl.boot(tval, conf.int = .99)[3]
    #Lower = mean(tval) - sd(tval),
    #Upper = mean(tval) + sd(tval)
    Lower=quantile(tval, .10),
    Upper=quantile(tval, .90)
  ) %>% ungroup()

allreps_agg <- allreps_agg %>% ungroup() %>% mutate(thresh_num = as.character(round(as.numeric(sub(".*_", "", thresh)), 2)),
                                      nodeType = factor(nodeType, levels=c("Positive Target", "Negative Target", "Not Targeted"),
                                                        labels=c("Positive", "Negative", "Comparator")))

allreps <- allreps %>% ungroup() %>% mutate(thresh_num = as.numeric(sub(".*_", "", thresh)))


ptheme <- theme_bw(base_size=16) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size = rel(0.9), hjust = 0.0, margin=margin(b=2)),
        plot.margin = unit(c(1,1,1,1), "lines"))

common <- list(scale_color_brewer("Node Type", palette="Set1"),
  geom_hline(yintercept=c(-1.985, 1.985), color="gray60"),
  geom_hline(yintercept=c(0), color="black"),
  #geom_pointrange(size=1.0, position=position_dodge(width=0.5)),
  geom_crossbar(width=0.4, position=position_dodge(width=0.5)),
  #geom_pointrange(width=0.4, position=position_dodge(width=0.4)),
  coord_cartesian(ylim=c(-8, 14.5)))

dthresh <- allreps_agg %>% filter(method=="density" & thresh_num %in% c(.05, .15, .25)) #c(.05, .1, .15, .2, .25))

dthresh_foranalysis <- allreps %>% filter(method=="density")
dthresh_foranalysis %>% group_by(nodeType) %>% summarize(mean(tval), sd(tval), mean(pval))

xtabs(~nodeType + thresh, dthresh_foranalysis)

library(ez)
library(emmeans)
library(lme4)
ezANOVA(data=dthresh_foranalysis, wid=replication, within=.(thresh, nodeType), dv = .(tval))
mm <- lmer(tval ~ thresh*nodeType + (1|replication), dthresh_foranalysis)
car::Anova(mm)

#so, there is a linear trend across densities... maybe should just treat it as continuous
df <- summary(emmeans(mm, ~thresh*nodeType))

df <- summary(emmeans(mm, ~thresh*nodeType))


df$nodeType <- factor(df$nodeType, levels = c("Positive Target", "Negative Target", "Not Targeted"),
                      labels=c("Positive", "Negative", "Comparator"))
pdf("tstat_by_density.pdf", width=7, height=5)
ggplot(df, aes(x=thresh, y=emmean, ymin=lower.CL, ymax=upper.CL)) +
  geom_line(aes(group=1), size=1.0, color="grey60") + 
  geom_pointrange(fatten=.9, size=1.1) + facet_wrap(~nodeType, ncol=1, scales="free_y") +
  #ggtitle("Average t statistic as a function of density threshold and node type") +
  scale_x_discrete("Density threshold", labels=sort(unique(dthresh_foranalysis$thresh_num))) +
  ylab(expression(paste("Average group difference ", italic("t")))) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  theme(axis.title.x = element_text(margin=margin(t=10)),
        panel.spacing = unit(15, "pt"))
dev.off()

mmlin <- lmer(tval ~ thresh_num*nodeType + (1|replication), dthresh_foranalysis)
summary(mmlin)
car::Anova(mmlin)
summary(emtrends(mmlin, ~nodeType, var="thresh_num"), infer=c(TRUE, TRUE), adjust="none")

dplot <- ggplot(dthresh, aes(x=thresh_num, y=m_t, color=nodeType, ymin=Lower, ymax=Upper)) + 
  #geom_bar(stat="identity", width=0.7, position=position_dodge(0.7)) +
  #geom_errorbar(width=0.5, position=position_dodge(0.7)) +
  #geom_crossbar(width=0.5) +
  labs(x="Density", y="Mean t") + ptheme + common +
  ggtitle("Density threshold")

pdf("dpanel.pdf", width=6, height=6)
dplot
dev.off()

#FC thresholding

##analysis
rthresh_foranalysis <- allreps %>% filter(method=="strength")
rthresh_foranalysis %>% group_by(nodeType) %>% summarize(mean(tval), sd(tval), mean(pval), sd(pval), std.error(pval))
rmm <- lmer(tval ~ thresh*nodeType + (1|replication), rthresh_foranalysis)
car::Anova(rmm)
summary(emmeans(rmm, ~thresh))
summary(emmeans(rmm, ~nodeType))

##plot

rthresh <- allreps_agg %>% filter(method=="strength" & thresh_num %in% c(.2, .3, .4))
rthresh$thresh_num <- factor(rthresh$thresh_num, levels=c("0.4", "0.3", "0.2"))
rplot <- ggplot(rthresh, aes(x=thresh_num, y=m_t, color=nodeType, ymin=Lower, ymax=Upper)) + 
  #geom_bar(stat="identity", width=0.7, position=position_dodge(0.7)) +
  #geom_errorbar(width=0.5, position=position_dodge(0.7)) +
  #geom_crossbar() +
  labs(x="Pearson r", y="Mean t") + ptheme + common + ggtitle("FC Threshold")

pdf("rpanel.pdf", width=6, height=6)
rplot
dev.off()

#weighted

#analysis
wthresh_foranalysis <- allreps %>% filter(method=="weighted")
wthresh_foranalysis %>% group_by(nodeType) %>% summarize(mean(tval), sd(tval), mean(pval), sd(pval), std.error(pval))


#plot
wthresh <- allreps_agg %>% filter(method=="weighted")
wplot <- ggplot(wthresh, aes(x=NA, y=m_t, color=nodeType, ymin=Lower, ymax=Upper)) + 
  #geom_pointrange() +
  #geom_bar(stat="identity", width=0.7, position=position_dodge(0.7)) +
  #geom_errorbar(width=0.5, position=position_dodge(0.7)) +
  #geom_crossbar() +
  labs(x="", y="Mean t") + ptheme + common + ggtitle("Weighted") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

pdf("wpanel.pdf", width=6, height=6)
wplot
dev.off()

#panel for including density as covariate
dcovthresh_foranalysis <- allreps %>% filter(method=="strength_dcov")
dcovthresh_foranalysis %>% group_by(nodeType) %>% summarize(mean(tval), sd(tval), mean(pval), sd(pval), std.error(pval))

dcov_foranalysis <- allreps %>% filter(method=="strength_dcov")
dcov_foranalysis %>% group_by(nodeType) %>% summarize(mean(tval), sd(tval), mean(pval), sd(pval), std.error(pval))
rmm <- lmer(tval ~ thresh*nodeType + (1|replication), dcov_foranalysis)
car::Anova(rmm)
summary(emmeans(rmm, ~thresh))
summary(emmeans(rmm, ~nodeType))

mmlin <- lmer(tval ~ thresh_num*nodeType + (1|replication), dcov_foranalysis)
summary(mmlin)
car::Anova(mmlin)
summary(emtrends(mmlin, ~nodeType, var="thresh_num"), infer=c(TRUE, TRUE), adjust="none")



#plot

dcovthresh <- allreps_agg %>% filter(method=="strength_dcov" & thresh_num %in% c(.2, .3, .4))
dcovthresh$thresh_num <- factor(dcovthresh$thresh_num, levels=c("0.4", "0.3", "0.2"))

dcovplot <- ggplot(dcovthresh, aes(x=thresh_num, y=m_t, color=nodeType, ymin=Lower, ymax=Upper)) + 
  #geom_bar(stat="identity", width=0.7, position=position_dodge(0.7)) +
  #geom_errorbar(width=0.5, position=position_dodge(0.7)) +
  #geom_pointrange() +
  labs(x="Pearson r", y="Mean t") + ggtitle("FC Threshold\nDens. Cov.") + 
  ptheme + common

pdf("dcovpanel.pdf", width=7, height=6)
dcovplot
dev.off()

library(cowplot)
legend <- get_legend(wplot)
pbase <- plot_grid(dplot + theme(legend.position="none"),
                   rplot + theme(legend.position="none") + ylab(""),
                   dcovplot + theme(legend.position="none") + ggtitle("FC threshold, covary density"),
                     #annotate(geom="text", x=2, y=6, label="Density covaried"),
                   wplot + ylab(""),
                   labels = c("a)", "b)", "c)", "d)"), nrow=2,
                   align = 'h')

pdf("whack_combined_increase14_largedecrease.pdf", width=7, height=5)
#plot_grid(pbase, legend, rel_widths = c(2, .5), rel_heights = c(3, .2))
plot_grid(pbase)
dev.off()

# ggplot(filter(allreps_agg, method !="strength_dcov"), aes(x=thresh_num, y=m_t, fill=nodeType, ymin=m_t - se_t, ymax=m_t + se_t)) + 
#   geom_bar(stat="identity", position=position_dodge(0.5)) +
#   geom_errorbar(position=position_dodge(0.5), width=0.5) +
#   facet_wrap(~method, ncol=1, scales="free_x") +
#   geom_hline(yintercept=c(-1.985, 1.985))

# ggplot(repsummaries, aes(x=thresh, y=tval, fill=nodeType)) + geom_bar(stat="identity", position="dodge") + facet_wrap(~method, ncol=1) +
#   geom_hline(yintercept=c(-1.985, 1.985))
