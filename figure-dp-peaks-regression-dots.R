works_with_R("3.1.2", reshape2="1.2.2", ggplot2="1.0",
             data.table="1.9.4",
             dplyr="0.4.0",
             xtable="1.7.3")

load("dp.peaks.baseline.RData")
load("dp.peaks.regression.RData")
load("regularized.all.RData")
load("unsupervised.error.RData")
load("oracle.regularized.RData")
load("multires.bins.RData")

bins <- 
multires.bins$test.error %>%
  group_by(set.name, set.i) %>%
  summarise(errors=sum(fp+fn),
            regions=n()) %>%
  mutate(percent=errors/regions*100,
         algorithm="multires.bins")
  
reg <- regularized.all$error %>%
  filter(model.name %in% c("L1.reg", "log.bases.log.max")) %>%
  group_by(set.name, set.i, model.name) %>%
  summarise(errors=sum(errors),
            regions=sum(regions)) %>%
  mutate(percent=errors/regions*100,
         algorithm=model.name)

oreg <- oracle.regularized$error %>%
  filter(model.name %in% c("L1.reg", "log.bases.log.max")) %>%
  group_by(set.name, set.i, model.name) %>%
  summarise(errors=sum(errors),
            regions=sum(regions)) %>%
  mutate(percent=errors/regions*100,
         algorithm=paste0("oracle.", model.name))

regression.set.i <- dp.peaks.regression %>%
  group_by(set.name, set.i) %>%
  summarise(errors=sum(errors),
            regions=sum(regions)) %>%
  mutate(algorithm="PeakSeg (cDPA)")
both.cols <- c("set.name", "set.i", "algorithm", "errors", "regions")

un <- unsupervised.error %>%
  group_by(set.name, set.i, algorithm) %>%
  summarise(errors=sum(errors),
            regions=sum(regions)) %>%
  mutate(percent=errors/regions*100)
bad.algos <- c("none", "AIC", "BIC", "mBIC")
un.bad <- un %>%
  filter(algorithm %in% bad.algos) %>%
  group_by(set.name, set.i) %>%
  summarise(min=min(errors),
            max=max(errors),
            algos=n())
## make sure AIC/BIC/mBIC/none are really all the same.
with(un.bad, stopifnot(min == max))

un.other <- un %>%
  filter(!algorithm %in% c(bad.algos, "grid loss"))
un.bad <- un %>%
  filter(algorithm == "none") %>%
  mutate(algorithm="AIC/BIC.0")
un.show <-
  rbind(un.other,
        un.bad)

ann <- function(df, parameters, supervision){
  parameters <- ifelse(df$algorithm == "oracle.trained", "1", parameters)
  data.table(df,
             parameters=factor(parameters,
                    c(">1", "1", "0")),
             supervision)
}
alg.rep <-
  c("oracle"="oracle.0",
    "grid lik"="AIC/BIC.1",
    "macs.default"="macs.0",
    "macs.trained"="macs.1",
    "hmcan.broad.default"="hmcan.broad.0",
    "hmcan.broad.trained"="hmcan.broad.1",
    L1.reg="AIC/BIC.41",
    log.bases.log.max="AIC/BIC.3",
    oracle.L1.reg="oracle.41",
    oracle.log.bases.log.max="oracle.3",
    "grid oracle"="oracle.1")
both <-
  rbind(##regression.set.i[, both.cols],
        ann(data.frame(un.show)[, both.cols],
            "0", "unsupervised"),
        ann(data.frame(reg)[, both.cols],
            ">1", "supervised"),
        ann(data.frame(oreg)[, both.cols],
            ">1", "supervised"),
        ann(dp.peaks.baseline[, both.cols], "1", "supervised"),
        ann(data.frame(bins)[, both.cols], ">1", "supervised")
        ) %>%
  group_by() %>%
  mutate(percent=errors/regions*100,
         set.name=as.character(set.name),
         algorithm=ifelse(algorithm %in% names(alg.rep),
           alg.rep[as.character(algorithm)], as.character(algorithm)),
         algo.type=ifelse(grepl("hmcan|macs", algorithm),
           "baselines", "PeakSeg (cDPA)"),
         algo.type=factor(algo.type, c("PeakSeg (cDPA)",
           "baselines")),
         learning=ifelse(grepl("[.]0$", algorithm), "unsupervised",
           ifelse(grepl("[.]1$", algorithm), "grid\nsearch",
                  ifelse(grepl("best", algorithm), "cheating",
                         "interval\nregression"))))

dot.counts <- with(both, table(algorithm, set.name))
stopifnot(dot.counts == 6)

region.sums <- both %>%
  group_by(set.name, set.i, algorithm) %>%
  summarise(regions=sum(regions)) %>%
  arrange(set.name, set.i, algorithm)
region.stats <- region.sums %>%
  group_by(set.name, set.i) %>%
  summarise(min.regions=min(regions),
            max.regions=max(regions))
stopifnot(with(region.stats, min.regions == max.regions))

wide <- dcast(both, set.name + set.i ~ algorithm, value.var="errors")
scatter <- 
ggplot()+
  geom_abline(aes(intercept=intercept, slope=slope),
              color="grey",
              data=data.frame(intercept=0, slope=1))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(. ~ set.name, labeller=function(var, val){
    gsub("_", "\n", val)
  })+
  coord_equal()

## Scatterplots show that L1.reg does not always have lower test error
## than competitors, but it does at least for a majority of train/test
## splits on each data set.
## scatter+
##   geom_point(aes(hmcan.broad.trained, L1.reg),
##              data=wide, pch=1)
## scatter+
##   geom_point(aes(macs.trained, L1.reg),
##              data=wide, pch=1)
## scatter+
##   geom_point(aes(log.bases.log.max, L1.reg),
##              data=wide, pch=1)

algo.colors <-
  c("cheating"="grey",
    "interval\nregression"="#D95F02",
    "grid\nsearch"="#1B9E77",
    "unsupervised"="#7570B3")
both$algo <- sub(".trained", "", both$algorithm)
both$algo <- both$algorithm

mean.both <- both %>%
  group_by(set.name, parameters, learning, supervision, algo, algo.type) %>%
  summarise(mean=mean(percent))
algo.stats <- mean.both %>%
  group_by(algo) %>%
  summarise(mean=mean(mean)) %>%
  arrange(mean)
algo.levs <- rev(algo.stats$algo)
algo.levs <-
  rev(c("best.DP",
        "multires.bins",
        ## "oracle.41", "AIC/BIC.41",
        ## "oracle.3", "AIC/BIC.3",
        ## "oracle.1", "AIC/BIC.1",
        ## "oracle.0", "AIC/BIC.0",
        "oracle.41", "oracle.3", "oracle.1", "oracle.0",
        "AIC/BIC.41", "AIC/BIC.3", "AIC/BIC.1", "AIC/BIC.0",
        "hmcan.broad.1", "hmcan.broad.0",
        "macs.1", "macs.0"))
both$algo <- factor(both$algo, algo.levs)
mean.both$algo <- factor(mean.both$algo, algo.levs)

## scatterplot comparing oracle.trained and L1.reg.
ref.diff <- both %>%
  filter(algo=="best.DP") %>%
  mutate(baseline=percent) %>%
  select(set.name, set.i, baseline) %>%
  inner_join(both, c("set.name", "set.i")) %>%
  mutate(diff=percent-baseline)

both.more <- both %>%
  mutate(algo.type=sub("[.].*$", "", algo),
         learned=sub("^.*[.]", "", algo))
oracle.bic <- both.more %>%
  filter(algo.type %in% c("oracle", "AIC/BIC"))
both.wide <-
  dcast(oracle.bic,
        set.name + set.i + learned ~ algo.type,
        value.var = "percent")
scatter.best <- 
scatter+
  facet_grid(learned ~ set.name, labeller=function(var, val){
    gsub("_", "\n", val)
  })+
  geom_point(aes(`AIC/BIC`, oracle),
             data=both.wide, pch=1)
pdf("figure-dp-peaks-regression-scatter.pdf", h=2)
print(scatter.best)
dev.off()

## scatter.oracle <- 
## scatter+
##   geom_point(aes(oracle.trained, oracle.default),
##              data=both.wide, pch=1)
## pdf("figure-dp-peaks-regression-scatter-oracle.pdf", h=3)
## print(scatter.oracle)
## dev.off()

dp.peaks.regression.dots <- both
save(dp.peaks.regression.dots, file="dp.peaks.regression.dots.RData")
boxes <- 
ggplot()+
  geom_boxplot(aes(algo, percent),
               data=both)+
  facet_grid(.~set.name, labeller=function(var, val){
    gsub("_", "\n", val)
  })+
  ylab("algorithm")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  coord_flip()+
  scale_y_continuous("percent test error", breaks=seq(0, 50, by=25))
pdf("figure-dp-peaks-regression-boxes.pdf", h=1.6)
print(boxes)
dev.off()

boxes.flip <- 
ggplot()+
  geom_boxplot(aes(algo, percent),
               data=both)+
  facet_grid(.~set.name, labeller=function(var, val){
    gsub("_", "\n", val)
  })+
  ylab("algorithm")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  scale_y_continuous("percent test error", breaks=seq(0, 50, by=25))
pdf("figure-dp-peaks-regression-boxes-flip.pdf", h=1.6)
print(boxes.flip)
dev.off()

best.diff <- ref.diff %>%
  group_by(set.name, algo) %>%
  summarise(mean=mean(diff)) %>%
  filter(algo != "best.DP") %>%
  group_by(set.name) %>%
  summarise(min=min(mean))

dots <-  #with differences.
ggplot()+
  ##geom_vline(xintercept=0, color="grey")+
  geom_vline(aes(xintercept=min), data=best.diff)+
  geom_point(aes(diff, algo, color=learning),
             data=ref.diff, pch=1)+
  facet_grid(algo.type ~ set.name, labeller=function(var, val){
    gsub("_", "\n", val)
  }, scales="free_y", space="free_y")+
  scale_y_discrete("algorithm")+
  theme_bw()+
  guides(color=guide_legend())+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position="top")+
  scale_color_manual("learning\nalgorithm", values=algo.colors,
                     breaks=names(algo.colors))+
  scale_x_continuous("test error difference from best possible DP model",
                     breaks=seq(0, 100, by=20))
##pdf("figure-dp-peaks-regression-dots-diff.pdf", h=3, w=8)
##print(dots)
##dev.off()

un.both <- both %>%
  filter(grepl("[.][0]$", algo)) %>%
  mutate(algo.type=sub(" ", "\n", algo.type),
         algo.fac=factor(algo.type, c("PeakSeg\n(cDPA)", "baselines")))
un.mean <- un.both %>%
  group_by(set.name, algo, algo.fac, learning) %>%
  summarise(mean=mean(percent))
un.best <- un.mean %>%
  group_by(set.name) %>%
  summarise(min=min(mean),
            algo=algo[which.min(mean)],
            learning=learning[which.min(mean)])

un.dots <-  #unsupervised
ggplot()+
  geom_vline(aes(xintercept=min), data=un.best)+
  geom_point(aes(mean, algo, color=learning),
             data=un.mean, alpha=0.25, size=4)+
  geom_point(aes(percent, algo, color=learning),
             data=un.both, pch=1)+
  facet_grid(algo.fac ~ set.name, labeller=function(var, val){
    gsub("_", "\n", val)
  }, scales="free_y", space="free_y")+
  scale_y_discrete("algorithm . parameters learned")+
  theme_bw()+
  guides(color=guide_legend())+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position="top")+
  scale_color_manual("learning\nalgorithm", values=algo.colors,
                     breaks=names(algo.colors))+
  scale_fill_manual("learning\nalgorithm", values=algo.colors,
                     breaks=names(algo.colors))+
  ##scale_x_continuous("percent test error", breaks=seq(0, 50, by=25))
  scale_x_continuous("percent incorrect peak region labels (test error)",
                     breaks=seq(0, 100, by=20),
                     limits=c(0, NA))
pdf("figure-dp-peaks-regression-dots-unsupervised.pdf", h=3.5, w=8)
print(un.dots)
dev.off()

grid.both <- both %>%
  filter(grepl("[.][01]$", algo))
grid.mean <- mean.both %>%
  filter(grepl("[.][01]$", algo))
grid.best <- grid.mean %>%
  group_by(set.name) %>%
  summarise(min=min(mean),
            algo=algo[which.min(mean)],
            learning=learning[which.min(mean)])

grid.dots <-  #with 1 set of facets.
ggplot()+
  geom_vline(aes(xintercept=min), data=grid.best)+
  geom_point(aes(mean, algo, color=learning),
             data=grid.mean, alpha=0.25, size=4)+
  geom_point(aes(percent, algo, color=learning),
             data=grid.both, pch=1)+
  facet_grid(algo.type ~ set.name, labeller=function(var, val){
    gsub("[ _]", "\n", val)
  }, scales="free_y", space="free_y")+
  scale_y_discrete("algorithm . parameters learned")+
  theme_bw()+
  guides(color=guide_legend())+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position="top")+
  scale_color_manual("learning\nalgorithm", values=algo.colors,
                     breaks=names(algo.colors))+
  scale_fill_manual("learning\nalgorithm", values=algo.colors,
                     breaks=names(algo.colors))+
  ##scale_x_continuous("percent test error", breaks=seq(0, 50, by=25))
  scale_x_continuous("percent incorrect peak region labels (test error)",
                     breaks=seq(0, 100, by=20))
pdf("figure-dp-peaks-regression-dots-grid.pdf", h=4.1, w=8)
print(grid.dots)
dev.off()

best.percent <- mean.both %>%
  filter(algo != "best.DP") %>%
  group_by(set.name) %>%
  summarise(min=min(mean),
            algo=algo[which.min(mean)],
            learning=learning[which.min(mean)])

dots <-  #with 1 set of facets.
ggplot()+
  geom_vline(aes(xintercept=min), data=best.percent)+
  geom_point(aes(mean, algo, color=learning),
             data=mean.both, alpha=0.25, size=4)+
  ## geom_point(aes(min, algo, fill=learning),
  ##            data=best.percent, pch=21, size=4, color="black")+
  geom_point(aes(percent, algo, color=learning),
             data=both, pch=1)+
  facet_grid(algo.type ~ set.name, labeller=function(var, val){
    gsub("_", "\n", val)
  }, scales="free_y", space="free_y")+
  scale_y_discrete("model . parameters learned")+
  theme_bw()+
  guides(color=guide_legend())+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position="top")+
  scale_color_manual("learning\nalgorithm", values=algo.colors,
                     breaks=names(algo.colors))+
  scale_fill_manual("learning\nalgorithm", values=algo.colors,
                     breaks=names(algo.colors))+
  ##scale_x_continuous("percent test error", breaks=seq(0, 50, by=25))
  scale_x_continuous("percent incorrect peak region labels (test error)",
                     breaks=seq(0, 100, by=20))
pdf("figure-dp-peaks-regression-dots.pdf", h=4.2, w=8)
print(dots)
dev.off()
