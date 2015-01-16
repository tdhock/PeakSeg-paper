works_with_R("3.1.2", reshape2="1.2.2", ggplot2="1.0",
             "Rdatatable/data.table@fd5464e4c587d4a96fa32c809eb94a45bb361133",
             dplyr="0.4.0",
             xtable="1.7.3")

load("dp.peaks.baseline.RData")
load("dp.peaks.regression.RData")
load("regularized.all.RData")

reg <- regularized.all$error %>%
  filter(model.name %in% c("L1.reg", "log.bases.log.max")) %>%
  group_by(set.name, set.i, model.name) %>%
  summarise(errors=sum(errors),
            regions=sum(regions)) %>%
  mutate(percent=errors/regions*100,
         algorithm=model.name)

regression.set.i <- dp.peaks.regression %>%
  group_by(set.name, set.i) %>%
  summarise(errors=sum(errors),
            regions=sum(regions)) %>%
  mutate(algorithm="PeakSeg")
both.cols <- c("set.name", "set.i", "algorithm", "errors", "regions")

both <- rbind(##regression.set.i[, both.cols],
              dp.peaks.baseline[, both.cols],
              data.frame(reg)[, both.cols]) %>%
  group_by() %>%
  mutate(percent=errors/regions*100,
         set.name=as.character(set.name))

wide <- dcast(both, set.name + set.i ~ algorithm, value.var="errors")

scatter <- 
ggplot()+
  geom_abline(aes(intercept=intercept, slope=slope),
              color="grey",
              data=data.frame(intercept=0, slope=1))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(. ~ set.name)+
  coord_equal()

## Scatterplots show that L1.reg does not always have lower test error
## than competitors, but it does at least for a majority of train/test
## splits on each data set.
scatter+
  geom_point(aes(hmcan.broad.trained, L1.reg),
             data=wide, pch=1)
scatter+
  geom_point(aes(macs.trained, L1.reg),
             data=wide, pch=1)

algo.colors <-
  c(macs.trained="#1B9E77", PeakSeg="#D95F02", hmcan.broad.trained="#7570B3")
both$algo <- sub(".trained", "", both$algorithm)

mean.both <- both %>%
  group_by(set.name, algo) %>%
  summarise(mean=mean(percent))
algo.stats <- mean.both %>%
  group_by(algo) %>%
  summarise(mean=mean(mean)) %>%
  arrange(mean)
algo.levs <- rev(algo.stats$algo)
both$algo <- factor(both$algo, algo.levs)
mean.both$algo <- factor(mean.both$algo, algo.levs)

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

dots <- 
ggplot()+
  geom_point(aes(algo, mean),
             data=mean.both, color="grey", size=5)+
  geom_point(aes(algo, percent),
             data=both, pch=1)+
  facet_grid(.~set.name, labeller=function(var, val){
    gsub("_", "\n", val)
  })+
  xlab("algorithm")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  coord_flip()+
  scale_y_continuous("percent test error", breaks=seq(0, 50, by=25))
pdf("figure-dp-peaks-regression-dots.pdf", h=3)
print(dots)
dev.off()
