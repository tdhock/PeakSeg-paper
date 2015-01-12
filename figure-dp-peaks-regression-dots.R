works_with_R("3.1.1", reshape2="1.2.2", ggplot2="1.0", dplyr="0.2",
             xtable="1.7.3")

load("dp.peaks.baseline.RData")
load("dp.peaks.regression.RData")

regression.set.i <- dp.peaks.regression %.%
  group_by(set.name, set.i) %.%
  summarise(errors=sum(errors),
            regions=sum(regions)) %.%
  mutate(algorithm="PeakSeg")
both.cols <- c("set.name", "set.i", "algorithm", "errors", "regions")

both <- rbind(regression.set.i[, both.cols], dp.peaks.baseline[, both.cols]) %.%
  group_by() %.%
  mutate(percent=errors/regions*100,
         set.name=as.character(set.name))
         
wide <- dcast(both, set.name + set.i ~ algorithm, value.var="errors")

competitors <-
  c(hmcan.broad.trained="H3K36me3",
    macs.trained="H3K4me3")

br <- c(0, 50, 150)
for(algorithm in names(competitors)){
  experiment <- competitors[[algorithm]]
  only <- subset(wide, grepl(experiment, set.name))
errplot <- 
ggplot()+
  geom_point(aes_string(algorithm, "PeakSeg"),
             data=only, pch=1)+
  coord_equal()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(.~set.name)+
  scale_y_continuous(breaks=br)+
  scale_x_continuous(breaks=br)+
  geom_abline(aes(slope=slope, intercept=intercept),
              data=data.frame(slope=1, intercept=0))
  short <- sub("[.].*", "", algorithm)
  pdf.name <- sprintf("figure-dp-peaks-regression-%s.pdf", short)
  print(pdf.name)
  pdf(pdf.name, h=3, w=8)
  print(errplot)
  dev.off()
}

algo.colors <-
  c(macs.trained="#1B9E77", PeakSeg="#D95F02", hmcan.broad.trained="#7570B3")
both$algo <- sub(".trained", "", both$algorithm)
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
pdf("figure-dp-peaks-regression-dots.pdf", h=1.6)
print(dots)
dev.off()
