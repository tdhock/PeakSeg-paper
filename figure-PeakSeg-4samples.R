works_with_R("3.1.1",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             dplyr="0.3.0.2",
             PeakSeg="2014.11.24")

load("PeakSeg4samples.RData")

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

error.list <- PeakSeg4samples$regions
errors <- NULL
for(peaks in 0:4){
  rname <- paste0("PeakSeg", peaks)
  e <- error.list[[rname]] %>%
    group_by(sample.id) %>%
    summarise(errors=sum(fp+fn))
  errors <- rbind(errors, data.frame(peaks, e))
}
errplot <- 
ggplot(errors, aes(peaks, errors))+
  geom_step()+
  geom_point()+
  theme_bw()+
  facet_grid(sample.id ~ .)+
  scale_y_continuous("incorrect annotated regions", minor_breaks=NULL)+
  scale_x_continuous("model complexity (peaks)", minor_breaks=NULL)+
  theme(panel.margin=grid::unit(0, "cm"))+
  coord_cartesian(ylim=c(-0.5, 3.5))
pdf("figure-PeakSeg-4samples-errors.pdf", h=5, w=7)
print(errplot)
dev.off()

zero.errors <- subset(errors, errors==0)
errtext <- errplot+
  geom_text(aes(2, 2, label=label), color="green",
            data=data.frame(sample.id=sprintf("McGill%04d", c(2, 4, 91, 322)),
              label=paste("best =", c("{2}", "{2}", "{1, 2, 3, 4}", "{1, 2}"))))+
  geom_point(color="green", data=zero.errors)
pdf("figure-PeakSeg-4samples-errors-text.pdf", h=5, w=7)
print(errtext)
dev.off()

model.regions <- error.list[[1]]

just.regions <- 
  ggplot()+
  scale_y_continuous("aligned read coverage",
                     labels=function(x){
                       sprintf("%.1f", x)
                     },
                     breaks=function(limits){
                       limits[2]
                     })+
  xlab("position on chr11 (kilo base pairs)")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ ., scales="free")+
  scale_fill_manual("annotation", values=ann.colors,
                    breaks=names(ann.colors))+
  coord_cartesian(xlim=c(118090, 118125))+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    fill=annotation),
                data=model.regions, alpha=1/2, color="grey")+
  guides(linetype=guide_legend(order=2,
           override.aes=list(fill="white")))+
  geom_step(aes(chromStart/1e3, coverage),
            data=PeakSeg4samples$signal, color="grey50")

png("figure-4samples-just-regions.png",
    units="in", res=200, width=7, height=4)
print(just.regions)
dev.off()

