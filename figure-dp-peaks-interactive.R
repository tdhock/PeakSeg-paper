works_with_R("3.1.1", reshape2="1.2.2", dplyr="0.2",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             "tdhock/animint@4b6a459b4e7425831ca4f3bb2b81f7305835252b")

load("dp.peaks.interactive.RData")
load("dp.peaks.regression.RData")

load("dp.peaks.baseline.RData")

roc <- do.call(rbind, dp.peaks.roc) %.%
  mutate(algorithm="PeakSegDP")
roc.chosen <- do.call(rbind, dp.peaks.roc.chosen) %.%
  mutate(algorithm="PeakSegDP")
broc <- do.call(rbind, dp.peaks.baseline.roc) %.%
  mutate(thresh=param.name,
         algorithm=sub(".trained", "", algorithm))
broc.chosen <- do.call(rbind, dp.peaks.baseline.chosen) %.%
  mutate(thresh=param.name,
         algorithm=sub(".trained", "", algorithm))

default.params <-
  c(macs.trained="1.30103", hmcan.broad.trained="2.30258509299405")

default.list <- list()
for(curve.name in names(dp.peaks.baseline.roc)){
  one.curve <- dp.peaks.baseline.roc[[curve.name]]
  default.param <- default.params[[as.character(one.curve$algorithm[1])]]
  one.param <- subset(one.curve, param.name == default.param) %.%
    mutate(algorithm=sub(".trained", "", algorithm),
           thresh=param.name)
  default.list[[curve.name]] <- one.param
}

default.df <- do.call(rbind, default.list)

str(roc)
str(roc.chosen)
str(broc)
str(broc.chosen)

roc.cols <- c("set.name", "set.i", "FPR", "TPR", "algorithm", "thresh")

all.roc <- rbind(roc[, roc.cols], broc[, roc.cols]) %.%
  mutate(set.name=as.character(set.name),
         testSet=paste(set.name, "split", set.i))
all.dots <-
  rbind(data.frame(roc.chosen[, roc.cols], model.type="trained"),
        data.frame(broc.chosen[, roc.cols], model.type="trained"),
        data.frame(default.df[, roc.cols], model.type="default")) %.%
    mutate(set.name=as.character(set.name),
           testSet=paste(set.name, "split", set.i))

longErrors <- dp.peaks.interactive$regionErrors %.%
  group_by(testSet, test.chunk, algorithm) %.%
  summarise(fp=sum(fp),
            fn=sum(fn)) %.%
  group_by(testSet) %.%
  mutate(chunk.x=as.numeric(factor(as.character(test.chunk))),
         errors=fp+fn)
longErrors
fp.fn <- melt(longErrors, measure.vars=c("fp", "fn"),
              variable.name="error.type",
              value.name="errors")
wideErrors <-
  dcast(longErrors,
        testSet + test.chunk ~ algorithm,
        value.var="errors")
comp <- function(algorithm){
  wideErrors$competitor <- algorithm
  wideErrors$competitor.errors <- wideErrors[,algorithm]
  wideErrors
}
competitors <- rbind(comp("macs.trained"), comp("hmcan.broad.trained"))

chunkErrors <- dp.peaks.interactive$regionErrors %.%
  group_by(testSet, test.chunk, algorithm, sample.id) %.%
  summarise(fp=sum(fp),
            fn=sum(fn)) %.%
  mutate(errors=fp+fn)
others <- chunkErrors %.%
  group_by(testSet, test.chunk, sample.id) %.%
  filter(algorithm != "PeakSeg") %.%
  summarise(min.other=min(errors),
            max.other=max(errors))
chunk.others <- inner_join(chunkErrors, others) 
chunk.diffs <- chunk.others %.%
  filter(algorithm=="PeakSeg") %.%
  mutate(diff1=errors-min.other,
         diff2=errors-max.other) %.%
  arrange(diff1, errors, diff2) %.%
  group_by(testSet, test.chunk) %.%
  mutate(rank=seq_along(diff1))
algo.order <- c("hmcan.broad.trained", "macs.trained", "PeakSeg")
afactor <- function(x)factor(x, algo.order)
algo.y <- -0.1*(3:1)
names(algo.y) <- algo.order
chunk.ordered <-
  inner_join(chunk.diffs, chunkErrors,
             c("testSet", "test.chunk", "sample.id")) %.%
  mutate(algorithm=afactor(algorithm.y),
         errors=errors.y)
dp.peaks.interactive$regionErrors %.%
  group_by(algorithm) %.%
  summarise(regions=n())

regions <- dp.peaks.interactive$regionErrors %.%
  filter(algorithm == algorithm[1]) 

regionErrors <- dp.peaks.interactive$regionErrors %.%
  mutate(algo.y=algo.y[as.character(algorithm)])

peaks <- dp.peaks.interactive$peaks %.%
  mutate(algo.y=algo.y[as.character(algorithm)])

## TODO: label parameter values as well.
peakLabels <- data.frame(algo.y, algorithm=names(algo.y))

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

baselines <- dp.peaks.interactive$baselineParams %.%
  mutate(competitor=algorithm)

test.pred <- dp.peaks.prediction %.%
  filter(set=="test") 

train.pred <- dp.peaks.prediction %.%
  filter(set=="train")

hgTracks <- "http://ucscbrowser.genap.ca/cgi-bin/hgTracks"

psplits <- dp.peaks.interactive$splitErrors %.%
  filter(algorithm=="PeakSeg")

cat(deparse(RColorBrewer::brewer.pal(3, "Dark2")))
algo.colors <-
  c(hmcan.broad="#1B9E77",
    macs="#D95F02",
    PeakSegDP="#7570B3",
    hmcan.broad.trained="#1B9E77",
    macs.trained="#D95F02",
    PeakSeg="#7570B3")
viz <-
  list(splits=ggplot()+
       scale_color_manual(values=algo.colors, guide="none")+
       geom_point(aes(percent, algorithm, clickSelects=testSet,
                      color=algorithm),
                  data=dp.peaks.interactive$splitErrors, size=5, alpha=0.6)+
       facet_grid(.~set.name)+
       ggtitle("select data set and train/test split")+
       xlab("percent test error")+
       theme_animint(width=1200, height=150),
       
       samples=ggplot()+
       scale_color_manual(values=algo.colors, guide="none")+
       theme(axis.line.x=element_blank(),
             axis.text.x=element_blank(),
             axis.ticks.x=element_blank())+
       theme_animint(height=150)+
       xlab("<- PeakSeg better ----- samples ----- PeakSeg worse ->")+
       ylab("algorithm")+
       ggtitle("error counts, select sample")+
       geom_text(aes(rank, algorithm, label=errors,
                     color=algorithm,
                     clickSelects=sample.id,
                     showSelected=testSet,
                     showSelected2=test.chunk),
                 data=chunk.ordered),

       competitors=ggplot()+
       geom_point(aes(competitor.errors, PeakSeg,
                      clickSelects=test.chunk,
                      showSelected=testSet),
                  data=competitors, size=5, alpha=0.6)+
       geom_text(aes(80, 100, color=algorithm,
                     label=paste0("selected parameter=", param.name),
                     showSelected=testSet),
                 data=baselines)+
       scale_color_manual(values=algo.colors, guide="none")+
       scale_x_continuous(paste("incorrectly predicted annotated regions",
                                "(competitor)"),
                          breaks=seq(0, 200, by=25))+
       scale_y_continuous("incorrectly predicted annotated regions (PeakSeg)",
                          breaks=seq(0, 200, by=25))+
       ggtitle("Test error comparison, select test chunk")+
       theme_animint(width=500, height=300)+
       facet_grid(.~competitor, labeller=label_both)+
       coord_equal()+
       geom_segment(aes(x, y, xend=xend, yend=yend),
                   data=data.frame(x=0, y=0, xend=100, yend=100)),

       roc=ggplot()+
       theme_animint(width=300, height=300)+
       geom_text(aes(0.5, 1.01, label=testSet,
                     showSelected=testSet),
                 data=psplits)+
       geom_path(aes(FPR, TPR, color=algorithm, group=algorithm,
                     showSelected=testSet),
                 data=all.roc)+
       geom_point(aes(FPR, TPR, color=algorithm, fill=model.type,
                      showSelected=testSet),
                  pch=21, data=all.dots)+
       scale_fill_manual(values=c(trained="white", default="black"))+
       coord_equal()+
       scale_color_manual(values=algo.colors)+
       scale_x_continuous("False positive rate",
                          breaks=seq(0, 1, by=1/2))+
       ggtitle("Test ROC curves")+
       scale_y_continuous("True positive rate",
                          breaks=seq(0, 1, by=1/2)),
       
       extrapolation=ggplot()+
       theme_animint(width=300, height=300)+
       geom_text(aes(5, 16.5, label=testSet,
                     showSelected=testSet),
                 data=psplits)+
       ggtitle("PeakSeg model feature space")+
       geom_contour(aes(log.max.coverage, log.total.weight, z=log.lambda,
                        color=..level..,
                        showSelected=testSet),
                    data=dp.peaks.grid)+
       geom_point(aes(log.max.coverage, log.total.weight,
                      fill=set,
                      showSelected=testSet),
                  data=train.pred, pch=21, alpha=1/10)+
       geom_polygon(aes(log.max.coverage, log.total.weight, group=test.chunk,
                        showSelected=testSet,
                        clickSelects=test.chunk),
                    data=dp.peaks.polygon, alpha=6/10, size=3)+
       ## geom_segment(aes(log.max.coverage, log.total.weight,
       ##                  xend=log.max.coverage.end, yend=log.total.weight.end,
       ##                  showSelected=testSet,
       ##                  clickSelects=test.chunk),
       ##              data=dp.peaks.segment, alpha=6/10, size=3)+
       xlab("log.max.coverage feature")+
       ylab("log.total.weight feature")+
       geom_text(aes(log.max.coverage, log.total.weight,
                     label=test.errors,
                     showSelected=testSet,
                     showSelected2=test.chunk,
                     clickSelects=sample.id),
                 data=test.pred)+
       geom_blank(aes(log.max.coverage, log.total.weight,
                      fill=set,
                      showSelected=testSet,
                      showSelected2=test.chunk),
                 data=test.pred)+
       guides(fill=guide_legend(
                override.aes=list(alpha=1/2)))+
       ##scale_color_gradient("log.penalty", low="white", high="blue")+
       scale_color_gradient2("log.penalty")+
       scale_fill_manual(values=c(train="red", test="black")),

       profile=ggplot()+
       geom_text(aes(0, 1, label=sprintf("%s max=%d", sample.id, max.coverage),
                     showSelected=test.chunk,
                     showSelected2=sample.id),
                 data=dp.peaks.interactive$profile.info, hjust=0)+
       geom_text(aes(0.5, 1,
                     label=sprintf("%.1fkb on %s", bases/1e3, chunkChrom),
                     href=sprintf("%s?db=hg19&position=%s:%d-%d",
                       hgTracks, chunkChrom, expandStart, expandEnd),
                     showSelected=test.chunk),
                 data=dp.peaks.interactive$chunks)+
       theme(axis.line=element_blank(),
             axis.text=element_blank(),
             axis.ticks=element_blank())+
       theme_animint(width=1300)+
       geom_rect(aes(xmin=chromStartNorm, xmax=chromEndNorm,
                     ymin=-0.4, ymax=1, 
                     fill=annotation,
                     showSelected=testSet,
                     showSelected2=test.chunk,
                     showSelected3=sample.id),
                 data=regions, linetype="none")+
       geom_rect(aes(xmin=chromStartNorm, xmax=chromEndNorm,
                     ymin=algo.y-0.025, ymax=algo.y+0.025, 
                     fill=annotation,
                     linetype=status,
                     showSelected=testSet,
                     showSelected2=test.chunk,
                     showSelected3=sample.id),
                 color="black", size=2,
                 data=regionErrors)+
       geom_segment(aes(ifelse(chromStartNorm < 0, 0, chromStartNorm),
                        algo.y,
                        xend=ifelse(chromEndNorm > 1, 1, chromEndNorm),
                        yend=algo.y,
                        showSelected=testSet,
                        showSelected2=test.chunk,
                        showSelected3=sample.id),
                    data=peaks,
                    color="deepskyblue")+
       geom_point(aes(ifelse(chromEndNorm > 1, 1, chromEndNorm),
                        algo.y,
                        showSelected=testSet,
                        showSelected2=test.chunk,
                        showSelected3=sample.id),
                    data=peaks,
                    color="deepskyblue")+
       geom_text(aes(0, algo.y, label=sub(".trained", "", algorithm),
                     color=algorithm),
                 data=peakLabels,
                 hjust=1)+
       scale_color_manual(values=algo.colors, guide="none")+
       guides(linetype=guide_legend(order=2,
                override.aes=list(fill="white")))+
       scale_linetype_manual("error type",
                             values=c(correct="none",
                               "false negative"="dashed",
                               "false positive"="solid"))+
       ylab("aligned read counts")+
       xlab("position on chromosome")+
       scale_fill_manual(values=ann.colors)+
       geom_line(aes(baseNorm, coverageNorm,
                     showSelected=test.chunk,
                     showSelected2=sample.id),
                 data=dp.peaks.interactive$profiles,
                 color="grey"),

       first=list(testSet="H3K36me3_TDH_immune split 5",
         test.chunk="H3K36me3_TDH_immune/3",
         sample.id="McGill0101"),

       title="Peak detection model comparison")

stopifnot(length(dp.peaks.segment) == 0)

with.facets <- viz$extrapolation+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(set.name ~ set.i, labeller=function(var, val){
    gsub("_", "\n", val)
  })

animint2dir(viz, "figure-dp-peaks-interactive")
##animint2gist(viz)
