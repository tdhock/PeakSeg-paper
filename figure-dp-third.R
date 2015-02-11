works_with_R("3.1.2",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             "tdhock/PeakSegDP@5bcee97f494dcbc01a69e0fe178863564e9985bc",
             "Rdatatable/data.table@200b5b40dd3b05112688c3a9ca2dd41319c2bbae",
             reshape2="1.2.2",
             dplyr="0.4.0")

chunk.name <- "H3K36me3_AM_immune/8"
chunk.name <- "H3K4me3_PGP_immune/7"

counts.file <- file.path("data", chunk.name, "counts.RData")
load(counts.file)

counts.list <- split(counts, counts$sample.id)
sample.id <- "McGill0026"
sample.counts <- counts.list[[sample.id]]
sample.counts$weight <- with(sample.counts, chromEnd-chromStart)
seg.info <- with(sample.counts, PoissonSeg(coverage, weight, maxSegments=3L))
loss.mat <- t(seg.info$res1)
end.mat <- t(seg.info$res2)
mean.mat <- t(seg.info$res3)

## need to compute best ends ourself, if we want to show
## feasible/infeasible.

## cumsum1 <- with(sample.counts, cumsum(coverage * weight))
## best.mean1 <- cumsum1/cumsum(sample.counts$weight)
## best.loss1 <- cumsum1 * (1-log(best.mean1))
## stopifnot(all.equal(best.mean1, mean.mat[,1]))
## to <- 1000
## to.data <- sample.counts[1:to,]
## to.mean <- with(to.data, sum(coverage * weight)/sum(weight))
## stopifnot(all.equal(to.mean, best.mean1[to]))
## to.loss <- with(to.data, PoissonLoss(coverage, to.mean, weight))
## stopifnot(all.equal(to.loss, best.loss1[to]))
## n <- nrow(sample.counts)
## seg2.ends <- 2:n
## best.seg2.start <- rep(NA, n)
## best.seg2.loss <- rep(NA, n)
## for(seg2.end in seg2.ends){
##   seg2.starts <- 2:seg2.end
##   seg1.ends <- seg2.starts-1L
##   seg1.loss <- best.loss1[seg1.ends]
##   seg1.mean <- best.mean1[seg1.ends]
##   rev.counts <- sample.counts[rev(seg2.starts), ]
##   rev.cumsum2 <- with(rev.counts, cumsum(coverage * weight))
##   rev.mean2 <- rev.cumsum2/cumsum(rev.counts$weight)
##   rev.loss2 <- ifelse(rev.mean2==0, 0, rev.cumsum2 * (1-log(rev.mean2)))
##   seg2.loss <- rev(rev.loss2)
##   seg2.mean <- rev(rev.mean2)
##   all.models <- 
##   data.table(is.feasible=seg2.mean > seg1.mean,
##              total.loss=seg1.loss + seg2.loss,
##              seg2.starts)
##   feasible.models <- all.models %>%
##     filter(is.feasible)
##   feasible.best <- feasible.models %>%
##     filter(total.loss == min(total.loss))
##   stopifnot(nrow(feasible.best) <= 1)
##   if(nrow(feasible.best) == 1){
##     best.seg2.start[seg2.end] <- feasible.best$seg2.starts
##     best.seg2.loss[seg2.end] <- feasible.best$total.loss
##   }
## }
n <- 3500
l <- 400
seg3.starts <- as.integer(seq(4, n, l=l))
loss.list <- list()
mean3.mat <- matrix(NA, l, 3)
for(model.i in seq_along(seg3.starts)){
  seg3.start <- seg3.starts[[model.i]]
  seg3.chromStart <- sample.counts$chromStart[seg3.start]
  ## index of last point on segment 2 is in  end.mat[,3] and index of
  ## last point on segment 1 is in end.mat[,2].

  ## end.mat[n, 3] gives last point on segment 2 for the best model up
  ## to and including data point n.
  seg2.end <- seg3.start-1L
  seg1.end <- end.mat[seg2.end, 2]
  seg.ends <- as.integer(c(seg1.end, seg2.end, n))
  seg.starts <- as.integer(c(1, seg1.end+1, seg3.start))
  
  for(seg.i in seq_along(seg.starts)){
    seg.start <- seg.starts[[seg.i]]
    seg.end <- seg.ends[[seg.i]]
    seg.data <- sample.counts[seg.start:seg.end, ]
    seg.mean <- with(seg.data, sum(coverage * weight)/sum(weight))
    seg.loss <- with(seg.data, PoissonLoss(coverage, seg.mean, weight))
    mean3.mat[model.i, seg.i] <- seg.mean
    loss.list[[paste(model.i, seg.i)]] <-
      data.table(model.i, seg.i,
                 seg3.chromStart, seg3.start,
                 seg.start, seg.end,
                 seg.chromStart=sample.counts$chromStart[seg.start],
                 seg.chromEnd=sample.counts$chromEnd[seg.end],
                 seg.mean, seg.loss)
  }
}
loss.dt <- do.call(rbind, loss.list)
model.dt <- loss.dt %>%
  group_by(model.i, seg3.chromStart) %>%
  summarise(loss=sum(seg.loss))
model.dt$feasible <-
  ifelse(mean3.mat[,3] < mean3.mat[,2], "yes", "no")
feasible <- model.dt %>%
  filter(feasible=="yes")

show.models <-
  sort(c(100, 250,
         which.min(model.dt$loss),
         350))
show.loss.list <- split(loss.dt, loss.dt$model.i)
show.model.list <- split(model.dt, model.dt$model.i)
png.list <- list()
last.base <- max(model.dt$seg3.chromStart/1e3)
best.loss <- min(model.dt$loss)
t.dt <- data.table(last.base, best.loss, what="loss")
for(show.model.i in seq_along(show.models)){
  model.i <- show.models[[show.model.i]]
  show.model <- show.model.list[[model.i]]
  show.loss <- show.loss.list[[model.i]]
  selectedPlot <- 
  ggplot()+
  geom_step(aes(chromStart/1e3, coverage),
            data=data.table(sample.counts, what="profile"),
            color="grey50")+
  geom_segment(aes(seg.chromStart/1e3, seg.mean,
                   xend=seg.chromEnd/1e3, yend=seg.mean),
               color="green",
               data=data.frame(show.loss, what="profile"))+
    geom_vline(aes(xintercept=last.base), data=t.dt,
               color="grey")+
    geom_text(aes(last.base, best.loss, label="t "),
              data=t.dt, hjust=1, vjust=0, color="grey")+
  ## geom_line(aes(seg3.chromStart/1e3, loss),
  ##           data=data.table(model.dt, what="loss"))+
  geom_point(aes(seg3.chromStart/1e3, loss, size=feasible),
             data=data.table(model.dt, what="loss"),
             pch=1)+
  geom_point(aes(seg3.chromStart/1e3, loss, size=feasible),
             data=data.table(show.model, what="loss"),
             pch=1,
             color="green")+
  geom_vline(aes(xintercept=seg3.chromStart/1e3),
             data=show.model,
             linetype="dotted",
             color="green")+
  geom_text(aes(seg3.chromStart/1e3,
                max(sample.counts$coverage),
                label="t' "),
             hjust=1,
             vjust=1,
             data=data.table(show.model, what="profile"),
             color="green")+
  scale_size_manual(values=c(yes=2, no=0.5))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(what ~ ., scales="free")+
  scale_y_continuous("",
                     breaks=c(seq(0, 30, by=10),
                       (13:17)*-10000))+
  xlab(paste("position on chromosome t' (kb = kilo bases)"))

  png(png.name <- sprintf("figure-dp-third-%d.png", show.model.i),
      units="in", res=200, width=6, height=3)
  print(selectedPlot)
  dev.off()

  png.list[[png.name]] <- png.name
}

pngs <- do.call(c, png.list)

png.tex <- sprintf("
\\begin{frame}
\\frametitle{Computation of optimal loss $\\mathcal L_{s, t}$
 for $s=3$ segments up to data point $t$}
  \\includegraphics[width=\\textwidth]{%s}

$$
\\mathcal L_{3, t} =
\\min_{
  t' < t
}
\\underbrace{
  \\mathcal L_{2, t'}
}_{
  \\text{optimal loss in 2 segments up to $t'$}
}
+
\\underbrace{
  c_{(t', t]}
}_{
  \\text{optimal loss of 3rd segment $(t', t]$}
}
$$

\\end{frame}
", pngs)

cat(png.tex, file="figure-dp-third.tex")
