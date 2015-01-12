works_with_R("3.1.1",
             "tdhock/PeakSegDP@5bcee97f494dcbc01a69e0fe178863564e9985bc",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             tikzDevice="0.7.0",
             quadmod="2013.8.23",
             dplyr="0.3.0.2",
             "tdhock/animint@1246883c2cd75a773f598a12fca4a2af4c9160e9")

load("PeakSeg4samples.RData")

error.list <- PeakSeg4samples$regions
errors <- NULL
for(peaks in 0:4){
  rname <- paste0("PeakSeg", peaks)
  e <- error.list[[rname]] %>%
    group_by(sample.id) %>%
    summarise(errors=sum(fp+fn))
  errors <- rbind(errors, data.frame(peaks, e))
}

loss.list <- split(PeakSeg4samples$loss,
                   PeakSeg4samples$loss$sample.id,
                   drop=TRUE)
err.list <- split(errors, errors$sample.id, drop=TRUE)
signal.list <- split(PeakSeg4samples$signal,
                     PeakSeg4samples$signal$sample.id, drop=TRUE)
exact.dfs <- NULL
intervals <- NULL
for(sample.id in names(loss.list)){
  one <- subset(loss.list[[sample.id]], peaks <= 4)
  exact.df <- with(one, exactModelSelection(error, peaks))
  err <- err.list[[sample.id]]
  rownames(err) <- err$peaks
  signal <- signal.list[[sample.id]]
  exact.df$errors <- err[as.character(exact.df$model.complexity), "errors"]
  indices <- with(exact.df, {
    largestContinuousMinimum(errors, max.log.lambda-min.log.lambda)
  })
  meta <- data.frame(sample.id, log.max.count=log(max(signal$coverage)))
  exact.dfs <- rbind(exact.dfs, {
    data.frame(meta, exact.df)
  })
  intervals <- rbind(intervals, {
    data.frame(meta,
               min.log.lambda=exact.df$min.log.lambda[indices$start],
               max.log.lambda=exact.df$max.log.lambda[indices$end])
  })
}

what.peaks <- "peaks $p_i^*(\\lambda)$"
what.error <- "error $E_i(\\lambda)$"
blank.df <- data.frame(x=10, what=what.peaks, sample.id="McGill0002")
exact.peaks <- data.frame(exact.dfs, what=what.peaks)
exact.error <- data.frame(exact.dfs, what=what.error)
funplot <- 
ggplot()+
  geom_segment(aes(min.log.lambda, model.complexity,
                   xend=max.log.lambda, yend=model.complexity),
               data=exact.peaks, size=2)+
  geom_segment(aes(min.log.lambda, errors,
                   xend=max.log.lambda, yend=errors),
               data=exact.error, size=2)+
  geom_blank(aes(x, 0), data=blank.df)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(what ~ sample.id, scales="free_y", space="free_y")+
  scale_x_continuous("$\\log \\lambda$", limits=c(7, 14),
                     breaks=seq(8, 12, by=2))+
  ylab("")
options(tikzDocumentDeclaration="\\documentclass{beamer}\\usepackage{amsmath,amssymb,amsthm}",
        tikzMetricsDictionary="slideMetrics")

zero.peaks <- subset(exact.peaks, errors==0)
zero.error <- subset(exact.error, errors==0)
funplot <- 
ggplot()+
  geom_segment(aes(min.log.lambda, model.complexity,
                   xend=max.log.lambda, yend=model.complexity),
               data=zero.peaks, size=3, color="green")+
  geom_segment(aes(min.log.lambda, errors,
                   xend=max.log.lambda, yend=errors),
               data=zero.error, size=3, color="green")+
  geom_segment(aes(min.log.lambda, model.complexity,
                   xend=max.log.lambda, yend=model.complexity),
               data=exact.peaks, size=2)+
  geom_segment(aes(min.log.lambda, errors,
                   xend=max.log.lambda, yend=errors),
               data=exact.error, size=2)+
  geom_blank(aes(x, 0), data=blank.df)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(what ~ sample.id, scales="free", space="free")+
  scale_x_continuous("$\\log \\lambda$", limits=c(7, 14),
                     breaks=seq(8, 12, by=2))+
  ylab("")

what.intervals <- data.frame(intervals, what=what.error)
what.intervals$log.lambda <- what.intervals$max.log.lambda
left.intervals <- subset(what.intervals, is.finite(min.log.lambda))
left.intervals$log.lambda <- left.intervals$min.log.lambda
target.text <-
  rbind(data.frame(label="$\\overline L_i$", what.intervals, hjust=1),
        data.frame(label="$\\underline L_i$", left.intervals, hjust=0))
tsize <- 3
fun2 <- funplot+
  geom_text(aes(log.lambda, 1.5, label=label),
            data=target.text)+
  geom_point(aes(min.log.lambda, 0), shape=21, fill="white",
             data=left.intervals,
             size=tsize)+
  geom_point(aes(max.log.lambda, 0), shape=21, fill="black",
             data=what.intervals, size=tsize)

max.margin <- function
### Support vector interval regression for separable data. The idea is
### that we first normalize the feature matrix, giving normalized
### features x_i in R^p. Then we use a linear function f(x_i) = w'x_i
### + b to predict a log(lambda) that falls between all the log limits
### L_i^left and L_i^right and maximizes the margin. So the
### optimization problem is: max_{M,f} subject to, for all finite
### limits, L_i^right - f(x_i) >= M and f(x_i) - L_i^left >= M. Since
### we assume f is linear the problem becomes min_{w,b,M} -M subject
### to -M - w'x - b >= -L_i^right and -M + w'x + b >= L_i^left. We
### call M margin, w weights, b intercept.
(features,
### Matrix n x p of inputs: n signals, each with p features. We will
### scale these internally.
 limits,
### Matrix n x 2 of output lambda. Each row corresponds to the lower
### and upper bound of an interval on the log(lambda) which is optimal
### with respect to annotation error. Lower bound can be -Inf and
### upper bound can be Inf, which correspond to zero asymptotic
### cost. 
 verbose=0,
 ...
### ignored.
 ){
  ## reality checks.
  stopifnot(nrow(features)==nrow(limits))
  if(ncol(limits)!=2){
    cat("str(limits)=\n")
    str(limits)
    stop("limits should be a 2-column matrix")
  }
  stopifnot(is.matrix(features))
  
  ## check if there are any flat error curves, which have no limits.
  has.limits <- apply(is.finite(limits),1,any)
  ## we train the model on this subset.
  some.limits <- limits[has.limits,]
  some.features <- features[has.limits,,drop=FALSE]

  scaled <- scale(some.features)
  mu <- attr(scaled,"scaled:center")
  sigma <- attr(scaled,"scaled:scale")

  n <- nrow(scaled)
  p <- ncol(scaled)
  vars <- make.ids(margin=1,intercept=1,weights=p)
  constraints <- list(vars$margin*1 >= 0)
  for(i in 1:n){
    if(verbose >= 1)cat(sprintf("example constraints %5d / %5d",i,n))

    left <- some.limits[i,1]
    if(is.finite(left)){
      ivars <- with(vars,{
        intercept * 1 + sum(weights)*scaled[i,] + margin*-1
      })
      constraints <- c(constraints,list(ivars >= left))
    }

    right <- some.limits[i,2]
    if(is.finite(right)){
      ivars <- with(vars,{
        intercept * -1 + sum(weights)*scaled[i,]*-1 +margin*-1
      })
      constraints <- c(constraints,list(ivars >=  - right))
    }

    if(verbose >= 1)cat("\n")

  }
  const.info <- standard.form.constraints(constraints,vars)
  n.vars <- length(unlist(vars))
  Dvec <- rep(1e-10,n.vars)
  D <- diag(Dvec)
  d <- rep(0,n.vars)
  d[vars$margin] <- 1
  if(verbose >= 1)cat(sprintf("solving for %d variables and %d constraints... ",
              n.vars,length(constraints)))
  sol <- solve.QP(D,d,const.info$A,const.info$b0)
  if(verbose >= 1)cat("solved!\n")
  sol$mu <- mu
  sol$sigma <- sigma
  sol$scaled <- scaled
  sol$log.limits <- some.limits
  sol$features <- some.features
  sol$weights <- sol$solution[vars$weights]
  sol$intercept <- sol$solution[vars$intercept]
  sol$margin <- sol$solution[vars$margin]
  ## this function will be applied to new data before applying the
  ## model.
  sol$normalize <- function(X){
    mu.mat <- matrix(mu,nrow(X),ncol(X),byrow=TRUE)
    s.mat <- matrix(sigma,nrow(X),ncol(X),byrow=TRUE)
    (X-mu.mat)/s.mat
  }
  sol$f <- function(x){
    sum(x*sol$weights)+sol$intercept
  }
  sol$predict <- function(X){
    stopifnot(is.matrix(X))
    X.norm <- sol$normalize(X)
    weights.mat <- matrix(sol$weights,nrow(X),ncol(X),byrow=TRUE)
    L.hat <- rowSums(X.norm * weights.mat) + sol$intercept
    L.hat
  }
  sol$L.pred <- apply(scaled,1,sol$f)
  sol$lambda.pred <- sol$predict(features)
  sol
### List of solver results. For a feature matrix X with p columns, you
### can use list$predict(X) to get model estimates of log(lambda).
}

fit <- with(intervals, {
  max.margin(cbind(log.max.count), cbind(min.log.lambda, max.log.lambda))
})

zero.error <- subset(exact.dfs, errors==0)

intervals$predicted <- fit$L.pred
intervals$mid.log.lambda <- with(intervals, (max.log.lambda+min.log.lambda)/2)
prediction.in.mid <- with(intervals, abs(mid.log.lambda-predicted) < 1e-6)
## max margin line should be in the middle of one interval.
stopifnot(any(prediction.in.mid)) 
min.x <- min(intervals$log.max.count)
max.x <- max(intervals$log.max.count)
small.slope <- (9-11)/(min.x-max.x)
intervals$bad.pred <- small.slope*(intervals$log.max.count - max.x) + 11
intervals$right.margin <- with(intervals, max.log.lambda-predicted)
intervals$left.margin <- with(intervals, predicted-min.log.lambda)
min.log.count <- min(intervals$log.max.count)
max.log.count <- max(intervals$log.max.count)
seg.df <-
  data.frame(min.log.lambda=9, min.log.count,
             max.log.lambda=c(11, 11.5), max.log.count) %>%
  mutate(slope=(min.log.lambda-max.log.lambda)/(min.log.count-max.log.count),
         intercept=min.log.lambda-slope*min.log.count)
reg.df <- rbind(seg.df[, c("slope", "intercept")],
                data.frame(slope=0, intercept=8.5))
count.grid <- c(3, 7)
penalty.grid <- NULL
for(reg.i in 1:nrow(reg.df)){
  r <- reg.df[reg.i, ]
  penalty.grid <- rbind(penalty.grid, {
    data.frame(count.grid, log.lambda=r$slope * count.grid + r$intercept, reg.i)
  })
}

extreme <- subset(intervals, log.max.count %in% range(log.max.count))
extreme$max.log.lambda[1] <- 11
extreme$min.log.lambda[2] <- extreme$predicted[2]

## Plot the model that was selected by the large margin model.
rownames(intervals) <- intervals$sample.id
selected <- exact.dfs %>%
  mutate(predicted=intervals[as.character(sample.id), "predicted"]) %>%
  filter(min.log.lambda < predicted & predicted < max.log.lambda)

profile.list <- list()
for(sample.i in 1:nrow(selected)){
  r <- selected[sample.i, ]
  model.name <- paste0("PeakSeg", r$model.complexity)
  print(model.name)
  for(data.type in c("regions", "peaks")){
    profile.list[[data.type]] <- rbind(profile.list[[data.type]], {
      subset(PeakSeg4samples[[data.type]][[model.name]],
             sample.id==as.character(r$sample.id))
    })
  }
}

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

sample.ids <- unique(profile.list$regions$sample.id)
compare.peak.list <-
  list(PeakSeg=profile.list$peaks[, c("sample.id", "chromStart", "chromEnd")])
sample.max.df <- PeakSeg4samples$signal %>%
  group_by(sample.id) %>%
  summarise(count=max(coverage))
sample.max <- sample.max.df$count
names(sample.max) <- as.character(sample.max.df$sample.id)

selectedPlot <- 
ggplot()+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    fill=annotation),
                data=profile.list$regions,
                color="grey",
                alpha=1/2)+
  geom_step(aes(chromStart/1e3, coverage),
            data=PeakSeg4samples$signal, color="grey50")+
  ## geom_point(aes(chromStart/1e3, 0),
  ##            data=profile.list$peaks,
  ##            pch=1, size=2, color="deepskyblue")+
  ## geom_text(aes(118120, y.mid, label=label),
  ##           data=compare.labels, hjust=1, size=3)+
  ## geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
  ##               ymin=y.min, ymax=y.max,
  ##               linetype=status),
  ##           data=compare.regions,
  ##           fill=NA, color="black", size=0.5)+
  ## scale_linetype_manual("error type",
  ##                       values=c(correct=0,
  ##                         "false negative"=3,
  ##                         "false positive"=1))+
  ## geom_segment(aes(chromStart/1e3, y.mid,
  ##                  xend=chromEnd/1e3, yend=y.mid),
  ##              data=compare.peaks, size=1.5, color="deepskyblue")+
  geom_segment(aes(chromStart/1e3, 0,
                   xend=chromEnd/1e3, yend=0),
               data=profile.list$peaks, size=1.5, color="#6A3D9A")+
  coord_cartesian(xlim=c(118090, 118125))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ ., scales="free")+
  scale_y_continuous("aligned read coverage",
                     labels=function(x){
                       sprintf("%.1f", x)
                     },
                     breaks=function(limits){
                       limits[2]
                     })+
  xlab("position on chr11 (kilo base pairs)")+
  scale_fill_manual("annotation", values=ann.colors,
                    breaks=names(ann.colors))

png("figure-PeakSeg-4samples-intervals-selected.png",
    units="in", res=200, width=7, height=4)
print(selectedPlot)
dev.off()
##system("display figure-PeakSeg-4samples-intervals-selected.png")
 
p <- 
ggplot()+
  geom_segment(aes(min.log.lambda, log.max.count,
                   xend=max.log.lambda, yend=log.max.count),
               data=intervals, size=1.5)+
  geom_segment(aes(min.log.lambda, log.max.count,
                   xend=max.log.lambda, yend=log.max.count),
               data=extreme, color="red")+
  ## geom_line(aes(predicted, log.max.count), data=intervals, color="blue")+
  ## geom_line(aes(bad.pred, log.max.count), data=intervals, color="blue")+
  ## geom_vline(aes(xintercept=8.3), color="blue")+
  geom_line(aes(log.lambda, count.grid, group=reg.i),
            data=penalty.grid, color="blue")+
  coord_cartesian(ylim=c(4.2, 5.7))+
  ## geom_segment(aes(min.log.lambda, min.log.count,
  ##                  xend=max.log.lambda, yend=max.log.count),
  ##              data=seg.df, color="blue")+
  geom_text(aes(8.6, 5.3, label="1 error\nconstant"))+
  geom_text(aes(10.25, 5.3, label="0 errors\nsmall margin"))+
  geom_text(aes(11.5, 5.3, label="0 errors\nlarge margin"))+
  geom_point(aes(min.log.lambda, log.max.count), data=zero.error, size=5,pch="|")+
  geom_point(aes(max.log.lambda, log.max.count), data=zero.error, size=5, pch="|")+
  geom_text(aes((min.log.lambda + max.log.lambda)/2, log.max.count,
                label=sprintf("%d peak%s", model.complexity,
                  ifelse(model.complexity==1, "", "s")),
                hjust=ifelse(is.finite(min.log.lambda), 0.5, 0)),
            data=zero.error, vjust=-0.5, size=3)+
  geom_text(aes(max.log.lambda, log.max.count, label=sample.id),
            data=intervals, hjust=0)+
  ggtitle("max margin interval regression, margin in red")

pdf("figure-PeakSeg-4samples-intervals.pdf", w=10)
print(p)
dev.off()
 
## Plot the max margin regression line.
text.df <- zero.error %>%
  ##filter(model.complexity != 3)
  mutate(label=sprintf("%d peak%s", model.complexity,
                   ifelse(model.complexity==1, "", "s")),
         label.penalty=(min.log.lambda + max.log.lambda)/2)
rownames(text.df) <- with(text.df, paste(sample.id, model.complexity))
text.df["McGill0004 2", "label.penalty"] <- 8.75
text.df["McGill0091 1", "label.penalty"] <- 10
text.df["McGill0002 2", "label.penalty"] <- 12
tsize <- 2.5
p <- 
ggplot()+
  theme_grey()+
  geom_segment(aes(min.log.lambda, log.max.count,
                   xend=max.log.lambda, yend=log.max.count),
               data=intervals, size=1.5)+
  geom_point(aes(min.log.lambda, log.max.count),
             data=subset(zero.error, is.finite(min.log.lambda)),
             size=tsize, pch=1)+
  geom_point(aes(max.log.lambda, log.max.count),
             data=zero.error,
             size=tsize, pch=1)+
  geom_point(aes(min.log.lambda, log.max.count),
             data=subset(intervals, is.finite(min.log.lambda)),
             size=tsize, shape=21, fill="white")+
  geom_point(aes(max.log.lambda, log.max.count),
             data=intervals,
             size=tsize, shape=21, fill="black")+
  geom_text(aes(label.penalty,
                log.max.count +
                ifelse(log.max.count==max(log.max.count), -1, 1)*0.03,
                label=label,
                vjust=ifelse(is.finite(min.log.lambda), 0.5, -0.5),
                hjust=ifelse(log.max.count==max(log.max.count), 1, 0)),
            data=text.df, size=3)+
  geom_text(aes(max.log.lambda, log.max.count, label=sample.id,
                hjust=ifelse(log.max.count==min(log.max.count), 0,
                  ifelse(log.max.count==max(log.max.count), 1, 0.5))),
            data=intervals, vjust=-0.5, size=3)+
  coord_cartesian(xlim=c(7.7, 14))+
  scale_x_continuous("penalty $\\log\\lambda_i$", 
                     minor_breaks=NULL)+
  scale_y_continuous("feature $x_i = \\log\\max\\mathbf y_i$",
                     minor_breaks=NULL)+
  coord_flip(ylim=c(3.6, 6.3), xlim=c(7, 13.2))

modelsPlot <-
  p+
  geom_segment(aes(min.log.lambda, log.max.count,
                   xend=predicted, yend=log.max.count),
               data=intervals["McGill0002",], color="red")+
  ## geom_vline(aes(xintercept=8.5), color="blue")+
  ## geom_line(aes(predicted, log.max.count), data=intervals, color="blue")+
  ## geom_line(aes(bad.pred, log.max.count), data=intervals, color="blue")+
  geom_line(aes(log.lambda, count.grid, group=reg.i),
            data=penalty.grid, color="blue")+
  geom_text(aes(12.5, 5.8, label="0 errors\nlarge margin"),
            hjust=0, vjust=0, color="blue", size=3)+
  geom_text(aes(11, 5.8, label="0 errors\nsmall margin"),
            hjust=0, vjust=1, color="blue", size=3)+
  geom_text(aes(9, 5.8, label="1 error\nconstant"),
            hjust=0, color="blue", size=3)

options(tikzDocumentDeclaration="\\documentclass{article}\\usepackage{amsmath,amssymb,amsthm}",
        tikzMetricsDictionary="tikzMetrics")
tikz("figure-interval-regression.tex", h=3, w=4.5)
print(modelsPlot)
dev.off()
