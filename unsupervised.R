works_with_R("3.1.2",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             directlabels="2014.6.13",
             dplyr="0.4.0",
             PeakSeg="2014.12.2",
             Segmentor3IsBack="1.8")

set.seed(1)
N <- 250 
x <- rpois(10*N, rep(c(8,1,5,3,16,33,2,12,7,1),each=N))
Kmax <- 40
res <- Segmentor(data=x, model=1, Kmax=Kmax)
Cr <- SelectModel(res, penalty='oracle', keep=TRUE)
Cr.mBIC <- SelectModel(res, penalty='mBIC', keep=TRUE)
plot(Cr$criterion, type="l")
best.k <- which.min(Cr$criterion)
points(best.k, Cr$criterion[best.k])
mu <- mean(x)
my.lik.1seg <- -sum(dpois(x, mu, log=TRUE))
stopifnot(all.equal(res@likelihood[1], my.lik.1seg))
## Segmentor likelihood computation is the same as dpois!
  
alice.oracle <- function
### An adapted version of Alice's Segmentor3IsBack::SelectModel code
### for penalty="oracle" -- this function can be used with any model
### (not just Segmentor S4 classes).
(n,
### number of base pairs to segment.
 Lik
### Numeric vector of Kmax -sum(dpois(x, mu, log=TRUE)) values.
 ){
  Kmax <- length(Lik)
  saut <- function(Lv, pen, Kseq, seuil = sqrt(n)/log(n), biggest = TRUE) {
    J = -Lv
    Kmax = length(J)
    k = 1
    kv = c()
    dv = c()
    pv = c()
    dmax = 1
    while (k < Kmax) {
      pk = (J[(k + 1):Kmax] - J[k])/(pen[k] - pen[(k + 1):Kmax])
      pm = max(pk)
      dm = which.max(pk)
      dv = c(dv, dm)
      kv = c(kv, k)
      pv = c(pv, pm)
      if (dm > dmax) {
        dmax = dm
        kmax = k
        pmax = pm
      }
      k = k + dm
    }
    if (biggest) {
      pv = c(pv, 0)
      kv = c(kv, Kmax)
      dv = diff(kv)
      dmax = max(dv)
      rt = max(dv)
      rt = which(dv == rt)
      pmax = pv[rt[length(rt)]]
      alpha = 2 * pmax
      km = kv[alpha >= pv]
      Kh = Kseq[km[1]]
      return(c(Kh, alpha))
    }
    else {
      paux <- pv[which(kv <= seuil)]
      alpha <- 2 * min(paux)
      km = kv[alpha >= pv]
      Kh = Kseq[km[1]]
      return(c(Kh, alpha))
    }
  }

  ## oracle penalty code.
  Kseq <- 1:Kmax
  in.square <- 1 + 4 * sqrt(1.1 + log(n/Kseq))
  pen <- Kseq * in.square * in.square
  from.saut <- saut(-Lik[Kseq], pen, Kseq, 
                    n/log(n), biggest = FALSE)
  list(crit=Lik + from.saut[2] * pen,
       segments=from.saut[1])
}

my.oracle <- function
### An adapted version of Alice's Segmentor3IsBack::SelectModel code
### for penalty="oracle" -- this function can be used with any model
### (not just Segmentor S4 classes).
(n,
### number of base pairs to segment.
 Lik,
### Numeric vector of Kmax -sum(dpois(x, mu, log=TRUE)) values.
 beta=10^seq(-2, 4, l=200)
### beta from equation (6) in arXiv:1301.2534.
 ){
  Kmax <- length(Lik)
  Kseq <- 1:Kmax
  in.square <- 1 + 4 * sqrt(1.1 + log(n/Kseq))
  pen <- Kseq * in.square * in.square
  sapply(beta, function(b){
    which.min(Lik + b * pen)
  })
}

### From cghseg.
getmBIC <- function(K,lv,mu,CGHo){  
  M   = length(names(mu))
  N   = M*mu[[1]]$end[dim(mu[[1]])[1]]
  Ent = sum(unlist(lapply(mu,FUN = function(x){log(x$end-x$begin+1)})))

  if (CGHo["calling"]==FALSE){
    mBIC = ((N-K+1)/2)*(lv*(2/N)+1+log(2*pi))-0.5*Ent -(K-M)*log(N)+lgamma((N-K+1)/2)-((N-K+1)/2)*log(N)
  } else {
    P  = CGHo["nblevels"]
    Np = sapply(mu,FUN=function(x){
	tabulate(rep(x$levels,x$end-x$begin+1),P)
    })
    Np = apply(Np,1,sum)
    Ent  = sum(log(Np))
    mBIC = ((N-P+1)/2)*(lv*(2/N)+1+log(2*pi))-0.5*Ent-(K-M)*log(N)+lgamma((N-P+1)/2)-((N-P+1)/2)*log(N)
  }
  return(mBIC)
}

alice.mBIC <- function
### An adapted version of Alice's Segmentor3IsBack::SelectModel code
### for penalty="mBIC" -- this function can be used with any model
### (not just Segmentor S4 classes).
(end.mat,
### Kmax x Kmax lower diagonal matrix. Row K has the last base indices
### of the model with K segments.
 Lik
### Numeric vector of Kmax -sum(dpois(x, mu, log=TRUE)) values.
 ){
  stopifnot(length(Lik) == nrow(end.mat))
  stopifnot(length(Lik) == ncol(end.mat))
  n <- end.mat[1, 1]
  sizenr <- function(k) {
    ## This is copied directly from Alice's code and I suspect it is
    ## incorrect -- the 1 should be 0.
    bases.vec <- diff(c(1, end.mat[k, 1:k]))
    sum(log(bases.vec))
    ## the number of base pairs is used in the penalty computation
    ## ... this will change with our weighted problem!
  }
  ## mBIC penalty code.
  entropy.term <- sapply(1:Kmax, sizenr)
  crit <- Lik + 0.5 * entropy.term + (1:Kmax - 0.5) * log(n)
  K <- which.min(crit)
  list(crit=crit,
       segments=K)
}

incorrect.mBIC <- function(end.mat, Lik, weight){
  stopifnot(length(Lik) == nrow(end.mat))
  stopifnot(length(Lik) == ncol(end.mat))
  n <- sum(weight)
  sizenr <- function(k) {
    model.ends <- end.mat[k, 1:k]
    model.starts <- c(1, 1+model.ends[-length(model.ends)])
    bases.vec <-
      sapply(seq_along(model.starts), function(segment.i){
        from <- model.starts[[segment.i]]
        to <- model.ends[[segment.i]]
        sum(weight[from:to])
      })
    bases.vec[1] <- bases.vec[1] - 1 # THIS LINE IS THE ONLY DIFFERENCE!
    sum(log(bases.vec))
  }
  ## mBIC penalty code.
  entropy.term <- sapply(1:Kmax, sizenr)
  crit <- Lik + 0.5 * entropy.term + (1:Kmax - 0.5) * log(n)
  K <- which.min(crit)
  list(crit=crit,
       segments=K)
}

correct.mBIC <- function(end.mat, Lik, weight){
  stopifnot(length(Lik) == nrow(end.mat))
  stopifnot(length(Lik) == ncol(end.mat))
  n <- sum(weight)
  Kmax <- length(Lik)
  sizenr <- function(k) {
    model.ends <- end.mat[k, 1:k]
    if(all(!is.na(model.ends))){
      model.starts <- c(1, 1+model.ends[-length(model.ends)])
      bases.vec <-
        sapply(seq_along(model.starts), function(segment.i){
          from <- model.starts[[segment.i]]
          to <- model.ends[[segment.i]]
          sum(weight[from:to])
        })
      sum(log(bases.vec))
    }else{
      NA
    }
  }
  ## mBIC penalty code.
  entropy.term <- sapply(1:Kmax, sizenr)
  crit <- Lik + 0.5 * entropy.term + (1:Kmax - 0.5) * log(n)
  K <- which.min(crit)
  list(crit=crit,
       segments=K)
}

## My copy of Alice's oracle criterion computation is the same as
## Alice's code in Segmentor3IsBack.
crit.info <- alice.oracle(length(x), res@likelihood)
stopifnot(all.equal(as.numeric(crit.info$crit), as.numeric(Cr$criterion)))

## My copy of Alice's mBIC criterion computation is the same as
## Alice's code in Segmentor3IsBack.
mBIC.info <- alice.mBIC(res@breaks, res@likelihood)
stopifnot(all.equal(as.numeric(mBIC.info$crit), as.numeric(Cr.mBIC$criterion)))

## Guillem's Unconstrained DP returns the same breaks as Alice's
## unconstrained Pruned DP.
w <- rep(1, length(x))
un.fit <- uPoissonSeg_(x, w, 40)
un.ends <- getPath(un.fit)
stopifnot(all.equal(as.integer(un.ends[40, ]),
                    as.integer(res@breaks[40, ])))
PoissonLik <- function(count, bases, end.mat){
  Kmax <- nrow(end.mat)
  lik <- rep(NA, Kmax)
  loss <- rep(NA, Kmax)
  for(segments in 1:Kmax){
    seg.lik <- rep(NA, segments)
    seg.loss <- rep(NA, segments)
    ends <- end.mat[segments, 1:segments]
    if(all(!is.na(ends))){
      breaks <- ends[-length(ends)]
      starts <- c(1, breaks+1)
      for(segment.i in 1:segments){
        first <- starts[segment.i]
        last <- ends[segment.i]
        seg.data <- count[first:last]
        seg.bases <- bases[first:last]
        seg.mean <- sum(seg.data * seg.bases)/sum(seg.bases)
        loglik.vec <- dpois(seg.data, seg.mean, log=TRUE)
        seg.lik[segment.i] <- -sum(loglik.vec * seg.bases)
        seg.loss[segment.i] <- PoissonLoss(seg.data, seg.mean, seg.bases)
      }
      lik[segments] <- sum(seg.lik)
      loss[segments] <- sum(seg.loss)
    }
  }
  attr(lik, "loss") <- loss
  lik
}

## My likelihood is the same as Alice's Segmentor computation.
un.lik <- PoissonLik(x, w, un.ends)
stopifnot(all.equal(as.numeric(un.lik),
                    as.numeric(res@likelihood)))
my.info <- alice.oracle(length(x), un.lik)
stopifnot(all.equal(as.numeric(my.info$crit), as.numeric(Cr$criterion)))

## Verify the same computation for weighted data.
x.rle <- rle(x)
comp.fit <- uPoissonSeg_(x.rle$values, x.rle$lengths, 40)
comp.ends <- getPath(comp.fit)
comp.lik <- PoissonLik(x.rle$values, x.rle$lengths, comp.ends)
stopifnot(all.equal(as.numeric(comp.lik),
                    as.numeric(res@likelihood)))
comp.info <- alice.oracle(length(x), comp.lik)
stopifnot(all.equal(as.numeric(comp.info$crit), as.numeric(Cr$criterion)))

## Verify mBIC for weighted data.
comp.incorrect <- incorrect.mBIC(comp.ends, comp.lik, x.rle$lengths)
comp.correct <- correct.mBIC(comp.ends, comp.lik, x.rle$lengths)
stopifnot(all.equal(as.numeric(comp.incorrect$crit),
                    as.numeric(Cr.mBIC$criterion)))
all.equal(as.numeric(comp.correct$crit),
          as.numeric(Cr.mBIC$criterion)) # not equal!

## Also compute AIC/BIC.
crit.BIC <- comp.lik + (1:Kmax) * log(length(x))
crit.AIC <- comp.lik + (1:Kmax) * 2

cr <- function(crit, penalty){
  data.frame(crit=as.numeric(crit),
             segments=seq_along(crit),
             penalty)
}
crit.df <- 
rbind(cr(comp.info$crit, "oracle"),
      cr(crit.BIC, "BIC"),
      cr(crit.AIC, "AIC"),
      cr(comp.lik, "none"),
      cr(comp.correct$crit, "mBIC"))
best.df <- 
crit.df %>%
  group_by(penalty) %>%
  filter(seq_along(crit) == which.min(crit))
with.legend <- 
ggplot()+
  geom_line(aes(segments, crit, color=penalty),
            data=crit.df)+
  geom_point(aes(segments, crit, color=penalty),
             data=best.df, pch=1)+
  scale_y_log10()
direct.label(with.legend)

unsupervised <- list()
oracle.segments <- list()
loss.segments <- list()
lik.segments <- list()
model.files <- Sys.glob("data/*/*/dp.model.RData")
for(model.file.i in seq_along(model.files)){
  model.file <- model.files[[model.file.i]]
  count.file <- sub("dp.model", "counts", model.file)
  chunk.name <- sub("data/", "", dirname(model.file))
  cat(sprintf("%4d / %4d %s\n", model.file.i, length(model.files), chunk.name))
  load(model.file)
  load(count.file)
  sample.list <- split(counts, counts$sample.id)
  segments.list <- list()
  oseg.list <- list()
  loss.list <- list()
  lik.list <- list()
  for(sample.id in names(sample.list)){
    sample.segments <- dp.model[[sample.id]]$segments
    sample.counts <- sample.list[[sample.id]]
    seg.list <- split(sample.segments, sample.segments$segments)
    end.mat <- matrix(NA, 19, 19)
    rownames(end.mat) <- 1:19
    for(seg.str in names(seg.list)){
      seg.num <- as.numeric(seg.str)
      segs <- seg.list[[seg.str]]
      end.mat[seg.str, 1:nrow(segs)] <- segs$last
    }
    weight <- with(sample.counts, chromEnd-chromStart)
    sample.lik <- 
      PoissonLik(sample.counts$coverage, weight, end.mat)
    sample.loss <- attr(sample.lik, "loss")
    segSeq <- seq(1, 19, by=2)
    lik.list[[sample.id]] <- sapply(10^seq(2, 6, l=200), function(pen){
      segSeq[which.min(sample.lik[segSeq] + pen * segSeq)]
    })
    loss.list[[sample.id]] <- sapply(10^seq(2, 6, l=200), function(pen){
      segSeq[which.min(sample.loss[segSeq] + pen * segSeq)]
    })
    qplot(sample.lik, sample.loss)+
      coord_equal()
    force.na <- is.na(sample.lik)
    sample.lik[seq(2, 18, by=2)] <- sample.lik[seq(1, 17, by=2)]
    while(any(isNA <- is.na(sample.lik))){
      na.i <- which(isNA)[1]
      sample.lik[na.i] <- sample.lik[na.i-1]
    }
    ##stopifnot(diff(sample.lik) <= 0)
    sample.lik[is.na(sample.lik)] <- sample.lik[1]
    sample.mBIC <- correct.mBIC(end.mat, sample.lik, weight)
    mbic <- sample.mBIC$crit[segSeq]
    mbic[is.na(mbic)] <- mbic[1]
    sample.oracle <- alice.oracle(sum(weight), sample.lik)
    oseg.list[[sample.id]] <- my.oracle(sum(weight), sample.lik)
    write.na <- force.na[segSeq]
    penalty.mat <-  #peaks x penalties
      cbind(oracle=sample.oracle$crit[segSeq],
            mBIC=mbic,
            BIC=sample.lik[segSeq] + segSeq * log(sum(weight)),
            AIC=sample.lik[segSeq] + segSeq * 2,
            none=sample.lik[segSeq])
    rownames(penalty.mat) <- segSeq
    penalty.mat[write.na, ] <- NA
    selected.rows <- apply(penalty.mat, 2, which.min)
    selected.segs <- segSeq[selected.rows]
    names(selected.segs) <- names(selected.rows)
    segments.list[[sample.id]] <- selected.segs
  }
  unsupervised[[chunk.name]] <- do.call(rbind, segments.list)
  oracle.segments[[chunk.name]] <- do.call(rbind, oseg.list)
  loss.segments[[chunk.name]] <- do.call(rbind, loss.list)
  lik.segments[[chunk.name]] <- do.call(rbind, lik.list)
}
test <- do.call(rbind, unsupervised)
stopifnot(!is.na(test))

save(unsupervised, oracle.segments, loss.segments, lik.segments,
     file="unsupervised.RData")
