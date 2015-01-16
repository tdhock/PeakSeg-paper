works_with_R("3.1.2",
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
  sizenr <- function(k) {
    sum(log(diff(c(1, end.mat[k, 1:k]))))
    ## the number of base pairs is used in the penalty computation
    ## ... this will change with our weighted problem!
  }
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
  pen <- Kseq * (1 + 4 * sqrt(1.1 + log(n/Kseq))) * 
    (1 + 4 * sqrt(1.1 + log(n/Kseq)))
  from.saut <- saut(-Lik[Kseq], pen, Kseq, 
                    n/log(n), biggest = FALSE)
  list(crit=Lik + from.saut[2] * pen,
       segments=from.saut[1])
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
    sum(log(diff(c(1, end.mat[k, 1:k]))))
    ## the number of base pairs is used in the penalty computation
    ## ... this will change with our weighted problem!
  }
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
  for(segments in 1:Kmax){
    seg.lik <- rep(NA, segments)
    ends <- end.mat[segments, 1:segments]
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
    }
    lik[segments] <- sum(seg.lik)
  }
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

## TODO: verify mBIC for weighted data.

if(FALSE){ # code from SelectModel
  if (penalty == "mBIC")
    entropy.term <- sapply(1:Kmax, sizenr)
  K <- which.min(crit <- Lik + 0.5 * entropy.term + (1:Kmax - 0.5) * log(n))
  if (penalty == "BIC") 
    K <- which.min(crit <- Lik + (1:K) * log(n))
  if (penalty == "AIC") 
    K <- which.min(crit <- Lik + (1:K) * 2)
  if (penalty == "oracle") {
    Kseq = 1:Kmax
    pen = Kseq * (1 + 4 * sqrt(1.1 + log(n/Kseq))) * 
      (1 + 4 * sqrt(1.1 + log(n/Kseq)))
    if (greatjump) {
      K = saut(-Lik[Kseq], pen, Kseq)
      crit <- Lik[Kseq] + K[2] * pen
      K <- K[1]
    }
    else {
      K = saut(-Lik[Kseq], pen, Kseq, 
        seuil, biggest = FALSE)
      crit <- Lik[Kseq] + K[2] * pen
      K <- K[1]
    }
  }
}

## TODO: For each (penalty, set, chunk, sample) compute the optimal
## number of peaks.
