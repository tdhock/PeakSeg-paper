works_with_R("3.1.2", Segmentor3IsBack="1.8")

set.seed(1)
N <- 250 
x <- rpois(10*N, rep(c(8,1,5,3,16,33,2,12,7,1),each=N))
Kmax <- 40
res <- Segmentor(data=x, model=3, Kmax=Kmax)
Cr <- SelectModel(res, penalty='oracle', keep=TRUE)
plot(Cr$criterion, type="l")
best.k <- which.min(Cr$criterion)
points(best.k, Cr$criterion[best.k])

break.mat <- res@breaks
n <- break.mat[1, 1]
sizenr <- function(k) {
  sum(log(diff(c(1, break.mat[k, 1:k]))))
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
    pk = (J[(k + 1):Kmax] - J[k])/(pen[k] - pen[(k + 
                                                 1):Kmax])
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

Lik <- res@likelihood

if (penalty == "mBIC") 
  K <- which.min(crit <- Lik + 0.5 * sapply(1:Kmax, 
                                                         sizenr) + (1:Kmax - 0.5) * log(n))
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

## oracle penalty code.
Kseq <- 1:Kmax
pen <- Kseq * (1 + 4 * sqrt(1.1 + log(n/Kseq))) * 
  (1 + 4 * sqrt(1.1 + log(n/Kseq)))
from.saut <- saut(-Lik[Kseq], pen, Kseq, 
          n/log(n), biggest = FALSE)
crit <- Lik + from.saut[2] * pen
K <- from.saut[1]

## My computation is the same as Alice's.
stopifnot(all.equal(as.numeric(crit), as.numeric(Cr$criterion)))
