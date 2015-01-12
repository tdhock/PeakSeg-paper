### List of functions, each a phi loss.
phi.list <- list(square=function(x){
  ifelse(x<1,(x-1)^2,0)
},huber=function(x){
  ifelse(x<0,1-x,ifelse(x<2,-x+x^2/4+1,0))
},log=function(x){
  log(1+exp(-x))
})

### if we have already calculated the linear predictor using
### fit$predict, this function can be useful.
calc.loss.from.lp.list <- lapply(phi.list,function(phi){
  force(phi)
  function(linear.predictor,lim){
    left.term <- phi(linear.predictor-lim[,1])
    right.term <- phi(lim[,2]-linear.predictor)
    sum(left.term+right.term)/nrow(lim)
  }
})

### List of interval regression loss functions: x, feat, lim =>
### numeric.
calc.loss.list <- lapply(calc.loss.from.lp.list,function(calc.from.lp){
  force(calc.from.lp)
  function(x,feat,lim){
    linear.predictor <- as.vector( cbind(1,feat) %*% x )
    calc.from.lp(linear.predictor,lim)
  }
})

### List of functions, each a derivative of a phi loss.
deriv.list <- list(huber=function(x){
  ifelse(x<0,-1,ifelse(x<2,x/2-1,0))
},square=function(x){
  ifelse(x<1,2*(x-1),0)
},log=function(x){
  -1/(exp(x)+1)
})

### List of calc.grad functions: x, features, limits -> gradient.
calc.grad.list <- lapply(deriv.list,function(phi.deriv){
  force(phi.deriv)
  function(x,feat,lim){
    linear.predictor <- as.vector( cbind(1,feat) %*% x )
    left.term <- phi.deriv(linear.predictor-lim[,1])
    right.term <- phi.deriv(lim[,2]-linear.predictor)
    full.grad <- cbind(1,feat) * (left.term-right.term)
    colSums(full.grad)/nrow(full.grad)
  }
})

### List of regression functions: features, limits -> list.
regression.funs <- lapply(names(calc.grad.list),function(loss.name){
  calc.grad <- calc.grad.list[[loss.name]]
  calc.loss <- calc.loss.list[[loss.name]]
  function(...){
    smooth.interval.regression(calc.grad=calc.grad,calc.loss=calc.loss,...)
  }
})
names(regression.funs) <- names(calc.grad.list)

positive.part <- function(x){
  ifelse(x<0, 0, x)
}

find.solution <- function
### Interval regression using a smooth loss to relax the annotation
### loss, and an accelerated proximal gradient descent solver. The
### idea is that we first normalize the feature matrix, giving
### normalized features x_i in R^p. Then we use a linear function
### f(x_i) = w'x_i + b to predict a log(lambda) that falls between the
### log limits L_i^left and L_i^right. So the optimization problem is:
### min_{w,b} gamma*||w||_1 + 1/n * sum_i phi[f(x_i)-L_i^left] +
### phi[L_i^right-f(x_i)], where phi(L) is a convex relaxation of the
### annotation loss, and should be specified in the calc.grad and
### calc.loss arguments. The optimization stops when we find an
### optimization variable for which each dimension is within a
### threshold of the subgradient optimality condition.
(features,
### feature matrix n x p which we assume has already been scaled.
 limits,
### limit matrix n x 2.
 gamma=0,
### regularization >= 0, by default 0.
 starting.iterate=rep(0,l=ncol(features)+1),
### Where to start the optimization, by default at the origin.
 threshold=1e-3,
### When the stopping criterion gets below this threshold, the
### solution is optimal.
 verbose=1,
### print optimization status?
 max.iterations=1e4,
### exit with an error if we haven't converged after this many
### iterations.
 calc.grad=stop("must specify calc.grad(x,features,limits)"),
### Function x,features,limits => gradient vector. This will be used
### for the optimization.
 calc.loss=stop("must specify calc.loss(x,features,limits)"),
### Function to display the calculate cost function, necessary for the
### backtracking line search.
 L0=ncol(features)+sqrt(ncol(features)),
### Lipshitz constant for step size.
 step="constant"
### constant or linesearch.
 ){
  stopifnot( length(starting.iterate) == (ncol(features)+1) )
  calc.penalty <- function(x){
    gamma * sum(abs(x[-1]))
  }
  calc.cost <- function(x){
    calc.loss(x,features,limits) + calc.penalty(x)
  }
  soft.threshold <- function(x,thresh){
    ifelse(abs(x) < thresh, 0, x-thresh*sign(x))
  }
  ## do not threshold the intercept.
  prox <- function(x,thresh){
    x[-1] <- soft.threshold(x[-1],thresh)
    x
  }
  ## p_L from the fista paper.
  pL <- function(x,L){
    grad <- calc.grad(x,features,limits)
    prox(x - grad/L, gamma/L)
  }
  ## Q_L from the fista paper, equation 2.5.
  QL <- function(x,y,L){
    fy <- calc.loss(y,features,limits)
    diff <- x-y
    inner.prod <- t(diff) %*% calc.grad(y,features,limits)
    quad.term <- L/2*t(diff) %*% diff
    gx <- calc.penalty(x)
    fy + inner.prod + quad.term + gx
  }

  iterate.count <- 1
  stopping.crit <- threshold
  last.iterate <- this.iterate <- y <- starting.iterate
  this.t <- 1
  last.L <- L0
  while(stopping.crit >= threshold){
    ## here we implement the FISTA method as described by in the Beck
    ## and Tebolle paper.
    last.iterate <- this.iterate
    this.iterate <- if(step=="linesearch" && 0==(iterate.count %% 2)){
      L.grid <- 2^seq(-2,20,l=5)
      pL.grid <- lapply(L.grid,function(L)pL(y,L))
      cost.grid <- sapply(pL.grid,calc.cost)
      best <- which.min(cost.grid)
      show.x <- log10(L.grid)
      show.y <- log10(cost.grid)
      plot(show.x,show.y,type="o")
      points(show.x[best],show.y[best],pch=20)
      last.L <- L.grid[best]
      pL.grid[[best]]
    }else{
      last.L <- L0
      pL(y, last.L)
    }
    last.t <- this.t
    this.t <- (1+sqrt(1+4*last.t^2))/2
    y <- this.iterate + (last.t - 1)/this.t*(this.iterate-last.iterate)
    ## here we calculate the subgradient optimality condition, which
    ## requires 1 more gradient evaluation per iteration.
    after.grad <- calc.grad(this.iterate,features,limits)
    dist2subgrad.opt <- function(w,g){
      ifelse(w==0,positive.part(abs(g)-gamma),
             ifelse(w<0,abs(-gamma+g),abs(gamma+g)))
    }
    w.dist <- dist2subgrad.opt(this.iterate[-1],after.grad[-1])
    zero.at.optimum <- c(abs(after.grad[1]),w.dist)
    stopping.crit <- max(zero.at.optimum)

    if(verbose >= 1){
      cost <- calc.cost(this.iterate)
      cat(sprintf("%10d cost %10f crit %10.7f L %f\n",
                  iterate.count,
                  cost,
                  stopping.crit,
                  last.L))
    }
    iterate.count <- iterate.count + 1
    if(iterate.count > max.iterations){
      stop(max.iterations," iterations, try increasing L0")
    }
  }
  this.iterate
}

smooth.interval.regression <- function
### Scale features and filter flat limits, then perform one interval
### regression and return a results list. The precise optimization
### problem is described in find.solution.
(features,
### Matrix n x p of inputs: n signals, each with p features that will
### be scaled.
 limits,
### Matrix n x 2 of output log(lambda). Each row corresponds to the
### lower and upper bound of an interval on the lambda which is
### optimal with respect to annotation error. Lower bound can be -Inf
### and upper bound can be Inf, which correspond to zero asymptotic
### cost.
 ...
### Passed to optimization code in find.solution. You must specify
### calc.grad and calc.loss.
 ){
  ## reality checks.
  stopifnot(nrow(features)==nrow(limits))
  ## dont know how to process missing data
  stopifnot(all(!is.na(features)))
  stopifnot(all(!is.na(limits)))
  if(ncol(limits)!=2){
    cat("str(limits)=\n")
    str(limits)
    stop("limits should be a 2-column matrix")
  }
  ## don't know how to process matrix with no colnames.
  stopifnot(is.matrix(features))
  stopifnot(!is.null(colnames(features)))
  
  ## check if there are any flat error curves, which have no limits.
  has.limits <- apply(is.finite(limits),1,any)
  ## filter zero-variance features.
  sigma <- apply(features[has.limits,,drop=FALSE],2,sd)
  zero.var <- sigma == 0
  if(any(zero.var)){
    cat("ignoring zero-variance variables:\n")
    print(names(sigma)[zero.var])
  }
  ## we train the model on this subset.
  some.limits <- limits[has.limits,,drop=FALSE]
  some.features <- features[has.limits,!zero.var,drop=FALSE]
  
  scaled <- scale(some.features)
  mu <- attr(scaled,"scaled:center")
  sigma <- attr(scaled,"scaled:scale")

  n <- nrow(scaled)
  p <- ncol(scaled)

  ## Solver:
  this.iterate <- find.solution(features=scaled, limits=some.limits, ...)
  sol <- list(intercept=this.iterate[1],
              weights=this.iterate[-1],
              mu=mu,sigma=sigma)
  sol$scaled <- scaled
  sol$log.limits <- some.limits
  ## this function will be applied to new data before applying the
  ## model.
  sol$normalize <- function(X){
    not.present <- !colnames(scaled) %in% colnames(X)
    if(any(not.present)){
      print(colnames(scaled)[not.present])
      stop("need all training features to predict")
    }
    X <- X[,colnames(scaled),drop=FALSE]
    mu.mat <- matrix(mu,nrow(X),ncol(X),byrow=TRUE)
    s.mat <- matrix(sigma,nrow(X),ncol(X),byrow=TRUE)
    (X-mu.mat)/s.mat
  }
  sol$f <- function(x){
    stopifnot(is.vector(x))
    stopifnot(length(x)==p)
    sum(x*sol$weights)+sol$intercept
  }
  sol$predict <- function(X){
    stopifnot(is.matrix(X))
    X.norm <- sol$normalize(X)
    L.hat <- (X.norm %*% sol$weights) + sol$intercept
    L.hat
  }
  sol$train.f <- apply(scaled,1,sol$f)
  sol$train.predict <- sol$predict(features)
  sol
### List of solver results. For a feature matrix X with p columns, you
### can use list$predict(X) to get model estimates of lambda.
}

regularized.interval.regression <- function
### Filter zero-variance features, scale features, filter flat limits,
### then perform a path of increasingly regularized interval
### regressions, returning a results list. We start at a small
### regularization parameter specified as gamma.initial. Then we use
### warm restarts, i.e. once we find the solution for one gamma, we
### use that as a starting point for the optimization problem with a
### larger gamma. We stop when gamma is so large that all the
### coefficients are 0 at the optimum.
(features,
### Matrix n x p of inputs: n signals, each with p features that will
### be scaled. Zero-variance features will be ignored.
 limits,
### Matrix n x 2 of output log(lambda). Each row corresponds to the
### lower and upper bound of an interval on the lambda which is
### optimal with respect to annotation error. Lower bound can be -Inf
### and upper bound can be Inf, which correspond to zero asymptotic
### cost.
 gamma.initial=1e-4,
### First regularization parameter in the path.
 gamma.factor=1.5,
### Multiplicative factor to increase gamma between steps in the path.
 ...
### Passed to optimization code in find.solution. You must specify
### calc.grad and calc.loss.
 ){
  ## reality checks.
  stopifnot(nrow(features)==nrow(limits))
  ## dont know how to process missing data
  stopifnot(all(!is.na(features)))
  stopifnot(all(!is.na(limits)))
  if(ncol(limits)!=2){
    cat("str(limits)=\n")
    str(limits)
    stop("limits should be a 2-column matrix")
  }
  stopifnot(is.matrix(features))
  
  ## check if there are any flat error curves, which have no limits.
  has.limits <- apply(is.finite(limits),1,any)
  ## filter zero-variance features.
  sigma <- apply(features[has.limits,,drop=FALSE],2,sd)
  zero.var <- sigma == 0
  if(any(zero.var)){
    cat("ignoring zero-variance variables:\n")
    print(names(sigma)[zero.var])
  }
  ## we train the model on this subset.
  some.limits <- limits[has.limits,,drop=FALSE]
  some.features <- features[has.limits,!zero.var,drop=FALSE]
  
  scaled <- scale(some.features)
  mu <- attr(scaled,"scaled:center")
  sigma <- attr(scaled,"scaled:scale")
  n <- nrow(scaled)
  p <- ncol(scaled)

  ## Pathwise solver:
  gamma <- gamma.initial
  coef.vec <- rep(1,p+1)
  coef.vec.list <- list()
  gamma.seq <- c()
  ## we stop the path when gamma is so large that the optimal solution
  ## is all zero coefficients.
  while(!all(coef.vec[-1] == 0)){
    cat(sprintf("gamma=%10f\n",gamma))
    gamma.seq <- c(gamma.seq,gamma)
    coef.vec <- find.solution(features=scaled,
                              limits=some.limits,
                              gamma=gamma,
                              starting.iterate=coef.vec,
                              ...)
    coef.vec.list[[length(gamma.seq)]] <- coef.vec
    gamma <- gamma * gamma.factor
  }
  coef.mat <- do.call(cbind,coef.vec.list)

  sol <- list(coefs=coef.mat,
              mu=mu,
              sigma=sigma,
              scaled=scaled,
              log.limits=some.limits,
              gamma.seq=gamma.seq)
  ## this function will be applied to new data before applying the
  ## model.
  sol$normalize <- function(X){
    not.present <- !colnames(scaled) %in% colnames(X)
    if(any(not.present)){
      print(colnames(scaled)[not.present])
      stop("need all training features to predict")
    }
    X <- X[,colnames(scaled)]
    mu.mat <- matrix(mu,nrow(X),ncol(X),byrow=TRUE)
    s.mat <- matrix(sigma,nrow(X),ncol(X),byrow=TRUE)
    (X-mu.mat)/s.mat
  }
  sol$predict <- function(X){
    stopifnot(is.matrix(X))
    X.norm <- sol$normalize(X)
    cbind(1,X.norm) %*% sol$coefs
  }
  sol$train.predict <- sol$predict(features)
  sol
### List of solver results. For a feature matrix X with p columns, you
### can use list$predict(X) to get model estimates of log(lambda).
}

library.install <- function(x,repos=getOption("repos"),type="source"){
  if(!require(x,character.only=TRUE)){
    install.packages(x,repos=repos,type=type)
    library(x,character.only=TRUE) ## library fails if pkg not found
  }
}
options(repos=c(#"http://cran.miroir-francais.fr/",
          "http://mirror.ibcp.fr/pub/CRAN/",
          #"http://cran.univ-lyon1.fr/",
          "http://cran.r-project.org"))

#library.install("quadmod",repos="http://r-forge.r-project.org")
library.install("quadprog")
library.install("ggplot2")
library.install("reshape2")
library.install("xtable")
library.install("tikzDevice")

hinge.interval.regression <- function
### Support vector interval regression using a quadratic programming
### (QP) solver. The idea is that we first normalize the feature
### matrix, giving normalized features x_i in R^p. Then we use a
### linear function f(x_i) = w'x_i + b to predict a log(lambda) that
### falls between the log limits L_i^left and L_i^right. So the
### optimization problem is: min_f ||f|| + sum_i C*hinge(L_i^left,
### L_i^right, f(x_i)). Since we assume f is linear the problem
### becomes min_{w,b,z} w'w + sum_i C*z_i, with these constraints for
### all relevant i: z_i >= 0, z_i >= 1 + L_i^left - b - w'x_i, z_i >=
### 1 - b - w'x + L_i^right. We call z_i slack, w weights, b intercept.
(features,
### Matrix n x p of inputs: n signals, each with p features. We will
### scale these internally.
 limits,
### Matrix n x 2 of output lambda. Each row corresponds to the lower
### and upper bound of an interval on the log(lambda) which is optimal
### with respect to annotation error. Lower bound can be -Inf and
### upper bound can be Inf, which correspond to zero asymptotic
### cost. 
 tune.C=1,
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
  vars <- make.ids(slack=n,intercept=1,weights=p)
  constraints <- list()
  for(i in 1:n){
    if(verbose >= 1)cat(sprintf("slack example constraints %5d / %5d",i,n))

    left <- some.limits[i,1]
    if(is.finite(left)){
      ivars <- with(vars,{
        intercept * 1 + sum(weights)*scaled[i,] + slack[i]
      })
      constraints <- c(constraints,list(ivars >= 1 + left))
    }

    right <- some.limits[i,2]
    if(is.finite(right)){
      ivars <- with(vars,{
        intercept * -1 + sum(weights)*scaled[i,]*-1 + slack[i]
      })
      constraints <- c(constraints,list(ivars >= 1 - right))
    }

    ## positivity.
    ## when we have 2 limits \__/ 3 constraints are necessary in this
    ## case, but not in this case \/. 
    gap <- right-left
    if(verbose >= 1)cat(sprintf(" gap=%4.2f",gap))
    if( (!is.finite(gap)) || (gap > 2) ){
      if(verbose >= 1)cat(" positivity constraint")
      constraints <- c(constraints,vars$slack[i] >= 0)
    }

    if(verbose >= 1)cat("\n")

  }
  const.info <- standard.form.constraints(constraints,vars)
  n.vars <- length(unlist(vars))
  Dvec <- rep(1e-6,n.vars)
  Dvec[vars$weights] <- 1
  D <- diag(Dvec)
  d <- rep(0,n.vars)
  d[vars$slack] <- -tune.C ## like C in svm
  if(verbose >= 1)cat(sprintf("solving for %d variables and %d constraints... ",
              n.vars,length(constraints)))
  sol <- solve.QP(D,d,const.info$A,const.info$b0)
  if(verbose >= 1)cat("solved!\n")
  sol$mu <- mu
  sol$sigma <- sigma
  sol$scaled <- scaled
  sol$log.limits <- some.limits
  sol$weights <- sol$solution[vars$weights]
  sol$intercept <- sol$solution[vars$intercept]
  sol$slack <- sol$solution[vars$slack]
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
