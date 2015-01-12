works_with_R("3.1.1",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732")

files <- Sys.glob("data/*/*/dp.model.RData")

## Parse the first occurance of pattern from each of several strings
## using (named) capturing regular expressions, returning a matrix
## (with column names).
str_match_perl <- function(string,pattern){
  stopifnot(is.character(string))
  stopifnot(is.character(pattern))
  stopifnot(length(pattern)==1)
  parsed <- regexpr(pattern,string,perl=TRUE)
  captured.text <- substr(string,parsed,parsed+attr(parsed,"match.length")-1)
  captured.text[captured.text==""] <- NA
  captured.groups <- do.call(rbind,lapply(seq_along(string),function(i){
    st <- attr(parsed,"capture.start")[i,]
    if(is.na(parsed[i]) || parsed[i]==-1)return(rep(NA,length(st)))
    substring(string[i],st,st+attr(parsed,"capture.length")[i,]-1)
  }))
  result <- cbind(captured.text,captured.groups)
  colnames(result) <- c("",attr(parsed,"capture.names"))
  result
}

pattern <-
  paste0("data/",
         "(?<set_name>.+?)",
         "/",
         "(?<chunk_id>[0-9]+)")
matched <- str_match_perl(files, pattern)
dp.peaks.features <- list()
for(file.i in seq_along(files)){
  r <- matched[file.i, ]
  set.name <- r[["set_name"]]
  chunk.id <- r[["chunk_id"]]
  regions.str <- paste0(set.name, "/", chunk.id)
  f <- files[[file.i]]
  cat(sprintf("%4d / %4d %s\n", file.i, length(files), f))
  load(f)
  count.f <- sub("dp.model", "counts", f)
  load(count.f)
  count.list <- split(counts, counts$sample.id)
  sample.ids <- names(dp.model)
  chunk.mat <-
    matrix(NA, length(sample.ids), 14*4,
           dimnames=list(sample.id=sample.ids, feature=NULL))
  for(sample.id in sample.ids){
    mod <- dp.model[[sample.id]]
    peaks.num <- 0:9
    p.err <- if(nrow(mod$error) == 1){
      rep(mod$error$error, 10)
    }else{
      poisson.err <- rep(NA, 10)
      names(poisson.err) <- peaks.num
      poisson.err[as.character(mod$error$peaks)] <- mod$error$error
      plot(peaks.num, poisson.err)
      p.approx <- approx(peaks.num, poisson.err, peaks.num)
      with(p.approx, lines(x, y, col="red", pch=20))
      p.approx$y
    }
    bases <- with(count.df, chromEnd-chromStart)
    long <- rep(count.df$coverage, bases)
    n.bases <- sum(bases)
    n.data <- nrow(count.df)
    n.segments <- peaks.num*2 + 1
    under.sqrt <- 1.1 + log(n.data/n.segments)
    in.square <- 1 + 4 * sqrt(under.sqrt)
    cleynen.factor <- in.square * in.square
    cleynen <- n.segments * cleynen.factor
    s <- ifelse(peaks.num==1, "", "s")
    peaks.str <- paste0(peaks.num, "peak", s)
    names(p.err) <- paste0("PoissonError", peaks.str)
    names(under.sqrt) <- names(in.square) <- names(cleynen) <- peaks.str
    count.df <- count.list[[sample.id]]
    feature.vec <-
      c(unweighted.quartile=quantile(count.df$coverage),
        weighted.quartile=quantile(long),
        unweighted.mean=mean(count.df$coverage),
        weighted.mean=mean(long),
        bases=n.bases,
        ## sqrt=under.sqrt,
        ## square=in.square,
        ## cleynen=cleynen,
        ## p.err,
        data=n.data)        
    log.features <-
      c(feature.vec,
        `log+1`=log(feature.vec+1),
        log=log(feature.vec),
        log.log=log(log(feature.vec)))
    chunk.mat[sample.id, ] <- log.features
  }
  colnames(chunk.mat) <- names(log.features)
  colnames(chunk.mat)[colMeans(is.finite(chunk.mat)) == 1]
  dp.peaks.features[[regions.str]] <- chunk.mat
}

save(dp.peaks.features, file="dp.peaks.features.RData")
