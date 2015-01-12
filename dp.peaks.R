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
dp.peaks <- list()
for(file.i in seq_along(files)){
  r <- matched[file.i, ]
  set.name <- r[["set_name"]]
  chunk.id <- r[["chunk_id"]]
  regions.str <- paste0(set.name, "/", chunk.id)
  f <- files[[file.i]]
  cat(sprintf("%4d / %4d %s\n", file.i, length(files), f))
  load(f)
  for(sample.id in names(dp.model)){
    peak.list <- dp.model[[sample.id]]$peaks
    dp.peaks[[regions.str]][[sample.id]] <- peak.list
  }
}

save(dp.peaks, file="dp.peaks.RData")
