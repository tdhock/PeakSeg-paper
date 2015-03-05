works_with_R("3.1.2",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             reshape2="1.2.2",
             data.table="1.9.4",
             directlabels="2014.6.13",
             dplyr="0.4.0")

load("dp.peaks.sets.RData")

all.best.list <- list()
for(set.name in names(dp.peaks.sets)){
  train.sets <- dp.peaks.sets[[set.name]]
  label.dir <- paste0("~/genomelabels/", set.name)
  RData.files <- Sys.glob(file.path(label.dir, "*.RData"))
  names(RData.files) <- sub(".RData$", "", basename(RData.files))
  RData.by.sample <- list()
  e <- new.env()
  fl.list <- list()
  error.list <- list()
  for(sample.id in names(RData.files)){
    RData.file <- RData.files[[sample.id]]
    objs <- load(RData.file, e)
    RData.by.sample[[sample.id]] <- as.list(e)
    error.list[[sample.id]] <- data.table(sample.id, e$errors)
    for(bases.per.bin in names(e$features.limits)){
      fl <- e$features.limits[[bases.per.bin]]
      chunk.list <- fl$limits
      for(chunk.id in names(chunk.list)){
        limits <- chunk.list[[chunk.id]]
        chunk.name <- paste0(set.name, "/", chunk.id)
        features <- fl$features[rownames(limits), , drop=FALSE ]
        rownames(features) <- rownames(limits) <-
          paste(chunk.name, bases.per.bin, rownames(limits))
        fl.list[[bases.per.bin]][[chunk.name]][[sample.id]] <-
          list(features=features, limits=limits)
      }
    }
  }
  all.errors <- do.call(rbind, error.list)
  all.errors[, chunk.name := paste0(set.name, "/", chunk.id)]
  for(set.i in seq_along(train.sets)){
    testSet <- paste(set.name, "split", set.i)
    train.validation <- train.sets[[set.i]]

    error.by.sample <- 
    all.errors %>%
      filter(chunk.name %in% train.validation) %>%
      group_by(bases.per.bin, sample.id) %>%
      summarise(errors=sum(errors),
                regions=sum(regions))

    each.sample <- 
    ggplot(error.by.sample, aes(bases.per.bin, errors, group=sample.id))+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free")+
    geom_line()+
    geom_point()+
    scale_x_log10()

    error.wide <- error.by.sample %>%
      group_by(bases.per.bin) %>%
      summarise(errors=sum(errors),
                regions=sum(regions)) %>%
      mutate(percent=errors/regions)
    print(error.wide)
    error.tall <-
      melt(error.wide, measure.vars=c("errors", "percent"))

    best.error <- error.tall %>%
      group_by(variable) %>%
      filter(seq_along(value) == which.min(value))

    ggplot(NULL, aes(bases.per.bin, value))+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(variable ~ ., scales="free")+
    geom_line(data=error.tall)+
    geom_point(data=best.error)+
    scale_x_log10()

    all.best.list[[testSet]] <- data.table(set.name, testSet, best.error)

    n.folds <- if(length(train.validation) == 2) 2 else 3
    fold.id <- rep(1:n.folds, l=length(train.validation))

    e.list <- list()
    w.list <- list()
    s.list <- list()
    for(validation.fold in 1:n.folds){
      is.validation <- fold.id == validation.fold
      ## Randomly divide the train set into 1/2 train, 1/2 validation.
      set.chunks <-
        list(train=train.validation[!is.validation],
             validation=train.validation[is.validation])

    }#validation.fold
  }#set.i
}#set.name

all.best <- do.call(rbind, all.best.list)

L <- split(all.best, all.best$set.name)
lapply(L, with, table(bases.per.bin, variable))
