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
cell.type <- as.character(sample.counts$cell.type[1])

sample.counts$weight <- with(sample.counts, chromEnd-chromStart)
seg1.ends <- as.integer(c(1000, 2000, 4000))
loss.list <- list()
for(model.i in seq_along(seg1.ends)){
  seg.end <- seg1.ends[[model.i]]
  seg.start <- 1
  seg.data <- sample.counts[seg.start:seg.end, ]
  seg.mean <- with(seg.data, sum(coverage * weight)/sum(weight))
  seg.loss <- with(seg.data, PoissonLoss(coverage, seg.mean, weight))
  loss.list[[paste(model.i)]] <-
    data.table(model.i, 
               seg.start, seg.end,
               seg.chromStart=sample.counts$chromStart[seg.start],
               seg.chromEnd=sample.counts$chromEnd[seg.end],
               seg.mean, seg.loss)
}
loss.dt <- do.call(rbind, loss.list)

show.loss.list <- split(loss.dt, loss.dt$model.i)
png.list <- list()
for(model.i in seq_along(show.loss.list)){
  show.loss <- show.loss.list[[model.i]]
  selectedPlot <- 
  ggplot()+
    geom_vline(aes(xintercept=seg.chromEnd/1e3), data=show.loss,
               color="grey")+
    geom_text(aes(seg.chromEnd/1e3, max(sample.counts$coverage), label="t "),
              data=show.loss, hjust=1, vjust=1, color="grey")+
  geom_step(aes(chromStart/1e3, coverage),
            data=sample.counts,
            color="grey50")+
  geom_segment(aes(seg.chromStart/1e3, seg.mean,
                   xend=seg.chromEnd/1e3, yend=seg.mean),
               color="green",
               data=show.loss)+
                 theme_bw()+
  ylab("count of aligned reads")+
  xlab(paste("position on chromosome (kb = kilo bases)"))

  png(png.name <- sprintf("figure-dp-first-%d.png", model.i),
      units="in", res=200, width=6, height=3)
  print(selectedPlot)
  dev.off()

  png.list[[png.name]] <- png.name
}

pngs <- do.call(c, png.list)

png.tex <- sprintf("
\\begin{frame}
\\frametitle{Computation of optimal loss $\\mathcal L_{s, t}$
  for $s=1$ segments up to data point $t$}
  \\includegraphics[width=\\textwidth]{%s}

$$
\\mathcal L_{1, t} =
\\underbrace{
  c_{(0, t]}
}_{
  \\text{optimal loss of 1st segment $(0, t]$}
}
$$

\\end{frame}
", pngs)

cat(png.tex, file="figure-dp-first.tex")
