works_with_R("3.1.1",
             dplyr="0.2",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             directlabels="2014.6.13")

load("dp.timings.RData")

dp.timings %.%
  arrange(seconds) %.%
  tail()

dp.timings %.%
  arrange(data) %.%
  tail()

refs <- data.frame(seconds=c(60, 60*60, 60*60*24, 60*60*24*7, 1),
                   unit=c("1 minute", "1 hour", "1 day", "1 week", "1 second"))
ref.color <- "grey50"
algo.colors <-
  c(DP="#1B9E77", Segmentor="#D95F02", cghseg="#7570B3")
  c(DP="#E41A1C", Segmentor="#377EB8", cghseg="#4DAF4A")

all.timings <-
  rbind(#data.frame(algorithm="Segmentor", Segmentor.timings),
        data.frame(algorithm="DP", dp.timings))

zoom.in <- 
  ggplot()+
  geom_text(aes(4000, seconds, label=unit),
            data=refs, vjust=1.5, hjust=0, color=ref.color)+
  geom_hline(aes(yintercept=seconds),
             data=refs, color=ref.color)+
  geom_point(aes(data, seconds, color=algorithm),
             data=all.timings, pch=1)+
  scale_x_continuous("data to segment")+
  scale_color_manual(values=algo.colors)+
  scale_fill_manual(values=algo.colors)+
  coord_cartesian(xlim=c(0, 5000), ylim=c(-5, 61))
pdf("figure-dp-timings-zoom.pdf", h=5)
print(zoom.in)
dev.off()

ggplot()+
  scale_color_manual(values=algo.colors)+
  scale_x_log10()+
  scale_y_log10()+
  geom_point(aes(bases, seconds, color=algorithm),
             data=all.timings, pch=1)

ggplot()+
  scale_color_manual(values=algo.colors)+
  scale_x_log10()+
  scale_y_log10()+
  geom_point(aes(bases, data, color=algorithm),
             data=all.timings, pch=1)

some.refs <- refs[1:2,]

p <- ggplot()+
  geom_text(aes(1e5, seconds, label=unit),
            data=some.refs, vjust=1.5, hjust=0, color=ref.color)+
  geom_hline(aes(yintercept=seconds),
             data=some.refs, color=ref.color)+
  scale_color_manual(values=algo.colors)+
  geom_point(aes(data, seconds, color=algorithm),
             data=all.timings, pch=1)+
  xlab("data points to segment")+
  ggtitle("DP algorithm time complexity quadratic in data size")

pdf("figure-dp-timings.pdf")
print(p)
dev.off()
