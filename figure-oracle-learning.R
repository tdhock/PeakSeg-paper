works_with_R("3.1.3",
             directlabels="2014.6.13",
             data.table="1.9.4",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             dplyr="0.4.0")

load("oracle.learning.RData")

error.lines <- oracle.learning %>%
  group_by(set.name, set.i, train.size, model.name) %>%
  summarise(errors=sum(errors),
            regions=sum(regions)) %>%
  mutate(percent=errors/regions*100,
         model.name=sub("default", "0", sub("trained", "1", model.name)),
         model.type=ifelse(grepl("oracle", model.name), "cDPA", "baselines"),
         algorithm=ifelse(grepl("[.]1$", model.name), "grid\nsearch",
           ifelse(grepl("[.]0$", model.name), "unsupervised",
                  "interval\nregression")))

data.frame(subset(error.lines, set.name==set.name[1] & set.i==set.i[1]))

## Make sure the same number of regions were used to compute test
## error for each of the models.
perror.check <- error.lines %>%
  group_by(set.name, set.i) %>%
  summarise(min.regions=min(regions),
            max.regions=max(regions))
with(error.check, stopifnot(min.regions == max.regions))

error.bands <- error.lines %>%
  group_by(set.name, train.size, algorithm, model.name, model.type) %>%
  summarise(mean=mean(percent),
            sd=sd(percent))

algo.colors <-
  c("cheating"="grey",
    "interval\nregression"="#D95F02",
    "grid\nsearch"="#1B9E77",
    "unsupervised"="#7570B3")

with.legend <- 
ggplot()+
  scale_color_manual(values=algo.colors)+
  scale_fill_manual(values=algo.colors)+
  geom_line(aes(train.size, mean,
                color=algorithm,
                group=model.name),
            data=error.bands,
            size=2)+
  geom_ribbon(aes(train.size, ymin=mean-sd, ymax=mean+sd,
                  group=model.name,
                  fill=algorithm),
              data=error.bands,
              alpha=0.25)+
  geom_line(aes(train.size, percent,
                color=algorithm,
                group=interaction(set.i, model.name)),
            data=error.lines)+
  geom_line(aes(train.size, mean, group=model.name),
            data=error.bands,
            size=1)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(set.name ~ model.type)+
  ylab("percent test error")

with.labels <-
  with.legend+
  coord_cartesian(ylim=c(0, 40))+
  scale_x_continuous("training set size (groups of labeled regions)",
                     breaks=seq(0, 10, by=2),
                     limits=c(-1, 13))+
  geom_dl(aes(train.size, mean,
                color=algorithm,
                label=model.name),
          data=error.bands,
          method=dl.combine("first.polygons", "last.polygons"))

pdf("figure-oracle-learning.pdf", w=10)
print(with.labels)
dev.off()
