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
  mutate(percent=errors/regions*100)

## Make sure the same number of regions were used to compute test
## error for each of the models.
error.check <- error.lines %>%
  group_by(set.name, set.i, train.size) %>%
  summarise(min.regions=min(regions),
            max.regions=max(regions))
with(error.check, stopifnot(min.regions == max.regions))

error.bands <- error.lines %>%
  group_by(set.name, train.size, model.name) %>%
  summarise(mean=mean(percent),
            sd=sd(percent))

no.legend <- 
ggplot()+
  geom_line(aes(train.size, mean, color=model.name),
            data=error.bands,
            size=2)+
  geom_ribbon(aes(train.size, ymin=mean-sd, ymax=mean+sd,
                  fill=model.name),
              data=error.bands,
              alpha=0.25)+
  geom_line(aes(train.size, percent, color=model.name,
                group=interaction(set.i, model.name)),
            data=error.lines)+
  geom_line(aes(train.size, mean, group=model.name),
            data=error.bands,
            size=1)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(set.name ~ .)

with.legend <-
  direct.label(no.legend, dl.combine("first.qp", "last.qp"))

## TODO: add oracle.0 and oracle.1 models?

pdf("figure-oracle-learning.pdf")
print(with.legend)
dev.off()
