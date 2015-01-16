works_with_R("3.1.2",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             "Rdatatable/data.table@fd5464e4c587d4a96fa32c809eb94a45bb361133",
             dplyr="0.4.0",
             directlabels="2014.6.13")

load("regularized.all.RData")

percent <- regularized.all$error %>%
  group_by(set.name, set.i, model.name) %>%
  summarise(errors=sum(errors),
            regions=sum(regions)) %>%
  mutate(percent=errors/regions*100)
wide <- dcast(percent, set.name + set.i ~ model.name, value.var="percent")

compare <-
  ggplot()+
  geom_point(aes(percent, model.name),
             data=percent, pch=1)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(. ~ set.name)

coef.list <- list()
for(set.name in names(regularized.all$models)){
  splits <- regularized.all$models[[set.name]]
  for(testSet in names(splits)){
    fit.list <- splits[[testSet]]
    for(model.name in names(fit.list)){
      fit <- fit.list[[model.name]]
      fit.dt <- 
      rbind(data.table(feature=names(fit$weights), weight=fit$weights),
            data.table(feature="intercept", weight=fit$intercept))
      coef.list[[paste(set.name, testSet, model.name)]] <-
        data.table(set.name, testSet, model.name, fit.dt)
    }
  }
}
coef.dt <- do.call(rbind, coef.list)
all.stats <- coef.dt %>%
  filter(feature != "intercept") %>%
  group_by(set.name, testSet, model.name) %>%
  summarise(nonzero=sum(weight != 0),
            mean=mean(weight),
            sd=sd(weight)) %>%
  arrange(-nonzero, -abs(mean))

coef.two <- coef.dt %>%
  filter(feature %in% c("log.bases", "log.unweighted.quartile.100%"),
         model.name != "L1.reg.log")
coef.wide <- dcast(coef.two, set.name + testSet + model.name ~ feature) %>%
  inner_join(all.stats)

br <- seq(-2, 2, by=0.5)
two.features.scatter <- 
ggplot()+
  guides(color=guide_legend(direction="horizontal", ))+
  geom_hline(yintercept=0, color="grey")+
  geom_vline(xintercept=0, color="grey")+
  ## geom_text(aes(`log.unweighted.quartile.100%`, log.bases,
  ##               label=nonzero,
  ##               color=model.name),
  ##            data=coef.wide, pch=1)+
  geom_point(aes(`log.unweighted.quartile.100%`, log.bases,
                color=model.name),
             data=coef.wide, pch=1)+
  scale_x_continuous("log(max(count)) normalized weight",
                     breaks=br)+
  scale_y_continuous("log(bases) normalized weight", breaks=br)+
  theme_bw()+
  coord_equal()+
  theme(panel.margin=grid::unit(0, "cm"),
        legend.position="top")+
  facet_grid(. ~ set.name, labeller=function(var, val){
    gsub("_", "\n", val)
  })

nonzero <- coef.dt %>%
  filter(weight != 0)

dots <- 
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(model.name ~ set.name)+
  geom_point(aes(weight, feature), data=nonzero, pch=1)

dots <- 
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(model.name ~ feature)+
  geom_point(aes(weight, set.name), data=nonzero, pch=1)

coef.l1 <- coef.dt %>%
  filter(model.name == "L1.reg",
         feature != "intercept")
coef.stats <- coef.l1 %>%
  group_by(feature) %>%
  summarise(nonzero=sum(weight != 0),
            mean=mean(weight),
            sd=sd(weight)) %>%
  arrange(-nonzero, -abs(mean))
coef.l1.nonzero <- coef.l1 %>%
  filter(weight != 0) %>%
  mutate(feature.fac=factor(feature, coef.stats$feature))
nonzero.plot <- 
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(. ~ set.name)+
  geom_point(aes(weight, feature), data=coef.l1.nonzero, pch=1)

pdf("figure-regularized-all.pdf", h=4, w=7)
print(two.features.scatter)
dev.off()
