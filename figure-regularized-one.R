works_with_R("3.1.2",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             "Rdatatable/data.table@31cf650e6fe54fe76b73086e71ad4815bc7eeb8c",
             dplyr="0.4.0",
             directlabels="2014.6.13")

load("regularized.one.RData")

selected <- regularized.one$errors %>%
  filter(set=="validation") %>%
  arrange(arclength) %>%
  group_by(seed) %>%
  filter(seq_along(errors) == which.min(errors)) %>%
  select(regularization, arclength)

best.weights <- inner_join(regularized.one$weights, selected)
nonzero.coefs <- best.weights %>%
  filter(value != 0)

feature.ranks <- best.weights %>%
  group_by(feature) %>%
  summarise(nonzero=sum(value != 0),
            mean=mean(value),
            sd=sd(value)) %>%
  arrange(-nonzero, -abs(mean))

best.stats <- best.weights %>%
  mutate(feature.fac=factor(feature, feature.ranks$feature))

dots <- 
ggplot()+
  ggtitle("features ranked by number of nonzero splits, mean weight")+
  xlab("weight")+
  ylab("feature")+
  geom_point(aes(value, feature.fac), data=best.stats, pch=1)

data.frame(feature.ranks)

lasso <-
  ggplot()+
  geom_vline(aes(xintercept=arclength), data=selected, color="grey")+
  xlab("model complexity (L1 norm of weights)")+
  geom_line(aes(arclength, value, group=feature, color=feature),
            data=regularized.one$weights)+
  geom_line(aes(arclength, percent, group=set, linetype=set),
            data=regularized.one$errors, show_guide=FALSE)+
  geom_dl(aes(arclength, percent, label=set),
          data=regularized.one$errors, method="lines2")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(what ~ seed, scales="free", labeller=function(var, val){
    if(var=="seed"){
      paste("split", val)
    }else{
      paste(val)
    }
  })

lasso.selected <- list(function(df, ...){
  subset(df, groups %in% nonzero.coefs$feature)
}, "lasso.labels")

lasso.dl <-
  direct.label(lasso, "lasso.selected")

pdf("figure-regularized-one.pdf", w=12)
print(lasso.dl)
print(dots)
dev.off()
