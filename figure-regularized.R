works_with_R("3.1.2",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             data.table="1.9.4",
             dplyr="0.4.0",
             directlabels="2014.6.13")

load("regularized.RData")

selected <- regularized$errors %>%
  filter(set=="validation") %>%
  group_by(seed) %>%
  filter(seq_along(errors) == which.min(errors)) %>%
  select(regularization, arclength)

best.weights <- inner_join(regularized$weights, selected)
nonzero.coefs <- best.weights %>%
  filter(value != 0)

feature.ranks <- best.weights %>%
  group_by(feature) %>%
  summarise(nonzero=sum(value != 0)) %>%
  arrange(-nonzero)

data.frame(feature.ranks)

lasso <-
  ggplot()+
  geom_vline(aes(xintercept=arclength), data=selected, color="grey")+
  xlab("model complexity (L1 norm of weights)")+
  geom_line(aes(arclength, value, group=feature, color=feature),
            data=regularized$weights)+
  geom_line(aes(arclength, errors, group=set, linetype=set),
            data=regularized$errors, show_guide=FALSE)+
  geom_dl(aes(arclength, errors, label=set),
          data=regularized$errors, method="lines2")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(what ~ seed, scales="free")

lasso.selected <- list(function(df, ...){
  subset(df, groups %in% nonzero.coefs$feature)
}, "lasso.labels")

lasso.dl <-
  direct.label(lasso, "lasso.selected")

pdf("figure-regularized.pdf", w=12)
print(lasso.dl)
dev.off()
