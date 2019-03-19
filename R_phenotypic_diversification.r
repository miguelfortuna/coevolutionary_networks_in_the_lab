##################################################################################
### DOES PHENOTYPIC DIVERSIFICATION INCREASE OVER TIME AND LEVEL OF RESOURCES? ###
##################################################################################
library(tidyverse)
library("lme4")

data <- read_csv('phenotypic_diversification.csv')
data <- data %>%
  unite(unique_ID, resources, replicate, remove=FALSE) %>%
  unite(indiv_ID, resources, replicate, organism, remove=FALSE) %>%
  mutate(replicate = as.factor(replicate),
         resources = as.factor(resources),
         organism = as.factor(organism))

lmer_phen <- lmer(nPhen ~
                # fixed effects
                resources * timeStep * organism
                # random effects
                + (1|unique_ID) # same as (1 | replicate:resources)
                + (1|indiv_ID), # same as (1 | replicate:resources:organism)
              data = data)
summary(lmer_phen)
anova(lmer_phen, type=1, ddf = "Kenward-Roger")

predict_effects <- ggeffect(lmer_phen, terms = c("timeStep","organism"))
data %>%
  ggplot(., aes(x=timeStep, y=nPhen, color=organism)) +
  geom_ribbon(data=rename(as.data.frame(predict_effects), organism=group), aes(x=x, ymin=conf.low, ymax=conf.high), inherit.aes = F, alpha=0.15) +
  geom_line(data=rename(as.data.frame(predict_effects), organism=group), aes(x=x, y=predicted), size=1.5) +
  #geom_point(aes(size=network_size, alpha=connectance)) +
  facet_wrap(~organism) +
  scale_color_manual(values=c("blue","red")) +
  ylab("Phenotypes") +
  scale_x_continuous(breaks=1:6, name="Time")
