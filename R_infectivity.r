library("tidyverse")
library("betapart")
library("effects")
library("car")
library("piecewiseSEM")
library("lme4")
library("sjstats")
library("visreg")
library("cowplot")
library("svglite")

#################################################################################################################
# function to test overdispersion (https://rdrr.io/github/markushuff/PsychHelperFunctions/src/R/overdisp_fun.R) #
#################################################################################################################
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

#################
### load data ###
#################
data <- read_csv('infectivity.csv') %>%
  distinct() %>% # from genoytpes (isolates) to phenotypes
  group_by(resources, replicate, bacteriaPhen, phagePhen) %>%
  filter(n() == 1 | (n() > 1 & timeStepBacteria > min(timeStepBacteria) | timeStepPhage > min(timeStepPhage))) %>%
  ungroup() %>%
  # 27,911 pairwise interactions among phenotypes sampled in only one point in time
  # 12,678 pairwise interactions among phenotypes sampled in more than one point in time
  mutate(timeShift = ifelse(timeStepPhage == timeStepBacteria, "contemporary",
                            ifelse(timeStepPhage > timeStepBacteria, "phage_more_evolved",
                                   "bacteria_more_evolved"))) %>%
  select(-timeStepBacteria, -timeStepPhage) %>%
  distinct() %>% # 7,022 interactions were removed because they were already considered the first time the phenotypes were sampled (i.e., duplicates)
  # 27,911 + 12,678 - 7,022 = 33,567 interactions used to infer coevolutionary dynamics
  mutate(resources = factor(resources)) %>%
  mutate(resources = relevel(resources, ref = "low")) %>%
  mutate(replicate = factor(replicate)) %>% 
  mutate(timeShift = factor(timeShift)) %>%
  mutate(timeShift = relevel(timeShift, ref = "contemporary"))

###################################################################################################
# generalized linear mixed model (accounts for non-independence through random effects: replicate #
###################################################################################################
glmer_interact <- glmer(infection ~  
                          timeShift * resources + # fixed effects
                          (1 | resources:replicate) +  # random effect nested within fixed effect (resources)
                          (1 | resources:replicate:timeShift),
                        data = data, family = binomial)
summary(glmer_interact)

###############################################################################
# test overdispersion of the response variable for this binomial distribution #
###############################################################################
overdisp_fun(glmer_interact)
# overdispersion of the distribution that we have used: ratio = 0.994; p-value=0.799 (there is no overdispersion)

#################################################################################
# get the predicted values for our response variable (probability of infection) #
#################################################################################
combos <- expand.grid(resources=c("low", "high"), timeShift=c("bacteria_more_evolved", "contemporary", "phage_more_evolved"))
combos_predict <- mutate(combos,
                         p_resistance = 1 - predict(glmer_interact, newdata = combos, re.form=NA, type="response"),
                         p_infection = predict(glmer_interact, newdata = combos, re.form=NA, type="response"))
combos_predict

#############################################################################################
# quantifying the effect of each predictor variable (main effects and the interaction term) #
#############################################################################################
# compare our full glmer model with one without the interaction term
anova(update(glmer_interact, .~. -resources:timeShift), glmer_interact, test = "LR")
# X^2,(df=2) = 10.15, p = 0.006 (the effect of the interaction term is significant)
# compare our glmer model without the interaction term with our glmer model without the interaction term plus resources
anova(update(glmer_interact, .~. -resources:timeShift -resources), update(glmer_interact, .~. -resources:timeShift), test = "LR")
# X^2,(df=1) = 4.38, p = 0.036 (the effect of the resources is statistically significant)
# compare our glmer model without the interaction term with our glmer model without the interaction term plus timeShift
anova(update(glmer_interact, .~. -resources:timeShift -timeShift), update(glmer_interact, .~. -resources:timeShift), test = "LR")
# X^2,(df=2) = 11.5, p = 0.003 (the effect of the timeShift is statistically very significant)

#################################################
# plot the interaction between the main effects #
#################################################
summary(effect(c("timeShift*resources"), glmer_interact))

### infectivity ###
infectivity_data <- as.data.frame(effect(c("timeShift*resources"), glmer_interact))
fig_2a <- ggplot(infectivity_data, aes(timeShift, fit, color=resources)) +
  geom_point() +
  ylim(0, 1) +
  geom_line(aes(colour=resources, group=resources)) +
  geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) +
  theme_classic(base_size=12)
save_plot(filename = "../figures/figure_2a.svg", plot = fig_2a, base_height = 7, base_width = 4)

### resistance ###
resistance_data <- as.data.frame(effect(c("timeShift*resources"), glmer_interact)) %>%
  mutate(fit=1-fit)
fig_2b <- ggplot(resistance_data, aes(timeShift, fit, color=resources)) +
  geom_point() +
  ylim(0, 1) +
  geom_line(aes(colour=resources, group=resources)) +
  geom_errorbar(aes(ymin=fit-se, ymax=fit+se), width=0.4) +
  theme_classic(base_size=12)
save_plot(filename = "../figures/figure_2b.svg", plot = fig_2b, base_height = 7, base_width = 4)