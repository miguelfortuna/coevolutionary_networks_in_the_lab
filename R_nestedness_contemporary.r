###############################
### LOAD REQUIRED LIBRARIES ###
###############################
library(tidyverse)
library(cowplot)
library(piecewiseSEM)
library(visreg)
library(lme4)
library(lmerTest) # for easy calculation of P-values
library(ggeffects)
library(boot) # for logit transformation
library(viridisLite)
library(svglite)

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
### LOAD DATA ###
#################
data <- read_csv("nestedness_contemporary.csv")
data <- data %>%
  unite(unique_ID, resources, replicate, remove = FALSE) %>%
  unite(indiv_ID, resources, replicate, timeStep, remove=FALSE) %>%
  mutate(network_size = n_phagesPhen*n_bacteriaPhen,
         replicate = as.factor(replicate),
         resources = as.factor(resources)) %>% # for ggpredict
  filter(!is.na(nestedness)) # remove NA values
data # 53 data points

#########################################################################
### NETWORK SIZE EXPLAINS THE PREVALENCE OF PERFECTLY NESTED NETWORKS ###
#########################################################################
###################################################################################################
# generalized linear mixed model (accounts for non-independence through random effects: replicate #
###################################################################################################
data <- mutate(data, nested = ifelse(nestedness < 1, 0, 1))
nested_glmer <- glmer(nested ~ 1 
                    + network_size+timeStep*resources
                    + (1|replicate:resources),
                    data=data,
                    family="binomial")
summary(nested_glmer)

###############################################################################
# test overdispersion of the response variable for this binomial distribution #
###############################################################################
overdisp_fun(nested_glmer)
# overdispersion of the distribution that we have used: 
# ratio = 0.639; p-value=0.97 (there is no overdispersion)

#####################
# fitting the model #
#####################
sem.model.fits(nested_glmer)
# fixed effects and interaction term explain 64% of the total variance (marginal)
# fixed effects + interaction term + replicate as random effect explain 77% of the total variance (conditional)
# compare our glmer model with a model without network size
anova(update(nested_glmer, .~. -network_size), nested_glmer, test = "LR")
# X^2,(df=1) = 22.93, p < 0.001
# compare our glmer model with a model without the interaction term
anova(update(nested_glmer, .~. -timeStep:resources), nested_glmer, test = "LR")
# X^2,(df=1) = 0.37, p = 0.544
# compare our glmer model without the interaction term with our glmer model without the interaction term plus resources
anova(update(nested_glmer, .~. -timeStep:resources - resources), update(nested_glmer, .~. -timeStep:resources), test = "LR")
# X^2,(df=1) = 3.69, p = 0.055
# compare our glmer model without the interaction term with our glmer model without the interaction term plus timeStep
anova(update(nested_glmer, .~. -timeStep:resources - timeStep), update(nested_glmer, .~. -timeStep:resources), test = "LR")
# X^2,(df=1) = 0.03, p = 0.544

########
# plot #
########
visreg(nested_glmer, xvar="network_size", by="resources", scale="response", overlay=T, gg=T)




########################################
### CHANGE IN NETWORK SIZE OVER TIME ###
########################################
###################################################################################################
# generalized linear mixed model (accounts for non-independence through random effects: replicate #
###################################################################################################

net_size <- glmer(network_size ~ timeStep+resources
                  + (1|replicate:resources)
                  + (1|replicate:timeStep:resources),
                  data=data, family=poisson)
summary(net_size)

##############################################################################
# test overdispersion of the response variable for this poisson distribution #
##############################################################################
overdisp_fun(net_size)
# overdispersion of the distribution that we have used: 
# ratio = 0.059; p-value=1 (there is no overdispersion)

#####################
# fitting the model #
#####################
sem.model.fits(net_size)
# fixed effects explain 14% of the total variance (marginal)
# fixed effects + replicate as random effect explain 59% of the total variance (conditional)
# compare our glmer model with a model without timeStep (i.e., effect of timeStep)
anova(update(net_size, .~. -timeStep), net_size, test = "LR")
# X^2,(df=1) = 10.30, p = 0.001

########
# plot #
########
visreg(net_size, xvar="timeStep", by="resources", scale="response", overlay=T)



######################################
### CHANGE IN NESTEDNESS OVER TIME ###
######################################
#############################
# explore nestedness values #
#############################
hist(filter(data, nestedness < 1)$nestedness)
hist(logit(filter(data, nestedness < 1)$nestedness)) # clear evidence that we should logit-transform the data

###################################################################################################
# linear mixed model (accounts for non-independence through random effects: replicate (unique_ID) #
###################################################################################################
contemporary_nestedness <- lmer(logit(nestedness) ~ 
                                  # fixed effects
                                  network_size + connectance + timeStep*resources +
                                  # random effects
                                  (1|unique_ID), # same as (1 | replicate:resources)
                                # subset data to nestedness values < 1
                                data=filter(data, nestedness < 1))
summary(contemporary_nestedness)
# both methods give the same qualitative results
anova(contemporary_nestedness, type=1, ddf = "Kenward-Roger") 
#anova(contemporary_nestedness, type=2, ddf = "Kenward-Roger")


###############
# plot Fig. 3 #
###############
predict_effects <- ggeffect(contemporary_nestedness, terms = c("timeStep","resources"))
my_breaks <- c(0,2,4,6)
my_labels <- round(inv.logit(my_breaks),2)

fig_3 <- filter(data, nestedness < 1) %>%
  ggplot(., aes(x=timeStep, y=logit(nestedness), color=resources)) +
  geom_ribbon(data=rename(as.data.frame(predict_effects), resources=group), aes(x=x, ymin=conf.low, ymax=conf.high), inherit.aes = F, alpha=0.15) +
  geom_line(data=rename(as.data.frame(predict_effects), resources=group), aes(x=x, y=predicted), size=1.5) +
  geom_point(aes(size=network_size, alpha=connectance)) +
  facet_wrap(~resources) +
  scale_color_manual(values=c("blue","red")) +
  ylab("logit(Nestedness)") +
  scale_x_continuous(breaks=1:6, name="Time")

save_plot(filename = "../figures/figure_3.svg", plot = fig_3, base_height = 7, base_width = 8)

