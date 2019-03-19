library(tidyverse)
library(cowplot)
library(visreg)
library(lme4)

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

#############################################################################
### COEVOLUTION NOR TIME EXPLAIN THE CONNECTANCE OF CONTEMPORARY NETWORKS ###
#############################################################################
data <- read_csv("nestedness_contemporary.csv")
data <- data %>%
  unite(unique_ID, resources, replicate, remove = FALSE) %>%
  unite(indiv_ID, resources, replicate, timeStep, remove=FALSE) %>%
  mutate(network_size = n_phagesPhen*n_bacteriaPhen) %>%
  filter(!is.na(nestedness))

# generalized linear mixed model
connectance_glmer <- glmer(connectance ~ timeStep*resources
                              + (1|replicate:resources)
                              + (1| replicate:resources:timeStep),
                              data=data, family="binomial", weights = network_size)
summary(connectance_glmer)

###############################################################################
# test overdispersion of the response variable for this binomial distribution #
###############################################################################
overdisp_fun(connectance_glmer)


# compare our glmer model with a model without the interaction term
anova(update(connectance_glmer, .~. -timeStep:resources), connectance_glmer, test = "LR")
#drop1(connectance_glmer, test="Chisq")
# X^2,(df=1) = 0.23, p = 0.632
# compare our glmer model without the interaction term with our glmer model without the interaction term plus resources
anova(update(connectance_glmer, .~. -timeStep:resources - resources), update(connectance_glmer, .~. -timeStep:resources), test = "LR")
# X^2,(df=1) = 1.31, p = 0.253
# compare our glmer model without the interaction term with our glmer model without the interaction term plus timeStep
anova(update(connectance_glmer, .~. -timeStep:resources - timeStep), update(connectance_glmer, .~. -timeStep:resources), test = "LR")
# X^2,(df=1) = 0.02, p = 0.879
#drop1(update(connectance_glmer, .~. -timeStep:resources), test="Chisq")
