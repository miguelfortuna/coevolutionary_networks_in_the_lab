################################################################
### DOES BETADIVERSITY INCREASE WITH THE LEVEL OF RESOURCES? ###
################################################################
data <- read_csv('betadiversity.csv')
data <- data %>%
  mutate(resources = as.factor(resources),
         organism = as.factor(organism))

# total beta-diversity
lm_betadiversity <- lm(betadiversity ~ resources * organism,
                    data = data)
summary(lm_betadiversity)
anova(lm_betadiversity)

# fraction of the total beta-diversity explained by phenotypic turnover
lm_betadiversity_turnover <- lm(betadiversity_turnover/betadiversity ~ resources * organism,
                       data = data)
summary(lm_betadiversity_turnover)
anova(lm_betadiversity_turnover)