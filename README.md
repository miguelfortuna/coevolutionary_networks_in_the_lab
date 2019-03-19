# <a href="https://github.com/miguelfortuna/coevolutionary_networks_in_the_lab/blob/master/Coevolutionary%20dynamics%20shape%20the%20structure%20of%20bacteria-phage%20infection%20networks.pdf" target="_blank">Coevolutionary dynamics shape the structure of bacteria-phage infection networks.</a>

Dataset and R code used in the following study: "Coevolutionary dynamics shape the structure of bacteria-phage infection networks", that will be published in Evolution (2019).

# Data were obtained from the file "database.csv"
## files provided
- database.csv
- phenotypic_diversification.csv
- betadiversity.csv
- infectivity.csv
- nestedness_contemporary.csv
- nestedness_global.csv

# R scripts for the statistical analysis
## files provided
- R_phenotypic_diversification.r
- R_betadiversity.r
- R_infectivity.r
- R_connectance.r
- R_nestedness_contemporary.r
- R_nestedness_global.r

# Statistical analysis.
## Phenotypic diversification
We used a linear mixed model to test the effect of resources on phenotypic diversification.
#### data
- phenotypic_diversification.csv
#### code
- R_phenotypic_diversification.r

## Beta-diversity
We used a linear model to analyze the effect of resources, type of organism, and their interaction on the total beta-diversity and on the fraction of the total beta-diversity explained by phenotypic turnover.
#### data
- betadiversity.csv
#### code
- R_betadiversity.r

## Phage infectivity to evolving and coevolving bacteria
The role of resources in explaining the probability for a phage to infect a coevolving bacterium compared to that of infecting a bacterium either from the past or from the future was analyzed using a generalized linear mixed model. We modeled the probability of infection with a binomial distribution (link function=logit). We specified the statistical interaction between the type of interaction and the resource level as fixed effects, and we included replicate as a random effect.
#### data
- infectivity.csv
#### code
- R_infectivity.r
## Nestedness
Since 43% of the contemporary networks were perfectly nested (i.e., N = 1), we first tested the role of network size in explaining the prevalence of perfect nestedness by using a generalized linear mixed model (binomial distribution; link function=logit). We specified network size and the interaction between time and resources as fixed effects, and included replicate as a random effect.
#### data
- nestedness_contemporary.csv
#### code
- R_nestedness_contemporary.r

Next, we explored changes in connectance over time for each resource level by using a generalized linear mixed model (binomial distribution; link function=logit). We specified resources and the interaction between time and resources as fixed effects, and included replicate as a random effect.
#### data
- nestedness_contemporary.csv
#### code
- R_connectance.r

After that, we focused on contemporary networks that were large enough to allow nestedness to vary (i.e., N < 1). We then used a linear mixed model to analyze the effect of the resource level in determining changes in nestedness (logit-transformed) over time. We specified connectance, network size, and the interaction between time and resources as fixed effects, and included replicate as a random effect.
#### data
- nestedness_contemporary.csv
#### code
- R_nestedness_contemporary.r

Finally, we used a linear model to analyze the effect of connectance, resources, and their interaction, on the nestedness of the global network.
#### data
- nestedness_global.csv
#### code
- R_nestedness_global.r
