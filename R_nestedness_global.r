##############################################################################
### DOES NESTEDNESS (GLOBAL NETWORK) INCREASE WITH THE LEVEL OF RESOURCES? ###
##############################################################################
data <- read_csv('nestedness_global.csv')
data$resources <- as.factor(data$resources)
data$resources <- factor(data$resources, levels=c("low", "high"))
data$resources <- relevel(data$resources, ref = "low")

lm_nestedness <- lm(nestedness_network ~ resources * connectance,
                  data = data)
plot(lm_nestedness)
# Q: are the residuals normally distributed?
hist(residuals(lm_nestedness))
shapiro.test(residuals(lm_nestedness)) # normality test
# A: yes they are, because p-value > 0.05
# Q: is the residual variation homogenous?
boxplot(residuals(lm_nestedness) ~ data$resources)
leveneTest(y = residuals(lm_nestedness), group = data$resources) # homogeneity of variance test
# A: yes it is, because p-value > 0.05
summary(lm_nestedness)
# F-statistic: 3.882 on 3 and 8 DF,  p-value: 0.0555
sem.model.fits(lm_nestedness)
# the model explains 59% of the variance
anova(lm_nestedness)
# F-statistic of the interaction term (resources*connectance): 10.8863 (p=0.011)
effect(c("resources*connectance"), lm_nestedness)
plot((effect(c("resources*connectance"), lm_nestedness, se=TRUE)), "connectance")
ggplot(data, aes(x=connectance,y=nestedness_network,color=resources))+geom_point() + geom_smooth(method="lm", level=0.95)