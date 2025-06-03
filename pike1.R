
library(ggplot2)
library(car)
library(mgcv)

#setwd("~/WCS/Manuscripts/Pike leech")

#healy.data <- read.csv("Data/pike.csv")
healy.data <- read.csv("Data/pike2.csv")


##Lipid content##
anova <- aov(leech_cm ~ lipid_level2, data = healy.data)
summary(anova)

tukey <- TukeyHSD(anova)
print(tukey)

p <- ggplot(healy.data, aes(lipid_level2, leech_cm))
p + geom_violin(alpha = 0.5, fill='burlywood', color="darkolivegreen")+
  geom_boxplot(width = 0.1, fill='burlywood', color="darkolivegreen")+
  ylab("Leeches Per cm Pike Length")+
  xlab("Percent Lipid Content")+
  theme_classic()+
  scale_y_continuous(limits=c(0,10),breaks=seq(0, 10, by = 1))+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1))

#log transform#
#levene test
result = leveneTest(leech_cm ~ lipid_level2, healy.data2)
print(result)

result = leveneTest(leech_cm ~ merc_level, healy.data2)
print(result)


anova3 <- aov(leech_cm ~ lipid_level2, data = healy.data2)
summary(anova3)

tukey <- TukeyHSD(anova3)
print(tukey)

p <- ggplot(healy.data2, aes(lipid_level2, leech_cm))
p + geom_violin(alpha = 0.5, fill='burlywood', color="darkolivegreen")+
  geom_boxplot(width = 0.1, fill='burlywood', color="darkolivegreen")+
  ylab("Log of Leeches Per cm Pike Length")+
  xlab("Percent Lipid Content")+
  theme_classic()+
  scale_y_continuous(limits=c(-1,1),breaks=seq(-1, 1, by = 0.1))+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1))





#
anova2 <- aov(leech_cm ~ lipid_level, data = healy.data)
summary(anova2)

tukey2 <- TukeyHSD(anova2)
print(tukey2)

ggplot(healy.data, aes(lipid_level, leech_cm)) + 
  geom_boxplot()

##Mercury##
anova <- aov(leech_cm ~ merc_level, data = healy.data2)
summary(anova)

tukey <- TukeyHSD(anova)
print(tukey)

ggplot(healy.data, aes(merc_level, leech_cm)) + 
  geom_boxplot()

p <- ggplot(healy.data2, aes(merc_level, leech_cm))
p + geom_violin(alpha = 0.5, fill='burlywood', color="darkolivegreen")+
  geom_boxplot(width = 0.1, fill='burlywood', color="darkolivegreen")+
  ylab("Log of Leeches Per cm Pike Length")+
  xlab("Mercury Concentration Fish Consumption Advisory Threshold")+
  theme_classic()+
  scale_y_continuous(limits=c(-1,1),breaks=seq(-1, 1, by = 0.1))+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1))

###GAM###
##Lipid##
#non-log data
gam_model <- gam(lipid_content ~ s(leech_cm), data = healy.data)
summary(gam_model)

#log data
gam_model2 <- gam(lipid_content ~ s(leech_cm), data = healy.data2)
summary(gam_model2)

##Mercury##
#non-log data
gam_model3 <- gam(mercury ~ s(leech_cm), data = healy.data)
summary(gam_model3)

#log data
gam_model4 <- gam(mercury ~ s(leech_cm), data = healy.data2)
summary(gam_model4)
