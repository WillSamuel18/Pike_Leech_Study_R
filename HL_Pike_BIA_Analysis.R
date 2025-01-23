################################################################################
#### Analysis of Northern Pike Body condition compared to leech infestation ####
######################## Project PI: Dr. Kevin Fraley ##########################
################ Primary author of this script: Will Samuel ####################
## Other project contributors: Dr. Morag Clinton, Taylor Cubbage, Joe Spencer ##
################################################################################

library(tidyverse)
library(ggplot2)



# Data manipulation --------------------------------------------------------
cap.data <- read.csv("Data/HL Capture data formatted.csv")
BIA.data <- read.csv("Data/HL Pike BIA data formatted.csv")
str(cap.data)
str(BIA.data)

#Remove burbot from the dataset
cap.data <- cap.data %>%   filter(!grepl("B", Fish_ID))
BIA.data <- BIA.data %>%   filter(!grepl("B", Fish_ID))


#Assign numerical
cap.data <- cap.data %>%
  mutate(across(c(Fork_length_mm), ~ as.numeric(as.integer(.)))) %>% 
  mutate(across(c(Girth_mm, Muscle_sample_weight_kg), ~ as.numeric(as.integer(.))))


BIA.data <- BIA.data %>%
  mutate(across(c(DML_res, DML_rea, DML_length, DV_res, DV_rea), 
                ~ as.numeric(as.integer(.)))) %>% 
 mutate(across(c(DML_length, DV_length, Vent_temp), 
                ~ as.numeric(as.character(.))))


#Add FL back to the BIA datasheet
FL.dat <- cap.data %>% select(Sample_num, Fork_length_mm, Weight_kg)

BIA.data <- merge(BIA.data, FL.dat, by = "Sample_num")
BIA.data <- BIA.data %>%
  mutate(Weight_kg = as.numeric(Weight_kg))  # Ensure it's numeric

#Remove the fish that don't have BIA data
BIA.data <- BIA.data %>%  
  filter(!is.na(DML_res)) %>% 
  filter(!is.na(Weight_kg)) #%>% #Have to filter out fish with incomplete data
  #mutate("num" = seq(1,n()))#, 
         #Population = "Native") #These help make it match up with Taylor's data and models




#An instrument error prevented us from collecting vent temp measurements in 2023
#Although, there was almost no variability in the vent temp measured in 2024
mean(BIA.data$Vent_temp, na.rm = T) #Mean = 1.08 C
sd(BIA.data$Vent_temp, na.rm = T) #SD = 0.08 C
#Since this variability is smaller than the measurement error of the thermometer,
#we feel it is appropriate to assume that the vent temp was 1.08 across  all fish 
#at both sampling events
BIA.data <- BIA.data %>% 
  mutate(Vent_temp = ifelse(is.na(Vent_temp), 1.08, Vent_temp),
         "whole_weight" = Weight_kg*1000) #need to conduct the BIA model in g not kg






# Run BIA model estimates -------------------------------------------------
attach(BIA.data)

######Temperature corrections
lrest = DML_res* ((-8.07*Vent_temp+438.31)/(-8.07*10+438.31))
lreat = DML_rea* ((-1.77*Vent_temp+129.57)/(-1.77*10+129.57))
vrest = DV_res* ((-5.67*Vent_temp+208.70)/(-5.67*10+208.70))
vreat = DV_rea* ((-2.41*Vent_temp+85.27)/(-2.41*10+85.27))

######Lateral surface
Rp = lrest+(lreat**2/lrest)
Xcp = lreat+(lrest**2/lreat)
Cpf = 3.1831E-18/Xcp
Z = (lrest**2+lreat**2)**0.5
Z2 = lrest*lreat/(lrest**2+lreat**2)**0.5

######Lateral covariates
LE1 = DML_length**2/lrest                            #Resistance in series (Rs)
LE2 = DML_length**2/Rp                               #Resistance in parallel (Rp)
LE3 = DML_length**2/lreat                            #Reactance in series (Xc)
LE4 = DML_length**2/Xcp                              #Reactance in parallel (Xcp)
LE5 = DML_length**2/Cpf                              #Capacitance (Cpf)
LE6 = DML_length**2/Z                                #Impedance in series (Zs)
LE7 = atan(lreat/lrest)*180/pi                         #Phase angle
LE8 = DML_length*LE7                                 #Standardized phase angle
LE9 = DML_length**2/Z2                               #Impedance in parallel


######Ventral surface
VRp = vrest+(vreat**2/vrest)
VXcp = vreat+(vrest**2/vreat)
VCpf = 3.1831E-18/VXcp
VZ = (vrest**2+vreat**2)**0.5
VZ2 = vrest*vreat/(vrest**2+vreat**2)**0.5

######Ventral covariates
VE1 = DV_length**2/vrest
VE2 = DV_length**2/VRp
VE3 = DV_length**2/vreat
VE4 = DV_length**2/VXcp
VE5 = DV_length**2/VCpf
VE6 = DV_length**2/VZ
VE7 = atan(vreat/vrest)*180/pi
VE8 = DV_length*VE7
VE9 = DV_length**2/VZ2

######Body condition index (Bentley & Schindler 2013)
BCI <- rstandard(glm(log10(Fork_length_mm)~log10(whole_weight)))                 

######Predict percent dry mass from top multiple linear regression model (Table 4 of manuscript) 
#Reading in model files
drymass_model <- readRDS("Example Pike BIA Model - Cubbage 2022/Pike_BIA_drymass_model.rds")
lipid_model <- readRDS("Example Pike BIA Model - Cubbage 2022/Pike_BIA_lipid_model.rds")
summary(drymass_model)
summary(lipid_model)

#Combining predictors
#lengths <- sapply(list(LE1, LE2, LE3, LE4, LE5, LE6, LE7, LE8, LE9, 
#                       VE1, VE2, VE3, VE4, VE5, VE6, VE7, VE8, VE9, 
#                       BCI, Weight_kg, Fork_length_mm), length)
#print(lengths)

BIA.data <- BIA.data %>% rename(Fork_Length = Fork_length_mm)

#lengths <- sapply(list(LE1, LE2, LE3, LE4, LE5, LE6, LE7, LE8, LE9, 
#                       VE1, VE2, VE3, VE4, VE5, VE6, VE7, VE8, VE9, 
#                       BCI, BIA.data$whole_weight, BIA.data$Fork_Length), length)
#print(lengths)

drymass_predictors <- data.frame(cbind(LE1,LE2,LE3,LE4,LE5,LE6,LE7,LE8,LE9,VE1,VE2,
                                       VE3,VE4,VE5,VE6,VE7,VE8,VE9,
                                       BCI, "whole_weight" = BIA.data$whole_weight)) #need to adjust from kg to g for the model
lipid_predictors <- data.frame(cbind(LE1,LE2,LE3,LE4,LE5,LE6,LE7,LE8,LE9,VE1,VE2,
                                     VE3,VE4,VE5,VE6,VE7,VE8,VE9,
                                     BCI, "whole_weight" = BIA.data$whole_weight, 
                                     "Fork_Length" = BIA.data$Fork_Length))



#Creating data frame with predicted dry mass and lipid values
pred.df <- data.frame("Sample_num" = BIA.data$Sample_num, "X" = BIA.data$X, "Fork_Length" = BIA.data$Fork_Length, 
                      "whole_weight" =  BIA.data$whole_weight,
                      predDM = predict(drymass_model, newdata = drymass_predictors),
                      predDL = predict(lipid_model, newdata = lipid_predictors))
#View(pred.df)
detach(BIA.data)



######Plot predicted dry mass and dry lipid ~ fork length
plot(pred.df$Fork_Length, pred.df$predDM,ylab="Dry Mass (%)",xlab = "Fork Length (mm)")
plot(pred.df$Fork_Length, pred.df$predDL,ylab="Dry Lipid (%)",xlab = "Fork Length (mm)")

######Plot histogram of predicted dry mass and dry lipid
hist(pred.df$predDM,xlab="Dry Mass (%)",main=NULL)
hist(pred.df$predDL,xlab="Dry Lipid (%)",main=NULL)


pred.df2 <- pred.df %>%  select(-Fork_Length, -whole_weight, Sample_num)

newBIA.data <- merge(BIA.data, pred.df2, by = "X")
newBIA.data <- newBIA.data %>% rename(Fork_length_mm = Fork_Length)

write.csv(newBIA.data, "Data/BIA.data_with_predictions.csv")



sum.pred.df <- pred.df %>%  
  select(-X) %>% 
  group_by(Sample_num) %>%  
  summarize(predDM = mean(predDM),    
            predDL = mean(predDL))

newcap.data <- merge(cap.data, sum.pred.df, by = "Sample_num")

write.csv(newcap.data, "Data/HL_Capture_data_formatted_with_predictions.csv")



#A few more model checks for the BIA results...
ggplot(newBIA.data, aes(whole_weight, predDM))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_text(aes(label = Date), hjust = 0.5, vjust = -0.5, size = 3) +  
  theme_bw()

ggplot(newBIA.data, aes(whole_weight, predDL))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_text(aes(label = Date), hjust = 0.5, vjust = -0.5, size = 3) +  
  theme_bw()


# Comparisons between body condition and leech load.  ---------------------


newcap.data <- read.csv("Data/HL_Capture_data_formatted_with_predictions.csv")
str(newcap.data)


#Calculate fultons condition factor
newcap.data <- newcap.data %>% 
  mutate("Ful_cond_fact" = (100*((Weight_kg*1000)/((Fork_length_mm/10)^3))), #calculating in mm and g. 
         #"BCI" = Weight_kg/Fork_length_mm,
         ######Body condition index (Bentley & Schindler 2013)
         "BCI" = rstandard(glm(log10(Fork_length_mm)~log10(Weight_kg*1000))), 
         "Leech_load" = Leech_count/Fork_length_mm) #standardizing to body length 


ggplot(newcap.data, aes(Ful_cond_fact, Leech_load))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

ggplot(newcap.data, aes(Ful_cond_fact, Leech_count))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

ggplot(newcap.data, aes(BCI, Leech_load))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

ggplot(newcap.data, aes(BCI, Leech_count))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()


ggplot(newcap.data, aes(Ful_cond_fact, Mercury))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

ggplot(newcap.data, aes(BCI, Mercury))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()


ggplot(newcap.data, aes(Ful_cond_fact, Age))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

ggplot(newcap.data, aes(BCI, Age))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()


ggplot(newcap.data, aes(Mercury, Leech_load))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

ggplot(newcap.data, aes(Mercury, Leech_count))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()


ggplot(newcap.data, aes(Sex, Leech_count))+
  geom_boxplot()+
  geom_violin()+
  geom_jitter()+
  theme_bw()

Female <- newcap.data %>% subset(Sex == "F")
Male <- newcap.data %>% subset(Sex == "M")
t.test(Female$Leech_count, Male$Leech_count) #not a significant difference, maybe with more data

ggplot(newcap.data, aes(Sex, Mercury))+
  geom_violin()+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()

t.test(Female$Mercury, Male$Mercury) 


ggplot(newcap.data, aes(Age, Mercury))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()


ggplot(newcap.data, aes(Fork_length_mm, Mercury))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

summary(lm(Mercury ~ Fork_length_mm, newcap.data)) #not a significant difference, maybe with more data


ggplot(newcap.data, aes(Age, Leech_load))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()














  