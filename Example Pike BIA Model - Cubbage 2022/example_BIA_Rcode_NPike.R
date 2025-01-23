######Pike BIA Example 

######Clear environment
rm(list=ls(all=TRUE))

######Set working directory to location of Pike_Example.csv and model.rds files
#setwd("/Users/tlcubbage/OneDrive - University of Alaska/Documents/UAF1/Thesis materials/pike physiology/BIA/BIA_supplement")

#####Import Pike BIA data
BIA_data <- data.frame(read.csv("Example Pike BIA Model - Cubbage 2022/Pike_Example_BIA.csv"))
str(BIA_data)
attach(BIA_data)

######Temperature corrections
lrest = DML_res* ((-8.07*Vent_temp1+438.31)/(-8.07*10+438.31))
lreat = DML_rea* ((-1.77*Vent_temp1+129.57)/(-1.77*10+129.57))
vrest = DV_res* ((-5.67*Vent_temp2+208.70)/(-5.67*10+208.70))
vreat = DV_rea* ((-2.41*Vent_temp2+85.27)/(-2.41*10+85.27))

######Lateral surface
Rp = lrest+(lreat**2/lrest)
Xcp = lreat+(lrest**2/lreat)
Cpf = 3.1831E-18/Xcp
Z = (lrest**2+lreat**2)**0.5
Z2 = lrest*lreat/(lrest**2+lreat**2)**0.5

######Lateral covariates
LE1 = DML_Detector**2/lrest                            #Resistance in series (Rs)
LE2 = DML_Detector**2/Rp                               #Resistance in parallel (Rp)
LE3 = DML_Detector**2/lreat                            #Reactance in series (Xc)
LE4 = DML_Detector**2/Xcp                              #Reactance in parallel (Xcp)
LE5 = DML_Detector**2/Cpf                              #Capacitance (Cpf)
LE6 = DML_Detector**2/Z                                #Impedance in series (Zs)
LE7 = atan(lreat/lrest)*180/pi                         #Phase angle
LE8 = DML_Detector*LE7                                 #Standardized phase angle
LE9 = DML_Detector**2/Z2                               #Impedance in parallel


######Ventral surface
VRp = vrest+(vreat**2/vrest)
VXcp = vreat+(vrest**2/vreat)
VCpf = 3.1831E-18/VXcp
VZ = (vrest**2+vreat**2)**0.5
VZ2 = vrest*vreat/(vrest**2+vreat**2)**0.5

######Ventral covariates
VE1 = DV_Detector**2/vrest
VE2 = DV_Detector**2/VRp
VE3 = DV_Detector**2/vreat
VE4 = DV_Detector**2/VXcp
VE5 = DV_Detector**2/VCpf
VE6 = DV_Detector**2/VZ
VE7 = atan(vreat/vrest)*180/pi
VE8 = DV_Detector*VE7
VE9 = DV_Detector**2/VZ2

######Body condition index (Bentley & Schindler 2013)
BCI <- rstandard(glm(log10(Fork_Length)~log10(whole_weight)))                 

######Predict percent dry mass from top multiple linear regression model (Table 4 of manuscript) 
#Reading in model files
drymass_model <- readRDS("Example Pike BIA Model - Cubbage 2022/Pike_BIA_drymass_model.rds")
lipid_model <- readRDS("Example Pike BIA Model - Cubbage 2022/Pike_BIA_lipid_model.rds")
#Combining predictors
lengths <- sapply(list(LE1, LE2, LE3, LE4, LE5, LE6, LE7, LE8, LE9, 
                       VE1, VE2, VE3, VE4, VE5, VE6, VE7, VE8, VE9, 
                       BCI, whole_weight, Fork_Length), length)
print(lengths)

drymass_predictors <- data.frame(cbind(LE1,LE2,LE3,LE4,LE5,LE6,LE7,LE8,LE9,VE1,VE2,VE3,VE4,VE5,VE6,VE7,VE8,VE9,BCI,whole_weight))
lipid_predictors <- data.frame(cbind(LE1,LE2,LE3,LE4,LE5,LE6,LE7,LE8,LE9,VE1,VE2,VE3,VE4,VE5,VE6,VE7,VE8,VE9,BCI,whole_weight, Fork_Length))


#Creating data frame with predicted dry mass and lipid values
pred.df <- data.frame(Pike_ID,Population,Fork_Length,whole_weight,
                      predDM = predict(drymass_model, newdata = drymass_predictors),
                      predDL = predict(lipid_model, newdata = lipid_predictors))
                     


######Plot predicted dry mass and dry lipid ~ fork length
plot(pred.df$Fork_Length, pred.df$predDM,ylab="Dry Mass (%)",xlab = "Fork Length (mm)")
plot(pred.df$Fork_Length, pred.df$predDL,ylab="Dry Lipid (%)",xlab = "Fork Length (mm)")

######Plot histogram of predicted dry mass and dry lipid
hist(pred.df$predDM,xlab="Dry Mass (%)",main=NULL)
hist(pred.df$predDL,xlab="Dry Lipid (%)",main=NULL)
boxplot(pred.df$predDM ~ pred.df$Population, xlab = "Population")
boxplot(pred.df$predDL ~ pred.df$Population, xlab = "Population")


