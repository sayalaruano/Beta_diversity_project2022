setwd("C:/Users/antonella.bernardi/OneDrive - Universidad de Las Américas/Francisco Cuesta/paper")

# variables explicativas parcelas
veg<-read.csv("20212704_MATRIZ_Thesis_final.csv")
head(veg)

### Verificar correlaciones entre variables explicativas, paso inicial para selección de estas variables

library(Hmisc)

test.data.soil = veg[, c(14, 15, 17, 20)]

rcorr(as.matrix(test.data.soil))



test.data.product = veg[, c(4, 5, 6, 7, 8, 9, 10, 21, 22, 23, 25)]

rcorr(as.matrix(test.data.product))



plot_data_fin = read.csv("AlphFD_Indexvalues_270122.csv", sep = ",")
head (plot_data_fin)
test.data.clima = veg[, c(33, 34, 38, 39, 40, 41, 42, 43, 44, 45)]

rcorr(as.matrix(test.data.clima))
cor.test(veg$PCA1_temp, veg$PC1.p)

# regresion lineales

# Modelo 1: Para funcional eveness 
# Primer paso revisar relación entre Feve y cada variable explicativa
FEveLM = lm(FEve ~ AGC, data = veg) 
summary (FEveLM)
coef(FEveLM)

FEveLM = lm(FEve ~ Cnratio, data = veg) 
summary (FEveLM)
coef(FEveLM)

FEveLM = lm(FEve ~ Pmgkg, data = veg) 
summary (FEveLM)

FEveLM = lm(FEve ~ PCA1_temp, data = veg) 
summary (FEveLM)

FEveLM = lm(FEve ~ minDHr, data = veg) 
summary (FEveLM)


# Revisar ajuste de modelo completo 
FEveLM = lm(FEve ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC, data = veg) 
summary (FEveLM)
coef(FEveLM)
car::vif(FEveLM)


residual = residuals(FEveLM)
preds = predict(FEveLM)


plot(residual ~  preds, xlab = "Predicted Values",  ylab = "Residuals")

library(nortest)
pearson.test(residual)
lillie.test(residual)


par(mfrow = c(2, 2))
plot(FEveLM)

op = par(no.readonly = TRUE)

par(op)

lmtest::bptest (FEveLM)

search <- step(FEveLM, ~.^2)

#Realizar seleccion de mejor modelo, stepwise, backwise y ambos con funcion stepAIC
library(MASS)


stepAIC(FEveLM, direction = c("both", "backward", "forward"))

#Step:  AIC=-67.8
#FEve ~ Cnratio + PCA1_temp + AGC

FEveLM = glm(FEve ~ Cnratio + PCA1_temp + AGC, data = veg, family = gaussian()) 
summary (FEveLM)

residual = residuals(FEveLM)
preds = predict(FEveLM)


plot(residual ~  preds, xlab = "Predicted Values",  ylab = "Residuals")

pearson.test(residual)
lillie.test(residual)


par(mfrow = c(2, 2))
plot(FEveLM)

op = par(no.readonly = TRUE)

par(op)

## Problemas con la normalidad

#comprobar modelo nulo 

FEveLMnull<-lm(FEve ~ 1, data = veg)
summary (FEveLMnull)
stepAIC(FEveLMnull)
# AIC = -57.76

### #Realizar seleccion de mejor modelo, stepwise, backwise y ambos con funcion glmulti

library(glmulti)

glmulti.FEve = glmulti(FEve ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC , data = veg, level = 1, method = "h", crit = "aic", confsetsize = 100, plotty = F, report = T, fitfunction = "glm")

# Best model:  FEve~1+Cnratio+PCA1_temp+AGC
# Crit= -25.7526884617164
# Mean crit= -21.7944654829656

coef(glmulti.FEve)
#### http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti

print(glmulti.FEve)

#Obtener los mejores modelos (aic + 2)
tmp = weightable(glmulti.FEve)
tmp = tmp[tmp$aic <= min(tmp$aic) + 2,]
tmp

#examinar al "mejor" modelo, ahi se puede ver beta, SE, y p
summary(glmulti.FEve@objects[[1]])






# Modelo 2:Para funcional divergence 
# Primer paso revisar relación entre Fdiv y cada variable explicativa
FDivLM = lm(FDiv ~ AGC, data = veg) 
summary (FDivLM)
coef(FDivLM)

FDivLM = lm(FDiv ~ Cnratio, data = veg) 
summary (FDivLM)
coef(FDivLM)

FDivLM = lm(FDiv ~ Pmgkg, data = veg) 
summary (FDivLM)
coef(FDivLM)

FDivLM = lm(FDiv ~ PCA1_temp, data = veg) 
summary (FDivLM)
coef(FDivLM)

FDivLM = lm(FDiv ~ minDHr, data = veg) 
summary (FDivLM)
coef(FDivLM)

# Revisar ajuste de modelo completo 

FDivLM = glm(FDiv ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC, data = veg) 
summary (FDivLM)
coef(FDivLM)
car::vif(FDivLM)


residual = residuals(FDivLM)
preds = predict(FDivLM)


plot(residual ~  preds, xlab = "Predicted Values",  ylab = "Residuals")

pearson.test(residual)
lillie.test(residual)


par(mfrow = c(2, 2))
plot(FDivLM)

op = par(no.readonly = TRUE)

par(op)

#Realizar seleccion de mejor modelo, stepwise, backwise y ambos con funcion stepAIC

stepAIC(FDivLM, direction = c("both", "backward", "forward"))

#Step:  AIC=-21.96
#FDiv ~ PCA1_temp + minDHr + AGC

FDivLM = lm(FDiv ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC, data = veg) 
stepAIC(FDivLM, direction = c("both", "backward", "forward"))

FDivLM = lm(FDiv ~ PCA1_temp + minDHr + AGC, data = veg) 
summary (FDivLM)

residual = residuals(FDivLM)
preds = predict(FDivLM)


plot(residual ~  preds, xlab = "Predicted Values",  ylab = "Residuals")

pearson.test(residual)
lillie.test(residual)


par(mfrow = c(2, 2))
plot(FDivLM)

op = par(no.readonly = TRUE)

par(op)

FDivLM = lm(FDiv ~ 1, data = veg) 
summary (FDivLM)

#comparar el modelo nulo 


FDivnull<-lm(FDiv ~ 1, data = veg)
summary (FDivnull)
stepAIC(FDivnull)
# AIC = -71.28

#Realizar seleccion de mejor modelo, stepwise, backwise y ambos con funcion glmulti

glmulti.FDiv = glmulti(FDiv ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC, data = veg, level = 1, method = "h", crit = "aic", confsetsize = 100, plotty = F, report = T, fitfunction = "lm")

# Best model:  FDiv~1+minDHr+AGC
# Crit= -24.0662732279094
# Mean crit= -21.0753470902433

coef(glmulti.FDiv)
#### http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti

print(glmulti.FDiv)

#Obtener los mejores modelos (aic + 2)
tmp = weightable(glmulti.FDiv)
tmp = tmp[tmp$aic <= min(tmp$aic) + 2,]
tmp

#examinar al "mejor" modelo, ahi se puede ver beta, SE, y p
summary(glmulti.FDiv@objects[[1]])




# Modelo 3:
FRicLM = lm(FRic ~ AGC, data = veg) 
summary (FRicLM)
coef(FRicLM)

FRicLM = lm(FRic ~ Cnratio, data = veg) 
summary (FRicLM)
coef(FRicLM)

FRicLM = lm(FRic ~ Pmgkg, data = veg) 
summary (FRicLM)
coef(FRicLM)

FRicLM = lm(FRic ~ PCA1_temp, data = veg) 
summary (FRicLM)
coef(FRicLM)

FRicLM = lm(FRic ~ minDHr, data = veg) 
summary (FRicLM)
coef(FRicLM)


FRicLM = lm(FRic ~ + PCA1_temp + Cnratio + AGC + minDHr + Pmgkg, data = veg) 
summary (FRicLM)
coef(FRicLM)
car::vif(FRicLM)


residual = residuals(FRicLM)
preds = predict(FRicLM)


plot(residual ~  preds, xlab = "Predicted Values",  ylab = "Residuals")

pearson.test(residual)
lillie.test(residual)


par(mfrow = c(2, 2))
plot(FRicLM)

op = par(no.readonly = TRUE)

par(op)

stepAIC(FRicLM, direction = c("both", "backward", "forward"))

#Step:  AIC=33.59
#FRic ~ PCA1_temp + Cnratio + minDHr + Pmgkg

FRicLM = lm(FRic ~ PCA1_temp + Cnratio + minDHr + Pmgkg, data = veg) 
summary (FRicLM)

residual = residuals(FRicLM)
preds = predict(FRicLM)


plot(residual ~  preds, xlab = "Predicted Values",  ylab = "Residuals")

pearson.test(residual)
lillie.test(residual)


par(mfrow = c(2, 2))
plot(FRicLM)

op = par(no.readonly = TRUE)

par(op)

FRicLM = lm(FRic ~ Cnratio + Pmgkg + minDHr , data = veg) 
summary (FRicLM)

#
#comprobar modelo nulo 

FRicLMnull<-lm(FRic ~ 1, data = veg)
stepAIC(FRicLMnull)
# AIC = 55.84


glmulti.FRic = glmulti(FRic ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC, data = veg, level = 1, method = "h", crit = "aic", confsetsize = 100, plotty = F, report = T, fitfunction = "lm")

# Best model: FRic~1+PCA1_temp+AGC
# Crit= 80.996788001781
# Mean crit= 92.2042259421701

coef(glmulti.FRic)

print(glmulti.FRic)

#Obtener los mejores modelos (aic + 2)
tmp = weightable(glmulti.FRic)
tmp = tmp[tmp$aic <= min(tmp$aic) + 2,]
tmp

#examinar al "mejor" modelo, ahi se puede ver beta, SE, y p
summary(glmulti.FRic@objects[[1]])



# Modelo 4:
FDisLM = lm(FDis ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC, data = veg) 
summary (FDisLM)
coef(FDisLM)
car::vif(FDisLM)


residual = residuals(FDisLM)
preds = predict(FDisLM)


plot(residual ~  preds, xlab = "Predicted Values",  ylab = "Residuals")

pearson.test(residual)
lillie.test(residual)


par(mfrow = c(2, 2))
plot(FDisLM)

op = par(no.readonly = TRUE)

par(op)

lmtest::bptest (FDisLM)

search <- step(FDisLM, ~.^2)


stepAIC(FDisLM, direction = c("both", "backward", "forward"))

#Step:  AIC=-73.16
#FEve ~ Pmgkg + minDHr

FDisLM = lm(FDis ~ PCA1_temp, data = veg) 
summary (FDisLM)

residual = residuals(FDisLM)
preds = predict(FDisLM)


plot(residual ~  preds, xlab = "Predicted Values",  ylab = "Residuals")

pearson.test(residual)
lillie.test(residual)


par(mfrow = c(2, 2))
plot(FDisLM)

op = par(no.readonly = TRUE)

par(op)

FDisLM = lm(FDis ~ PCA1_temp + AGC, data = veg) 
summary (FDisLM)

#comprobar modelo nulo 

FDisLMnull<-lm(FDis ~ 1, data = veg)

stepAIC(FDisLMnull)
# AIC = -70.85

### glmulti
glmulti.FDis = glmulti(FDis ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC, data = veg, level = 1, method = "h", crit = "aic", confsetsize = 100, plotty = F, report = T, fitfunction = "lm")

# Best model:  FEve~1+Pmgkg+minDHr
# Crit= -25.7526884617164
# Mean crit= -21.7944654829656

coef(glmulti.FDis)
#### http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti

print(glmulti.FDis)

#Obtener los mejores modelos (aic + 2)
tmp = weightable(glmulti.FDis)
tmp = tmp[tmp$aic <= min(tmp$aic) + 2,]
tmp

#examinar al "mejor" modelo, ahi se puede ver beta, SE, y p
summary(glmulti.FDis@objects[[1]])




## con taxonomicos

#Hilldiversity

HillDivLM = lm(Hill_Diversity ~ PCA1_temp + AGC + Cnratio + minDHr + Pmgkg, data = veg) 
summary (HillDivLM)
coef(HillDivLM)
car::vif(HillDivLM)


residual = residuals(HillDivLM)
preds = predict(HillDivLM)


plot(residual ~  preds, xlab = "Predicted Values",  ylab = "Residuals")

pearson.test(residual)
lillie.test(residual)


par(mfrow = c(2, 2))
plot(HillDivLM)

op = par(no.readonly = TRUE)

par(op)
lmtest::bptest (HillDivLM)

search <- step(HillDivLM, ~.^2)


#Stepwise regression
HillDivLM = lm(Hill_Diversity ~ PCA1_temp + AGC + Cnratio + minDHr + Pmgkg, data = veg)
HillDivLM = stepAIC(HillDivLM, direction = c("both", "backward", "forward"))


HillDivLM = lm(Hill_Diversity ~ PCA1_temp + AGC + Cnratio + Pmgkg, data = veg) 
summary (HillDivLM)

residual = residuals(HillDivLM)
preds = predict(HillDivLM)


plot(residual ~  preds, xlab = "Predicted Values",  ylab = "Residuals")

pearson.test(residual)
lillie.test(residual)


par(mfrow = c(2, 2))
plot(HillDivLM)

op = par(no.readonly = TRUE)

par(op)

HillDivLM = lm(Hill_Diversity ~ PCA1_temp + AGC + Cnratio + Pmgkg, data = veg) 
summary (HillDivLM)


#comprobar modelo nulo 

Hillnull<-lm(Hill_Diversity ~ 1, data = veg)
summary (Hillnull)
stepAIC(Hillnull)


# AIC = 89.33


glmulti.HillDivLM = glmulti(Hill_Diversity ~ PCA1_temp + AGC + Cnratio + minDHr + Pmgkg, 
                            data = veg, level = 1, method = "h", crit = "aic", confsetsize = 100, plotty = F, report = T, fitfunction = "lm")


coef(glmulti.HillDivLM)

#Obtener los mejores modelos (aic + 2)
tmp = weightable(glmulti.HillDivLM)
tmp = tmp[tmp$aic <= min(tmp$aic) + 2,]
tmp

#examinar al "mejor" modelo, ahi se puede ver beta, SE, y p
summary(glmulti.HillDivLM@objects[[1]])



#### FOURTH CORNER


library(ade4)

library(mvabund)
library(Rcpp)

library(lattice)
library(gllvm)

######## loading data ######

## table with traits of species, for species without data these were completed usin gthe genus level
traits<-read.csv("20211108_MeanTraits_FourthCorner.csv",stringsAsFactors = TRUE)
head(traits1)
rownames(traits)<-traits$Scientific_name
traits1 = traits[c("Ht2019","LA","SLA")]
colnames(traits1)<-c("Height","LA","SLA")

sp_summit<-read.csv("20211108_BA_FuncTraits_Repres sp.csv")
rownames(sp_summit)<-sp_summit$X
head(sp_summit1)
del1<-c(sp_summit,"X")
sp_summit1<-sp_summit[,-which(names(sp_summit) %in% del1)]

env1<-read.csv("20211117_Environmental_FourthCorner.csv")
rownames(env1)<-env1$PLOT
head(env1)
env1$PLOT = NULL

#### fourth corner with ade4 ###

### model where permutation 
four<-fourthcorner(env1,sp_summit1,traits1, modeltype = 3, nrepet = 3000)
summary(four)


## figura
setEPS()
svg("fourth_corner.svg", width = 3.6, height = 3.6, pointsize = 8)
plot(four, stat = "G")
dev.off()

plot(four, stat = "D2")

### model 6 (default): combination of the outputs of models 2 and 4. Dray and Legendre (2008) and ter Braak et al. (20012) 
## showed that all models (except model 6) have inflated type I error. 
four1<-fourthcorner(env1,sp_summit1,traits1, modeltype = 1, nrepet = 10000)
summary(four1)
plot(four1, stat = "G")
plot(four1, stat = "D2")
