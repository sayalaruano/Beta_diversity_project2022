# Script to run linear models - GLMs

# Import libraries
library(Hmisc)
library(nortest)
library(MASS)
library(glmulti) # this package requires rJava

# Load data
df_var <- read.csv("Data/20212704_MATRIZ_Thesis_final.csv")
df_fd_alphadiv <- read.csv("Data/AlphFD_Indexvalues_270122.csv")

#######################################################
################ CORRELATION ANALYSIS #################
#######################################################

# Verify correlations between variables to select the best ones
data_soil <- df_var[, c(14, 15, 17, 20)]
corr_data_soil <- rcorr(as.matrix(data_soil))

data_prod <- df_var[, c(4, 5, 6, 7, 8, 9, 10, 21, 22, 23, 25)]
corr_data_prod <- rcorr(as.matrix(data_prod))

data_climate = df_var[, c(33, 34, 38, 39, 40, 41, 42, 43, 44, 45)]
corr_data_climate <- rcorr(as.matrix(data_climate))

#######################################################
################ LINEAR REGRESSIONS  ##################
#######################################################

# Function for simple linear regression
simple_lr <- function(x, y){
  if(length(y)==1){
    
    lr <- lm(x ~ y, data = df_var)
  }
  x_test <- "FEve"
  y_test <- "Cnratio + Pmgkg + PCA1_temp + minDHr + AGC"
  formula_test <- paste(x_test, y_test, sep= " ~ ")
  FEveLM_all_test <- lm(as.formula(formula_test), data = df_var)
  summary(FEveLM_all_test)
  
}



#######################################################
############## MODEL 1: EVENNESS ######################
#######################################################

# Check relationships between Feve and other variables
FEveLM_AGC <- lm(FEve ~ AGC, data = df_var)
summary(FEveLM_AGC)

FEveLM_Cnratio <- lm(FEve ~ Cnratio, data = df_var) 
summary(FEveLM_Cnratio)

FEveLM_Pmgkg <- lm(FEve ~ Pmgkg, data = df_var) 
summary(FEveLM_Pmgkg)

FEveLM_PCA1_temp <- lm(FEve ~ PCA1_temp, data = df_var) 
summary(FEveLM_PCA1_temp)

FEveLM_minDHr <- lm(FEve ~ minDHr, data = df_var) 
summary(FEveLM_minDHr)

# Run the complete model
FEveLM_all <- lm(FEve ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC, data = df_var) 

# Check the adjustment of the complete model
summary(FEveLM_all)
coef(FEveLM_all)
car::vif(FEveLM_all)

# Plot results of the complete linear regression model
par(mfrow = c(2, 2))
plot(FEveLM_all)

# Obtain predicted values and residuals, and some tests
FEveLM_all_residual <- residuals(FEveLM_all)
FEveLM_all_preds <- predict(FEveLM_all)

pearson.test(FEveLM_all_residual)
lillie.test(FEveLM_all_residual)
lmtest::bptest(FEveLM_all)
search <- step(FEveLM_all, ~.^2)

# Plot of predicted values vs residuals 
par(mfrow = c(1, 1))
plot(FEveLM_all_residual ~  FEveLM_all_preds, xlab = "Predicted Values",  ylab = "Residuals")

# Select the best model: stepwise, backwise and both of them with stepAIC function
FEveLM_all_best_model <- stepAIC(FEveLM_all, direction = c("both", "backward", "forward"))

# Test null model of linear model
FEve_nullLM <- lm(FEve ~ 1, data = df_var)
summary(FEve_nullLM)
stepAIC(FEve_nullLM)

# Run the generalized linear model with FEve ~ Cnratio + PCA1_temp + AGC
FEveGLM_all <- glm(FEve ~ Cnratio + PCA1_temp + AGC, data = df_var, 
                  family = gaussian()) 
summary (FEveGLM_all)

# Plot results of the complete general linear regression model
par(mfrow = c(2, 2))
plot(FEveGLM_all)

# Obtain predicted values and residuals, and some tests
residuals_FEveglm <- residuals(FEveGLM_all)
preds_FEveglm <- predict(FEveGLM_all)

pearson.test(residuals_FEveglm)
lillie.test(residuals_FEveglm)

# Plot of predicted values vs residuals 
par(mfrow = c(1, 1))
plot(residuals_FEveglm ~ preds_FEveglm, xlab = "Predicted Values",  ylab = "Residuals")

#Select the best model: stepwise, backwise and both of them with glmulti function
FEve_glmulti <- glmulti(FEve ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC , 
                       data = df_var, level = 1, method = "h", crit = "aic", 
                       confsetsize = 100, plotty = F, report = T, 
                       fitfunction = "glm")

# Best model: FEve~1+Cnratio+PCA1_temp+AGC
# Crit= -25.7526884617164
# Mean crit= -21.7944654829656

# More info: http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti

# Obtain the best GML models (aic + 2)
best_models1_glm <- weightable(FEve_glmulti)
best_models1_glm <- best_models1_glm[best_models1_glm$aic <= min(best_models1_glm$aic) + 2,]

# Check the best model to obtain beta, SE, and p value
summary(FEve_glmulti@objects[[1]])

#######################################################
####### MODEL 2: FUNCTIONAL DIVERGENCE ################
#######################################################

# Check relationships between FDiv and other variables
FDivLM_AGC <- lm(FDiv ~ AGC, data = df_var)
summary(FDivLM_AGC)

FDivLM_Cnratio <- lm(FDiv ~ Cnratio, data = df_var) 
summary(FDivLM_Cnratio)

FDivLM_Pmgkg <- lm(FDiv ~ Pmgkg, data = df_var) 
summary(FDivLM_Pmgkg)

FDivLM_PCA1_temp <- lm(FDiv ~ PCA1_temp, data = df_var) 
summary(FDivLM_PCA1_temp)

FDivLM_minDHr <- lm(FDiv ~ minDHr, data = df_var) 
summary(FDivLM_minDHr)

# Run the complete model
FDivLM_all <- lm(FDiv ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC, data = df_var) 

# Check the adjustment of the complete model
summary(FDivLM_all)
coef(FDivLM_all)
car::vif(FDivLM_all)

# Plot results of the complete linear regression model
par(mfrow = c(2, 2))
plot(FDivLM_all)

# Obtain predicted values and residuals, and some tests
FDivLM_all_residual <- residuals(FDivLM_all)
FDivLM_all_preds <- predict(FDivLM_all)

pearson.test(FDivLM_all_residual)
lillie.test(FDivLM_all_residual)
lmtest::bptest(FDivLM_all)
search <- step(FDivLM_all, ~.^2)

# Plot of predicted values vs residuals 
par(mfrow = c(1, 1))
plot(FDivLM_all_residual ~  FDivLM_all_preds, xlab = "Predicted Values",  ylab = "Residuals")

# Select the best model: stepwise, backwise and both of them with stepAIC function
FDivLM_all_best_model <- stepAIC(FDivLM_all, direction = c("both", "backward", "forward"))

# Test null model of linear model
FDiv_nullLM <- lm(FDiv ~ 1, data = df_var)
summary(FDiv_nullLM)
stepAIC(FDiv_nullLM)

# Run the generalized linear model with FDiv ~ Cnratio + PCA1_temp + AGC
FDivGLM_all <- glm(FDiv ~ Cnratio + PCA1_temp + AGC, data = df_var, 
                   family = gaussian()) 
summary(FDivGLM_all)

# Plot results of the complete general linear regression model
par(mfrow = c(2, 2))
plot(FDivGLM_all)

# Obtain predicted values and residuals, and some tests
residuals_FDivglm <- residuals(FDivGLM_all)
preds_FDivglm <- predict(FDivGLM_all)

pearson.test(residuals_FDivglm)
lillie.test(residuals_FDivglm)

# Plot of predicted values vs residuals 
par(mfrow = c(1, 1))
plot(residuals_FDivglm ~ preds_FDivglm, xlab = "Predicted Values",  ylab = "Residuals")

#Select the best model: stepwise, backwise and both of them with glmulti function
FDiv_glmulti = glmulti(FDiv ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC , 
                       data = df_var, level = 1, method = "h", crit = "aic", 
                       confsetsize = 100, plotty = F, report = T, 
                       fitfunction = "glm")

# More info: http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti

# Obtain the best models (aic + 2)
best_models2_glm = weightable(FDiv_glmulti)
best_models2_glm = best_models2_glm[best_models2_glm$aic <= min(best_models2_glm$aic) + 2,]

# Check the best model to obtain beta, SE, and p value
summary(FDiv_glmulti@objects[[1]])

#######################################################
####### MODEL 3: FUNCTIONAL RICHNESS ##################
#######################################################

# Check relationships between FRic and other variables
FRicLM_AGC <- lm(FRic ~ AGC, data = df_var)
summary(FRicLM_AGC)

FRicLM_Cnratio <- lm(FRic ~ Cnratio, data = df_var) 
summary(FRicLM_Cnratio)

FRicLM_Pmgkg <- lm(FRic ~ Pmgkg, data = df_var) 
summary(FRicLM_Pmgkg)

FRicLM_PCA1_temp <- lm(FRic ~ PCA1_temp, data = df_var) 
summary(FRicLM_PCA1_temp)

FRicLM_minDHr <- lm(FRic ~ minDHr, data = df_var) 
summary(FRicLM_minDHr)

# Run the complete model
FRicLM_all <- lm(FRic ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC, data = df_var) 

# Check the adjustment of the complete model
summary(FRicLM_all)
coef(FRicLM_all)
car::vif(FRicLM_all)

# Plot results of the complete linear regression model
par(mfrow = c(2, 2))
plot(FRicLM_all)

# Obtain predicted values and residuals, and some tests
FRicLM_all_residual <- residuals(FRicLM_all)
FRicLM_all_preds <- predict(FRicLM_all)

pearson.test(FRicLM_all_residual)
lillie.test(FRicLM_all_residual)
lmtest::bptest(FRicLM_all)
search <- step(FRicLM_all, ~.^2)

# Plot of predicted values vs residuals 
par(mfrow = c(1, 1))
plot(FRicLM_all_residual ~ FRicLM_all_preds, xlab = "Predicted Values",  ylab = "Residuals")

# Select the best model: stepwise, backwise and both of them with stepAIC function
FRicLM_all_best_model <- stepAIC(FRicLM_all, direction = c("both", "backward", "forward"))

# Test null model of linear model
FRic_nullLM <- lm(FRic ~ 1, data = df_var)
summary(FRic_nullLM)
stepAIC(FRic_nullLM)

# Run the generalized linear model with FRic ~ Cnratio + PCA1_temp + AGC
FRicGLM_all <- glm(FRic ~ Cnratio + PCA1_temp + AGC, data = df_var, 
                   family = gaussian()) 
summary(FRicGLM_all)

# Plot results of the complete general linear regression model
par(mfrow = c(2, 2))
plot(FRicGLM_all)

# Obtain predicted values and residuals, and some tests
residuals_FRicglm <- residuals(FRicGLM_all)
preds_FRicglm <- predict(FRicGLM_all)

pearson.test(residuals_FRicglm)
lillie.test(residuals_FRicglm)

# Plot of predicted values vs residuals 
par(mfrow = c(1, 1))
plot(residuals_FRicglm ~ preds_FRicglm, xlab = "Predicted Values",  ylab = "Residuals")

#Select the best model: stepwise, backwise and both of them with glmulti function
FRic_glmulti = glmulti(FRic ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC , 
                       data = df_var, level = 1, method = "h", crit = "aic", 
                       confsetsize = 100, plotty = F, report = T, 
                       fitfunction = "glm")

# More info: http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti

# Obtain the best models (aic + 2)
best_models3_glm = weightable(FRic_glmulti)
best_models3_glm = best_models3_glm[best_models3_glm$aic <= min(best_models3_glm$aic) + 2,]

# Check the best model to obtain beta, SE, and p value
summary(FRic_glmulti@objects[[1]])

#######################################################
####### MODEL 4: FUNCTIONAL DIS #######################
#######################################################

# Run the complete model
FDisLM_all <- lm(FDis ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC, data = df_var) 

# Check the adjustment of the complete model
summary(FDisLM_all)
coef(FDisLM_all)
car::vif(FDisLM_all)

# Plot results of the complete linear regression model
par(mfrow = c(2, 2))
plot(FDisLM_all)

# Obtain predicted values and residuals, and some tests
FDisLM_all_residual <- residuals(FDisLM_all)
FDisLM_all_preds <- predict(FDisLM_all)

pearson.test(FDisLM_all_residual)
lillie.test(FDisLM_all_residual)
lmtest::bptest(FDisLM_all)
search <- step(FDisLM_all, ~.^2)

# Plot of predicted values vs residuals 
par(mfrow = c(1, 1))
plot(FDisLM_all_residual ~ FDisLM_all_preds, xlab = "Predicted Values",  ylab = "Residuals")

# Select the best model: stepwise, backwise and both of them with stepAIC function
FDisLM_all_best_model <- stepAIC(FDisLM_all, direction = c("both", "backward", "forward"))

# Test null model of linear model
FDis_nullLM <- lm(FDis ~ 1, data = df_var)
summary(FDis_nullLM)
stepAIC(FDis_nullLM)

# Run the generalized linear model with FDis ~ Cnratio + PCA1_temp + AGC
FDisGLM_all <- glm(FDis ~ Cnratio + PCA1_temp + AGC, data = df_var, 
                   family = gaussian()) 
summary(FDisGLM_all)

# Plot results of the complete general linear regression model
par(mfrow = c(2, 2))
plot(FDisGLM_all)

# Obtain predicted values and residuals, and some tests
residuals_FDisglm <- residuals(FDisGLM_all)
preds_FDisglm <- predict(FDisGLM_all)

pearson.test(residuals_FDisglm)
lillie.test(residuals_FDisglm)

# Plot of predicted values vs residuals 
par(mfrow = c(1, 1))
plot(residuals_FDisglm ~ preds_FDisglm, xlab = "Predicted Values",  ylab = "Residuals")

#Select the best model: stepwise, backwise and both of them with glmulti function
FDis_glmulti = glmulti(FDis ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC , 
                       data = df_var, level = 1, method = "h", crit = "aic", 
                       confsetsize = 100, plotty = F, report = T, 
                       fitfunction = "glm")

# More info: http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti

# Obtain the best models (aic + 2)
best_models4_glm = weightable(FDis_glmulti)
best_models4_glm = best_models4_glm[best_models4_glm$aic <= min(best_models4_glm$aic) + 2,]

# Check the best model to obtain beta, SE, and p value
summary(FDis_glmulti@objects[[1]])

#######################################################
####### MODEL 5: TAXONOMIC DIV ########################
#######################################################

# Run the complete model
HillDivLM_all <- lm(Hill_Diversity ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC, data = df_var) 

# Check the adjustment of the complete model
summary(HillDivLM_all)
coef(HillDivLM_all)
car::vif(HillDivLM_all)

# Plot results of the complete linear regression model
par(mfrow = c(2, 2))
plot(HillDivLM_all)

# Obtain predicted values and residuals, and some tests
HillDivLM_all_residual <- residuals(HillDivLM_all)
HillDivLM_all_preds <- predict(HillDivLM_all)

pearson.test(HillDivLM_all_residual)
lillie.test(HillDivLM_all_residual)
lmtest::bptest(HillDivLM_all)
search <- step(HillDivLM_all, ~.^2)

# Plot of predicted values vs residuals 
par(mfrow = c(1, 1))
plot(HillDivLM_all_residual ~ HillDivLM_all_preds, xlab = "Predicted Values",  
     ylab = "Residuals")

# Select the best model: stepwise, backwise and both of them with stepAIC function
HillDivLM_all_best_model <- stepAIC(HillDivLM_all, direction = c("both", 
                                                                 "backward", 
                                                                 "forward"))

# Test null model of linear model
HillDiv_nullLM <- lm(Hill_Diversity ~ 1, data = df_var)
summary(HillDiv_nullLM)
stepAIC(HillDiv_nullLM)

# Run the generalized linear model with FDis ~ Cnratio + PCA1_temp + AGC
HillDivGLM_all <- glm(Hill_Diversity ~ Cnratio + PCA1_temp + AGC, data = df_var, 
                      family = gaussian()) 
summary(HillDivGLM_all)

# Plot results of the complete general linear regression model
par(mfrow = c(2, 2))
plot(HillDivGLM_all)

# Obtain predicted values and residuals, and some tests
residuals_HillDivglm <- residuals(HillDivGLM_all)
preds_HillDivglm <- predict(HillDivGLM_all)

pearson.test(residuals_HillDivglm)
lillie.test(residuals_HillDivglm)

# Plot of predicted values vs residuals 
par(mfrow = c(1, 1))
plot(residuals_HillDivglm ~ preds_HillDivglm, xlab = "Predicted Values",  ylab = "Residuals")

#Select the best model: stepwise, backwise and both of them with glmulti function
HillDiv_glmulti = glmulti(Hill_Diversity ~ Cnratio + Pmgkg + PCA1_temp + minDHr + AGC, 
                          data = df_var, level = 1, method = "h", crit = "aic", 
                          confsetsize = 100, plotty = F, report = T, 
                          fitfunction = "glm")

# More info: http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti

# Obtain the best models (aic + 2)
best_models5_glm = weightable(HillDiv_glmulti)
best_models5_glm = best_models5_glm[best_models5_glm$aic <= min(best_models5_glm$aic) + 2,]

# Check the best model to obtain beta, SE, and p value
summary(HillDiv_glmulti@objects[[1]])
