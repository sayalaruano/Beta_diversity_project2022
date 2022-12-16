# Script to run linear models - GLMs

# Import libraries
library(Hmisc)
library(nortest)
library(MASS)
library(glmulti) # this package requires rJava
# More info: http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti

# Load data
df_amb_var <- read.csv("Data/Environmental_variables.csv")
df_fucnt_alphadiv <- read.csv("Data/Funct_alpha_div_indices2022.csv")
df_tax_alphadiv <- read.csv("Data/Taxon_div_indices2022.csv")

# Create unique dataset of env, funct, and tax data
# Select columns and add plot names as indices and delete column with names
df_amb_var_filt <- df_amb_var[,colnames(df_amb_var) %in% c('Plot','C_N',
                                                        'AGCp', 'PCA1_temp', 
                                                        'PC2.p')]

row.names(df_amb_var_filt) <- df_amb_var_filt$Plot
df_amb_var_filt[1] <- NULL

df_fucnt_alphadiv_filt <- df_fucnt_alphadiv[,colnames(df_fucnt_alphadiv) %in% c('Plot','fric',
                                                           'feve')]

row.names(df_fucnt_alphadiv_filt) <- df_fucnt_alphadiv_filt$Plot
df_fucnt_alphadiv_filt[1] <- NULL

df_tax_alphadiv_filt <- df_tax_alphadiv[,colnames(df_tax_alphadiv) %in% c('Plot','Hill_Diversity',
                                                           'Richness')]

row.names(df_tax_alphadiv_filt) <- df_tax_alphadiv_filt$Plot
df_tax_alphadiv_filt[1] <- NULL

# Merge datasets
amb_funct_df <- merge(df_amb_var_filt, df_fucnt_alphadiv_filt, by = "row.names")

row.names(amb_funct_df) <- amb_funct_df$Row.names
amb_funct_df[1] <- NULL

amb_funct_tax_df <- merge(amb_funct_df, df_tax_alphadiv_filt, by = "row.names")

row.names(amb_funct_tax_df) <- amb_funct_tax_df$Row.names
amb_funct_tax_df[1] <- NULL

#######################################################
################ TAXONOMIC DIVERSITY  #################
#######################################################

#######################################################
############## MODEL 1: HILL NUMBERS ##################
#######################################################

# Run the generalized linear model with Hill ~ Cnratio + PCA1_temp + AGCp + Pc2p
HillGLM_all <- glm(Hill_Diversity ~ C_N + PCA1_temp + AGCp + PC2.p, data = amb_funct_tax_df, 
                  family = gaussian()) 

summary(HillGLM_all)

# Plot results of the complete general linear regression model
par(mfrow = c(2, 2))
plot(HillGLM_all)

# Obtain predicted values and residuals, and some tests
residuals_HillGLM_all <- residuals(HillGLM_all)
preds_HillGLM_all <- predict(HillGLM_all)

pearson.test(residuals_HillGLM_all)
lillie.test(residuals_HillGLM_all)

# Plot of predicted values vs residuals 
par(mfrow = c(1, 1))
plot(residuals_HillGLM_all ~ preds_HillGLM_all, xlab = "Predicted Values",  ylab = "Residuals")

#Select the best model: stepwise, backwise and both of them with glmulti function
HillGLM_glmulti <- glmulti(Hill_Diversity ~ C_N + PCA1_temp + AGCp + PC2.p, data = amb_funct_tax_df, 
                       level = 1, method = "h", crit = "aic", 
                       confsetsize = 100, plotty = F, report = T, 
                       fitfunction = "glm")

# Obtain the best GML models (aic + 2)
best_models1_glm <- weightable(HillGLM_glmulti)
best_models1_glm <- best_models1_glm[best_models1_glm$aic <= min(best_models1_glm$aic),]

# Check the best model to obtain beta, SE, and p value
HillGLM_glmulti_summary <- summary(HillGLM_glmulti@objects[[1]])

# Export outputs
sink("Outputs/Linear_models/Summary_best_HillGLM.txt")
print(best_models1_glm)
print(HillGLM_glmulti_summary)
sink()
#unlink("Outputs/Linear_models/Summary_best_HillGLM.txt")

#######################################################
############## MODEL 2: RICHNESS ######################
#######################################################

# Run the generalized linear model with Hill ~ Cnratio + PCA1_temp + AGCp + Pc2p
RichGLM_all <- glm(Richness ~ C_N + PCA1_temp + AGCp + PC2.p, data = amb_funct_tax_df, 
                   family = gaussian()) 

summary(RichGLM_all)

# Plot results of the complete general linear regression model
par(mfrow = c(2, 2))
plot(RichGLM_all)

# Obtain predicted values and residuals, and some tests
residuals_RichGLM_all <- residuals(RichGLM_all)
preds_RichGLM_all <- predict(RichGLM_all)

pearson.test(residuals_RichGLM_all)
lillie.test(residuals_RichGLM_all)

# Plot of predicted values vs residuals 
par(mfrow = c(1, 1))
plot(residuals_RichGLM_all ~ preds_RichGLM_all, xlab = "Predicted Values",  ylab = "Residuals")

#Select the best model: stepwise, backwise and both of them with glmulti function
RichGLM_glmulti <- glmulti(Richness ~ C_N + PCA1_temp + AGCp + PC2.p, data = amb_funct_tax_df, 
                           level = 1, method = "h", crit = "aic", 
                           confsetsize = 100, plotty = F, report = T, 
                           fitfunction = "glm")

# Obtain the best GML models (aic + 2)
best_models2_glm <- weightable(RichGLM_glmulti)
best_models2_glm <- best_models2_glm[best_models2_glm$aic <= min(best_models2_glm$aic), ]

# Check the best model to obtain beta, SE, and p value
RichGLM_glmulti_summary <- summary(RichGLM_glmulti@objects[[1]])

# Export outputs
sink("Outputs/Linear_models/Summary_best_RichGLM.txt")
print(best_models2_glm)
print(RichGLM_glmulti_summary)
sink()
#unlink("Outputs/Linear_models/Summary_best_RichGLM.txt")

#######################################################
################ FUNCTIONAL DIVERSITY  ################
#######################################################

#######################################################
############## MODEL 3: FRICHNESS #####################
#######################################################

# Run the generalized linear model with Hill ~ Cnratio + PCA1_temp + AGCp + Pc2p
FricGLM_all <- glm(fric ~ C_N + PCA1_temp + AGCp + PC2.p, data = amb_funct_tax_df, 
                   family = gaussian()) 

summary(FricGLM_all)

# Plot results of the complete general linear regression model
par(mfrow = c(2, 2))
plot(FricGLM_all)

# Obtain predicted values and residuals, and some tests
residuals_FricGLM_all <- residuals(FricGLM_all)
preds_FricGLM_all <- predict(FricGLM_all)

pearson.test(residuals_FricGLM_all)
lillie.test(residuals_FricGLM_all)

# Plot of predicted values vs residuals 
par(mfrow = c(1, 1))
plot(residuals_FricGLM_all ~ preds_FricGLM_all, xlab = "Predicted Values",  ylab = "Residuals")

#Select the best model: stepwise, backwise and both of them with glmulti function
FricGLM_glmulti <- glmulti(fric ~ C_N + PCA1_temp + AGCp + PC2.p, data = amb_funct_tax_df, 
                           level = 1, method = "h", crit = "aic", 
                           confsetsize = 100, plotty = F, report = T, 
                           fitfunction = "glm")

# Obtain the best GML models (aic + 2)
best_models3_glm <- weightable(FricGLM_glmulti)
best_models3_glm <- best_models3_glm[best_models3_glm$aic <= min(best_models3_glm$aic),]

# Check the best model to obtain beta, SE, and p value
FricGLM_glmulti_summary <- summary(FricGLM_glmulti@objects[[1]])

# Export outputs
sink("Outputs/Linear_models/Summary_best_FricGLM.txt")
print(best_models3_glm)
print(FricGLM_glmulti_summary)
sink()
# unlink("Outputs/Linear_models/Summary_best_FricGLM.txt")

#######################################################
############## MODEL 4: FEVENESS ######################
#######################################################

# Run the generalized linear model with Hill ~ Cnratio + PCA1_temp + AGCp + Pc2p
FeveGLM_all <- glm(feve ~ C_N + PCA1_temp + AGCp + PC2.p, data = amb_funct_tax_df, 
                   family = gaussian()) 

summary(FeveGLM_all)

# Plot results of the complete general linear regression model
par(mfrow = c(2, 2))
plot(FeveGLM_all)

# Obtain predicted values and residuals, and some tests
residuals_FeveGLM_all <- residuals(FeveGLM_all)
preds_FeveGLM_all <- predict(FeveGLM_all)

pearson.test(residuals_FeveGLM_all)
lillie.test(residuals_FeveGLM_all)

# Plot of predicted values vs residuals 
par(mfrow = c(1, 1))
plot(residuals_FeveGLM_all ~ preds_FeveGLM_all, xlab = "Predicted Values",  ylab = "Residuals")

#Select the best model: stepwise, backwise and both of them with glmulti function
FeveGLM_glmulti <- glmulti(feve ~ C_N + PCA1_temp + AGCp + PC2.p, data = amb_funct_tax_df, 
                           level = 1, method = "h", crit = "aic", 
                           confsetsize = 100, plotty = F, report = T, 
                           fitfunction = "glm")

# Obtain the best GML models (aic + 2)
best_models4_glm <- weightable(FeveGLM_glmulti)
best_models4_glm <- best_models4_glm[best_models4_glm$aic <= min(best_models4_glm$aic), ]

# Check the best model to obtain beta, SE, and p value
FeveGLM_glmulti_summary <- summary(FeveGLM_glmulti@objects[[1]])

# Export outputs
sink("Outputs/Linear_models/Summary_best_FeveGLM.txt")
print(best_models4_glm)
print(FeveGLM_glmulti_summary)
sink()
# unlink("Outputs/Linear_models/Summary_best_FeveGLM.txt")

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
