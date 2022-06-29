############# CHARCOAL ANALYSIS ###############

#############   DIFFERENCES BETWEEN FOREST TYPES: FIRE INTENSITY ###############
######################       FITTING LINEAR MODEL         #####################           IV)

num200<- read.table("d:/oxford/congo/analysis/estadistica_phd/charcoal/richard/av_ch_200yr.txt", header=T)
#num250<- read.table("d:/oxford/congo/analysis/estadistica_phd/richard/av_ch_250yr.txt", header=T)
num200<-num200[num200$bin<2600,]
#num250<-num250[num250$bin<2600,]

# from data exploration there is hetereogeneity,
# first try a linear model: fixed effects and interactions

####### Fit the Model with GLS                                                            IV.1)
# In this step we fit the model using the gls function. It allows us to compare the
# linear regression model with the mixed effects model that we will calculate using
# the lme function

library(nlme)

form200 <- formula(ch_flow ~ forest_name * bin)
form200.gls <- gls(form200, na.action=na.omit, data = num200)
form200.lme <- lme(form200, random = ~ 1 | core, method = "REML", data = num200,
                   na.action=na.omit)
anova(form200.gls, form200.lme)

###### Checking assumptions using plots                                                   IV.2)

# Best model was lme based on AIC
# Check for residuals (normality, heterogeneity and independence)

res200.lme<-resid(form200.lme,type="normalized")
fit200.lme<- fitted(form200.lme)

plot(x=fit200.lme,y=res200.lme,xlab="fitted values",ylab="normalized residuals")              ## heterogeneity
plot(x=num200$bin[!is.na(num200$ch_flow)],y=res200.lme, ylab="residuals", xlab="time")        ## heterogeneity by time

hist(res200.lme)                                                                              ## normality
boxplot(res200.lme, main="boxplot of residuals")                                              ## normality

boxplot(res200.lme~core, data=num200[!is.na(num200$ch_flow),],na.action=na.omit,
        ylim=c(-3,3),main="Residuals by core")                                                ## independence by core
abline(0,0); axis(2)

boxplot(res200.lme~forest_name, data=num200[!is.na(num200$ch_flow),],na.action=na.omit,
        ylim=c(-3,3),main="Residuals by forest type")                                         ## independence by forest type
abline(0,0); axis(2)

boxplot(res200.lme~bin, data=num200[!is.na(num200$ch_flow),],na.action=na.omit,         ## independence by time (bloxplot)
        ylim=c(-3,3),main="Residuals by bin")
abline(0,0); axis(2)

### Boxplots by core: residuals should lie in a cloud around the line without any
### patterns, for each x value, all residuals should not be above or below
### the zero line, otherwise the x variable has to be included in the model

### Autocorrelation
### it is important to include the NA values

a<- !is.na(num200$ch_flow)                                                          ## independence by time (autocorrelation)
res200.lme_full<- vector(length = length(num200$ch_flow))
res200.lme_full<- NA
res200.lme_full[a]<-res200.lme
acf(res200.lme_full, na.action = na.pass, main = "Auto-correlation plot for residuals")

library(lattice)                                                                    ## independence by time (adding a smoother by forest type)
xyplot(res200.lme ~ bin | forest_name,
       data = num200, ylab = "Residuals",
       xlab = "bins",
       ylim = c(-2,2),             ######## to constraint values due to large values in core 13
       panel = function(x,y){
         panel.grid(h = -1, v = 2)
         panel.points(x, y, col = 1)
         panel.loess(x, y, span = 0.7, col = 1, lwd=2)})

library(lattice)                                                                    ## independence by time (adding a smoother by core)
xyplot(res200.lme ~ bin | core,
       data = num200, ylab = "Residuals",
       xlab = "bins",
       ylim = c(-2,2),
       panel = function(x,y){
         panel.grid(h = -1, v = 2)
         panel.points(x, y, col = 1)
         panel.loess(x, y, span = 0.5, col = 1, lwd=2)})

#### Plots with smoother: We have to be sure about independence between samples.
### the multipanel shows residuals of generalized linear model vs bins for each
### forest type. The added LOESS smoother should not show any pattern if
### linear model ok, if we observe a trend therefore we fit an additive model


##### FIXING SOME OF THE PROBLEMS

# Heterogeneity

library(nlme)

m200_lme2<-lme(ch_flow ~ forest_name + bin, weights=varIdent(form=~ 1|core),
               random = ~ 1 | core, method = "REML", data = num200,na.action=na.omit)

m200_lme3<-lme(ch_flow ~ forest_name + bin, weights=varIdent(form=~ 1|forest_name),
               random = ~ 1 | core, method = "REML",data = num200,na.action=na.omit)

anova(m200_lme2,m200_lme3)

### anova shows lower values of AIC for lme2

plot(m200_lme2)                                                                      ## heterogeneity

a<- !is.na(num200$ch_flow)                                                          ## independence by time (autocorrelation)
res200_lme2_full<- vector(length = length(num200$ch_flow))
res200_lme2_full<- NA
res200_lme2_full[a]<-resid(m200_lme2,type="normalized")
acf(res200_lme2_full, na.action = na.pass, main = "Auto-correlation plot for residuals")



### Better spread of residuals and less AIC values for lme3 but temporal autocorrelation
### Explore linearity of bins

library(lattice)
xyplot(resid(m200_lme2,type="normalized") ~ num200$bin[!is.na(num200$ch_flow)],
       data = num200, ylab = "Residuals",
       xlab = "bins",
       panel = function(x,y){
         panel.grid(h = -1, v = 2)
         panel.points(x, y, col = 1)
         panel.loess(x, y, span = 0.5, col = 1, lwd=2)})


### Still temporal autocorrelation, try correlation component corCompSymm

m200_lme4<-lme(ch_flow ~ forest_name + bin, weights=varIdent(form=~ 1|forest_name),
               random = ~ 1 | core, method = "REML", correlation = corCompSymm(form=~ bin),
               data = num200,na.action=na.omit)

plot(m200_lme4)

a<- !is.na(num200$ch_flow)                                                          ## independence by time (autocorrelation)
res200_lme4_full<- vector(length = length(num200$ch_flow))
res200_lme4_full<- NA
res200_lme4_full[a]<-resid(m200_lme4,type="normalized")
acf(res200_lme4_full, na.action = na.pass, main = "Auto-correlation plot for residuals")

library(lattice)
xyplot(resid(m200_lme4,type="normalized") ~ num200$bin[!is.na(num200$ch_flow)],
       data = num200, ylab = "Residuals",
       xlab = "bins",
       panel = function(x,y){
         panel.grid(h = -1, v = 2)
         panel.points(x, y, col = 1)
         panel.loess(x, y, span = 0.5, col = 1, lwd=2)})

### still autocorrelation, try another correlation structure ARIMA

m200_lme5<-lme(ch_flow ~ forest_name + bin, weights=varIdent(form=~ 1|forest_name),
               random = ~ 1 | core, method = "REML", correlation = corAR1(form=~ bin),
               data = num200,na.action=na.omit)

plot(m200_lme5)

a<- !is.na(num200$ch_flow)                                                          ## independence by time (autocorrelation)
res200_lme5_full<- vector(length = length(num200$ch_flow))
res200_lme5_full<- NA
res200_lme5_full[a]<-resid(m200_lme5,type="normalized")
acf(res200_lme5_full, na.action = na.pass, main = "Auto-correlation plot for residuals")

library(lattice)
xyplot(resid(m200_lme5,type="normalized") ~ num200$bin[!is.na(num200$ch_flow)],
       data = num200, ylab = "Residuals",
       xlab = "bins",
       panel = function(x,y){
         panel.grid(h = -1, v = 2)
         panel.points(x, y, col = 1)
         panel.loess(x, y, span = 0.5, col = 1, lwd=2)})

AIC(form200.lme,m200_lme2, m200_lme4, m200_lme5)

### other correlation bin/core
m200_lme6<-lme(ch_flow ~ forest_name + bin, weights=varIdent(form=~ 1|forest_name),
               random = ~ 1 | core, method = "REML", correlation = corAR1(form=~ bin),
               data = num200,na.action=na.omit)

library(lattice)
xyplot(resid(m200_lme6,type="normalized") ~ num200$bin[!is.na(num200$ch_flow)],
       data = num200, ylab = "Residuals",
       xlab = "bins",
       panel = function(x,y){
         panel.grid(h = -1, v = 2)
         panel.points(x, y, col = 1)
         panel.loess(x, y, span = 0.5, col = 1, lwd=2)})
### Model does not improve and there is evidence of non-linearity of bins,
### so try additive models

