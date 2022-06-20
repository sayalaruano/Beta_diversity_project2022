
# Set Working directory
wDpaper = "C:/Users/antonella.bernardi/OneDrive - Universidad de Las Américas/Francisco Cuesta/paper"

setwd(wDpaper)

### Running correlations between all climatical, soil and hydrological variables ######

library(corrplot)
library(Hmisc)

###### LOADING DATA ########

dat<-read.csv("20210810_MeanTraits.csv")
head(dat)


dat2 = veg[, c(14, 15, 17, 20)]

dat2 = veg[, c(14, 15, 17, 20, 4, 5, 6, 7, 8, 9, 10, 21, 22, 23, 25)]

dat2= veg[, c(14, 15, 17, 20, 4, 5, 6, 7, 8, 9, 10, 33, 34, 39, 40, 41, 42, 43, 44, 45)]

### removing one plot with no data and selecting variables for the analysis
dat2<-dat[c("Mean.LA", "Mean.SLA", "Mean.N.por", "Mean.P.por")]
head(dat2)

# "Mean.LBT", "Mean.LA", "Mean.SLA", "Mean.LDMC", 

###### RUNNING CORRELATION ANALYSIS ####

#cor0<-cor(as.matrix(dat2), use = "complete.obs", method = "spearman")

corr<-rcorr(as.matrix(dat2), type = "spearman")

## correlation matrix
corr_matrix<-corr$r

## p-values matrix
corr_p.values<-corr$P

##### GRAPH and saving results

corrplot(corr_matrix, type="upper", order="hclust", 
         p.mat = corr_p.values, sig.level = 0.05, na.rm = TRUE)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

setEPS()
svg("../analysis/variables_correlations.svg", width = 4.39, height = 4.28, pointsize = 8)
corrplot(corr_matrix, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and 
         tl.cex = 0.9, # text size
         tl.offset = 1, # text label
         number.cex = 0.91, # size of the correlation coefficient
         # Combine with significance
         p.mat = corr_p.values, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

dev.off()

#### writing results

write.csv(corr_matrix,file="../analysis/correlation_matrix.csv",row.names = FALSE)
write.csv(corr_p.values,file = "../analysis/correlation_p.values.csv",row.names = FALSE)



dat.pca <- prcomp(dat[,c(2,5:8,13,14)], center = TRUE, scale. = TRUE, na.rm = TRUE)

summary(dat.pca)

