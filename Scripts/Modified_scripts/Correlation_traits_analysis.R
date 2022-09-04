# Script to calculate correlations between all climactic, soil 
# and hydrological traits

# Import libraries
library(corrplot)
library(corrmorant)
library(Hmisc)

# Load data
df_traits_orig <- read.csv("Data/Functional_traits/20210811_MeanTraits.csv")
df_env_var <- read.csv("Data/Environmental_variables.csv")
# Add plot names as indices and delete column with names
row.names(df_env_var) <- df_env_var$PLOT
df_env_var[1] <- NULL

#######################################################
################ FUNCTIONAL TRAITS ####################
#######################################################

# Select traits for the analysis
df_funct_traits <- df_traits_orig[c("WD", "Ht2019", "Sum_BA_2019.m2.",
                              "Mean.LA", "Mean.SLA", 
                              "Mean.LBT", "Mean.LDMC")]

# Run correlation analysis 
corr_analysis_funct <-rcorr(as.matrix(df_funct_traits), type = "spearman")

# Obtain correlation matrix
corr_matrix_funct <- corr_analysis_funct$r

# Obtain p-values matrix
corr_pvalues_funct <- corr_analysis_funct$P

# Create correlation plot with corrplot
# Define colors for plot
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# Function to save plot
setEPS()
png("Outputs/Correlation_traits/corr_plot_funct1.png", width = 700, height = 350)

corrplot(corr_matrix_funct, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and 
         tl.cex = 0.9, # text size
         tl.offset = 1, # text label
         number.cex = 0.91, # size of the correlation coefficient
         # Combine with significance
         p.mat = corr_pvalues_funct, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

dev.off()

# Create correlation plot with corrmorant
corr_plot_funct2 <- ggcorrm(data = df_funct_traits) +
  theme_corrm(base_size = 6)+
  theme(axis.text.x = element_text(angle = 90, size = 4.5), 
        axis.text.y = element_text(size = 4.5),
        legend.text = element_text(size = 4.5),
        legend.title = element_text(size = 6)) +
  lotri(geom_point(alpha = 0.5)) + 
  lotri(geom_smooth(colour = "red4")) +
  utri_heatmap(alpha = 0.5, corr_method = "spearman") +
  utri_corrtext(corr_method = "spearman", size = 4) +
  dia_histogram(lower = 0.1, fill = "grey80", color = 1)+
  dia_density(lower = 0.1, alpha = .1, colour = "red4") +
  scale_fill_gradient2(low = "white", mid = "red3", high ="red4", 
                       midpoint = 0.5, space = "rgb", guide = guide_colorbar(title = "Correlation coefficient"),
                       limits = c(0,1))

# Export plot figure
ggsave(filename = "Outputs/Correlation_traits/corr_plot_funct2.png",
       plot = corr_plot_funct2, 
       width = 10,
       height = 8,
       units = "in")

# Export correlation matrices
write.csv(corr_matrix_funct,file="Outputs/Correlation_traits/corr_matrix_funct.csv",row.names = FALSE)
write.csv(corr_pvalues_funct,file = "Outputs/Correlation_traits/corr_pvalues_funct.csv",row.names = FALSE)

#######################################################
########### PRODUCTIVITY VARIABLES ####################
#######################################################

# Select traits for the analysis
df_prod_var <- df_env_var[c("AGC", "AGCnc", "AGCp")]

# Run correlation analysis 
corr_prodvar <- rcorr(as.matrix(df_prod_var), type = "spearman")

# Obtain correlation matrix
corr_mat_prodvar <- corr_prodvar$r

# Obtain p-values matrix
corr_pval_prodvar <- corr_prodvar$P

# Create correlation plot with corrplot
# Define colors for plot
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# Function to save plot
setEPS()
png("Outputs/Correlation_traits/corr_plot_prodvar1.png", width = 700, height = 350)

corrplot(corr_mat_prodvar, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and 
         tl.cex = 0.9, # text size
         tl.offset = 1, # text label
         number.cex = 0.91, # size of the correlation coefficient
         # Combine with significance
         p.mat = corr_pval_prodvar, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

dev.off()

# Create correlation plot with corrmorant
corr_plot_prodvar2 <- ggcorrm(data = df_prod_var) +
  theme_corrm(base_size = 6)+
  theme(axis.text.x = element_text(angle = 90, size = 4.5), 
        axis.text.y = element_text(size = 4.5),
        legend.text = element_text(size = 4.5),
        legend.title = element_text(size = 6)) +
  lotri(geom_point(alpha = 0.5)) + 
  lotri(geom_smooth(colour = "red4")) +
  utri_heatmap(alpha = 0.5, corr_method = "spearman") +
  utri_corrtext(corr_method = "spearman", size = 4) +
  dia_histogram(lower = 0.1, fill = "grey80", color = 1)+
  dia_density(lower = 0.1, alpha = .1, colour = "red4") +
  scale_fill_gradient2(low = "white", mid = "red3", high ="red4", 
                       midpoint = 0.5, space = "rgb", guide = guide_colorbar(title = "Correlation coefficient"),
                       limits = c(0,1))

# Export plot figure
ggsave(filename = "Outputs/Correlation_traits/corr_plot_prodvar2.png",
       plot = corr_plot_prodvar2, 
       width = 7,
       height = 5,
       units = "in")

# Export correlation matrices
write.csv(corr_mat_prodvar, file="Outputs/Correlation_traits/corr_mat_prodvar.csv",row.names = FALSE)
write.csv(corr_pval_prodvar,file = "Outputs/Correlation_traits/corr_pval_prodvar.csv",row.names = FALSE)

#######################################################
########### CLIMATIC VARIABLES ####################
#######################################################

# Select traits for the analysis
df_climvar <- df_env_var[c("PCA1_temp", "AP", "PC1.p", "PC2.p")]

# Run correlation analysis 
corr_climvar <- rcorr(as.matrix(df_clim_var), type = "spearman")

# Obtain correlation matrix
corr_mat_climvar <- corr_climvar$r

# Obtain p-values matrix
corr_pval_climvar <- corr_climvar$P

# Create correlation plot with corrplot
# Define colors for plot
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# Function to save plot
setEPS()
png("Outputs/Correlation_traits/corr_plot_climvar1.png", width = 700, height = 350)

corrplot(corr_mat_climvar, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and 
         tl.cex = 0.9, # text size
         tl.offset = 1, # text label
         number.cex = 0.91, # size of the correlation coefficient
         # Combine with significance
         p.mat = corr_pval_climvar, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

dev.off()

# Create correlation plot with corrmorant
corr_plot_climvar2 <- ggcorrm(data = df_climvar) +
  theme_corrm(base_size = 6)+
  theme(axis.text.x = element_text(angle = 90, size = 4.5), 
        axis.text.y = element_text(size = 4.5),
        legend.text = element_text(size = 4.5),
        legend.title = element_text(size = 6)) +
  lotri(geom_point(alpha = 0.5)) + 
  lotri(geom_smooth(colour = "red4")) +
  utri_heatmap(alpha = 0.5, corr_method = "spearman") +
  utri_corrtext(corr_method = "spearman", size = 4) +
  dia_histogram(lower = 0.1, fill = "grey80", color = 1)+
  dia_density(lower = 0.1, alpha = .1, colour = "red4") +
  scale_fill_gradient2(low = "white", mid = "red3", high ="red4", 
                       midpoint = 0.5, space = "rgb", guide = guide_colorbar(title = "Correlation coefficient"),
                       limits = c(0,1))

# Export plot figure
ggsave(filename = "Outputs/Correlation_traits/corr_plot_climvar2.png",
       plot = corr_plot_climvar2, 
       width = 7,
       height = 5,
       units = "in")

# Export correlation matrices
write.csv(corr_mat_climvar, file="Outputs/Correlation_traits/corr_mat_climvar.csv",row.names = FALSE)
write.csv(corr_pval_climvar,file = "Outputs/Correlation_traits/corr_pval_climvar.csv",row.names = FALSE)

#######################################################
########### SOIL VARIABLES ############################
#######################################################

# Select traits for the analysis
df_soilvar <- df_env_var[c("C_N", "C", "N", "S")]

# Run correlation analysis 
corr_soilvar <- rcorr(as.matrix(df_soilvar), type = "spearman")

# Obtain correlation matrix
corr_mat_soilvar <- corr_soilvar$r

# Obtain p-values matrix
corr_pval_soilvar <- corr_soilvar$P

# Create correlation plot with corrplot
# Define colors for plot
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# Function to save plot
setEPS()
png("Outputs/Correlation_traits/corr_plot_soilvar1.png", width = 700, height = 350)

corrplot(corr_mat_soilvar, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and 
         tl.cex = 0.9, # text size
         tl.offset = 1, # text label
         number.cex = 0.91, # size of the correlation coefficient
         # Combine with significance
         p.mat = corr_pval_soilvar, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

dev.off()

# Create correlation plot with corrmorant
corr_plot_soilvar2 <- ggcorrm(data = df_soilvar) +
  theme_corrm(base_size = 6)+
  theme(axis.text.x = element_text(angle = 90, size = 4.5), 
        axis.text.y = element_text(size = 4.5),
        legend.text = element_text(size = 4.5),
        legend.title = element_text(size = 6)) +
  lotri(geom_point(alpha = 0.5)) + 
  lotri(geom_smooth(colour = "red4")) +
  utri_heatmap(alpha = 0.5, corr_method = "spearman") +
  utri_corrtext(corr_method = "spearman", size = 4) +
  dia_histogram(lower = 0.1, fill = "grey80", color = 1)+
  dia_density(lower = 0.1, alpha = .1, colour = "red4") +
  scale_fill_gradient2(low = "white", mid = "red3", high ="red4", 
                       midpoint = 0.5, space = "rgb", guide = guide_colorbar(title = "Correlation coefficient"),
                       limits = c(0,1))

# Export plot figure
ggsave(filename = "Outputs/Correlation_traits/corr_plot_soilvar2.png",
       plot = corr_plot_soilvar2, 
       width = 7,
       height = 5,
       units = "in")

# Export correlation matrices
write.csv(corr_mat_soilvar, file="Outputs/Correlation_traits/corr_mat_soilvar.csv",row.names = FALSE)
write.csv(corr_pval_soilvar,file = "Outputs/Correlation_traits/corr_pval_soilvar.csv",row.names = FALSE)

#######################################################
########### ALL VARIABLES #############################
#######################################################

# Run correlation analysis 
corr_allvar <- rcorr(as.matrix(df_env_var), type = "spearman")

# Obtain correlation matrix
corr_mat_allvar <- corr_allvar$r

# Obtain p-values matrix
corr_pval_allvar <- corr_allvar$P

# Create correlation plot with corrplot
# Define colors for plot
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# Function to save plot
setEPS()
png("Outputs/Correlation_traits/corr_plot_allvar1.png", width = 700, height = 350)

corrplot(corr_mat_allvar, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and 
         tl.cex = 0.9, # text size
         tl.offset = 1, # text label
         number.cex = 0.91, # size of the correlation coefficient
         # Combine with significance
         p.mat = corr_pval_allvar, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

dev.off()

# Create correlation plot with corrmorant
corr_plot_allvar2 <- ggcorrm(data = df_env_var) +
  theme_corrm(base_size = 6)+
  theme(axis.text.x = element_text(angle = 90, size = 4.5), 
        axis.text.y = element_text(size = 4.5),
        legend.text = element_text(size = 4.5),
        legend.title = element_text(size = 6)) +
  lotri(geom_point(alpha = 0.5)) + 
  lotri(geom_smooth(colour = "red4")) +
  utri_heatmap(alpha = 0.5, corr_method = "spearman") +
  utri_corrtext(corr_method = "spearman", size = 4) +
  dia_histogram(lower = 0.1, fill = "grey80", color = 1)+
  dia_density(lower = 0.1, alpha = .1, colour = "red4") +
  scale_fill_gradient2(low = "white", mid = "red3", high ="red4", 
                       midpoint = 0.5, space = "rgb", guide = guide_colorbar(title = "Correlation coefficient"),
                       limits = c(0,1))

# Export plot figure
ggsave(filename = "Outputs/Correlation_traits/corr_plot_allvar2.png",
       plot = corr_plot_allvar2, 
       width = 9,
       height = 7,
       units = "in")

# Export correlation matrices
write.csv(corr_mat_allvar, file="Outputs/Correlation_traits/corr_mat_allvar.csv",row.names = FALSE)
write.csv(corr_pval_allvar,file = "Outputs/Correlation_traits/corr_pval_allvar.csv",row.names = FALSE)
