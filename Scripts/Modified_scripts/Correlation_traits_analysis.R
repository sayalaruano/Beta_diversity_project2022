# Script to calculate correlations between all climactic, soil 
# and hydrological traits

# Import libraries
library(corrplot)
library(corrmorant)
library(Hmisc)

# Load data
df_traits_orig <- read.csv("Data/20210811_MeanTraits.csv")

# Select traits for the analysis
df_traits <- df_traits_orig[c("Mean.LA", "Mean.SLA", "Mean.N.por", "Mean.P.por")]

# Run correlation analysis 
corr_analysis <-rcorr(as.matrix(df_traits), type = "spearman")

# Obtain correlation matrix
corr_matrix <- corr_analysis$r

# Obtain p-values matrix
corr_pvalues <- corr_analysis$P

# Create correlation plot with corrplot
# Define colors for plot
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# Function to save plot
setEPS()
svg("Outputs/Correlation_traits/correlation_traits_plot.svg", width = 4.39, height = 4.28, pointsize = 8)

corrplot(corr_matrix, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and 
         tl.cex = 0.9, # text size
         tl.offset = 1, # text label
         number.cex = 0.91, # size of the correlation coefficient
         # Combine with significance
         p.mat = corr_pvalues, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

dev.off()

# Create correlation plot with corrmorant
corr_traits2 <- ggcorrm(data = df_traits) +
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
ggsave(filename = "Outputs/Correlation_traits/correlation_traits_plot2.svg",
       plot = corr_traits2, 
       width = 6.61,
       height = 4.18,
       units = "in")

# Export correlation matrices
write.csv(corr_matrix,file="Outputs/Correlation_traits/correlation_matrix.csv",row.names = FALSE)
write.csv(corr_p.values,file = "Outputs/Correlation_traits/correlation_p.values.csv",row.names = FALSE)
