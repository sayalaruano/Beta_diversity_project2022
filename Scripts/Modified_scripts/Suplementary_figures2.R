# Imports
library(ggpubr)
library(ggplot2)

# Load data
df_indices = read.csv("Data/AlphFD_Indexvalues_270122.csv", sep = ",")
df_betafd = read.csv("Data/Beta_funct_div_matrix.csv", sep = ",")
df_betatd = read.csv("Data/Beta_tax_div_ba_rs_matrix.csv", sep = ",")
elev1 <- read.csv("Data/Altitudes_13052021.csv", sep = ";")

# Change name of the plot column and duplicate it to merge later with the pairwise distances df
colnames(elev1) <- c('Plot1','Elevation1')
elev2 = elev1
colnames(elev2) <- c('Plot2','Elevation2')

# 1)  Taxonomic dissimilarities among plots based on their altitudinal distance 
# All vs all
# Add plot names as indices and delete column with names
row.names(df_betatd) <- df_betatd$X
df_betatd$X <- NULL
# Convert dataframe into a matrix
df_betatd = as.matrix(df_betatd)

# Create dataframe with pairwise distances
df_betatd = data.frame(col = colnames(df_betatd)[col(df_betatd)], row = rownames(df_betatd)[row(df_betatd)], 
                       dist = c(df_betatd))

colnames(df_betatd) = c("Plot1", "Plot2", "Distances")

# Add elevation columns of plot1 and 2
df_betatd = merge(df_betatd, elev1, by = "Plot1")
df_betatd = merge(df_betatd, elev2, by = "Plot2")

# Delete records with 0 and NAs
df_betatd = subset(df_betatd, Distances != 0.0000000)
df_betatd = subset(df_betatd, Distances != "NA")

# Convert elevations columns into numeric types
df_betatd$Elevation <- as.numeric(as.character(df_betatd$Elevation1))
df_betatd$Elevation2 <- as.numeric(as.character(df_betatd$Elevation2))

# Calculate elevation change
df_betatd$Elev_change = abs(df_betatd$Elevation1 - df_betatd$Elevation2)

# Create the plot
tax_dis_elev = ggscatter(df_betatd, x = "Elev_change", y = "Distances",
                         ylab= "Taxonomic dissimilarities",
                         xlab = "Elevation change",
                         #ylim = c(0.3,1),
                         color = "black", shape = 20, size = 3,
                         add = "reg.line",
                         add.params = list(color = "firebrick3", linetype = "solid"))+ # Customize reg. line +
  stat_cor(label.x = 0, label.y = 1.1) +
  stat_regline_equation(label.x = 0, label.y = 1.05)

# 2) Functional dissimilarities among plots based on their altitudinal distance
# Add plot names as indices and delete column with names
row.names(df_betafd) <- df_betafd$X
df_betafd$X <- NULL
# Convert dataframe into a matrix
df_betafd = as.matrix(df_betafd)

# Create dataframe with pairwise distances
df_betafd = data.frame(col = colnames(df_betafd)[col(df_betafd)], row = rownames(df_betafd)[row(df_betafd)], 
                       dist = c(df_betafd))

colnames(df_betafd) = c("Plot1", "Plot2", "Distances")

# Add elevation columns of plot1 and 2
df_betafd = merge(df_betafd, elev1, by = "Plot1")
df_betafd = merge(df_betafd, elev2, by = "Plot2")

# Delete records with 0 and NAs
df_betafd = subset(df_betafd, Distances != 0.0000000)
df_betafd = subset(df_betafd, Distances != "NA")

# Convert elevations columns into numeric types
df_betafd$Elevation <- as.numeric(as.character(df_betafd$Elevation1))
df_betafd$Elevation2 <- as.numeric(as.character(df_betafd$Elevation2))

# Calculate elevation change
df_betafd$Elev_change = abs(df_betafd$Elevation1 - df_betafd$Elevation2)

# Create the plot
funct_dis_elev = ggscatter(df_betafd, x = "Elev_change", y = "Distances",
                           ylab= "Functional dissimilarities",
                           xlab = "Elevation change",
                           color = "black", shape = 20, size = 3,
                           add = "reg.line",
                           add.params = list(color = "firebrick3", linetype = "solid"))+ # Customize reg. line +
  stat_cor(label.x = 0, label.y = 1.1) +
  stat_regline_equation(label.x = 0, label.y = 1.05)

# Put the four plots in a single figure
final_plot = ggarrange(tax_dis_elev, funct_dis_elev,
                       ncol = 1, nrow = 2)

# Export final plot
ggsave(filename = "Outputs/Figure2_SI.png",
       plot = final_plot, 
       width = 13,
       height = 11,
       units = "in")



