# Imports
library(ggpubr)
library(ggplot2)

# Load data
df_indices = read.csv("Data/AlphFD_Indexvalues_270122.csv", sep = ",")
df_betafd = read.csv("Data/20012022_BetaFD_Stand_NoRC.csv", sep = ",")
df_betatd = read.csv("Data/Beta_funct_div_matrix.csv", sep = ",")
elev1 <- read.csv("Data/Altitudes_13052021.csv", sep = ";")

# Add plot names as indices and delete column with names
row.names(df_betatd) <- df_betatd$X
df_betatd[1] <- NULL

row.names(df_betafd) <- df_betafd$X
df_betafd[1] <- NULL

# a) Regression between alpha diversity (expressed with Hill numbers) and elevation
alp_vs_elev = ggscatter(df_indices, x = "Elevation", y = "TD_Hill",
                        ylab= "Taxonomic Hill diversity",
                      color = "black", shape = 20, size = 3,
                      conf.int = TRUE,
                      add = "reg.line",
                      add.params = list(color = "firebrick3", linetype = "solid"))+ # Customize reg. line +
            stat_cor(label.x = 2750, label.y = 75) +
            stat_regline_equation(label.x = 2800, label.y = 80)


# b) Regression between functional diversity (expressed with Hill numbers) and elevation
funct_vs_elev = ggscatter(df_indices, x = "Elevation", y = "FD_Hills_q2",
                        ylim  = c(0,6),
                        ylab= "Functional Hill diversity",
                        color = "black", shape = 20, size = 3,
                        conf.int = TRUE,
                        add = "reg.line",
                        add.params = list(color = "firebrick3", linetype = "solid"))+ # Customize reg. line +
                stat_cor(label.x = 2550, label.y = 5.6) +
                stat_regline_equation(label.x = 2600, label.y = 6)

# c) Taxonomic dissimilarities among plots based on their altitudinal distance
# Create the summary table
# Create a dataframe with elevation ranges
elev_change_num <- c(21, 174, 191, 259, 363, 
                 189, 50, 324, 9, 101, 
                 179, 440, 177, 312, 86)

elev_change_cat <- c('632-653', '653-827', '827-1018', '1018-1277', '1277-1640', 
                     '1640-1829', '1829-1879', '1879-2203', '2203-2212', '2212-2313', 
                     '2313-2492', '2492-2932', '2932-3109', '3109-3421', '3421-3507')

df_elevchangetd <- data.frame(elev_change_num, stringsAsFactors = TRUE)

# Obtain data from the observed values and put it into the previous dataframe
obs_val <- c()
idx1 <- 1
idx2 <- 2

for(i in 1:length(df_elevchangetd$elev_change_num)){
  obs_val <- c(obs_val, df_betatd[idx1,idx2])
  idx1 <- idx1 + 1
  idx2 <- idx2 + 1
}

df_elevchangetd$'Distances' <- obs_val

# Create the plot
tax_dis_elev <- ggscatter(df_elevchangetd, x = "elev_change_num", y = "Distances",
                         ylab= "Taxonomic dissimilarities",
                         xlab = "Elevation change",
                         #ylim = c(0.3,1),
                         color = "black", shape = 20, size = 3,
                         conf.int = TRUE,
                         add = "reg.line",
                         add.params = list(color = "firebrick3", linetype = "solid"))+ # Customize reg. line +
stat_cor(label.x = 325, label.y = 1.1) +
stat_regline_equation(label.x = 325, label.y = 1.05)

# d) Functional dissimilarities among plots based on their altitudinal distance
df_elevchangefd <- data.frame(elev_change_num, stringsAsFactors = TRUE)

# Obtain data from the observed values and put it into the previous dataframe
obs_val <- c()
idx1 <- 1
idx2 <- 2

for(i in 1:length(df_elevchangefd$elev_change_num)){
  obs_val <- c(obs_val, df_betafd[idx1,idx2])
  idx1 <- idx1 + 1
  idx2 <- idx2 + 1
}

df_elevchangefd$'Distances' <- obs_val

# Create the plot
funct_dis_elev = ggscatter(df_elevchangefd, x = "elev_change_num", y = "Distances",
                         ylab= "Functional dissimilarities",
                         xlab = "Elevation change",
                         #ylim = c(0.3,1),
                         color = "black", shape = 20, size = 3,
                         conf.int = TRUE,
                         add = "reg.line",
                         add.params = list(color = "firebrick3", linetype = "solid"))+ # Customize reg. line +
  stat_cor(label.x = 325, label.y = 1.1) +
  stat_regline_equation(label.x = 325, label.y = 1.05)

# Put the four plots in a single figure
final_plot = ggarrange(alp_vs_elev, funct_vs_elev, 
          tax_dis_elev, funct_dis_elev,
          ncol = 2, nrow = 2)

# Export final plot
ggsave(filename = "Outputs/Figure1.png",
       plot = final_plot, 
       width = 13,
       height = 11,
       units = "in")

