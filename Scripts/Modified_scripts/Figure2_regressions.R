# Imports
library(ggpubr)
library(ggplot2)

 # Load data
df_indices = read.csv("Data/AlphFD_Indexvalues_270122.csv", sep = ",")
df_betafd = read.csv("Data/Beta_funct_div_matrix.csv", sep = ",")
df_betatd = read.csv("Data/Beta_tax_div_ba_rs_matrix.csv", sep = ",")

# Add plot names as indices and delete column with names
row.names(df_betatd) <- df_betatd$X
df_betatd[1] <- NULL

row.names(df_betafd) <- df_betafd$X
df_betafd[1] <- NULL

# a) Regression between alpha functional and taxonomic diversity (expressed with Hill numbers)
alp_func_tax = ggscatter(df_indices, x = "TD_Hill", y = "FD_Hills_q2",
                        color = "Elevation",
                        ylab = "Taxonomic diversity",
                        xlab = "Functional diversity",
                        label = "Plot",
                        shape = 20, size = 4.5,
                        add = "reg.line",
                        add.params = list(color = "#AD1010", linetype = "solid"))+ # Customize reg. line +
              theme(legend.position="right")+
              scale_colour_gradient2(low = "#662506", mid= "#CC4C02", high = "#FEE391")+
              stat_cor(label.x = 5, label.y = 5) +
              stat_regline_equation(label.x = 5, label.y = 4.8)

# b) Regression between taxonomic and functional dissimilarities

# Create the summary table
# Create a dataframe with elevation ranges
elev_change_num <- c(21, 174, 191, 259, 363, 
                     189, 50, 324, 9, 101, 
                     179, 440, 177, 312, 86)

elev_change_cat <- c('632-653', '653-827', '827-1018', '1018-1277', '1277-1640', 
                     '1640-1829', '1829-1879', '1879-2203', '2203-2212', '2212-2313', 
                     '2313-2492', '2492-2932', '2932-3109', '3109-3421', '3421-3507')

df_elevchangetd <- data.frame(elev_change_num, stringsAsFactors = TRUE)

# Obtain data from the observed values of tax div and put it into the previous dataframe
obs_val <- c()
idx1 <- 1
idx2 <- 2

for(i in 1:length(df_elevchangetd$elev_change_num)){
  obs_val <- c(obs_val, df_betatd[idx1,idx2])
  idx1 <- idx1 + 1
  idx2 <- idx2 + 1
}

df_elevchangetd$'Distances' <- obs_val

# Dafraframe for functional div elev change
df_elevchangefd <- data.frame(elev_change_num, stringsAsFactors = TRUE)

# Obtain data from the observed values of funct div and put it into the previous dataframe
obs_val <- c()
idx1 <- 1
idx2 <- 2

for(i in 1:length(df_elevchangefd$elev_change_num)){
  obs_val <- c(obs_val, df_betafd[idx1,idx2])
  idx1 <- idx1 + 1
  idx2 <- idx2 + 1
}

df_elevchangefd$'Distances' <- obs_val

# Merge datasets of taxonomic and functional div
# Add elevation columns of plot1 and 2
df_beta_tax_div = merge(df_elevchangefd, df_elevchangetd, by = "elev_change_num")
colnames(df_beta_tax_div) = c('Elevation', 'Funct', 'Tax')

# Create the plot
funct_dis_tax_dis = ggscatter(df_beta_tax_div, x = "Funct", y = "Tax",
                           color = "Elevation",
                           ylab= "Taxonomic dissimilarities",
                           xlab = "Functional dissimilarities",
                           #ylim = c(0.3,1),
                           shape = 20, size = 4.5,
                           #conf.int = TRUE,
                           add = "reg.line",
                           add.params = list(color = "#AD1010", linetype = "solid"))+ # Customize reg. line +
                  theme(legend.position="right")+
                  scale_colour_gradient2("Elevation change", low = "#662506", mid= "#CC4C02", high = "#FEE391")+
                  stat_cor(label.x = 0.23, label.y = 0.8) +
                  stat_regline_equation(label.x = 0.23, label.y = 0.78)


# Put the four plots in a single figure
final_plot = ggarrange(alp_func_tax, funct_dis_tax_dis,
          ncol = 2, nrow = 1)

# Export final plot
ggsave(filename = "Outputs/Figure2.png",
       plot = final_plot, 
       width = 22,
       height = 8,
       units = "in")

