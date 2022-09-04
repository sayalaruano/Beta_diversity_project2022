# Imports
library(ggpubr)
library(ggplot2)

 # Load data
df_indices = read.csv("Data/AlphFD_Indexvalues_270122.csv", sep = ",")

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
              scale_color_gradient(low = "#A2DCEB", high = "#415E63") +
              stat_cor(label.x = 5, label.y = 5) +
              stat_regline_equation(label.x = 5, label.y = 4.8)

# b) Regression between taxonomic and functional dissimilarities


# Put the four plots in a single figure
final_plot = ggarrange(alp_vs_elev, funct_vs_elev, 
          tax_dis_elev, funct_dis_elev,
          ncol = 2, nrow = 2)

# Export final plot
ggsave(filename = "Outputs/Figure2.png",
       plot = alp_func_tax, 
       width = 12,
       height = 8,
       units = "in")

