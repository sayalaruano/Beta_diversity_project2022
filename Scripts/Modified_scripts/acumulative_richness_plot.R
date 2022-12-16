# Import libraies 
library(BiodiversityR) # also loads vegan
library(ggplot2)
library(ggsci)
library(reshape2)
library(dplyr)

# ggplot theme settings
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 10),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 12, colour = "gray25"),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 12),
  legend.key = element_blank(),
  legend.position="bottom")

# Load data
df_div <- read.csv("Data/2019_Diversidad.csv", sep = ",")
loc_names <- read.csv("Data/Location_names.csv", sep = ",")
elev <- read.csv("Data/Altitudes_13052021.csv", sep = ";")

# Add plot names as indices and delete column with names
row.names(loc_names) <- loc_names$Plot
loc_names[1] <- NULL

# Alpha diversity count by species abundance
# Create the frequency matrix
ad_ab_freq <- rename(count(df_div, Plot, Scientific_name), Freq = n)
ad_ab_freqmat <- recast(ad_ab_freq, Plot ~ Scientific_name, id.var = c("Scientific_name", "Plot"))

# Delete frequency column without scientific names
ad_ab_freqmat[2] <- NULL

# Replace NAs by 0s
ad_ab_freqmat[is.na(ad_ab_freqmat)] <- 0

# Add plot names as indices and delete column with names
row.names(ad_ab_freqmat) <- ad_ab_freqmat$Plot
ad_ab_freqmat[1] <- NULL

# Calculate sums of species by plot
loc_names$site.totals <- apply(ad_ab_freqmat,1,sum)

# test
accum_sp <- accumcomp(ad_ab_freqmat, y=loc_names, factor='Location', method='exact', 
                      legend=FALSE, conditioned=TRUE, scale='site.totals')

accumplot(accum_sp)

# Calculate accumulated species richness
#accum_sp <- accumcomp(ad_ab_freqmat, y=loc_names, factor='site.totals', 
 #                    method='exact', conditioned=FALSE, plotit=FALSE)

# Change format of accumlation data
accum_sp_long <- accumcomp.long(accum_sp, ci=NA, label.freq=5)

# Plot
# Create pallete
col_pallete <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
                 "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6")


accum_plot <- ggplot(data=accum_sp_long, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=2) +
  geom_point(data=subset(accum_sp_long, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=3)+
  xlim(100, 1600)+
  #geom_ribbon(aes(colour=Grouping, fill=after_scale(alpha(colour, 0.1))), 
   #           show.legend=FALSE) + 
  BioR.theme +
  scale_colour_npg()+
  labs(x = "Number of individuals", y = "Richness", colour = "Location", shape="Location")+
  scale_shape_manual(values = c(4,8,15,16,17,18,21,22,3))

# Export final plot
ggsave(filename = "Outputs/accumulative_richness.png",
       plot = accum_plot, 
       width = 12,
       height = 8,
       units = "in")
