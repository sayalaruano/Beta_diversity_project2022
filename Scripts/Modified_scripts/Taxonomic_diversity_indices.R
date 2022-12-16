# Script to calculate alpha and beta taxonomic diversity

# Import libraries
library(reshape2)
library(hillR)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(vegan)

# Load data
df_div <- read.csv("Data/2019_Div.csv", sep = ",")
elev <- read.csv("Data/Altitudes_13052021.csv", sep = ";")

# Add plot names as indices and delete column with names
row.names(elev) <- elev$Plot
elev[1] <- NULL

#######################################################
################ ALPHA DIVERSITY ######################
#######################################################

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

# Calculate hill numbers using the frequency matrix
ad_ab_hill <- hill_taxa(ad_ab_freqmat, q = 1, MARGIN = 1, base = exp(1))
ad_ab_hill <- as.data.frame(ad_ab_hill)
colnames(ad_ab_hill) <- c("Hill_Diversity")

# Add elevation column 
ad_ab_hill = merge(ad_ab_hill, elev, by = "row.names")

# Alpha diversity by basal area
ad_ba = subset(df_div, select = c("Plot", "Scientific_name", "BA"))

# Create a dataframe with sums of basal areas by plots
ad_ba_sum = aggregate(x = ad_ba$BA,
                  by = list(ad_ba$Plot, ad_ba$Scientific_name),
                  FUN = sum, na.rm = TRUE)

# Change column names
colnames(ad_ba_sum) <- c("Plot", "Scientific_name", "BAm2")

# Convert df into a matrix
ad_ba_summat = recast(ad_ba_sum, Plot ~ Scientific_name, id.var = c("Scientific_name", "Plot"))

# Delete frequency column without scientific names
ad_ba_summat[2] <- NULL

# Replace NAs by 0s
ad_ba_summat[is.na(ad_ba_summat)] <- 0

# Add plot names as indices and delete column with names
row.names(ad_ba_summat) <- ad_ba_summat$Plot
ad_ba_summat[1] <- NULL

# Calculate hill numbers using the frequency matrix
ad_ba_hill = hill_taxa(ad_ba_summat, q = 1, MARGIN = 1, base = exp(1))
ad_ba_hill <- as.data.frame(ad_ba_hill)
colnames(ad_ba_hill) <- c("Hill_Diversity")

# Obtain normalized values (calculating their square root) of basal area
ad_ba_sum$Raiz_BA <- sqrt(ad_ba_sum$BAm2)

ad_ba_sum_rs <- subset(ad_ba_sum, select = c("Scientific_name", "Plot", "Raiz_BA"))

# Convert df into a matrix
ad_ba_sum_rsmat <- recast(ad_ba_sum_rs, Plot ~ Scientific_name, id.var = c("Scientific_name", "Plot"))

# Delete frequency column without scientific names
ad_ba_sum_rsmat[2] <- NULL

# Replace NAs by 0s
ad_ba_sum_rsmat[is.na(ad_ba_sum_rsmat)] <- 0

# Add plot names as indices and delete column with names
row.names(ad_ba_sum_rsmat) <- ad_ba_sum_rsmat$Plot
ad_ba_sum_rsmat[1] <- NULL

# Calculate hill numbers using the basal area sum matrix
ad_ba_hill_rs <- hill_taxa(ad_ba_sum_rsmat, q = 1, MARGIN = 1, base = exp(1))
ad_ba_hill_rs <- as.data.frame(ad_ba_hill_rs)
colnames(ad_ba_hill_rs) <- c("Hill_Diversity")

# Add elevation column 
ad_ba_hill_rs <- merge(elev, ad_ba_hill_rs, by = "row.names")

# Add plot names as indices and delete column with names
row.names(ad_ba_hill_rs) <- ad_ba_hill_rs$Row.names
ad_ba_hill_rs[1] <- NULL

# Order by elevation 
ad_ba_hill_rs <- ad_ba_hill_rs[order(ad_ba_hill_rs$Elevation),]

#######################################################
################ BETA DIVERSITY ######################
#######################################################

bd_ba_rs <- betadiver(ad_ba_sum_rsmat, "w")
bd_ba_rs <- as.data.frame(as.matrix(bd_ba_rs))

# Add elevation column 
bd_ba_rs <- merge(elev, bd_ba_rs, by = "row.names")

# Add plot names as indices and delete column with names
row.names(bd_ba_rs) <- bd_ba_rs$Row.names
bd_ba_rs[1] <- NULL

# Order by elevation 
bd_ba_rs <- bd_ba_rs[order(bd_ba_rs$Elevation),]

# Reorder columns by elevation
bd_ba_rs <- bd_ba_rs[,c(1,11,10, 9, 8, 12, 13, 7, 6, 3, 5, 2, 4, 15, 16, 14, 17)]

# Delete elevation column
bd_ba_rs[1] <- NULL

# Save beta diversity matrix
write.csv(bd_ba_rs,"Data/Beta_tax_div_ba_rs_matrix.csv", row.names = TRUE)

#######################################################
################ RAREFACTION ##########################
#######################################################

raref <- as.data.frame(rarefy(ad_ab_freqmat, min(rowSums(ad_ab_freqmat))))
colnames(raref) <- c("Rarefaction")

# Add rarefaction column with previous dataset
final_df <- merge(ad_ba_hill_rs, raref, by = "row.names")

# Add plot names as indices and delete column with names
row.names(final_df) <- final_df$Row.names
final_df[1] <- NULL

# Create a dataframe with raref information
rarefcurve_ba <- specaccum(ad_ba_sum_rsmat)
df_raref = data.frame(rarefcurve_ba$sites, rarefcurve_ba$richness)
colnames(df_raref) = c('sites', 'richness')

# Create plot
raref_plot <- ggscatter(as.data.frame(df_raref), x = "sites", y = "richness",
                     color = "black", shape = 20, size = 2.5)+
  ylab("Estimated richness")+
  xlab("nplots")+
  geom_line(size = 0.6)

# Export final plot
ggsave(filename = "Outputs/Taxonomic_diversity/rarefaction_curve.png",
       plot = raref_plot, 
       width = 12,
       height = 8,
       units = "in")

#######################################################
################ EVENESS ##############################
#######################################################

# Calculate Shannon Weiner diversity 
H <- diversity(ad_ba_sum_rsmat)

# Calculate Pielou evenness
evenness <- as.data.frame(H/log(specnumber(ad_ba_sum_rsmat)))
colnames(evenness) <- c("Evenness")

# Add eveness column with previous dataset
final_df <- merge(final_df, evenness, by = "row.names")

# Add plot names as indices and delete column with names
row.names(final_df) <- final_df$Row.names
final_df[1] <- NULL

#######################################################
################ CHAO1 ################################
#######################################################
chao1 <- as.data.frame(t(estimateR(ad_ab_freqmat)))

chao1 <- subset(chao1, select = c("S.chao1", "se.chao1"))
colnames(chao1) <- c("Chao1", "SE_Chao1")

# Add chao1 column with previous dataset
final_df <- merge(final_df, chao1, by = "row.names")

# Add plot names as indices and delete column with names
row.names(final_df) <- final_df$Row.names
final_df[1] <- NULL

#######################################################
################ RICHNESS #############################
#######################################################

richness <- as.data.frame(specnumber(ad_ab_freqmat))
colnames(richness) <- c("Richness")

# Add richness column with previous dataset
final_df <- merge(final_df, richness, by = "row.names")

# Add plot names as indices and delete column with names
row.names(final_df) <- final_df$Row.names
final_df[1] <- NULL

# Express richness with rarefaction
Srar <- as.data.frame(rarefy(ad_ab_freqmat, min(rowSums(ad_ab_freqmat))))
Srar$Plot <- rownames(Srar)

#######################################################
############ EXPORT TABULAR DATA ######################
#######################################################

write.csv(final_df,"Outputs/Taxonomic_diversity/Taxon_div_indices.csv", row.names = TRUE)

#######################################################
################ NMDS ANALYSIS ########################
#######################################################

# Calculate the Nonmetric Multidimensional Scaling (NMDS) of the frequency 
# matrix to simplify multivariate data into a few important axes
nmds1 <- metaMDS(ad_ab_freqmat, distance = "bray", k = 2)

# Create dataset of NMDS from plots
nmds_plots = as.data.frame(nmds1$points)

# Create plot
nmds_sp <- ggscatter(as.data.frame(nmds1$species), x = "MDS1", y = "MDS2",
          color = "black", shape = 20, size = 3)+
          ylab("NMDS2")+
          xlab("NMDS1")+
          geom_point(data=nmds_plots, color="#AD1010")

# Export final plot
ggsave(filename = "Outputs/Taxonomic_diversity/nmds_species.png",
       plot = nmds_sp, 
       width = 12,
       height = 8,
       units = "in")

# Calculate distance matrix without species - considering plots
ad_ab_freqmat2 <- vegdist(ad_ab_freqmat, method = "bray")

# Convert vegan object into a matrix
ad_ab_freqmat2 <- as.matrix(ad_ab_freqmat2, labels = T)

# Calculate NMDS
nmds2 <- metaMDS(ad_ab_freqmat2, distance = "bray", k = 2, maxit = 999, trymax = 500,
                wascores = TRUE)

# Create dataset of NMDS from plots
nmds2_plots = as.data.frame(nmds2$points)

# Add elevation column 
nmds2_plots = merge(nmds2_plots, elev, by = "row.names")

colnames(nmds2_plots) = c("Plot", "NMDS1", "NMDS2", "Elevation")

# Create plot
palette_yelbrown = brewer.pal(n = 9, name = "YlOrBr")
nmds_plot <- ggscatter(nmds2_plots, x = "NMDS1", y = "NMDS2",
                     color = "Elevation", shape = 20, size = 3, 
                     label = "Plot")+
            theme(legend.position="right")+
            scale_colour_gradient2(low = "#662506", mid= "#CC4C02", high = "#FEE391")
            #scale_color_distiller(palette = "YlOrBr")

# Export final plot
ggsave(filename = "Outputs/Taxonomic_diversity/nmds_plots.png",
       plot = nmds_plot, 
       width = 12,
       height = 8,
       units = "in")
