# Script to calculate alpha and beta taxonomic diversity

# Import libraries
library(reshape2)
library(hillR)
library(dplyr)
library(ggplot2)
library(MBI)
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
ad_ba_hill_rs <- merge(ad_ba_hill_rs, elev, by = "row.names")

#######################################################
################ BETA DIVERSITY ######################
#######################################################

bd_ba_rs <- betadiver(ad_ba_sum_rsmat, "w")
bd_ba_rs <- as.data.frame(bd_ba_rs)

#######################################################
################ RAREFACTION ##########################
#######################################################

raref <- as.data.frame(rarefy(ad_ab_freqmat, min(rowSums(ad_ab_freqmat))))
colnames(raref) <- c("Rarefaction")

rarefcurve_ba <- specaccum(ad_ba_sum_rsmat)

plot(rarefcurve_ba)

S <- as.data.frame(specnumber(ad_ab_freqmat))

#######################################################
################ EVENESS ##############################
#######################################################

# Calculate Shannon Weiner diversity 
H <- diversity(ad_ba_sum_rsmat)

# Calculate Pielou evenness
evenness <- as.data.frame(H/log(specnumber(ad_ba_sum_rsmat)))
colnames(evenness) <- c("Evenness")

#######################################################
################ CHAO1 ################################
#######################################################
chao1 <- as.data.frame(estimateR(ad_ab_freqmat))

#######################################################
################ RICHNESS #############################
#######################################################

# Express richness with rarefaction
Srar <- as.data.frame(rarefy(ad_ab_freqmat, min(rowSums(ad_ab_freqmat))))
Srar$Plot <- rownames(Srar)
rownames(Srar) <- NULL
names(Srar)[1] = "Riqueza"
