# Import libraries
library(mFD)

# Load data
df_betafd = read.csv("Data/Functional_traits/20210811_MeanTraits.csv", sep = ",")
df_traits_cat = read.csv("Data/Functional_traits/MeanTraits_CAT_matrix.csv", sep = ",")
elev <- read.csv("Data/Altitudes_13052021.csv", sep = ";")

# Select columns 
df_betadf_filt <- df_betafd[,colnames(df_betafd) %in% c('Plot','Scientific_name',
                                                              'Ht2019_m', 'Mean.LBT_mm', 
                                                              'Mean.LA_cm2', 'Mean.SLA_cm2_g',
                                                              'Mean.LDMC_mg_g')]
# Normalize values 
df_betadf_filt <- df_betadf_filt %>% mutate_at(c("Ht2019_m","Mean.LBT_mm","Mean.LA_cm2",
                                       "Mean.SLA_cm2_g", "Mean.LDMC_mg_g"), 
                                     ~(scale(.) %>% as.vector))

# Create a compound key with scientific name and plot
df_betadf_filt$Scientific_name_Plot <- paste(df_betadf_filt$Scientific_name, df_betadf_filt$Plot)

# Create the frequency matrix
bd_freq <- rename(count(df_betadf_filt, Plot, Scientific_name_Plot), Freq = n)
bd_freqmat <- recast(bd_freq, Plot ~ Scientific_name_Plot, id.var = c("Scientific_name_Plot", "Plot"))

# Replace NAs by 0s
bd_freqmat[is.na(bd_freqmat)] <- 0

# Add plot names as indices and delete column with names
row.names(bd_freqmat) <- bd_freqmat$Plot
bd_freqmat[1] <- NULL

# Convert df into a matrix
bd_freqmat <- data.matrix(bd_freqmat)

# Create df to calculate distances
df_traits <- df_betadf_filt

# Add scientific and plot names as indices and delete column with names
row.names(df_traits) <- df_traits$Scientific_name_Plot

# Delete additional columns
df_traits <- df_traits[,!colnames(df_betadf_filt) %in% c('Plot', 'Scientific_name',
                                                              'Scientific_name_Plot')]
# Replace NAs by 0s
df_traits[is.na(df_traits)] <- 0

# Compute distances between species based on functional traits
# This function returns a dist object with traits-based distances 
# between all pairs of species
sp_dist_traits <- mFD::funct.dist(
  sp_tr         = df_traits,
  tr_cat        = df_traits_cat,
  metric        = "euclidean",
  scale_euclid  = "range",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Compute beta functional hill indices:
bd_plots <- beta.fd.hill(
  asb_sp_w         = bd_freqmat, 
  sp_dist          = sp_dist_traits, 
  q                = c(0,1,2), 
  tau              = 'mean',
  beta_type        = 'Jaccard', 
  check_input      = TRUE, 
  details_returned = TRUE)

# Convert matrix into a dataframe
beta_funct_div_plots <- data.frame(as.matrix(bd_plots$beta_fd_q$q0))

# Save beta diversity matrix
write.csv(beta_funct_div_plots,"Data/Beta_funct_div_matrix.csv", row.names = TRUE)
