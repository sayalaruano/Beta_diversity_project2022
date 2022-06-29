# Script to calculate the functional diversity indices

# You can find an extensive documentation of the mFD package in the following links:
# https://cmlmagneville.github.io/mFD/articles/mFD_general_workflow.html#know-your-data
# https://cran.r-project.org/web/packages/mFD/mFD.pdf

# Import libraries
library(mFD)
library(reshape2)
library(dplyr)

# Add a seed to obtain reproducible results. 
# If you want the results to vary, turn off set.seed()
set.seed(1234)

# Load data
df_traits_orig <- read.csv("Data/20210810_MeanTraits_FD.csv",stringsAsFactors = TRUE)
traits_cat <- read.csv("Data/20211108_MeanTraits_CAT.csv",stringsAsFactors = TRUE)

#######################################################
################ DATASET PREPARATION ##################
#######################################################

# Add scientific names and plots as indices
rownames(df_traits_orig) <- df_traits_orig$Scientific_name_plot

# Delete unnecessary columns 
df_traits <- df_traits_orig[,!colnames(df_traits_orig) %in% c('Plot','Scientific_name',
                                                   'Scientific_name_plot', 
                                                   'BA', 'RC_BA')]
# Normalize values 
df_traits <- df_traits %>% mutate_at(c("Ht2019","LA","SLA"), 
                                     ~(scale(.) %>% as.vector))

# Reshape the original dataframe of basal area by plot and scientific names columns
sp_BA_plot <- dcast(df_traits_orig, Plot ~ Scientific_name_plot, 
                    value.var = "RC_BA", sum)

# Add plots as indices of the matrix and delete the plots column
rownames(sp_BA_plot) <- sp_BA_plot$Plot
sp_BA_plot[1] <- NULL

# Convert dataframe into a matrix
sp_BA_plot <- as.matrix(sp_BA_plot)

#######################################################
########### FUNCTIONAL DIVERSITY ANALYSIS #############
#######################################################

# Summarizing data
# Species traits summary
traits_summ <- mFD::sp.tr.summary(
  tr_cat     = traits_cat,   
  sp_tr      = df_traits, 
  stop_if_NA = FALSE)

# Summary of the plot by species dataframe
asb_sp_traits_summ <- mFD::asb.sp.summary(asb_sp_w = sp_BA_plot)
asb_sp_traits_occ <- asb_sp_traits_summ$"asb_sp_occ"

# Compute distances between species based on functional traits
# This function returns a dist object with traits-based distances 
# between all pairs of species
sp_dist_traits <- mFD::funct.dist(
  sp_tr         = df_traits,
  tr_cat        = traits_cat,
  metric        = "euclidean",
  scale_euclid  = "range",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Compute multi-dimensional functional spaces and assess their quality
fspaces_quality_traits <- mFD::quality.fspaces(
  sp_dist             = sp_dist_traits,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")

df_fspaces_quality_traits <- fspaces_quality_traits$quality_fspaces

# Show the quality metrics of space a data frame gathering for each space 
# (in rows), values of quality metric(s) (in columns) 
round(fspaces_quality_traits$"quality_fspaces", 3)

# Create a plot about the quality of the selected functional spaces
# This function generates a figure with three panels (in rows) for each selected 
# functional space (in columns). Each column represents a functional space, the value 
# of the quality metric is written on the top of each column. The x-axis of all panels 
# represents trait-based distances. The y-axis is different for each row
# More info: https://cmlmagneville.github.io/mFD/articles/mFD_general_workflow.html#illustrating-the-quality-of-the-selected-functional-spaces
funct_sp_q_plot <- mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_traits,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d"), 
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")

# Test correlation between functional axes and traits
sp_faxes_coord_traits <- fspaces_quality_traits$"details_fspaces"$"sp_pc_coord"

Traits_faxes <- mFD::traits.faxes.cor(
  sp_tr          = df_traits, 
  sp_faxes_coord = sp_faxes_coord_traits[ , c("PC1", "PC2", "PC3")], 
  plot           = TRUE)

df_corr_traits_axes <- Traits_faxes$tr_faxes_stat

# Store traits with significant effect
signif_traits <- Traits_faxes$"tr_faxes_stat"[which(Traits_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]

# Create the functional space plot - PCoA axis
big_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord_traits[ , c("PC1", "PC2", "PC3")],
  faxes           = c("PC1", "PC2", "PC3"),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "darkgreen",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "blueviolet",
  fill_vert       = "blueviolet",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)

# Visualize the graph with all pairs of axes
big_plot$patchwork

#######################################################
########### ALPHA DIVERSITY ###########################
#######################################################

# Calculate functional alpha diversity indices in a multidimensional space
alpha_fd_indices_traits <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_traits[ , c("PC1", "PC2", "PC3")],
  asb_sp_w         = sp_BA_plot,
  ind_vect         = c("fdis", "feve", "fric", "fdiv"),
  scaling          = FALSE,
  check_input      = TRUE,
  details_returned = TRUE)

df_fd_ind_values_traits <- alpha_fd_indices_traits$"functional_diversity_indices"

details_list_traits <- alpha_fd_indices_traits$"details"

# Create a plot of the alpha functional indices for 2 assemblages (the highest 
# and lowest plots by elevation) in a multidimensional space 
plots_alpha <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_traits,
  plot_asb_nm              = c("MAPI_02", "YANA_01"),
  ind_nm                   = c("fdis", "feve", "fric", "fdiv"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  color_vert               = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_sp                  = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_vert                = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  color_ch                 = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_ch                  = c(pool = "white", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 

# Visualize plots
plots_alpha$"fric"$"patchwork"
plots_alpha$"fdis"$"patchwork"
plots_alpha$"feve"$"patchwork"
plots_alpha$"fdiv"$"patchwork"

# Generalization of Hill numbers for alpha functional diversity
Traits_gower <- mFD::funct.dist(
  sp_tr         = df_traits,
  tr_cat        = traits_cat,
  metric        = "euclidean",
  scale_euclid  = "noscale",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

traits_FD2mean <- mFD::alpha.fd.hill(
  asb_sp_w = sp_BA_plot, 
  sp_dist  = Traits_gower, 
  tau      = "mean", 
  q        = 1)

df_fd_ind_values_traits["Hill_numbers"] <- traits_FD2mean$asb_FD_Hill

#######################################################
########### BETA DIVERSITY ############################
#######################################################

# Calculate functional beta diversity indices in a multidimensional space
beta_fd_indices_traits <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_traits[ , c("PC1", "PC2", "PC3")],
  asb_sp_occ       = asb_sp_traits_occ,
  check_input      = TRUE,
  beta_family      = c("Sorensen"),
  details_returned = TRUE)

# Obtain beta div indices values for each pair of assemblages as dataframe
df_beta_div_indices <- dist.to.df(beta_fd_indices_traits$pairasb_fbd_indices)

# Obtain a list containing details such as inputs, vertices of the global 
# pool and of each assemblage and FRic values for each assemblage
beta_div_details <- beta_fd_indices_traits$"details"

# Obtain a vector containing the FRic value for each assemblage 
# retrieved through the details_beta list
beta_div_fric <- as_data_frame(beta_fd_indices_traits$"details"$"asb_FRic")
colnames(beta_div_fric) <- "FRic_beta_div"

# Obtain a list of vectors containing names of species being vertices of the 
# convex hull for each assemblage retrieved through the details_beta list
beta_div_scient_names <- beta_fd_indices_traits$"details"$"asb_vertices"

# Create a plot of the beta functional indices for 2 assemblages
# in a multidimensional space 
beta_plot_traits <- mFD::beta.multidim.plot(
  output_beta_fd_multidim = beta_fd_indices_traits,
  plot_asb_nm             = c("CEDR_01", "RIBR_01"),
  beta_family             = c("Sorensen"),
  faxes                   = paste0("PC", 1:3),
  name_file               = NULL,
  faxes_nm                = NULL,
  range_faxes             = c(NA, NA),
  color_bg                = "grey95",
  shape_sp                = c("pool" = 3.0, asb1 = 22, asb2 = 21),
  size_sp                 = c("pool" = 0.8, asb1 =  1, asb2 =  1),
  color_sp                = c("pool" = "grey50", asb1 = "blue", asb2 = "red"),
  fill_sp                 = c("pool" = NA, asb1 = "white", asb2 = "white"),
  fill_vert               = c("pool" = NA, asb1 = "blue", asb2 = "red"),
  color_ch                = c("pool" = NA, asb1 = "blue", asb2 = "red"),
  fill_ch                 = c("pool" = "white", asb1 = "blue", asb2 = "red"),
  alpha_ch                = c("pool" = 1, asb1 = 0.3, asb2 = 0.3),
  nm_size                 = 3,
  nm_color                = "black",
  nm_fontface             = "plain",
  check_input             = TRUE)

# Visualize plots
beta_plot_traits$"patchwork"

#######################################################
############ EXPORT TABULAR DATA ######################
#######################################################

# Create a list with the results
results <- list(funct_spaces_quality = df_fspaces_quality_traits, 
              corr_traits_axes = df_corr_traits_axes, 
              alpha_div_indices =  df_fd_ind_values_traits,
              beta_div_indices = df_beta_div_indices, 
              beta_div_frich = beta_div_fric)

# Create a folder to store the results. You should replace the name of the new 
# folder_path to store the results
newfolder_path <- "Outputs/Functional_diversity_analysis/"
dir.create(newfolder_path)

# Export the results as csv files
for (j in 1:length(results)){
  # Replace the name of the 
  write.csv(results[[j]],paste0(newfolder_path, "/", names(results)[j],'.csv'))
  
}



