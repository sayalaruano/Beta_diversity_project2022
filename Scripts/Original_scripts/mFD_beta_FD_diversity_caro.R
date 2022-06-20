setwd("C:/Users/antonella.bernardi/OneDrive - Universidad de Las Américas/Francisco Cuesta/paper")
library(mFD)
library(reshape2)
library(dplyr)

fruits_traits
fruits_traits_cat
baskets_fruits_weights

#https://cmlmagneville.github.io/mFD/articles/mFD_general_workflow.html#know-your-data
#https://cran.r-project.org/web/packages/mFD/mFD.pdf


set.seed(1234)


####### Load data ###############

#### traits
  traits0 <- read.csv("20210810_MeanTraits_FD.csv",stringsAsFactors = TRUE)
  head(traits0)
  rownames(traits0) <- traits0$Scientific_name_plot
  traits<-traits0[,c("Scientific_name_plot","Ht2019","LA","SLA")]
  traits1 <- c(traits,"Scientific_name_plot")
  traits <- traits[,-which(names(traits) %in% traits1)]
  
  traits2 <- traits %>% mutate_at(c("Ht2019","LA","SLA"), ~(scale(.) %>% as.vector))
  head(traits2)
  
# Display the table:
  knitr::kable(head(traits2),
               caption = "Species x traits data frame")

#### site x sp data
  sp_BA_plot<-dcast(traits0, Plot ~ Scientific_name_plot, value.var = "RC_BA", sum)
  sp_BA_plot1<-as.matrix(sp_BA_plot)
  rownames(sp_BA_plot)<-sp_BA_plot$Plot
  head(sp_BA_plot)
  del1 <- c(sp_BA_plot,"Plot")
  sp_BA_plot1 <- as.matrix(sp_BA_plot[,-which(names(sp_BA_plot) %in% del1)])
  
  # Display the table:
  knitr::kable(as.data.frame(sp_BA_plot1[1:16, 1:16]), 
               centering = TRUE,
               caption = "Species x plot matrix based on basal area")

#### trait data type  
  traits_cat <- read.csv("20211108_MeanTraits_CAT.csv",stringsAsFactors = TRUE)
  head(traits_cat)
  str(traits_cat)
  
  # Display the table:
  knitr::kable(head(traits_cat), 
               caption = "Traits types")

######## Summarising data #######

  # Species traits summary:
traits_summ <- mFD::sp.tr.summary(
  tr_cat     = traits_cat,   
  sp_tr      = traits2, 
  stop_if_NA = FALSE)

# Summary of the plot * species dataframe:
asb_sp_traits_summ <- mFD::asb.sp.summary(asb_sp_w = sp_BA_plot1)
asb_sp_traits_occ <- asb_sp_traits_summ$"asb_sp_occ"

##### Computing distances between species based on functional traits ####
sp_dist_traits <- mFD::funct.dist(
  sp_tr         = traits2,
  tr_cat        = traits_cat,
  metric        = "euclidean",
  scale_euclid  = "range",
  weight_type   = "equal",
  stop_if_NA    = TRUE)
#This function returns a dist object with traits-based distances between all pairs of species:

round(sp_dist_traits, 3)                 # Output of the function mFD::funct.dist()


# Compute multimensional functional spaces and assess their quality
fspaces_quality_traits <- mFD::quality.fspaces(
  sp_dist             = sp_dist_traits,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")

round(fspaces_quality_traits$"quality_fspaces", 3)  # Quality metrics of space a data frame gathering for each space (in rows), values of quality metric(s) (in columns) 

#Illustrating the quality of the selected functional spaces
mFD::quality.fspaces.plot(
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

#This function generates a figure with three panels (in rows) for each selected 
#functional space (in columns). Each column represents a functional space, the value 
#of the quality metric is written on the top of each column. The x-axis of all panels 
#represents trait-based distances. The y-axis is different for each row:

##on the first (top) row, the y-axis represents species functional distances in 
#the multidimensional space. Thus, the closer species are to the 1:1 line, the 
#better distances in the functional space fit trait-based ones.
##on the second row, the y-axis shows the raw deviation of species distances in 
#the functional space compared to trait-based distances. Thus, the raw deviation 
#reflects the distance to the 1:1 line.
##on the third row, the y-axis shows the absolute or squared deviation of the 
#("scaled") distance in the functional space. It is the deviation that is taken 
#into account for computing the quality metric.

##### Test correlation between functional axes and traits ###########

sp_faxes_coord_traits <- fspaces_quality_traits$"details_fspaces"$"sp_pc_coord"

Traits_faxes <- mFD::traits.faxes.cor(
  sp_tr          = traits2, 
  sp_faxes_coord = sp_faxes_coord_traits[ , c("PC1", "PC2", "PC3")], 
  plot           = TRUE)

# Print traits with significant effect:
Traits_faxes$"tr_faxes_stat"[which(Traits_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]
## Return plots:
Traits_faxes$"tr_faxes_plot"

#Plot functional space
#Here are the plots the PCoA axis:
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


# Plot the graph with all pairs of axes:
big_plot$patchwork

#Functional alpha diversity indices in a multidimensional space
alpha_fd_indices_traits <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_traits[ , c("PC1", "PC2", "PC3")],
  asb_sp_w         = sp_BA_plot1,
  ind_vect         = c("fdis", "feve", "fric", "fdiv"),
  scaling          = FALSE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values_traits <- alpha_fd_indices_traits$"functional_diversity_indices"
fd_ind_values_traits

details_list_traits <- alpha_fd_indices_traits$"details"

#plot functional indices for 2 assemblages
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

plots_alpha$"fric"$"patchwork"
plots_alpha$"fdis"$"patchwork"

# Generalisation of Hill numbers for alpha functional diversity
Traits_gower <- mFD::funct.dist(
  sp_tr         = traits2,
  tr_cat        = traits_cat,
  metric        = "euclidean",
  scale_euclid  = "noscale",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

traits_FD2mean <- mFD::alpha.fd.hill(
  asb_sp_w = sp_BA_plot1, 
  sp_dist  = Traits_gower, 
  tau      = "mean", 
  q        = 1)


# Functional beta diversity indices based on multidimensional space

beta_fd_indices_traits <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_traits[ , c("PC1", "PC2", "PC3")],
  asb_sp_occ       = asb_sp_traits_occ,
  check_input      = TRUE,
  beta_family      = c("Sorensen"),
  details_returned = TRUE)

# returns a dist object with beta indices values for each pair of assemblages:
head(beta_fd_indices_traits$"pairasb_fbd_indices", 10)

#returns a list containing details such as inputs, vertices of the global pool and of each assemblage and FRic values for each assemblage
beta_fd_indices_traits$"details"

#returns a vector containing the FRic value for each assemblage retrieved through the details_beta list:
beta_fd_indices_traits$"details"$"asb_FRic"

# returns a list of vectors containing names of species being vertices of the convex hull for each assemblage retrieved through the details_beta list:
beta_fd_indices_traits$"details"$"asb_vertices"

# illustrate functional beta-diversity indices for a pair of assemblages in a multidimensional space 
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

beta_plot_traits$"patchwork"


