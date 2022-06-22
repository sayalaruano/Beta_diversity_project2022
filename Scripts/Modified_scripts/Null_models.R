# Script to run the null models

# Import libraries
library(vegan)
library(tidyr)
library(ggplot2)
library(purrr)
library(boot)
library(Rmisc)

# Define a function to run the null models
calc_zscore <- function(file_to_read,file_elevation,spr.t.f,method.null,
                        meth.diss,n.simulations,n.year){
  
    # Add a seed to obtain reproducible results. 
    # If you want the results to vary, turn off set.seed()
    set.seed(953)

    # Load data 
    ba_ind <- read.csv(file_to_read)
    elev <- read.csv(file_elevation, sep = ";")


    # Validate if the dataset includes information about years
    if (isTRUE(spr.t.f)){
      year<-n.year
      ba_year <- as.data.frame(ba_ind[,c('Plot',"Scientific_name",year)])
      ba_spread <- spread(ba_year,Scientific_name,year)
    } else {
      ba_spread <- ba_ind
    }
  
    # Replace NAs by 0s
    ba_spread[is.na(ba_spread)] <- 0
  
    # Merge elevation data
    ba_spread <- merge(ba_spread,elev,by=c('Plot'))
    
    # Add plot names as indices
    row.names(ba_spread) <- ba_spread$Elevation
    
    # Delete plot and elevation columns
    ba_spread <- ba_spread[,!colnames(ba_spread) %in% c('Plot','Elevation')]
    
    # Order the dataframe by row names (elevation)
    ba_spread <- ba_spread[order(as.numeric(row.names(ba_spread))),]
  
    # Calculate the dissimilarity index of the data (jaccard, bray curtis, or other)
    dist.o <- vegdist(ba_spread, method=meth.diss, binary=FALSE)
  
    # List with the numbers of beta diversity indices to be calculated
    lindex <- c(1,2,8,10,11,13,15,22)
  
    # Create empty list to store beta diversity indices results 
    list_index <- list()
  
    # Calculate beta diversity indices
    for (i in 1:length(lindex)){
      b1 <- betadiver(ba_spread,method=i)
      list_index[[i]]<-as.matrix(b1)
    }
  
    # Null models
    # Define a signature for the commsim function 
    f <- function(x, n, ...)  array(replicate(n, sample(x)), c(dim(x), n))
    
    # Create an object of the commsim class to run a null model
    cso <- commsim(method.null, fun=f, binary=FALSE, isSeq=FALSE, mode="double")
    
    # Creates an object of the nullmodel class that can serve as a basis for 
    # Null Model simulation via the simulate method
    qw1a <- nullmodel(ba_spread, cso)
    
    #qd1a <- nullmodel(dist.o, cso)
  
    #qw1a <- nullmodel(ba_spread, cso)
    
    # Run the null models with the simulate function 
    sqw1a <- simulate(qw1a, nsim=n.simulations)
    #sqd1a <- simulate(qd1a, nsim=nsimm)
  
    # Calculate z-scores of betadiversity
    # Create empty lists 
    lbig1 <- list()
    ldist <- list()
    
    # Run null models considering the number of simulations 
    # and previous parameters
    for (j in 1:n.simulations){
      
      nm <- as.data.frame(t(sqw1a[1,,j]))
      #nd <- as.data.frame(t(sqd1a[1,,j]))
      
      for(i in 2:nrow(ba_spread)){
        nm[i,] <- as.data.frame(t(sqw1a[i,,j]))
        #nd[i,] <- as.data.frame(t(sqd1a[i,,j]))
      }
      row.names(nm) <- row.names(ba_spread)
      #row.names(nd) <- row.names(dist.o)
      
      bt1 <- betadiver(nm,method=1) ##1 whittaker
      dist.a <- vegdist(nm, method=meth.diss, binary=FALSE)
    
      lbig1[[j]] <- as.matrix(bt1) ##$statistic[3]
      ldist[[j]] <- as.matrix(dist.a)
      
    }
    
    # Reduce lists into single values by iteratively summing the numbers of lists
    v <- reduce(lbig1, `+`)
    dd <- reduce(ldist, `+`)
    
    # Normalize the values of the matrices
    v <- as.matrix(v)/length(lbig1)
    dd <- as.matrix(dd)/length(ldist)
  
  ##inlcuir loop el modelo nulo se hace en base a la matriz de disimilitudes bray curtis.
    
    # Create empty matrices to store parameters of the null models
    sds <- matrix(0,n.simulations,1)
    sds1 <- matrix(0,nrow(ba_spread),nrow(ba_spread))
    
    sdsd <- matrix(0,n.simulations,1)
    sdsd1 <- matrix(0,nrow(ba_spread),nrow(ba_spread))
    
    ci1 <- matrix(0,nrow(ba_spread),nrow(ba_spread))
    ci2 <- matrix(0,nrow(ba_spread),nrow(ba_spread))
    
    cid1 <- matrix(0,nrow(ba_spread),nrow(ba_spread))
    cid2 <- matrix(0,nrow(ba_spread),nrow(ba_spread))
    
    # Fill empty matrices with parameters of null models
    for (j in 1:nrow(as.matrix(lbig1[[1]]))){
      for (k in 1: ncol(as.matrix(lbig1[[1]]))){
        for(i in 1:n.simulations){
          # Fill matrices with null models coefficients 
          sds[i]<-as.matrix(lbig1[[i]])[j,k]
          sdsd[i]<-as.matrix(ldist[[i]])[j,k]
        }
        # Fill matrices with null models standard deviations
        sds1[j,k]<-sd(sds)
        sdsd1[j,k]<-sd(sdsd)
        
        # Fill matrices with null models upper and lower confidence intervals 
        ci1[j,k]<-CI(sds,ci=0.95)[[1]]
        ci2[j,k]<-CI(sds,ci=0.95)[[3]]
        
        cid1[j,k]<-CI(sdsd,ci=0.95)[[1]]
        cid2[j,k]<-CI(sdsd,ci=0.95)[[3]]
      }
    }
  
    # Calculate z-scores of the Whittaker index
    zsc <- as.matrix(list_index[[1]]) - v
    zscsws500 <- zsc/sds1
    
    # Calculate z-scores of the Bray-curtis index
    zsd <- as.matrix(dist.o) - dd
    zsds1000 <- zsd/sdsd1
    
    # Create a list with the results
    lexit <- list(obs_whitaker = list_index[[1]], obs_diss = ldist[[1]], mean_null_whitaker = v, mean_null_diss = dd,
                 zscore_whitaker = zscsws500, zscore_diss =  zsds1000, ci_low_whitaker = ci2, ci_high_whitaker = ci1,
                 ci_low_diss = cid2, ci_high_diss = cid1)
    
    return (lexit)
}

# Run the function with a dataset (replace the name of the file in the 
# file_to_read argument)
results <- calc_zscore(file_to_read = 'Data/17082021_BaseModelosNulos_BasalArea.csv', file_elevation = 'Data/Altitudes_13052021.csv', spr.t.f = FALSE,
                    method.null = 'c0_both', meth.diss = 'bray',n.simulations = 1000, n.year = NA)

# Create a folder to store the results. You should replace the name of the new 
# folder_path to store the results
newfolder_path <- "Outputs/Null_models/Species_Bray"
dir.create(newfolder_path)

# Export the results as csv files
for (j in 1:length(results)){
  # Replace the name of the 
  write.csv(results[[j]],paste0(newfolder_path, "/", names(results)[j],'_species_bray.csv'))
  
}
