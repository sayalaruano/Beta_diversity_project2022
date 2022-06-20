
library(vegan)
library(tidyr)
library(ggplot2)
library(purrr)
library(boot)
library(Rmisc)

setwd("C:/Users/anton/Documents/NOROCCIDENTE/paper/Outputs_NullModels_AreaBasal/Bases")

calc_zscore<-function(file_to_read,file_elevation,spr.t.f,method.null,meth.diss,n.simulations,n.year){
  
#turn off set.seed() if you want the results to vary
set.seed(953)

#ba_plot <- read.csv('ba_plot.csv')
#ba_ind <- read.csv('ba_ind.csv')
#ba_ind <- read.csv('trich.csv')
#ba_ind <- read.csv('trich_species.csv')

  ba_ind <- read.csv(file_to_read)
  #elev <- read.csv('elevation.csv')
  #elev <- read.csv('elev_blanca.csv')
  elev <- read.csv(file_elevation)
  
#row.names(ba_plot)<-ba_plot$Plot
#ba_plot<- ba_plot[,!colnames(ba_plot) %in% 'Plot']
#ba_ba <- spread(ba_plot, key=SN, value=Ba)

#library(betapart)

#plot.s.core <- betapart.core(ba_plot)


############ind
  if (isTRUE(spr.t.f)){
    #year<-'nind_2015'
    year<-n.year
    ba_year <- as.data.frame(ba_ind[,c('Plot',"Scientific_name",year)])
    
    #ba_year$Plot <- ba_ind$Plot
    
    ba_spread <- spread(ba_year,Scientific_name,year)
  } else {
    ba_spread <- ba_ind
    View(ba_spread)
  }
  
  

  ba_spread[is.na(ba_spread)] <- 0
  
  
  ba_spread <- merge(ba_spread,elev,by=c('Plot'))
  row.names(ba_spread)<-ba_spread$Altitud
  
  ba_spread <- ba_spread[order(ba_spread$Altitud),]
  
  ba_spread<- ba_spread[,!colnames(ba_spread) %in% c('Plot','Altitud')]
  
  dist.o <- vegdist(ba_spread, method=meth.diss, binary=FALSE)
  
  lindex <- c(1,2,8,10,11,13,15,22)
  
  list_index <- list()
  
  for (i in 1:length(lindex)){
    b1 <- betadiver(ba_spread,method=i)
    list_index[[i]]<-as.matrix(b1)
  }
  
  
  f <- function(x, n, ...)  array(replicate(n, sample(x)), c(dim(x), n))
  
  #cso <- commsim("c0_both", fun=f, binary=FALSE,
  cso <- commsim(method.null, fun=f, binary=FALSE,
                isSeq=FALSE, mode="integer")
  
  qw1a <- nullmodel(ba_spread, cso)
  #qd1a <- nullmodel(dist.o, cso)
  
  #qw1a <- nullmodel(ba_spread, cso)
  
  nsimm <- n.simulations
  
  sqw1a <- simulate(qw1a, nsim=nsimm)
  #sqd1a <- simulate(qd1a, nsim=nsimm)
  
  ##volver a calcular la diversidad beta
  
  # tnodf <- nestednodf(ba_spread)
  
  #
  # 
  # nm <- data.frame()
  # lbig <- matrix(0,nsimm,1)
  # 
  # for (j in 1:nsimm){
  #   
  #   nm <- as.data.frame(t(sr1a[1,,j]))
  #   
  #     for(i in 2:nrow(ba_spread)){
  #     nm[i,] <- as.data.frame(t(sr1a[i,,j]))
  #   }
  #   row.names(nm) <- row.names(ba_spread)
  #   tn1 <- nestednodf(nm)
  #   
  #   lbig[j] <- tn1$statistic[3]
  # }
  
  ##generar z-scores de betadiversidad
  
  
  lbig1 <- list()
  ldist <- list()
  
  for (j in 1:nsimm){
    
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
  
    
    lbig1[[j]] <- as.matrix(bt1)##$statistic[3]
    ldist[[j]] <- as.matrix(dist.a)
    
  }
  
  v <- reduce(lbig1, `+`)
  v <- as.matrix(v)/length(lbig1)
  
  dd <- reduce(ldist, `+`)
  dd <- as.matrix(dd)/length(ldist)
  
  
  #sss<-unlist(lbig1,TRUE,TRUE)
  #aaa<- array (sss,c(length(lbig1), nrow(as.matrix(lbig1[[1]])), ncol(as.matrix(lbig1[[1]]))))
  #ccc<- apply(aaa,c(1,2),sd)
  
  ##inlcuir loop el modelo nulo se hace en base a la matriz de disimilitudes bray curtis.
  
  sds <- matrix(0,nsimm,1)
  sds1 <- matrix(0,nrow(ba_spread),nrow(ba_spread))
  
  sdsd <- matrix(0,nsimm,1)
  sdsd1 <- matrix(0,nrow(ba_spread),nrow(ba_spread))
  
  ci1 <- matrix(0,nrow(ba_spread),nrow(ba_spread))
  ci2 <- matrix(0,nrow(ba_spread),nrow(ba_spread))
  
  cid1 <- matrix(0,nrow(ba_spread),nrow(ba_spread))
  cid2 <- matrix(0,nrow(ba_spread),nrow(ba_spread))
  
  for (j in 1:nrow(as.matrix(lbig1[[1]]))){
    for (k in 1: ncol(as.matrix(lbig1[[1]]))){
      for(i in 1:nsimm){
        sds[i]<-as.matrix(lbig1[[i]])[j,k]
        sdsd[i]<-as.matrix(ldist[[i]])[j,k]
      }
      sds1[j,k]<-sd(sds)
      sdsd1[j,k]<-sd(sdsd)
      
      ci1[j,k]<-CI(sds,ci=0.95)[[1]]
      ci2[j,k]<-CI(sds,ci=0.95)[[3]]
      
      cid1[j,k]<-CI(sdsd,ci=0.95)[[1]]
      cid2[j,k]<-CI(sdsd,ci=0.95)[[3]]
    }
  }
  
  #Whittaker
  zsc<-as.matrix(list_index[[1]]) - v
  zscsws500 <- zsc/sds1
  
  #Bray-curtis
  zsd<-as.matrix(dist.o) - dd
  zsds1000 <- zsd/sdsd1
  
  lexit<- list(obs_whitaker = list_index[[1]], obs_diss = ldist[[1]], mean_null_whitaker = v, mean_null_diss = dd,
               zscore_whitaker =zscsws500, zscore_diss =  zsds1000, ci_low_whitaker = ci2, ci_high_whitaker = ci1,
               ci_low_diss = cid2, ci_high_diss = cid1)
  return (lexit)
}

#mett<-'jaccard'



blanca_species_jaccard <- calc_zscore(file_to_read = '17082021_BaseModelosNulos_rcBasalArea.csv', file_elevation = 'elevation2.csv', spr.t.f = FALSE,
                              method.null = 'c0_both', meth.diss = 'jaccard',n.simulations = 1000, n.year = NA)


blanca_species_bray <- calc_zscore(file_to_read = '17082021_BaseModelosNulos_rcBasalArea.csv', file_elevation = 'elevation2.csv', spr.t.f = FALSE,
                                      method.null = 'c0_both', meth.diss = 'bray',n.simulations = 1000, n.year = NA) 

blanca_genera_bray <- calc_zscore(file_to_read = '17082021_BaseModelosNulos_Gen_rcBasalArea.csv', file_elevation = 'elevation2.csv', spr.t.f = FALSE,
                                   method.null = 'c0_both', meth.diss = 'bray',n.simulations = 1000, n.year = NA)

blanca_genera_jaccard <- calc_zscore(file_to_read = '17082021_BaseModelosNulos_Gen_rcBasalArea.csv', file_elevation = 'elevation2.csv', spr.t.f = FALSE,
                                  method.null = 'c0_both', meth.diss = 'jaccard',n.simulations = 1000, n.year = NA)

# blanca_species_jaccard <- calc_zscore(file_to_read = 'trich_species.csv', file_elevation = 'elev_blanca_hml.csv', spr.t.f = FALSE,
#                                       method.null = 'c0_both', meth.diss = 'jaccard',n.simulations = 1000, n.year = NA) 


for (i in 1:length(blanca_species_bray)){
  write.csv(blanca_species_bray[[i]],paste0(names(blanca_species_bray)[i],'_bray_species.csv'))
  write.csv(blanca_genera_bray[[i]],paste0(names(blanca_genera_bray)[i],'_bray_genera.csv'))
  
  write.csv(blanca_species_jaccard[[i]],paste0(names(blanca_species_jaccard)[i],'_jaccard_species.csv'))
  write.csv(blanca_genera_jaccard[[i]],paste0(names(blanca_genera_jaccard)[i],'_jaccard_genera.csv'))
  
}

# write.csv(blanca_species_bray[['obs_whitaker']],'obs_whitaker_species.csv')
# write.csv(blanca_genera_bray[['obs_whitaker']],'obs_whitaker_gen.csv')
# 
# write.csv(blanca_species_bray[["zscore_whitaker"]],'zscore_whitaker_species.csv')
# write.csv(blanca_genera_bray[['zscore_whitaker']],'zscore_whitaker_gen.csv')
# 
# write.csv(blanca_species_bray[[]],'obs_whitaker_species.csv')
# write.csv(blanca_genera_bray[['obs_whitaker']],'obs_whitaker_gen.csv')
# 
# write.csv(blanca_species_bray[["zscore_whitaker"]],'zscore_whitaker_species.csv')
# write.csv(blanca_genera_bray[['zscore_whitaker']],'zscore_whitaker_gen.csv')
# 

#write.csv(zscsws500,'zscore_species.csv')

#write.csv(ldist[[1]],'obs_dist_jaccard_species.csv')
#write.csv(zsds1000,'zscore_dist_jaccard_species.csv')


#write.csv(ci1,'upper_ci_genera.csv')
#write.csv(ci2,'lower_ci_genera.csv')

#positive means Betadiversity was higher than expected by chance
## negative by random/chance

##dejar bonito

#diversidad alfa

#calcular chao1 (spader - ChaoSpecies /fossil chao1)

# fisher.alpha, diversity shannon, hill=TRUE (vegan)

## Pielou's evenness (J):
# data(BCI)
# H <- diversity(BCI)
# simp <- diversity(BCI, "simpson")
# invsimp <- diversity(BCI, "inv")
# r.2 <- rarefy(BCI, 2)
# alpha <- fisher.alpha(BCI)
# pairs(cbind(H, simp, invsimp, r.2, alpha), pch="+", col="blue")
# ## Species richness (S) and Pielou's evenness (J):
# S <- specnumber(BCI) ## rowSums(BCI > 0) does the same...
# J <- H/log(S)  ##eveness

## glm dep=beta/alfa, independent= hr, temp (suelo y aire), biomasa (nAGC), 
## productividad medida como el delta entre censos, fertilidad suelo, prec(potencialmente), 
## diversidad funcional (e.g., shannon (calcular para traits), )









