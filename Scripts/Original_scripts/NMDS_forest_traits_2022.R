# Set Working directory
wDpaper = "C:/Users/antonella.bernardi/OneDrive - Universidad de Las Américas/Francisco Cuesta/paper"

setwd(wDpaper)

bio<-read.csv("Matriz_abundanciaFREC_plots.csv", h=T, sep = ";")
bio[is.na(bio)] <- 0
row.names(bio)<-bio$PLOT
bio[1] <- NULL
bio = as.matrix(bio)
head(bio)
class(bio)


library(vegan)

set.seed(2)
nmds <- metaMDS(bio, distance = "bray", k = 2)


plot(nmds)
nmds

treat = c(rep("Treatment1", 4), rep("Treatment2", 4), rep("Treatment3", 4), rep("Treatment4", 4))
orditorp(nmds, display = "sites", )
ordiellipse(nmds, groups = treat, draw = "polygon", col = "grey90", label = F)

## Con matriz de distancias sin especies

bio_distmat <- vegdist(bio, method = "bray")

bio_distmat <- as.matrix(bio_distmat, labels = T)


nmds <- metaMDS(bio_distmat, distance = "bray", k = 2, maxit = 999, trymax = 500,
        wascores = TRUE)

plot(nmds)
nmds

treat = c(rep("Treatment1", 4), rep("Treatment2", 4), rep("Treatment3", 4), rep("Treatment4", 4))
orditorp(nmds, display = "sites", )
ordiellipse(nmds, groups = treat, draw = "polygon", col = "grey90", label = F)


# Para agregar las variables ambientales

env0 <- read.csv("20212704_envi_Plots.csv", sep = ";")
env = subset(env0, select = c("PLOT", "Elevationmasl", "AGC_2019Mgha1", "PCA1_temp", "PC1.p", "PC2.p"))
head(env)

#extract NMDS scores (x and y coordinates)
data.scores1 = as.data.frame(scores(nmds))
data.scores1$PLOT = rownames(data.scores1)

#add columns to data frame 
data.scores = merge(data.scores1, env, by = "PLOT")

head(data.scores)
