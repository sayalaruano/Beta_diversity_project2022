library(reshape2)
library(hillR)
library(dplyr)
library(ggplot2)
library(MBI)
library(vegan)


# Set Working directory
wDpaper = "C:/Users/antonella.bernardi/OneDrive - Universidad de Las Américas/Francisco Cuesta/paper"

setwd(wDpaper)

# alpha diver Count by species abundance


DIV = read.csv("2019_Div.csv", sep = ",")

head(DIV)

DIV2019 = rename(count(DIV, Plot, Scientific_name), Freq = n)
DIV2019Mat = recast(DIV2019, Plot ~ Scientific_name, id.var = c("Scientific_name", "Plot"))
DIV2019Mat[2] <- NULL

DIV2019Mat[is.na(DIV2019Mat)] <- 0
row.names(DIV2019Mat) <- DIV2019Mat$Plot
DIV2019Mat[1] <- NULL

BD2019 = hill_taxa(DIV2019Mat, q = 1, MARGIN = 1, base = exp(1))
BD2019 <- as.data.frame(BD2019)
BD2019$Plot <- rownames(BD2019)
rownames(BD2019) <- NULL
names(BD2019)[1] = "Hill_Diversity"
print(BD2019)

BD2019$Elevation = BD2019$Plot


BD2019["Elevation"][BD2019["Elevation"] == "BECL_01"] <- "2313"
BD2019["Elevation"][BD2019["Elevation"] == "BECL_03"] <- "2203"
BD2019["Elevation"][BD2019["Elevation"] == "CEDR_01"] <- "2492"
BD2019["Elevation"][BD2019["Elevation"] == "CEDR_03"] <- "2212"
BD2019["Elevation"][BD2019["Elevation"] == "INTI_01"] <- "1879"
BD2019["Elevation"][BD2019["Elevation"] == "INTI_02"] <- "1829"
BD2019["Elevation"][BD2019["Elevation"] == "MALO_01"] <- "1018"
BD2019["Elevation"][BD2019["Elevation"] == "MALO_02"] <- "827"
BD2019["Elevation"][BD2019["Elevation"] == "MAPI_01"] <- "653"
BD2019["Elevation"][BD2019["Elevation"] == "MAPI_02"] <- "632"
BD2019["Elevation"][BD2019["Elevation"] == "MIND_01"] <- "1277"
BD2019["Elevation"][BD2019["Elevation"] == "RIBR_01"] <- "1640"
BD2019["Elevation"][BD2019["Elevation"] == "VERD_01"] <- "3421"
BD2019["Elevation"][BD2019["Elevation"] == "VERD_02"] <- "2932"
BD2019["Elevation"][BD2019["Elevation"] == "VERD_03"] <- "3109"
BD2019["Elevation"][BD2019["Elevation"] == "YANA_01"] <- "3507"


## alfa div by basal area


DIV2019 = read.csv("2019_Div.csv", sep = ",")

head(DIV2019)

DIV2019. = subset(DIV2019, select = c("Plot", "Scientific_name", "BA"))


sumBA2019 = aggregate(x = DIV2019.$BA,
                  by = list(DIV2019.$Plot, DIV2019$Scientific_name),
                  FUN = sum, na.rm = TRUE)

names(sumBA2019)[3] = "BAm2_2019"

names(sumBA2019)[2] = "Scientific_name"

names(sumBA2019)[1] = "Plot"


Sum2019Mat = recast(sumBA2019, Plot ~ Scientific_name, id.var = c("Scientific_name", "Plot"))
Sum2019Mat[2] <- NULL

Sum2019Mat[is.na(Sum2019Mat)] <- 0
row.names(Sum2019Mat) <- Sum2019Mat$Plot
Sum2019Mat[1] <- NULL

BD2019BA = hill_taxa(Sum2019Mat, q = 1, MARGIN = 1, base = exp(1))
BD2019BA <- as.data.frame(BD2019BA)
BD2019BA$Plot <- rownames(BD2019BA)
rownames(BD2019BA) <- NULL
names(BD2019BA)[1] = "Hill_Diversity"
print(BD2019BA)

### con raiz cuadrada de BA para alpha diver. con Hill numbers

sumBA2019$Raiz_BA ='^'(sumBA2019$BAm2_2019,1/2)

sumBA2019 = subset(sumBA2019, select = c("Scientific_name", "Plot", "Raiz_BA"))

Sum2019Mat = recast(sumBA2019, Plot ~ Scientific_name, id.var = c("Scientific_name", "Plot"))
Sum2019Mat[2] <- NULL

Sum2019Mat[is.na(Sum2019Mat)] <- 0
row.names(Sum2019Mat) <- Sum2019Mat$Plot
Sum2019Mat[1] <- NULL

BD2019BA = hill_taxa(Sum2019Mat, q = 1, MARGIN = 1, base = exp(1))
BD2019BA <- as.data.frame(BD2019BA)
BD2019BA$Plot <- rownames(BD2019BA)
rownames(BD2019BA) <- NULL
names(BD2019BA)[1] = "Hill_Diversity"
print(BD2019BA)

# agregar elevacion a plots para graficar

### rarefaccion, eveness y riqueza

Rare = rarefy(DIV2019Mat, min(rowSums(DIV2019Mat)))


curvaBA <- specaccum(Sum2019Mat)

plot(curvaBA)


S <- specnumber(DIV2019Mat)



# evenness Pielou, con indice de Shannon

H <- diversity(Sum2019Mat)

J <- as.data.frame(H/log(specnumber(Sum2019Mat)))

J$Plot <- rownames(J)
rownames(J) <- NULL
names(J)[1] = "Pielou_even"


# Chao1/Chao2


specpool(DIV2019Mat)

estimateR(DIV2019Mat)


chao2(DIV2019Mat, taxa.row = FALSE)


# expresar riqueza con rarefacción

Srar <- as.data.frame(rarefy(DIV2019Mat, min(rowSums(DIV2019Mat))))

Srar$Plot <- rownames(Srar)
rownames(Srar) <- NULL
names(Srar)[1] = "Riqueza"


####

#### Diversidad Beta
DIV2019 = read.csv("2019_Div.csv", sep = ",")

head(DIV2019)
DIV2019 = aggregate(x = DIV2019$BA,
                           by = list(DIV2019$Plot, DIV2019$Scientific_name),
                           FUN = sum, na.rm = TRUE)

names(DIV2019)[2] = "Scientific_name"

names(DIV2019)[1] = "Plot"

names(DIV2019)[3] = "BA"

# para usar con raiz cuadrada
DIV2019$Raiz_BA ='^'(DIV2019$BA,1/2)
BDCWMGrad = read.csv("Altitudes_13052021.csv", sep = ",")
head(BDCWMGrad)


BDCWMGrad = subset(BDCWMGrad, select = c("Plot", "Elevation"))
DIV2019 = merge(DIV2019, BDCWMGrad, by = "Plot")

DIV2019 = subset(DIV2019, select = c("Scientific_name", "Elevation", "Raiz_BA"))


BetaBA2019Mat_SR = recast(DIV2019, Elevation ~ Scientific_name, id.var = c("Scientific_name", "Elevation"))
BetaBA2019Mat_SR[2] <- NULL

BetaBA2019Mat_SR[is.na(BetaBA2019Mat_SR)] <- 0
row.names(BetaBA2019Mat_SR)<-BetaBA2019Mat_SR$Elevation
BetaBA2019Mat_SR[1] <- NULL

BetaBA2019Mat_SR = as.matrix(BetaBA2019Mat_SR)

BetaDiv_SR = betadiver(BetaBA2019Mat_SR, "w")

BetaDiv_SR = as.matrix(BetaDiv_SR)






