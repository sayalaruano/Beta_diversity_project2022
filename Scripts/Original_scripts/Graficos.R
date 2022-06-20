

FDHills = read.csv("AlphFD_Indexvalues_270122.csv", sep = ",")
head(FDHills)

FDHillsdiv0 = ggplot(FDHills, aes(x = TD_Hill, y = FDHill_q1_stand, colour = Elevation)) + 
  geom_point() + scale_colour_gradient(low="red", high="blue") + 
  geom_smooth(method = lm, se = TRUE)

FDHillsdiv1 = ggscatter(FDHills, x = "TD_Hill", y = "FDHill_q1_stand",
                        add = "reg.line",conf.int = TRUE,
                        palette = "npg") +stat_cor(label.x = 20, label.y = 5) +
  stat_regline_equation(label.x = 40, label.y = 5)

FDHillsdiv = ggplot(FDHills, aes(x = Elevation, y = FDHill_q1_stand)) + 
  geom_point() +
  geom_smooth(method = lm, se = TRUE)

TDHillsdiv = ggplot(FDHills, aes(x = Elevation, y = TD_Hill)) + 
  geom_point() +
  geom_smooth(method = lm, se = TRUE)

finalFEve = ggscatter(FDHills, x = "Elevation", y = "Hill_Diversity",
                      add = "reg.line",conf.int = TRUE,
                      palette = "npg") +stat_cor(label.x = 1000, label.y = 1) +
  stat_regline_equation(label.x = 2000, label.y = 1)

### Usar dataframes creados en taxonomic_diversity


### 

BetaDiv_SR = read.csv("20012022_BetaFD_Stand_NoRC.csv", sep = ",")
head(BetaDiv_SR)
str(BetaDiv_SR)
rownames(BetaDiv_SR)<-BetaDiv_SR$X
head(BetaDiv_SR)
del1 <- c(BetaDiv_SR,"X")
BetaDiv_SR <- as.matrix(BetaDiv_SR[,-which(names(BetaDiv_SR) %in% del1)])

BetaDiv_SR1 = data.frame(Elevation = colnames(BetaDiv_SR)[col(BetaDiv_SR)], Elevation2 = rownames(BetaDiv_SR)[row(BetaDiv_SR)], SR_BA_dist = c(BetaDiv_SR))


BetaDiv_SR1 = subset(BetaDiv_SR1, SR_BA_dist != 0.0000000)

BetaDiv_SR1$Elevation <- as.numeric(as.character(BetaDiv_SR1$Elevation))

BetaDiv_SR1$Elevation2 <- as.numeric(as.character(BetaDiv_SR1$Elevation2))

BetaDiv_SR1$Elev_change = abs(BetaDiv_SR1$Elevation - BetaDiv_SR1$Elevation2)



BETA_SRBA_allGrad2019 = ggplot(BetaDiv_SR1, aes(x = Elev_change, y = SR_BA_dist)) + geom_point()

BETA_SRBA_1vsresto = subset(BetaDiv_SR1, Elevation2 == "632")

BETA_SRBA_1vsrestoGrad = ggplot(BETA_SRBA_1vsresto, aes(x = Elev_change, y = SR_BA_dist)) + geom_point() + geom_smooth(method = lm, se = FALSE)


BetaDiv2 = read.csv("BetaDiver_periodicalchange.csv", sep = ",")
head(BetaDiv2)
str(BetaDiv2)
BETAallGrad2019 = ggplot(BetaDiv2, aes(x = Elevation, y = FD_Turnover)) + geom_point() +
  stat_summary(fun.data = BetaDiv2) + 
  geom_smooth(method='lm') 


BETAallGrad2019 = ggplot(BetaDiv2, aes(x = TD_Whittaker_Dissimilarity, y = FD_Dissimilarity, colour = Elevation)) + geom_point() +
  scale_colour_gradient(low="red", high="blue") + 
  geom_smooth(method='lm', se = FALSE) 

FD_TD_diss = ggscatter(BetaDiv2, x = "TD_Whittaker_Dissimilarity", y = "FD_Dissimilarity",
                       add = "reg.line",conf.int = TRUE,
                       palette = "npg") +stat_cor(label.x = 0.5, label.y = 1) +
  stat_regline_equation(label.x = 1, label.y = 1)


residual = residuals(BetaAll)
preds = predict(BetaAll)


plot(residual ~  preds, xlab = "Predicted Values",  ylab = "Residuals")

library(nortest)
pearson.test(residual)
lillie.test(residual)


par(mfrow = c(2, 2))
plot(BetaAll)

op = par(no.readonly = TRUE)

par(op)

lmtest::bptest (BetaAll)



BetaDiv2 = read.csv("20012022_BetaFD_Stand_NoRC.csv", sep = ",")
rownames(BetaDiv2) <- BetaDiv2$X
BetaDiv1 <- c(BetaDiv2,"X")
BetaDiv2 <- BetaDiv2[,-which(names(BetaDiv2) %in% BetaDiv1)]
BetaDiv2 = as.matrix(BetaDiv2)

BetaDiv3 = data.frame(col = colnames(BetaDiv2)[col(BetaDiv2)], row = rownames(BetaDiv2)[row(BetaDiv2)], dist = c(BetaDiv2))

BetaDiv3$Elevation = BetaDiv3$col


BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "BECL_01"] <- "2313"
BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "BECL_03"] <- "2203"
BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "CEDR_01"] <- "2492"
BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "CEDR_03"] <- "2212"
BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "INTI_01"] <- "1879"
BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "INTI_02"] <- "1829"
BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "MALO_01"] <- "1018"
BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "MALO_02"] <- "827"
BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "MAPI_01"] <- "653"
BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "MAPI_02"] <- "632"
BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "MIND_01"] <- "1277"
BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "RIBR_01"] <- "1640"
BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "VERD_01"] <- "3421"
BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "VERD_02"] <- "2932"
BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "VERD_03"] <- "3109"
BetaDiv3["Elevation"][BetaDiv3["Elevation"] == "YANA_01"] <- "3507"

BetaDiv3$Elevation2 = BetaDiv3$row


BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "BECL_01"] <- "2313"
BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "BECL_03"] <- "2203"
BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "CEDR_01"] <- "2492"
BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "CEDR_03"] <- "2212"
BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "INTI_01"] <- "1879"
BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "INTI_02"] <- "1829"
BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "MALO_01"] <- "1018"
BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "MALO_02"] <- "827"
BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "MAPI_01"] <- "653"
BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "MAPI_02"] <- "632"
BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "MIND_01"] <- "1277"
BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "RIBR_01"] <- "1640"
BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "VERD_01"] <- "3421"
BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "VERD_02"] <- "2932"
BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "VERD_03"] <- "3109"
BetaDiv3["Elevation2"][BetaDiv3["Elevation2"] == "YANA_01"] <- "3507"

BetaDiv3 = subset(BetaDiv3, dist != 0.0000000)
BetaDiv3 = subset(BetaDiv3, dist != "NA")

BetaDiv3$Elevation <- as.numeric(as.character(BetaDiv3$Elevation))

BetaDiv3$Elevation2 <- as.numeric(as.character(BetaDiv3$Elevation2))

BetaDiv3$Elev_change = abs(BetaDiv3$Elevation - BetaDiv3$Elevation2)



BETAallGrad2019 = ggplot(BetaDiv3, aes(x = Elev_change, y = dist)) + geom_point() +
  stat_summary(fun.data = BetaDiv3) + 
  geom_smooth(method='lm') 

BetaAll = lm(dist ~ Elev_change, data = BetaDiv3) 
summary (BetaAll)

BETA1vsresto = subset(BetaDiv3, row == "MAPI_02"| col == "MAPI_02")

BETA1vsrestoGrad = ggplot(BETA1vsresto, aes(x = Elev_change, y = dist)) + geom_point() +
  stat_summary(fun.data=BETA1vsresto$dist) + 
  geom_smooth(method='lm') 

BetaAllMapi = lm(dist ~ Elev_change, data = BETA1vsresto) 
summary (BetaAllMapi)




