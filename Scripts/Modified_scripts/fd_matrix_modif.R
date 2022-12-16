# Load data
df_betafd = read.csv("Data/20012022_BetaFD_Stand_NoRC.csv", sep = ",")
elev <- read.csv("Data/Altitudes_13052021.csv", sep = ";")

# Add plot names as indices and delete column with names
row.names(df_betafd) <- df_betafd$X
df_betafd[1] <- NULL

# Add plot names as indices and delete column with names
row.names(elev) <- elev$Plot
elev[1] <- NULL

# Add elevation column 
df_betafd <- merge(elev, df_betafd, by = "row.names")

# Add plot names as indices and delete column with names
row.names(df_betafd) <- df_betafd$Row.names
df_betafd[1] <- NULL

# Order by elevation 
df_betafd <- df_betafd[order(df_betafd$Elevation),]

# Reorder columns by elevation
df_betafd <- df_betafd[,c(1,11,10, 9, 8, 12, 13, 7, 6, 3, 5, 2, 4, 15, 16, 14, 17)]

# Delete elevation column
df_betafd[1] <- NULL

# Save beta diversity matrix
write.csv(df_betafd,"Data/20012022_BetaFD_Stand_NoRC.csv", row.names = TRUE)
