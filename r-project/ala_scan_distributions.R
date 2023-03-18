setwd("C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab")

par(mfrow=c(1,1))
df_energies1 <- as.data.frame(read.csv("db\\csv\\energy_distribution_XXXXXXXXX_CDR-H3.csv", header=FALSE, sep=","))
mean_values1 <- sort(colMeans(df_energies1))
barplot(mean_values1, main = "XXXXXXXXX CDR-H3", col='red')
boxplot(mean_values1, col='red', main="XXXXXXXXX CDR-H3")

par(mfrow=c(1,9))
for (i in 1:ncol(df_energies1)) {
  boxplot(df_energies1[,i], main = names(df_energies1)[i], col = 'red')
}

par(mfrow=c(1,1))
df_energies2 <- as.data.frame(read.csv("db\\csv\\energy_distribution_XXXX_CDR-H3.csv", header=FALSE, sep=","))
mean_values2 <- sort(colMeans(df_energies2))
barplot(mean_values2, main = "XXXX CDR-H3", col='red')
boxplot(mean_values2, col='red', main = "XXXX CDR-H3")

par(mfrow=c(1,4))
for (i in 1:ncol(df_energies2)) {
  boxplot(df_energies2[,i], main = names(df_energies2)[i], col = 'red')
}

par(mfrow=c(1,1))
df_energies3 <- as.data.frame(read.csv("db\\csv\\energy_distribution_XXXXX_CDR-H3.csv", header=FALSE, sep=","))
df_energies3 <- na.omit(df_energies3)
mean_values3 <- sort(colMeans(df_energies3))
barplot(mean_values3, main = "XXXXX CDR-H3", col='red')
boxplot(mean_values3, col='red', main = "XXXXX CDR-H3")

par(mfrow=c(1,5))
for (i in 1:ncol(df_energies3)) {
  boxplot(df_energies3[,i], main = names(df_energies3)[i], col = 'red')
}


par(mfrow=c(1,1))
df_energies4 <- as.data.frame(read.csv("db\\csv\\energy_distribution_XXXXXX_CDR-H3.csv", header=FALSE, sep=","))
df_energies4 <- na.omit(df_energies4)
mean_values4 <- sort(colMeans(df_energies4))
barplot(mean_values4, main = "XXXXXX CDR-H3", col='red')
boxplot(mean_values4, col='red', main = "XXXXXX CDR-H3")

par(mfrow=c(1,6))
for (i in 1:ncol(df_energies4)) {
  boxplot(df_energies4[,i], main = names(df_energies4)[i], col = 'red')
}

par(mfrow=c(1,1))
df_energies5 <- as.data.frame(read.csv("db\\csv\\energy_distribution_XXXXXXX_CDR-H3.csv", header=FALSE, sep=","))
df_energies5 <- na.omit(df_energies5)
mean_values5 <- sort(colMeans(df_energies5))
barplot(mean_values5, main = "XXXXXXX CDR-H3", col='red')
boxplot(mean_values5, col='red', main = "XXXXXXX CDR-H3")

par(mfrow=c(1,7))
for (i in 1:ncol(df_energies5)) {
  boxplot(df_energies5[,i], main = names(df_energies5)[i], col = 'red')
}


#### CDR - L1 

par(mfrow=c(1,1))
df_energies6 <- as.data.frame(read.csv("db\\csv\\energy_distribution_XXX_CDR-L1.csv", header=FALSE, sep=","))
df_energies6 <- na.omit(df_energies6)
mean_values6 <- sort(colMeans(df_energies6))
barplot(mean_values6, main = "XXX CDR-L1", col='blue')
boxplot(mean_values6, col='blue', main = "XXX CDR-L1")

par(mfrow=c(1,3))
for (i in 1:ncol(df_energies6)) {
  boxplot(df_energies6[,i], main = names(df_energies6)[i], col = 'blue')
}

### LFR - 3

par(mfrow=c(1,1))
df_energies7 <- as.data.frame(read.csv("db\\csv\\energy_distribution_XXXX_LFR3.csv", header=FALSE, sep=","))
df_energies7 <- na.omit(df_energies7)
mean_values7 <- sort(colMeans(df_energies7))
barplot(mean_values7, main = "XXXX LFR3", col='blue')
boxplot(mean_values7, col='blue', main = "XXXX LFR3")

par(mfrow=c(1,4))
for (i in 1:ncol(df_energies7)) {
  boxplot(df_energies7[,i], main = names(df_energies7)[i], col = 'blue')
}


## CDR - H3 with gap

par(mfrow=c(1,1))
df_energies8 <- as.data.frame(read.csv("db\\csv\\energy_distribution_XXXX1X_CDR-H3.csv", header=FALSE, sep=","))
df_energies8 <- na.omit(df_energies8)
names(df_energies8)[5] <- "G5"
mean_values8 <- sort(colMeans(df_energies8))
barplot(mean_values8, main = "XXXX1X CDR-H3", col='red')
boxplot(mean_values8, col='red', main = "XXXX1X CDR-H3")

par(mfrow=c(1,6))
for (i in 1:ncol(df_energies8)) {
  boxplot(df_energies8[,i], main = names(df_energies8)[i], col = 'red')
}


par(mfrow=c(1,1))
df_energies9 <- as.data.frame(read.csv("db\\csv\\energy_distribution_XXXX1XX_CDR-H3.csv", header=FALSE, sep=","))
df_energies9 <- na.omit(df_energies9)
names(df_energies9)[5] <- "G5"
mean_values9 <- sort(colMeans(df_energies9))
barplot(mean_values9, main = "XXXX1XX CDR-H3", col='red')
boxplot(mean_values9, col='red', main = "XXXX1XX CDR-H3")

par(mfrow=c(1,7))
for (i in 1:ncol(df_energies9)) {
  boxplot(df_energies9[,i], main = names(df_energies9)[i], col = 'red')
}

par(mfrow=c(1,1))
df_energies10 <- as.data.frame(read.csv("db\\csv\\energy_distribution_XXXXX1X_CDR-H3.csv", header=FALSE, sep=","))
df_energies10 <- na.omit(df_energies10)
names(df_energies10)[6] <- "G6"
mean_values10 <- sort(colMeans(df_energies10))
barplot(mean_values10, main = "XXXX1X CDR-H3", col='red')
boxplot(mean_values10, col='red', main = "XXXX1X CDR-H3")

par(mfrow=c(1,7))
for (i in 1:ncol(df_energies10)) {
  boxplot(df_energies10[,i], main = names(df_energies10)[i], col = 'red')
}



















