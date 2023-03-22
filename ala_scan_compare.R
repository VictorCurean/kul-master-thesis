setwd("C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab")

par(mfrow=c(1,1))

df <- as.data.frame(read.csv("db\\csv\\motif_ala_energies.csv", header=TRUE, sep=","))
df <- na.omit(df)

#remove outliers
outliers_ie <- boxplot(df$interaction_energy, plot=FALSE)$out
outliers_hc_m <- boxplot(df$hc_motif_energy, plot=FALSE)$out
outliers_lc_m <- boxplot(df$lc_motif_energy, plot=FALSE)$out
outliers_hc_nm <- boxplot(df$hc_non_motif_energy, plot=FALSE)$out
outliers_lc_nm <- boxplot(df$lc_non_motif_energy, plot=FALSE)$out

df <- df[-which(df$interaction_energy %in% outliers_ie),]
df <- df[-which(df$hc_motif_energy %in% outliers_hc_m),]
df <- df[-which(df$hc_non_motif_energy %in% outliers_hc_nm),]
df <- df[-which(df$lc_motif_energy %in% outliers_lc_m),]
df <- df[-which(df$lc_non_motif_energy %in% outliers_lc_nm),]

df_num <- subset(df, select = c('interaction_energy', 'hc_motif_energy', 'hc_non_motif_energy', 'lc_motif_energy', 'lc_non_motif_energy'))

colmeans_ala <- colMeans(df_num)
barplot(colmeans_ala)

cor(df_num)

plot(df_num$interaction_energy, df_num$hc_motif_energy)

plot(df_num$interaction_energy, df_num$lc_motif_energy)

hist(df_num$hc_motif_energy)
hist(df_num$lc_motif_energy)
hist(df_num$hc_non_motif_energy)
hist(df_num$lc_non_motif_energy)

#NO INTERACTION TERM
lm1 <- lm(interaction_energy~hc_motif_energy+lc_motif_energy, data=df_num)
lm2 <- lm(interaction_energy~hc_motif_energy+lc_motif_energy+hc_non_motif_energy+lc_non_motif_energy, data=df_num)

#INTERACTION TERMS
lm3 <-lm(interaction_energy~hc_motif_energy+lc_motif_energy+hc_motif_energy:lc_motif_energy, data=df_num)

summary(lm1)
summary(lm2)
summary(lm3)


