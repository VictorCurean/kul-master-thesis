setwd("C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab")

par(mfrow=c(1,1))

df <- as.data.frame(read.csv("db\\csv\\motif_ala_energies.csv", header=TRUE, sep=","))
df <- na.omit(df)
df_num <- subset(df, select = c('interaction_energy', 'hc_motif_energy', 'hc_non_motif_energy', 'lc_motif_energy', 'lc_non_motif_energy'))
colmeans_ala <- colMeans(df_num)
barplot(colmeans_ala)

cor(df_num)