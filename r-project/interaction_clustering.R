library(cluster)

setwd("C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab")

df_interaction <- as.data.frame(read.csv("db\\csv\\interaction_data.csv", header=TRUE, sep=","))
df_interaction <- df_interaction[!(df_interaction$interaction_energy > 50),]

#some plots
hist(df_interaction$interaction_energy)
plot(df_interaction$interface_residues, df_interaction$interaction_energy)



#eliminate outliers


#prepare data for kmeans (scale, remove NA)
df1 <- subset(df_interaction, select = -c(id, pdb_id,sloop_entropy,mloop_entropy,water_bridge,partial_covalent_bonds,entropy_complex))
df1 <- na.omit(df1)
df1 <- scale(df1)

wss <- (nrow(df1)-1)*sum(apply(df1,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(df1, centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

fit <- kmeans(df1, 3)

clusplot(df1, fit$cluster, color=TRUE, shade=TRUE, lines=0)

#correlations
cor(df1)

