library(glmnet)


setwd("C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab")

dm1 <- as.data.frame(read.csv("foldx\\datamatrix\\datamatrix.txt", header=TRUE, sep="\t"))
dm1 <- dm1[!(dm1$energy > 50),]

dm2 <- na.omit(dm1)

#remove outliers - 5ANY, 3J8W, 3IYW

pairs(dm1[,2:21])
hist(dm1$energy)
hist(dm1$no_motifs)
hist(dm1$conn_res_paratope)
hist(dm1$conn_res_epitope)


x <- as.matrix(subset(dm2, select = -c(pdb, energy)))
y <- dm2$energy

cv_model <- cv.glmnet(x, y, alpha=1)
best_lambda <- cv_model$lambda.min
plot(cv_model)

best_model <- glmnet(x, y, alpha=1, lambda=best_lambda)
coef(best_model)


plot(dm1$energy, dm1$conn_res_paratope)
abline(lm(dm1$energy~dm1$conn_res_paratope), col='red')
plot(dm1$energy, dm1$conn_res_epitope)
