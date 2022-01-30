library(snowfall)
library(e1071)
library(optimbase)
library(Matrix)
library(optimsimplex)
library(neldermead)        
library(classyfire)
library(R.methodsS3)
library(R.oo)
library(R.utils)
library(R.matlab)

#importing feacal data
bgw <- readMat("/Users/goncaloleiria/Desktop/gcms_data/faecal/BWG_FA_CDvCTRL.mat")

#getting essential elements for machine learning
X <-bgw$XTIC #GC-MS data for each sample
y <-bgw$CLASS #Disease state (1 = Control, 2 = Crohn's Disease (CD))

# copy sample names from SAM to row names of X
rownames(X) <-as.character(unlist(bgw$SAM))

# get the retention times(RT)
RT <-bgw$RT
i = 19# set index number to sample you want to plot
plot(RT, X[i,], type="l", xlab="RT(min)", ylab="Intensity",main=sprintf("TIC profile for sample: %s", rownames(X)[i]))


##PCA###
X_PCA <- prcomp(X, center = TRUE, scale = FALSE)#)
summary(X_PCA)

par(mfrow=c(1,1))
## make a scree plot
X_PCA.var <- X_PCA$sdev^2
X_PCA.var.per <- round(X_PCA.var/sum(X_PCA.var)*100, 1)
barplot(X_PCA.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

#plot 2 principal components
col <- c(rep("green",12), rep("blue",11))
plot(X_PCA$x, bg=col, pch=21)

#PCA results in df and subset
mat <- as.data.frame(X_PCA$x)
subset_PCA <- mat[c(1:5)]



########################################
#MACHINE LEARNING
#######################################
test1_feacal <- cfBuild(subset_PCA, inputClass = y, bootNum = 50, ensNum = 50, parallel = TRUE, cpus = 6)


#########################################
#PERMUTATION
#########################################
test1_Permut <- cfPermute(subset_PCA, inputClass = y, bootNum = 50, ensNum = 50, permNum = 50, parallel = TRUE, cpus = 6)
#histogram of permutation results
ggPermHist(test1_Permut, density=TRUE)

#permutation descriptive stats
test1_Permut_Stats <- getPerm5Num(test1_Permut)

#Fused histogram from cfBuild and cfPermute
ggFusedHist(test1, test1_Permut)


#creating dataframe with Acc test and Acc permutation
Ensemble <- test1_feacal[["testAcc"]]
Acc_table <- as.data.frame(Ensemble)
Acc_table$Permutation <- test1_Permut[["avgAcc"]]
par(mfrow=c(1,2))
#####
dens <- apply(Acc_table, 2, density)
plot(NA, xlab="Overall Test Accuracies (%CC)", xlim=range(sapply(dens, "[", "x")),ylab="Density", ylim=range(sapply(dens, "[", "y")))
mapply(lines, dens, col=1:length(dens))
legend("topright", legend=names(dens), fill=1:length(dens))

#boxplot
boxplot(Acc_table, xlab="Type", ylab="Accuracy")


#statistical test
Stattest1 <-t.test(Acc_table$Ensemble,Acc_table$Permutation, var.equal = TRUE)
Stattest1







