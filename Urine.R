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
urineD <- readMat("/Users/goncaloleiria/Desktop/gcms_data/urine/BWG_UR_CDvCTRL.mat")

#getting essential elements for machine learning
X <-urineD$XTIC #GC-MS data for each sample
y <-urineD$CLASS #Disease state (1 = Control, 2 = Crohn's Disease (CD))

# copy sample names from SAM to row names of X
rownames(X) <-as.character(unlist(urineD$SAM))

# get the retention times(RT)
RT <-urineD$RT
i = 19# set index number to sample you want to plot
plot(RT, X[i,], type="l", xlab="RT(min)", ylab="Intensity",main=sprintf("TIC profile for sample: %s", rownames(X)[i]))

###################
##PCA###
##################
X_PCA <- prcomp(X, center = TRUE, scale = FALSE)#)
summary(X_PCA)

## make a scree plot
X_PCA.var <- X_PCA$sdev^2
X_PCA.var.per <- round(X_PCA.var/sum(X_PCA.var)*100, 1)
barplot(X_PCA.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

#plot 2 principal components
col <- c(rep("green",12), rep("blue",11))
plot(X_PCA$x, bg=col, pch=21)

#PCA results in df and subset
mat <- as.data.frame(X_PCA$x)
subset_PCA <- mat[c(1:3)] #3 first PC account for 94% of the variation

###########################
#MACHINE LEARNING#     SVMs
############################
test1_urine <- cfBuild(subset_PCA, inputClass = y, bootNum = 50, ensNum = 50, parallel = TRUE, cpus = 6)

############################
#PERMUTATION#
############################
test1_Permut_urine <- cfPermute(subset_PCA, inputClass = y, bootNum = 50, ensNum = 50, permNum = 50, parallel = TRUE, cpus = 6)

#creating dataframe with Acc test and Acc permutation
Ensemble <- test1_urine[["testAcc"]]
Acc_table <- as.data.frame(Ensemble)
Acc_table$Permutation <- test1_Permut_urine[["avgAcc"]]
par(mfrow=c(1,2))
#####
dens <- apply(Acc_table, 2, density)
plot(NA, xlab="Overall Test Accuracies (%CC)", xlim=range(sapply(dens, "[", "x")),ylab="Density", ylim=range(sapply(dens, "[", "y")))
mapply(lines, dens, col=1:length(dens))
legend("topright", legend=names(dens), fill=1:length(dens),bty="n")
#boxplot
boxplot(Acc_table, xlab="Type", ylab="Accuracy")

#statistical test
Stattest4 <-t.test(Acc_table$Ensemble,Acc_table$Permutation, var.equal = TRUE)
Stattest4


