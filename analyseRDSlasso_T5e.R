#setwd("~")

setwd("/Users/isobelridler/Desktop/RosalindHPC/outputT5e")

library(psych);library(xtable); library(tidyverse)

# read chain00_**.rds files in
temp <- list.files(pattern="*.rds")
res5 <- lapply(temp, readRDS)

# check lists are correct length/contents ok?
res6 <- res5[lapply(res5,length)>2]

for (i in length(res6)){
  print(res6[[i]][[]])
}

# get all results into one list called res11 (1:1000 in this case as there are 1000 repeats of each cell so 1000 outputs for each of 30 cells)
count=0;res11 <- list()
for(i in 1:length(res6)){
  for(j in 1:1000){
    count=count+1
    res11[[count]] <- res6[[i]][[j]]
  }
}


# make a numeric vector 'test' of all trueRep values, dropping nulls
test <- lapply(res11, function(x) {as.numeric(paste(x[[6]]))})
test1 <- unlist(test)

# make a numeric vector 'samp' of all sample sizes used in each run, dropping nulls
samp <- lapply(res11, function(x) {as.numeric(paste(x[[5]]))})
samp <- unlist(samp)

# check how many failed to converge (out of 900,000 because 30 lambda and 30,000 total runs) = 22134 i.e. 2.46%
conv1 <- lapply(res11, function(x) {as.numeric(paste(x[[1]]$fits[,"conv"]))})
conv1 <- unlist(conv1)
conv.check <- conv1 != 0
sum(conv1 != 0)

# check that all sample sizes were not = 80, and then use ind to select only those that were not and make a new object called res12
ind <- samp != 80
nunique <- sum(ind)
res12 <- res11[ind]

# specify empty objects to fill with incoming for loop output
location <- rep(NA,nunique)
coll <- rep(NA,nunique)
samp2<- rep(NA,nunique)
# find BIC and convergence of regsem models (30 lambda so 30 BIC for each cell), then find minimum BIC and make sure it converged, and use location to give an index of which line that BIC is so we can pull out corresponding values 
for(i in 1:length(res12)){
  if(is.na(res12[[i]][[1]])){
    print(NA)
  }else{
    BIC = res12[[i]][[1]]$fits[,"BIC"]
    conv = res12[[i]][[1]]$fits[,"conv"]
    loc = which(abs(BIC)==min(abs(BIC[conv == 0 & is.nan(BIC)==FALSE & is.na(conv)==FALSE])))
    location[i] = loc
  }
}

for(i in 1:length(res12)){
  
}

# pull out the collinearity and sample size for each cell into coll and samp2 respectively
for(i in 1:nunique){
  coll[i] = res12[[i]][[4]]
  samp2[i] = res12[[i]][[5]]
}

# set up objects to place the regsem and lavaan final coefficients and pvalues (edit numbers to reflect number of estimated paramteters)
reg.coef = matrix(NA,nunique,80)
lav.coef <- matrix(NA,nunique,80)
lav.pvalue <- matrix(NA,nunique,80)

for(i in 1:nunique){
  if(is.na(location[i])){
    reg.coef[i,] = NA
    lav.coef[i,] = NA
    lav.pvalue[i,]= NA
  }else{
    reg.coef[i,] <- res12[[i]][[1]]$parameters[location[i],]
    lav.coef[i,] <- res12[[i]]$lav.coef
    lav.pvalue[i,] <- res12[[i]]$lav.pvalue
  }
}



# parameters actual simulated values (45 specified f1/f2=~, y~f1/f2,x#, f1~f2 - x1~~x1 etc, and y~~y unspecified)
pars2 <- c(rep(0.3,4),rep(0.7,4),rep(0.5,7),rep(0.3,4),rep(0.7,4),rep(0.5,7),0.4,0.6,0,0.2,0.4,0.6,0,0.2,0.4,0.6,0,0.2,0.4,0.6,0,0.2,0.4,0.6,0.3,rep(1,31))


# create an empty matrix to fill in with all unique combinations (12 in this case as 4 x 3 reps), for each of 22 total estimated params
sqe.diff <- matrix(NA,nunique,80)
# difference between lavaan estimate and true (squared because want root mean square error)
for(i in 1:nunique){
  sqe.diff[i,] = (lav.coef[i,] - pars) * (lav.coef[i,] - pars)
}

###################################
########****** RMSE *******########
###################################

################### MLE SEM ######################

################### COLLINEARITY #####################

# mean, sd etc for each collinearity
mse.lav.coll <- describeBy(sqe.diff,trim=0.02,group=coll)

# make empty object ready for...
mse.lav.coll2 <- matrix(NA,length(unique(coll)),80)
# (average? YES, avergage trimmed to 2dp) sqrt of mse (RMSE) of parameter estimates for each collinearity - why sqrt??
for(i in 1:length(unique(coll))){
  mse.lav.coll2[i,] <- sqrt(mse.lav.coll[[i]]$trimmed)
}

# averages of RMSE for type of parameter
mse.lav.coll3 <- cbind(rowMeans(mse.lav.coll2[, c(1:4, 16:19)]),
                       rowMeans(mse.lav.coll2[, c(9:15, 24:30)]),
                       rowMeans(mse.lav.coll2[, c(5:8, 20:23)]),
                       rowMeans(mse.lav.coll2[, 31:32]),
                       rowMeans(mse.lav.coll2[, c(33, 37, 41, 45)]),
                       rowMeans(mse.lav.coll2[, c(34, 38, 42, 46)]),
                       rowMeans(mse.lav.coll2[, c(35, 39, 43, 47)]),
                       rowMeans(mse.lav.coll2[, c(36, 40, 44, 48)]),
                       mse.lav.coll2[, 49],
                       rowMeans(mse.lav.coll2[, 50:79]),
                       mse.lav.coll2[, 80])

# set nsReps to how many MC reps were used in lasso.R simulation

n <- 1000

# Monte Carlo SE of RMSE - MCSE = sigma/sqrt(nsReps) ?
mse.lav.coll3.se <- cbind(apply(mse.lav.coll2[, c(1:4, 16:19)],1,sd)/sqrt(n),
                       apply(mse.lav.coll2[, c(9:15, 24:30)],1,sd)/sqrt(n),
                       apply(mse.lav.coll2[, c(5:8, 20:23)],1,sd)/sqrt(n),
                       apply(mse.lav.coll2[, 31:32],1,sd)/sqrt(n),
                       apply(mse.lav.coll2[, c(33, 37, 41, 45)],1,sd)/sqrt(n),
                       apply(mse.lav.coll2[, c(34, 38, 42, 46)],1,sd)/sqrt(n),
                       apply(mse.lav.coll2[, c(35, 39, 43, 47)],1,sd)/sqrt(n),
                       apply(mse.lav.coll2[, c(36, 40, 44, 48)],1,sd)/sqrt(n),
                       mse.lav.coll2[, 49]/sqrt(n),
                       apply(mse.lav.coll2[, 50:79],1,sd)/sqrt(n),
                       mse.lav.coll2[, 80]/sqrt(n))

(mat1 <- round(mse.lav.coll3,3))
(mat1.se <- round(mse.lav.coll3.se,3))
row.names(mat1) <- c(.20,.30,.50,.70,.85)
colnames(mat1) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
row.names(mat1.se) <- c(.20,.30,.50,.70,.85)
colnames(mat1.se) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
xtable(mat1,"Lavaan Collinearity RMSE")

###################### SAMPLE SIZE #########################

# mean, sd etc for each N
mse.lav.samp <- describeBy(sqe.diff,trim=0.02,group=samp2)
# make empty object ready for output
mse.lav.samp2 <- matrix(NA,length(unique(samp2)),80)
# populate matrix with RMSE
for(i in 1:length(unique(samp2))){
  mse.lav.samp2[i,] <- sqrt(mse.lav.samp[[i]]$trimmed)
}

# averages of RMSE for type of parameter
mse.lav.samp3 <- cbind(rowMeans(mse.lav.samp2[, c(1:4, 16:19)]),
                       rowMeans(mse.lav.samp2[, c(9:15, 24:30)]),
                       rowMeans(mse.lav.samp2[, c(5:8, 20:23)]),
                       rowMeans(mse.lav.samp2[, 31:32]),
                       rowMeans(mse.lav.samp2[, c(33, 37, 41, 45)]),
                       rowMeans(mse.lav.samp2[, c(34, 38, 42, 46)]),
                       rowMeans(mse.lav.samp2[, c(35, 39, 43, 47)]),
                       rowMeans(mse.lav.samp2[, c(36, 40, 44, 48)]),
                       mse.lav.samp2[, 49],
                       rowMeans(mse.lav.samp2[, 50:79]),
                       mse.lav.samp2[, 80])

# Monte Carlo SE of RMSE - MCSE = sigma/sqrt(nsReps) ?
mse.lav.samp3.se <- cbind(apply(mse.lav.samp2[, c(1:4, 16:19)],1,sd)/sqrt(n),
                          apply(mse.lav.samp2[, c(9:15, 24:30)],1,sd)/sqrt(n),
                          apply(mse.lav.samp2[, c(5:8, 20:23)],1,sd)/sqrt(n),
                          apply(mse.lav.samp2[, 31:32],1,sd)/sqrt(n),
                          apply(mse.lav.samp2[, c(33, 37, 41, 45)],1,sd)/sqrt(n),
                          apply(mse.lav.samp2[, c(34, 38, 42, 46)],1,sd)/sqrt(n),
                          apply(mse.lav.samp2[, c(35, 39, 43, 47)],1,sd)/sqrt(n),
                          apply(mse.lav.samp2[, c(36, 40, 44, 48)],1,sd)/sqrt(n),
                          mse.lav.samp2[, 49]/sqrt(n),
                          apply(mse.lav.samp2[, 50:79],1,sd)/sqrt(n),
                          mse.lav.samp2[, 80]/sqrt(n))

round(mse.lav.samp3,3)
(mat2 <- round(mse.lav.samp3,3))
row.names(mat2) <- c(100, 150, 200, 500, 1000, 2000)
colnames(mat2) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
(mat2.se <- round(mse.lav.samp3.se,3))
row.names(mat2.se) <- c(100, 150, 200, 500, 1000, 2000)
colnames(mat2.se) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
xtable(mat2,"Lavaan Sample RMSE")

################# regSEM #######################

# empty matrix for regsem differences
sqe.diff2 <- matrix(NA,nunique,80)
# regsem and real par diff
for(i in 1:nunique){
  sqe.diff2[i,] = (reg.coef[i,] - pars) * (reg.coef[i,] - pars)
}

################# COLLINEARITY ##################

mse.reg.coll <- describeBy(sqe.diff2,trim=0.02,group=coll)
mse.reg.coll2 <- matrix(NA,length(unique(coll)),80)
for(i in 1:length(unique(coll))){
  mse.reg.coll2[i,] <- sqrt(mse.reg.coll[[i]]$trimmed)
}

mse.reg.coll3 <- cbind(rowMeans(mse.reg.coll2[, c(1:4, 16:19)]),
                       rowMeans(mse.reg.coll2[, c(9:15, 24:30)]),
                       rowMeans(mse.reg.coll2[, c(5:8, 20:23)]),
                       rowMeans(mse.reg.coll2[, 31:32]),
                       rowMeans(mse.reg.coll2[, c(33, 37, 41, 45)]),
                       rowMeans(mse.reg.coll2[, c(34, 38, 42, 46)]),
                       rowMeans(mse.reg.coll2[, c(35, 39, 43, 47)]),
                       rowMeans(mse.reg.coll2[, c(36, 40, 44, 48)]),
                       mse.reg.coll2[, 49],
                       rowMeans(mse.reg.coll2[, 50:79]),
                       mse.reg.coll2[, 80])
                      
mse.reg.coll3.se <- cbind(apply(mse.reg.coll2[, c(1:4, 16:19)],1,sd)/sqrt(n),
                          apply(mse.reg.coll2[, c(9:15, 24:30)],1,sd)/sqrt(n),
                          apply(mse.reg.coll2[, c(5:8, 20:23)],1,sd)/sqrt(n),
                          apply(mse.reg.coll2[, 31:32],1,sd)/sqrt(n),
                          apply(mse.reg.coll2[, c(33, 37, 41, 45)],1,sd)/sqrt(n),
                          apply(mse.reg.coll2[, c(34, 38, 42, 46)],1,sd)/sqrt(n),
                          apply(mse.reg.coll2[, c(35, 39, 43, 47)],1,sd)/sqrt(n),
                          apply(mse.reg.coll2[, c(36, 40, 44, 48)],1,sd)/sqrt(n),
                          mse.reg.coll2[, 49]/sqrt(n),
                          apply(mse.reg.coll2[, 50:79],1,sd)/sqrt(n),
                          mse.reg.coll2[, 80]/sqrt(n))

(mat3 <- round(mse.reg.coll3,3))
row.names(mat3) <- c(.20,.30,.50,.70,.85)
colnames(mat3) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
(mat3.se <- round(mse.reg.coll3.se,3))
row.names(mat3.se) <- c(.20,.30,.50,.70,.85)
colnames(mat3.se) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
xtable(mat3,"Regsem lasso Collinearity RMSE")

############ SAMPLE SIZE #################

mse.reg.samp <- describeBy(sqe.diff2,trim=0.02,group=samp2)
mse.reg.samp2 <- matrix(NA,length(unique(samp2)),80)
for(i in 1:length(unique(samp2))){
  mse.reg.samp2[i,] <- sqrt(mse.reg.samp[[i]]$trimmed)
}

mse.reg.samp3 <- cbind(rowMeans(mse.reg.samp2[, c(1:4, 16:19)]),
                       rowMeans(mse.reg.samp2[, c(9:15, 24:30)]),
                       rowMeans(mse.reg.samp2[, c(5:8, 20:23)]),
                       rowMeans(mse.reg.samp2[, 31:32]),
                       rowMeans(mse.reg.samp2[, c(33, 37, 41, 45)]),
                       rowMeans(mse.reg.samp2[, c(34, 38, 42, 46)]),
                       rowMeans(mse.reg.samp2[, c(35, 39, 43, 47)]),
                       rowMeans(mse.reg.samp2[, c(36, 40, 44, 48)]),
                       mse.reg.samp2[, 49],
                       rowMeans(mse.reg.samp2[, 50:79]),
                       mse.reg.samp2[, 80])

# Monte Carlo SE of RMSE - MCSE = sigma/sqrt(nsReps) ?
mse.reg.samp3.se <- cbind(apply(mse.reg.samp2[, c(1:4, 16:19)],1,sd)/sqrt(n),
                          apply(mse.reg.samp2[, c(9:15, 24:30)],1,sd)/sqrt(n),
                          apply(mse.reg.samp2[, c(5:8, 20:23)],1,sd)/sqrt(n),
                          apply(mse.reg.samp2[, 31:32],1,sd)/sqrt(n),
                          apply(mse.reg.samp2[, c(33, 37, 41, 45)],1,sd)/sqrt(n),
                          apply(mse.reg.samp2[, c(34, 38, 42, 46)],1,sd)/sqrt(n),
                          apply(mse.reg.samp2[, c(35, 39, 43, 47)],1,sd)/sqrt(n),
                          apply(mse.reg.samp2[, c(36, 40, 44, 48)],1,sd)/sqrt(n),
                          mse.reg.samp2[, 49]/sqrt(n),
                          apply(mse.reg.samp2[, 50:79],1,sd)/sqrt(n),
                          mse.reg.samp2[, 80]/sqrt(n))

(mat4 <- round(mse.reg.samp3,3))
row.names(mat4) <- c(100, 150, 200, 500, 1000, 2000)
colnames(mat4) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
(mat4.se <- round(mse.reg.samp3.se,3))
row.names(mat4.se) <- c(100, 150, 200, 500, 1000, 2000)
colnames(mat4.se) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
xtable(mat4,"Regsem lasso Sample RMSE")


###################################
########****** BIAS *******########
###################################


############ LAVAAN ##################

############ COLLINEARITY ###################

bias.lav.coll <- matrix(NA,length(unique(coll)),80)
means.lav.coll <- describeBy(lav.coef,trim=0.02,group=coll)
for(i in 1:length(unique(coll))){
  for(j in 1: length(pars)){
    # if theta = 0, use average bias, if not, use average relative bias
    bias.lav.coll[i,j] <- 100*ifelse(pars[j]==0,means.lav.coll[[i]]$trimmed[j]-pars[j],(means.lav.coll[[i]]$trimmed[j]-pars[j])/pars[j])
  }
}
# Average Absolute Bias/Average Absolute Relative Bias
bias.lav.coll2 <- cbind(rowMeans(abs(bias.lav.coll[, c(1:4, 16:19)])),
                        rowMeans(abs(bias.lav.coll[, c(9:15, 24:30)])),
                        rowMeans(abs(bias.lav.coll[, c(5:8, 20:23)])),
                        rowMeans(abs(bias.lav.coll[, 31:32])),
                        rowMeans(abs(bias.lav.coll[, c(33, 37, 41, 45)])),
                        rowMeans(abs(bias.lav.coll[, c(34, 38, 42, 46)])),
                        rowMeans(abs(bias.lav.coll[, c(35, 39, 43, 47)])),
                        rowMeans(abs(bias.lav.coll[, c(36, 40, 44, 48)])),
                        abs(bias.lav.coll[, 49]),
                        rowMeans(abs(bias.lav.coll[, 50:79])),
                        abs(bias.lav.coll[, 80]))

bias.lav.coll2.se <- cbind(apply(abs(bias.lav.coll[, c(1:4, 16:19)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.coll[, c(9:15, 24:30)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.coll[, c(5:8, 20:23)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.coll[, 31:32]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.coll[, c(33, 37, 41, 45)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.coll[, c(34, 38, 42, 46)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.coll[, c(35, 39, 43, 47)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.coll[, c(36, 40, 44, 48)]),1,sd)/sqrt(n),
                           abs(bias.lav.coll[, 49])/sqrt(n),
                           apply(abs(bias.lav.coll[, 50:79]),1,sd)/sqrt(n),
                           abs(bias.lav.coll[, 80])/sqrt(n))

(mat11 <- round(bias.lav.coll2,3))
row.names(mat11) <- c(.20,.30,.50,.70,.85)
colnames(mat11) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
(mat11.se <- round(bias.lav.coll2.se,3))
row.names(mat11.se) <- c(.20,.30,.50,.70,.85)
colnames(mat11.se) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
xtable(mat11,"Lavaan Collinearity Relative Bias")

############### SAMPLE SIZE ##################

bias.lav.samp <- matrix(NA,length(unique(samp2)),80)
means.lav.samp <- describeBy(lav.coef,trim=0.02,group=samp2)
for(i in 1:length(unique(samp2))){
  for(j in 1: length(pars)){
    bias.lav.samp[i,j] <- 100*ifelse(pars[j]==0,means.lav.samp[[i]]$trimmed[j]-pars[j],(means.lav.samp[[i]]$trimmed[j]-pars[j])/pars[j])
  }
}

bias.lav.samp2 <- cbind(rowMeans(abs(bias.lav.samp[, c(1:4, 16:19)])),
                        rowMeans(abs(bias.lav.samp[, c(9:15, 24:30)])),
                        rowMeans(abs(bias.lav.samp[, c(5:8, 20:23)])),
                        rowMeans(abs(bias.lav.samp[, 31:32])),
                        rowMeans(abs(bias.lav.samp[, c(33, 37, 41, 45)])),
                        rowMeans(abs(bias.lav.samp[, c(34, 38, 42, 46)])),
                        rowMeans(abs(bias.lav.samp[, c(35, 39, 43, 47)])),
                        rowMeans(abs(bias.lav.samp[, c(36, 40, 44, 48)])),
                        abs(bias.lav.samp[, 49]),
                        rowMeans(abs(bias.lav.samp[, 50:79])),
                        abs(bias.lav.samp[, 80]))

bias.lav.samp2.se <- cbind(apply(abs(bias.lav.samp[, c(1:4, 16:19)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.samp[, c(9:15, 24:30)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.samp[, c(5:8, 20:23)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.samp[, 31:32]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.samp[, c(33, 37, 41, 45)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.samp[, c(34, 38, 42, 46)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.samp[, c(35, 39, 43, 47)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.samp[, c(36, 40, 44, 48)]),1,sd)/sqrt(n),
                           abs(bias.lav.samp[, 49])/sqrt(n),
                           apply(abs(bias.lav.samp[, 50:79]),1,sd)/sqrt(n),
                           abs(bias.lav.samp[, 80])/sqrt(n))

(mat22 <- round(bias.lav.samp2,3))
row.names(mat22) <- c(100, 150, 200, 500, 1000, 2000)
colnames(mat22) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
(mat22.se <- round(bias.lav.samp2.se,3))
row.names(mat22.se) <- c(100, 150, 200, 500, 1000, 2000)
colnames(mat22.se) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
xtable(mat22,"Lavaan Sample Relative Bias")

############### regSEM ###################

############## COLLINEARTY #################

bias.reg.coll <- matrix(NA,length(unique(coll)),80)
means.reg.coll <- describeBy(reg.coef,trim=0.02,group=coll)
for(i in 1:length(unique(coll))){
  for(j in 1: length(pars)){
    bias.reg.coll[i,j] <- 100*ifelse(pars[j]==0,means.reg.coll[[i]]$trimmed[j]-pars[j],(means.reg.coll[[i]]$trimmed[j]-pars[j])/pars[j])
  }
}

bias.reg.coll2 <- cbind(rowMeans(abs(bias.reg.coll[, c(1:4, 16:19)])),
                        rowMeans(abs(bias.reg.coll[, c(9:15, 24:30)])),
                        rowMeans(abs(bias.reg.coll[, c(5:8, 20:23)])),
                        rowMeans(abs(bias.reg.coll[, 31:32])),
                        rowMeans(abs(bias.reg.coll[, c(33, 37, 41, 45)])),
                        rowMeans(abs(bias.reg.coll[, c(34, 38, 42, 46)])),
                        rowMeans(abs(bias.reg.coll[, c(35, 39, 43, 47)])),
                        rowMeans(abs(bias.reg.coll[, c(36, 40, 44, 48)])),
                        abs(bias.reg.coll[, 49]),
                        rowMeans(abs(bias.reg.coll[, 50:79])),
                        abs(bias.reg.coll[, 80]))

bias.reg.coll2.se <- cbind(apply(abs(bias.reg.coll[, c(1:4, 16:19)]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.coll[, c(9:15, 24:30)]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.coll[, c(5:8, 20:23)]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.coll[, 31:32]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.coll[, c(33, 37, 41, 45)]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.coll[, c(34, 38, 42, 46)]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.coll[, c(35, 39, 43, 47)]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.coll[, c(36, 40, 44, 48)]),1,sd)/sqrt(n),
                           abs(bias.reg.coll[, 49])/sqrt(n),
                           apply(abs(bias.reg.coll[, 50:79]),1,sd)/sqrt(n),
                           abs(bias.reg.coll[, 80])/sqrt(n))

(mat33 <- round(bias.reg.coll2,3))
row.names(mat33) <- c(.20,.30,.50,.70,.85)
colnames(mat33) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
(mat33.se <- round(bias.reg.coll2.se,3))
row.names(mat33.se) <- c(.20,.30,.50,.70,.85)
colnames(mat33.se) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
xtable(mat33,"Regsem lasso Collinearity Relative Bias")


############### SAMPLE SIZE ####################

bias.reg.samp <- matrix(NA,length(unique(samp2)),80)
means.reg.samp <- describeBy(reg.coef,trim=0.02,group=samp2)
for(i in 1:length(unique(samp2))){
  for(j in 1: length(pars)){
    bias.reg.samp[i,j] <- 100*ifelse(pars[j]==0,means.reg.samp[[i]]$trimmed[j]-pars[j],(means.reg.samp[[i]]$trimmed[j]-pars[j])/pars[j])
  }
}

bias.reg.samp2 <- cbind(rowMeans(abs(bias.reg.samp[, c(1:4, 16:19)])),
                        rowMeans(abs(bias.reg.samp[, c(9:15, 24:30)])),
                        rowMeans(abs(bias.reg.samp[, c(5:8, 20:23)])),
                        rowMeans(abs(bias.reg.samp[, 31:32])),
                        rowMeans(abs(bias.reg.samp[, c(33, 37, 41, 45)])),
                        rowMeans(abs(bias.reg.samp[, c(34, 38, 42, 46)])),
                        rowMeans(abs(bias.reg.samp[, c(35, 39, 43, 47)])),
                        rowMeans(abs(bias.reg.samp[, c(36, 40, 44, 48)])),
                        abs(bias.reg.samp[, 49]),
                        rowMeans(abs(bias.reg.samp[, 50:79])),
                        abs(bias.reg.samp[, 80]))

bias.reg.samp2.se <- cbind(apply(abs(bias.reg.samp[, c(1:4, 16:19)]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.samp[, c(9:15, 24:30)]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.samp[, c(5:8, 20:23)]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.samp[, 31:32]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.samp[, c(33, 37, 41, 45)]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.samp[, c(34, 38, 42, 46)]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.samp[, c(35, 39, 43, 47)]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.samp[, c(36, 40, 44, 48)]),1,sd)/sqrt(n),
                           abs(bias.reg.samp[, 49])/sqrt(n),
                           apply(abs(bias.reg.samp[, 50:79]),1,sd)/sqrt(n),
                           abs(bias.reg.samp[, 80])/sqrt(n))

(mat44 <- round(bias.reg.samp2,3))
row.names(mat44) <- c(100, 150, 200, 500, 1000, 2000)
colnames(mat44) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
(mat44.se <- round(bias.reg.samp2.se,3))
row.names(mat44.se) <- c(100, 150, 200, 500, 1000, 2000)
colnames(mat44.se) <- c("small loadings","medium loadings","large loadings","y ~ LV","noise","small","medium","large","LV co-variance","residual variance","out variance")
xtable(mat44,"Regsem lasso Sample Size Relative Bias")


#################################################
########****** variable selection *******########
#################################################

###### REGSEM ######

# REGRESSIONS

null <- 1*as.matrix(abs(reg.coef[,c(33, 37, 41, 45)]) > .001)
null.samp.reg <- cbind(c(100, 150, 200, 500, 1000, 2000),
                       matrix(unlist(describeBy(rowMeans(null),group=samp2)),6,13,byrow=T)[,3])
null.samp.reg.se <- cbind(c(100, 150, 200, 500, 1000, 2000),
                          matrix(unlist(describeBy(apply(null,1,sd)/sqrt(n),group=samp2)),6,13,byrow=T)[,3])
null.coll.reg <- cbind(c(.20,.30,.50,.70,.85),
                       matrix(unlist(describeBy(rowMeans(null),group=coll)),5,13,byrow=T)[,3])
null.coll.reg.se <- cbind(c(.20,.30,.50,.70,.85),
                          matrix(unlist(describeBy(apply(null,1,sd)/sqrt(n),group=coll)),5,13,byrow=T)[,3])

small <- 1*as.matrix(abs(reg.coef[,c(34, 38, 42, 46)]) < .001)
small.samp.reg <- cbind(c(100, 150, 200, 500, 1000, 2000),
                      matrix(unlist(describeBy(rowMeans(small),group=samp2)),6,13,byrow=T)[,3])
small.samp.reg.se <- cbind(c(100, 150, 200, 500, 1000, 2000),
                         matrix(unlist(describeBy(apply(small,1,sd)/sqrt(n),group=samp2)),6,13,byrow=T)[,3])
small.coll.reg <- cbind(c(.20,.30,.50,.70,.85),
                      matrix(unlist(describeBy(rowMeans(small),group=coll)),5,13,byrow=T)[,3])
small.coll.reg.se <- cbind(c(.20,.30,.50,.70,.85),
                         matrix(unlist(describeBy(apply(small,1,sd)/sqrt(n),group=coll)),5,13,byrow=T)[,3])

med <- 1*as.matrix(abs(reg.coef[,c(35, 39, 43, 47)]) < .001)
med.samp.reg <- cbind(c(100, 150, 200, 500, 1000, 2000),
                      matrix(unlist(describeBy(rowMeans(med),group=samp2)),6,13,byrow=T)[,3])
med.samp.reg.se <- cbind(c(100, 150, 200, 500, 1000, 2000),
                         matrix(unlist(describeBy(apply(med,1,sd)/sqrt(n),group=samp2)),6,13,byrow=T)[,3])
med.coll.reg <- cbind(c(.20,.30,.50,.70,.85),
                      matrix(unlist(describeBy(rowMeans(med),group=coll)),5,13,byrow=T)[,3])
med.coll.reg.se <- cbind(c(.20,.30,.50,.70,.85),
                         matrix(unlist(describeBy(apply(med,1,sd)/sqrt(n),group=coll)),5,13,byrow=T)[,3])

large <- 1*as.matrix(abs(reg.coef[,c(36, 40, 44, 48)]) < .001)
large.samp.reg <- cbind(c(100, 150, 200, 500, 1000, 2000),
                      matrix(unlist(describeBy(rowMeans(large),group=samp2)),6,13,byrow=T)[,3])
large.samp.reg.se <- cbind(c(100, 150, 200, 500, 1000, 2000),
                         matrix(unlist(describeBy(apply(large,1,sd)/sqrt(n),group=samp2)),6,13,byrow=T)[,3])
large.coll.reg <- cbind(c(.20,.30,.50,.70,.85),
                      matrix(unlist(describeBy(rowMeans(large),group=coll)),5,13,byrow=T)[,3])
large.coll.reg.se <- cbind(c(.20,.30,.50,.70,.85),
                         matrix(unlist(describeBy(apply(large,1,sd)/sqrt(n),group=coll)),5,13,byrow=T)[,3])

LV <-  1*as.matrix(abs(reg.coef[,c(31:32)]) < .001)
LV.samp.reg <- cbind(c(100, 150, 200, 500, 1000, 2000),
                        matrix(unlist(describeBy(rowMeans(LV),group=samp2)),6,13,byrow=T)[,3])
LV.samp.reg.se <- cbind(c(100, 150, 200, 500, 1000, 2000),
                           matrix(unlist(describeBy(apply(LV,1,sd)/sqrt(n),group=samp2)),6,13,byrow=T)[,3])
LV.coll.reg <- cbind(c(.20,.30,.50,.70,.85),
                        matrix(unlist(describeBy(rowMeans(LV),group=coll)),5,13,byrow=T)[,3])
LV.coll.reg.se <- cbind(c(.20,.30,.50,.70,.85),
                           matrix(unlist(describeBy(apply(LV,1,sd)/sqrt(n),group=coll)),5,13,byrow=T)[,3])


#################################################
##############****** ERRORS *******##############
#################################################

#REGRESSIONS

reg.samp.err <- cbind(null.samp.reg, small.samp.reg, med.samp.reg, large.samp.reg, LV.samp.reg)[,c(1,2,4,6,8,10)]
reg.coll.err <- cbind(null.coll.reg, small.coll.reg, med.coll.reg, large.coll.reg, LV.coll.reg)[,c(1,2,4,6,8,10)]

colnames(reg.samp.err) <- c("Sample Size","Noise","Small","Medium","Large", "LV")
colnames(reg.coll.err) <- c("Collinearity","Noise","Small","Medium", "Large", "LV")

reg.samp.err.se <- cbind(null.samp.reg.se, small.samp.reg.se, med.samp.reg.se, large.samp.reg.se, LV.samp.reg.se)[,c(1,2,4,6,8,10)]
reg.coll.err.se <- cbind(null.coll.reg.se, small.coll.reg.se, med.coll.reg.se, large.coll.reg.se, LV.coll.reg.se)[,c(1,2,4,6,8,10)]

colnames(reg.samp.err.se) <- c("Sample Size","Noise","Small","Medium","Large", "LV")
colnames(reg.coll.err.se) <- c("Collinearity","Noise","Small","Medium", "Large", "LV")

xtable(reg.samp.err)
xtable(reg.coll.err)


############## LAVAAN #############


null.lav <- 1*as.matrix(lav.pvalue[,c(33, 37, 41, 45)] < .05)
null.samp.lav <- cbind(c(100, 150, 200, 500, 1000, 2000),
                       matrix(unlist(describeBy(rowMeans(null.lav),group=samp2)),6,13,byrow=T)[,3])
null.samp.lav.se <- cbind(c(100, 150, 200, 500, 1000, 2000),
                          matrix(unlist(describeBy(apply(null.lav,1,sd)/sqrt(n),group=samp2)),6,13,byrow=T)[,3])
null.coll.lav <- cbind(c(.20,.30,.50,.70,.85),
                       matrix(unlist(describeBy(rowMeans(null.lav),group=coll)),5,13,byrow=T)[,3])
null.coll.lav.se <- cbind(c(.20,.30,.50,.70,.85),
                          matrix(unlist(describeBy(apply(null.lav,1,sd)/sqrt(n),group=coll)),5,13,byrow=T)[,3])

small.lav <- 1*as.matrix(lav.pvalue[,c(34, 38, 42, 46)] > .05)
small.samp.lav <- cbind(c(100, 150, 200, 500, 1000, 2000),
                      matrix(unlist(describeBy(rowMeans(small.lav),group=samp2)),6,13,byrow=T)[,3])
small.samp.lav.se <- cbind(c(100, 150, 200, 500, 1000, 2000),
                         matrix(unlist(describeBy(apply(small.lav,1,sd)/sqrt(n),group=samp2)),6,13,byrow=T)[,3])
small.coll.lav <- cbind(c(.20,.30,.50,.70,.85),
                      matrix(unlist(describeBy(rowMeans(small.lav),group=coll)),5,13,byrow=T)[,3])
small.coll.lav.se <- cbind(c(.20,.30,.50,.70,.85),
                         matrix(unlist(describeBy(apply(small.lav,1,sd)/sqrt(n),group=coll)),5,13,byrow=T)[,3])

med.lav <- 1*as.matrix(lav.pvalue[,c(35, 39, 43, 47)] > .05)
med.samp.lav <- cbind(c(100, 150, 200, 500, 1000, 2000),
                      matrix(unlist(describeBy(rowMeans(med.lav),group=samp2)),6,13,byrow=T)[,3])
med.samp.lav.se <- cbind(c(100, 150, 200, 500, 1000, 2000),
                         matrix(unlist(describeBy(apply(med.lav,1,sd)/sqrt(n),group=samp2)),6,13,byrow=T)[,3])
med.coll.lav <- cbind(c(.20,.30,.50,.70,.85),
                      matrix(unlist(describeBy(rowMeans(med.lav),group=coll)),5,13,byrow=T)[,3])
med.coll.lav.se <- cbind(c(.20,.30,.50,.70,.85),
                         matrix(unlist(describeBy(apply(med.lav,1,sd)/sqrt(n),group=coll)),5,13,byrow=T)[,3])

large.lav <- 1*as.matrix(lav.pvalue[,c(36, 40, 44, 48)] > .05)
large.samp.lav <- cbind(c(100, 150, 200, 500, 1000, 2000),
                        matrix(unlist(describeBy(rowMeans(large.lav),group=samp2)),6,13,byrow=T)[,3])
large.samp.lav.se <- cbind(c(100, 150, 200, 500, 1000, 2000),
                           matrix(unlist(describeBy(apply(large.lav,1,sd)/sqrt(n),group=samp2)),6,13,byrow=T)[,3])
large.coll.lav <- cbind(c(.20,.30,.50,.70,.85),
                        matrix(unlist(describeBy(rowMeans(large.lav),group=coll)),5,13,byrow=T)[,3])
large.coll.lav.se <- cbind(c(.20,.30,.50,.70,.85),
                           matrix(unlist(describeBy(apply(large.lav,1,sd)/sqrt(n),group=coll)),5,13,byrow=T)[,3])

LV.lav <- 1*as.matrix(lav.pvalue[,c(31:32)] > .05)
LV.samp.lav <- cbind(c(100, 150, 200, 500, 1000, 2000),
                        matrix(unlist(describeBy(rowMeans(LV.lav),group=samp2)),6,13,byrow=T)[,3])
LV.samp.lav.se <- cbind(c(100, 150, 200, 500, 1000, 2000),
                           matrix(unlist(describeBy(apply(LV.lav,1,sd)/sqrt(n),group=samp2)),6,13,byrow=T)[,3])
LV.coll.lav <- cbind(c(.20,.30,.50,.70,.85),
                        matrix(unlist(describeBy(rowMeans(LV.lav),group=coll)),5,13,byrow=T)[,3])
LV.coll.lav.se <- cbind(c(.20,.30,.50,.70,.85),
                           matrix(unlist(describeBy(apply(LV.lav,1,sd)/sqrt(n),group=coll)),5,13,byrow=T)[,3])
##### ERRORS ######

lav.samp.err <- cbind(null.samp.lav, small.samp.lav, med.samp.lav, large.samp.lav, LV.samp.lav)[,c(1,2,4,6,8,10)]
lav.coll.err <- cbind(null.coll.lav, small.coll.lav, med.coll.lav, large.coll.lav, LV.coll.lav)[,c(1,2,4,6,8,10)]

colnames(lav.samp.err) <- c("Sample Size","Noise","Small","Medium","Large", "LV")
colnames(lav.coll.err) <- c("Collinearity","Noise","Small","Medium", "Large", "LV")

lav.samp.err.se <- cbind(null.samp.lav.se, small.samp.lav.se, med.samp.lav.se, large.samp.lav.se, LV.samp.lav.se)[,c(1,2,4,6,8,10)]
lav.coll.err.se <- cbind(null.coll.lav.se, small.coll.lav.se, med.coll.lav.se, large.coll.lav.se, LV.coll.lav.se)[,c(1,2,4,6,8,10)]

colnames(lav.samp.err.se) <- c("Sample Size","Noise","Small","Medium","Large", "LV")
colnames(lav.coll.err.se) <- c("Collinearity","Noise","Small","Medium", "Large", "LV")

xtable(lav.samp.err)
xtable(lav.coll.err)

# Combine all

all.objects=list(mat1,mat2,mat3,mat4,mat11,mat22,mat33,mat44,
                 lav.samp.err,lav.coll.err,reg.samp.err,reg.coll.err,
                 mat1.se,mat2.se,mat3.se,mat4.se,mat11.se,mat22.se,mat33.se,mat44.se,
                 lav.samp.err.se,lav.coll.err.se,reg.samp.err.se,reg.coll.err.se)


save(all.objects,file="/Users/isobelridler/Desktop/RosalindHPC/outputT5e/tables2T5e210301.RData")
