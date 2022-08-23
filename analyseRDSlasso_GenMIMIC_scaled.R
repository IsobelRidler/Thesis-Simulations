

setwd("/Users/isobelridler/Desktop/RosalindHPC/outputGenMIMIC/RegsOnlyScaled/")

library(psych);library(xtable); library(tidyverse)

# read chain00_**.rds files in
temp <- list.files(pattern="RegsOnly*")
res5 <- lapply(temp, readRDS)

# check lists are correct length/contents ok?
res6 <- res5[lapply(res5,length)>2]

# get all results into one list called res7 (1:100 in this case as there are 100 repeats of each cell so 100 outputs for each of 36 cells)
#count=0;res7 <- list()
#for(i in 1:length(res6)){
#  for(j in 1:100){
#    count=count+1
#    res7[[count]] <- res6[[i]][[j]]
#  }
#}

res7 <- unlist(res6, recursive = FALSE)

#re-order to get reps in line
res8 <- res7[order(sapply(res7,'[[',7))]

# to get nsnps = 
nsnps <- c(30, 60, 90, 120)

for (i in nsnps) {
  nsnps.v <- i
  res9 <- res8[which(sapply(res8,'[[',6) == nsnps.v)]


# make a numeric vector 'test' of all trueRep values, dropping nulls
# setdiff(3601:4800, test1)
test <- lapply(res9, function(x) {as.numeric(paste(x[[7]]))})
test1 <- unlist(test)

# check how many failed to converge: 0.029% (30); 0.047% (60); 0.091% (90); 0.154% (120)
conv1 <- lapply(res9, function(x) {as.numeric(paste(x[[1]]$fits[,"conv"]))})
conv1 <- unlist(conv1)
conv2 <- na.omit(conv1)
conv.check <- conv2 != 0
sum(conv2 != 0)
# number of combinations
n <- 12
# number of lambda
l <- 20

converged.perc <- ((sum(conv2 != 0))/((n*100*l)-(length(conv1)-length(conv2))))*100

print(converged.perc)

# make a numeric vector 'samp' of all sample sizes used in each run, dropping nulls
samp <- lapply(res9, function(x) {as.numeric(paste(x[[5]]))})
samp <- unlist(samp)
# check that all sample sizes were not = 80, and then use ind to select only those that were not and make a new object called res12
ind <- samp != 200
nunique <- sum(ind)
res10 <- res9[ind]

# specify empty objects to fill with incoming for loop output
location <- rep(NA,nunique)
coll <- rep(NA,nunique)
samp1<- rep(NA,nunique)
nsnps1 <- rep(NA,nunique)
# find BIC and convergence of regsem models (20 lambda so 20 BIC for each cell), then find minimum BIC and make sure it converged, and use location to give an index of which line that BIC is so we can pull out corresponding values 
for(i in 1:length(res10)){
  if(is.na(res10[[i]][[1]]$final_pars)){
    print(NA)
  }else{
    BIC = res10[[i]][[1]]$fits[,"BIC"]
    conv = res10[[i]][[1]]$fits[,"conv"]
    loc = which(abs(BIC)==min(abs(BIC[conv == 0 & is.nan(BIC)==FALSE & is.na(conv)==FALSE])))
    location[i] = loc
  }
}

# pull out the collinearity and sample size for each cell into coll, samp1 and nsnps respectively
for(i in 1:nunique){
  coll[i] = res10[[i]][[4]]
  samp1[i] = res10[[i]][[5]]
  nsnps1[i] = res10[[i]][[6]]
}

# set up objects to place the regsem and lavaan final coefficients and pvalues (edit numbers to reflect number of estimated paramteters)

reg.coef = matrix(NA,nunique,(nsnps.v+6))
lav.coef <- matrix(NA,nunique,(nsnps.v+6))
lav.pvalue <- matrix(NA,nunique,(nsnps.v+6))

for(i in 1:nunique){
  if(is.na(location[i])){
    reg.coef[i,] = NA
    lav.coef[i,] = NA
    lav.pvalue[i,]= NA
  }else{
    reg.coef[i,] <- res10[[i]][[1]]$parameters[location[i],]
    lav.coef[i,] <- res10[[i]]$lav.coef
    lav.pvalue[i,] <- res10[[i]]$lav.pvalue
  }
}

# parameters actual simulated values (specified snp->f1, , and  unspecified)
pars <- c(rep(0, nsnps.v/3), rep(0.001, nsnps.v/6), rep(0.005, nsnps.v/6), rep(0.01, nsnps.v/6), rep(0.05, nsnps.v/6), 0.6, 0.4, 0.3, rep(1, 3))


# set nsReps to how many MC reps were used in lasso.R simulation
n <- 100
# variables for selecting which parameters to group and keep for tables and graphs
a <- nsnps.v/3
b <- nsnps.v/6

# create an empty matrix to fill in with all unique combinations (12 in this case as 4 x 3 reps), for each of nsnps + 6 total estimated params
sqe.diff <- matrix(NA,nunique,(nsnps.v+6))
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
mse.lav.coll2 <- matrix(NA,length(unique(coll)),(nsnps.v+6))
# (average? YES, avergage trimmed to 2dp) sqrt of mse (RMSE) of parameter estimates for each collinearity - why sqrt??
for(i in 1:length(unique(coll))){
  mse.lav.coll2[i,] <- sqrt(mse.lav.coll[[i]]$trimmed)
}

# averages of RMSE for type of parameter
mse.lav.coll3 <- cbind(rowMeans(mse.lav.coll2[, 1:a]),
                       rowMeans(mse.lav.coll2[, (a+1):(2*a)]),
                       rowMeans(mse.lav.coll2[, ((2*a)+1):nsnps.v]),
                       rowMeans(mse.lav.coll2[, (nsnps.v+1):(nsnps.v+3)]),
                       rowMeans(mse.lav.coll2[, (nsnps.v+4):(nsnps.v+6)]))

# Monte Carlo SE of RMSE - MCSE = sigma/sqrt(nsReps) ?
mse.lav.coll3.se <- cbind(apply(mse.lav.coll2[, 1:a],1,sd)/sqrt(n),
                          apply(mse.lav.coll2[, (a+1):(2*a)],1,sd)/sqrt(n),
                          apply(mse.lav.coll2[, ((2*a)+1):nsnps.v],1,sd)/sqrt(n),
                          apply(mse.lav.coll2[, (nsnps.v+1):(nsnps.v+3)],1,sd)/sqrt(n),
                          apply(mse.lav.coll2[, (nsnps.v+4):(nsnps.v+6)],1,sd)/sqrt(n))

(mat1 <- round(mse.lav.coll3,5))
(mat1.se <- round(mse.lav.coll3.se,5))
row.names(mat1) <- c(0, 0.1, 0.2)
colnames(mat1) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
row.names(mat1.se) <- c(0, 0.1, 0.2)
colnames(mat1.se) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
xtable(mat1,"Lavaan Collinearity RMSE")

###################### SAMPLE SIZE #########################

# mean, sd etc for each N
mse.lav.samp <- describeBy(sqe.diff,trim=0.02,group=samp1)
# make empty object ready for output
mse.lav.samp2 <- matrix(NA,length(unique(samp1)),(nsnps.v+6))
# populate matrix with RMSE
for(i in 1:length(unique(samp1))){
  mse.lav.samp2[i,] <- sqrt(mse.lav.samp[[i]]$trimmed)
}

# averages of RMSE for type of parameter
mse.lav.samp3 <- cbind(rowMeans(mse.lav.samp2[, 1:a]),
                       rowMeans(mse.lav.samp2[, (a+1):(2*a)]),
                       rowMeans(mse.lav.samp2[, ((2*a)+1):nsnps.v]),
                       rowMeans(mse.lav.samp2[, (nsnps.v+1):(nsnps.v+3)]),
                       rowMeans(mse.lav.samp2[, (nsnps.v+4):(nsnps.v+6)]))

# Monte Carlo SE of RMSE - MCSE = sigma/sqrt(nsReps) ?
mse.lav.samp3.se <- cbind(apply(mse.lav.samp2[, 1:a],1,sd)/sqrt(n),
                          apply(mse.lav.samp2[, (a+1):(2*a)],1,sd)/sqrt(n),
                          apply(mse.lav.samp2[, ((2*a)+1):nsnps.v],1,sd)/sqrt(n),
                          apply(mse.lav.samp2[, (nsnps.v+1):(nsnps.v+3)],1,sd)/sqrt(n),
                          apply(mse.lav.samp2[, (nsnps.v+4):(nsnps.v+6)],1,sd)/sqrt(n))


(mat2 <- round(mse.lav.samp3,5))
row.names(mat2) <- c(1000, 4000, 7000, 10000)
colnames(mat2) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
(mat2.se <- round(mse.lav.samp3.se,5))
row.names(mat2.se) <- c(1000, 4000, 7000, 10000)
colnames(mat2.se) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
xtable(mat2,"Lavaan Sample RMSE")


################# regSEM #######################

# empty matrix for regsem differences
sqe.diff2 <- matrix(NA,nunique,(nsnps.v+6))
# regsem and real par diff
for(i in 1:nunique){
  sqe.diff2[i,] = (reg.coef[i,] - pars) * (reg.coef[i,] - pars)
}

################# COLLINEARITY ##################

mse.reg.coll <- describeBy(sqe.diff2,trim=0.02,group=coll)
mse.reg.coll2 <- matrix(NA,length(unique(coll)),(nsnps.v+6))
for(i in 1:length(unique(coll))){
  mse.reg.coll2[i,] <- sqrt(mse.reg.coll[[i]]$trimmed)
}

mse.reg.coll3 <- cbind(rowMeans(mse.reg.coll2[, 1:a]),
                       rowMeans(mse.reg.coll2[, (a+1):(2*a)]),
                       rowMeans(mse.reg.coll2[, ((2*a)+1):nsnps.v]),
                       rowMeans(mse.reg.coll2[, (nsnps.v+1):(nsnps.v+3)]),
                       rowMeans(mse.reg.coll2[, (nsnps.v+4):(nsnps.v+6)]))

mse.reg.coll3.se <- cbind(apply(mse.reg.coll2[, 1:a],1,sd)/sqrt(n),
                          apply(mse.reg.coll2[, (a+1):(2*a)],1,sd)/sqrt(n),
                          apply(mse.reg.coll2[, ((2*a)+1):nsnps.v],1,sd)/sqrt(n),
                          apply(mse.reg.coll2[, (nsnps.v+1):(nsnps.v+3)],1,sd)/sqrt(n),
                          apply(mse.reg.coll2[, (nsnps.v+4):(nsnps.v+6)],1,sd)/sqrt(n))

(mat3 <- round(mse.reg.coll3,5))
row.names(mat3) <- c(0, 0.1, 0.2)
colnames(mat3) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
(mat3.se <- round(mse.reg.coll3.se,5))
row.names(mat3.se) <- c(0, 0.1, 0.2)
colnames(mat3.se) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
xtable(mat3,"Regsem lasso Collinearity RMSE")

############ SAMPLE SIZE #################

mse.reg.samp <- describeBy(sqe.diff2,trim=0.02,group=samp1)
mse.reg.samp2 <- matrix(NA,length(unique(samp1)),(nsnps.v+6))
for(i in 1:length(unique(samp1))){
  mse.reg.samp2[i,] <- sqrt(mse.reg.samp[[i]]$trimmed)
}

mse.reg.samp3 <- cbind(rowMeans(mse.reg.samp2[, 1:a]),
                       rowMeans(mse.reg.samp2[, (a+1):(2*a)]),
                       rowMeans(mse.reg.samp2[, ((2*a)+1):nsnps.v]),
                       rowMeans(mse.reg.samp2[, (nsnps.v+1):(nsnps.v+3)]),
                       rowMeans(mse.reg.samp2[, (nsnps.v+4):(nsnps.v+6)]))

# Monte Carlo SE of RMSE - MCSE = sigma/sqrt(nsReps) ?
mse.reg.samp3.se <- cbind(apply(mse.reg.samp2[, 1:a],1,sd)/sqrt(n),
                          apply(mse.reg.samp2[, (a+1):(2*a)],1,sd)/sqrt(n),
                          apply(mse.reg.samp2[, ((2*a)+1):nsnps.v],1,sd)/sqrt(n),
                          apply(mse.reg.samp2[, (nsnps.v+1):(nsnps.v+3)],1,sd)/sqrt(n),
                          apply(mse.reg.samp2[, (nsnps.v+4):(nsnps.v+6)],1,sd)/sqrt(n))

(mat4 <- round(mse.reg.samp3,5))
row.names(mat4) <- c(1000, 4000, 7000, 10000)
colnames(mat4) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
(mat4.se <- round(mse.reg.samp3.se,5))
row.names(mat4.se) <- c(1000, 4000, 7000, 10000)
colnames(mat4.se) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
xtable(mat4,"Regsem lasso Sample RMSE")


###################################
########****** BIAS *******########
###################################


############ LAVAAN ##################

############ COLLINEARITY ###################

bias.lav.coll <- matrix(NA,length(unique(coll)),(nsnps.v+6))
means.lav.coll <- describeBy(lav.coef,trim=0.02,group=coll)
for(i in 1:length(unique(coll))){
  for(j in 1: length(pars)){
    # if theta = 0, use average bias, if not, use average relative bias
    bias.lav.coll[i,j] <- 100*ifelse(pars[j]==0,means.lav.coll[[i]]$trimmed[j]-pars[j],(means.lav.coll[[i]]$trimmed[j]-pars[j])/pars[j])
  }
}
# Average Absolute Bias/Average Absolute Relative Bias
bias.lav.coll2 <- cbind(rowMeans(abs(bias.lav.coll[, 1:a])),
                        rowMeans(abs(bias.lav.coll[, (a+1):(2*a)])),
                        rowMeans(abs(bias.lav.coll[, ((2*a)+1):nsnps.v])),
                        rowMeans(abs(bias.lav.coll[, (nsnps.v+1):(nsnps.v+3)])),
                        rowMeans(abs(bias.lav.coll[, (nsnps.v+4):(nsnps.v+6)])))

bias.lav.coll2.se <- cbind(apply(abs(bias.lav.coll[, 1:a]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.coll[, (a+1):(2*a)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.coll[, ((2*a)+1):nsnps.v]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.coll[, (nsnps.v+1):(nsnps.v+3)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.coll[, (nsnps.v+4):(nsnps.v+6)]),1,sd)/sqrt(n))

(mat11 <- round(bias.lav.coll2,3))
row.names(mat11) <- c(0, 0.1, 0.2)
colnames(mat11) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
(mat11.se <- round(bias.lav.coll2.se,3))
row.names(mat11.se) <- c(0, 0.1, 0.2)
colnames(mat11.se) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
xtable(mat11,"Lavaan Collinearity Relative Bias")

############### SAMPLE SIZE ##################

bias.lav.samp <- matrix(NA,length(unique(samp1)),(nsnps.v+6))
means.lav.samp <- describeBy(lav.coef,trim=0.02,group=samp1)
for(i in 1:length(unique(samp1))){
  for(j in 1: length(pars)){
    bias.lav.samp[i,j] <- 100*ifelse(pars[j]==0,means.lav.samp[[i]]$trimmed[j]-pars[j],(means.lav.samp[[i]]$trimmed[j]-pars[j])/pars[j])
  }
}

bias.lav.samp2 <- cbind(rowMeans(abs(bias.lav.samp[, 1:a])),
                        rowMeans(abs(bias.lav.samp[, (a+1):(2*a)])),
                        rowMeans(abs(bias.lav.samp[, ((2*a)+1):nsnps.v])),
                        rowMeans(abs(bias.lav.samp[, (nsnps.v+1):(nsnps.v+3)])),
                        rowMeans(abs(bias.lav.samp[, (nsnps.v+4):(nsnps.v+6)])))

bias.lav.samp2.se <- cbind(apply(abs(bias.lav.samp[, 1:a]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.samp[, (a+1):(2*a)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.samp[, ((2*a)+1):nsnps.v]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.samp[, (nsnps.v+1):(nsnps.v+3)]),1,sd)/sqrt(n),
                           apply(abs(bias.lav.samp[, (nsnps.v+4):(nsnps.v+6)]),1,sd)/sqrt(n))

(mat22 <- round(bias.lav.samp2,3))
row.names(mat22) <- c(1000, 4000, 7000, 10000)
colnames(mat22) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
(mat22.se <- round(bias.lav.samp2.se,3))
row.names(mat22.se) <- c(1000, 4000, 7000, 10000)
colnames(mat22.se) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
xtable(mat22,"Lavaan Sample Relative Bias")

############### regSEM ###################

############## COLLINEARTY #################

bias.reg.coll <- matrix(NA,length(unique(coll)),(nsnps.v+6))
means.reg.coll <- describeBy(reg.coef,trim=0.02,group=coll)
for(i in 1:length(unique(coll))){
  for(j in 1: length(pars)){
    bias.reg.coll[i,j] <- 100*ifelse(pars[j]==0,means.reg.coll[[i]]$trimmed[j]-pars[j],(means.reg.coll[[i]]$trimmed[j]-pars[j])/pars[j])
  }
}

bias.reg.coll2 <- cbind(rowMeans(abs(bias.reg.coll[, 1:a])),
                        rowMeans(abs(bias.reg.coll[, (a+1):(2*a)])),
                        rowMeans(abs(bias.reg.coll[, ((2*a)+1):nsnps.v])),
                        rowMeans(abs(bias.reg.coll[, (nsnps.v+1):(nsnps.v+3)])),
                        rowMeans(abs(bias.reg.coll[, (nsnps.v+4):(nsnps.v+6)])))

bias.reg.coll2.se <- cbind(apply(abs(bias.reg.coll[, 1:a]),1,sd)/sqrt(n),
                          apply(abs(bias.reg.coll[, (a+1):(2*a)]),1,sd)/sqrt(n),
                          apply(abs(bias.reg.coll[, ((2*a)+1):nsnps.v]),1,sd)/sqrt(n),
                          apply(abs(bias.reg.coll[, (nsnps.v+1):(nsnps.v+3)]),1,sd)/sqrt(n),
                          apply(abs(bias.reg.coll[, (nsnps.v+4):(nsnps.v+6)]),1,sd)/sqrt(n))

(mat33 <- round(bias.reg.coll2,3))
row.names(mat33) <- c(0, 0.1, 0.2)
colnames(mat33) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
(mat33.se <- round(bias.reg.coll2.se,3))
row.names(mat33.se) <- c(0, 0.1, 0.2)
colnames(mat33.se) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
xtable(mat33,"Regsem lasso Collinearity Relative Bias")


############### SAMPLE SIZE ####################

bias.reg.samp <- matrix(NA,length(unique(samp1)),(nsnps.v+6))
means.reg.samp <- describeBy(reg.coef,trim=0.02,group=samp1)
for(i in 1:length(unique(samp1))){
  for(j in 1: length(pars)){
    bias.reg.samp[i,j] <- 100*ifelse(pars[j]==0,means.reg.samp[[i]]$trimmed[j]-pars[j],(means.reg.samp[[i]]$trimmed[j]-pars[j])/pars[j])
  }
}

bias.reg.samp2 <- cbind(rowMeans(abs(bias.reg.samp[, 1:a])),
                        rowMeans(abs(bias.reg.samp[, (a+1):(2*a)])),
                        rowMeans(abs(bias.reg.samp[, ((2*a)+1):nsnps.v])),
                        rowMeans(abs(bias.reg.samp[, (nsnps.v+1):(nsnps.v+3)])),
                        rowMeans(abs(bias.reg.samp[, (nsnps.v+4):(nsnps.v+6)])))

bias.reg.samp2.se <- cbind(apply(abs(bias.reg.samp[, 1:a]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.samp[, (a+1):(2*a)]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.samp[, ((2*a)+1):nsnps.v]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.samp[, (nsnps.v+1):(nsnps.v+3)]),1,sd)/sqrt(n),
                           apply(abs(bias.reg.samp[, (nsnps.v+4):(nsnps.v+6)]),1,sd)/sqrt(n))

(mat44 <- round(bias.reg.samp2,3))
row.names(mat44) <- c(1000, 4000, 7000, 10000)
colnames(mat44) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
(mat44.se <- round(bias.reg.samp2.se,3))
row.names(mat44.se) <- c(1000, 4000, 7000, 10000)
colnames(mat44.se) <- c("noise","0.001/0.005","0.01/0.05","loadings","Item RV")
xtable(mat44,"Regsem lasso Sample Size Relative Bias")


#################################################
########****** variable selection *******########
#################################################

###### REGSEM ######

# REGRESSIONS

null <- 1*as.matrix(abs(reg.coef[,1:a]) >= .001)
null.samp.reg <- cbind(c(1000, 4000, 7000, 10000),
                       matrix(unlist(describeBy(rowMeans(null),group=samp1)),4,13,byrow=T)[,3])
null.samp.reg.se <- cbind(c(1000, 4000, 7000, 10000),
                          matrix(unlist(describeBy(apply(null,1,sd)/sqrt(n),group=samp1)),4,13,byrow=T)[,3])
null.coll.reg <- cbind(c(0, 0.1, 0.2),
                       matrix(unlist(describeBy(rowMeans(null),group=coll)),3,13,byrow=T)[,3])
null.coll.reg.se <- cbind(c(0, 0.1, 0.2),
                          matrix(unlist(describeBy(apply(null,1,sd)/sqrt(n),group=coll)),3,13,byrow=T)[,3])

small <- 1*as.matrix(abs(reg.coef[,(a+1):(2*a)]) < .001)
small.samp.reg <- cbind(c(1000, 4000, 7000, 10000),
                        matrix(unlist(describeBy(rowMeans(small),group=samp1)),4,13,byrow=T)[,3])
small.samp.reg.se <- cbind(c(1000, 4000, 7000, 10000),
                           matrix(unlist(describeBy(apply(small,1,sd)/sqrt(n),group=samp1)),4,13,byrow=T)[,3])
small.coll.reg <- cbind(c(0, 0.1, 0.2),
                        matrix(unlist(describeBy(rowMeans(small),group=coll)),3,13,byrow=T)[,3])
small.coll.reg.se <- cbind(c(0, 0.1, 0.2),
                           matrix(unlist(describeBy(apply(small,1,sd)/sqrt(n),group=coll)),3,13,byrow=T)[,3])

med <- 1*as.matrix(abs(reg.coef[,((2*a)+1):nsnps.v]) < .001)
med.samp.reg <- cbind(c(1000, 4000, 7000, 10000),
                      matrix(unlist(describeBy(rowMeans(med),group=samp1)),4,13,byrow=T)[,3])
med.samp.reg.se <- cbind(c(1000, 4000, 7000, 10000),
                         matrix(unlist(describeBy(apply(med,1,sd)/sqrt(n),group=samp1)),4,13,byrow=T)[,3])
med.coll.reg <- cbind(c(0, 0.1, 0.2),
                      matrix(unlist(describeBy(rowMeans(med),group=coll)),3,13,byrow=T)[,3])
med.coll.reg.se <- cbind(c(0, 0.1, 0.2),
                         matrix(unlist(describeBy(apply(med,1,sd)/sqrt(n),group=coll)),3,13,byrow=T)[,3])


#################################################
##############****** ERRORS *******##############
#################################################

#REGRESSIONS

reg.samp.err <- cbind(null.samp.reg, small.samp.reg, med.samp.reg)[,c(1,2,4,6)]
reg.coll.err <- cbind(null.coll.reg, small.coll.reg, med.coll.reg)[,c(1,2,4,6)]

colnames(reg.samp.err) <- c("Sample Size","Noise",".001/.005",".01/.05")
colnames(reg.coll.err) <- c("Collinearity","Noise",".001/.005",".01/.05")

reg.samp.err.se <- cbind(null.samp.reg.se, small.samp.reg.se, med.samp.reg.se)[,c(1,2,4,6)]
reg.coll.err.se <- cbind(null.coll.reg.se, small.coll.reg.se, med.coll.reg.se)[,c(1,2,4,6)]

colnames(reg.samp.err.se) <- c("Sample Size","Noise",".001/.005",".01/.05")
colnames(reg.coll.err.se) <- c("Collinearity","Noise",".001/.005",".01/.05")

xtable(reg.samp.err)
xtable(reg.coll.err)


############## LAVAAN #############


null.lav <- 1*as.matrix(lav.pvalue[,1:a] < .05)
null.samp.lav <- cbind(c(1000, 4000, 7000, 10000),
                       matrix(unlist(describeBy(rowMeans(null.lav),group=samp1)),4,13,byrow=T)[,3])
null.samp.lav.se <- cbind(c(1000, 4000, 7000, 10000),
                          matrix(unlist(describeBy(apply(null.lav,1,sd)/sqrt(n),group=samp1)),4,13,byrow=T)[,3])
null.coll.lav <- cbind(c(0, 0.1, 0.2),
                       matrix(unlist(describeBy(rowMeans(null.lav),group=coll)),3,13,byrow=T)[,3])
null.coll.lav.se <- cbind(c(0, 0.1, 0.2),
                          matrix(unlist(describeBy(apply(null.lav,1,sd)/sqrt(n),group=coll)),3,13,byrow=T)[,3])

small.lav <- 1*as.matrix(lav.pvalue[,(a+1):(2*a)] > .05)
small.samp.lav <- cbind(c(1000, 4000, 7000, 10000),
                        matrix(unlist(describeBy(rowMeans(small.lav),group=samp1)),4,13,byrow=T)[,3])
small.samp.lav.se <- cbind(c(1000, 4000, 7000, 10000),
                           matrix(unlist(describeBy(apply(small.lav,1,sd)/sqrt(n),group=samp1)),4,13,byrow=T)[,3])
small.coll.lav <- cbind(c(0, 0.1, 0.2),
                        matrix(unlist(describeBy(rowMeans(small.lav),group=coll)),3,13,byrow=T)[,3])
small.coll.lav.se <- cbind(c(0, 0.1, 0.2),
                           matrix(unlist(describeBy(apply(small.lav,1,sd)/sqrt(n),group=coll)),3,13,byrow=T)[,3])

med.lav <- 1*as.matrix(lav.pvalue[,((2*a)+1):nsnps.v] > .05)
med.samp.lav <- cbind(c(1000, 4000, 7000, 10000),
                      matrix(unlist(describeBy(rowMeans(med.lav),group=samp1)),4,13,byrow=T)[,3])
med.samp.lav.se <- cbind(c(1000, 4000, 7000, 10000),
                         matrix(unlist(describeBy(apply(med.lav,1,sd)/sqrt(n),group=samp1)),4,13,byrow=T)[,3])
med.coll.lav <- cbind(c(0, 0.1, 0.2),
                      matrix(unlist(describeBy(rowMeans(med.lav),group=coll)),3,13,byrow=T)[,3])
med.coll.lav.se <- cbind(c(0, 0.1, 0.2),
                         matrix(unlist(describeBy(apply(med.lav,1,sd)/sqrt(n),group=coll)),3,13,byrow=T)[,3])


##### ERRORS ######

lav.samp.err <- cbind(null.samp.lav, small.samp.lav, med.samp.lav)[,c(1,2,4,6)]
lav.coll.err <- cbind(null.coll.lav, small.coll.lav, med.coll.lav)[,c(1,2,4,6)]

colnames(lav.samp.err) <- c("Sample Size","Noise",".001/.005",".01/.05")
colnames(lav.coll.err) <- c("Collinearity","Noise",".001/.005",".01/.05")

lav.samp.err.se <- cbind(null.samp.lav.se, small.samp.lav.se, med.samp.lav.se)[,c(1,2,4,6)]
lav.coll.err.se <- cbind(null.coll.lav.se, small.coll.lav.se, med.coll.lav.se)[,c(1,2,4,6)]

colnames(lav.samp.err.se) <- c("Sample Size","Noise",".001/.005",".01/.05")
colnames(lav.coll.err.se) <- c("Collinearity","Noise",".001/.005",".01/.05")

xtable(lav.samp.err)
xtable(lav.coll.err)

# Combine all

all.objects=list(mat1,mat2,mat3,mat4,mat11,mat22,mat33,mat44,
                 lav.samp.err,lav.coll.err,reg.samp.err,reg.coll.err,
                 mat1.se,mat2.se,mat3.se,mat4.se,mat11.se,mat22.se,mat33.se,mat44.se,
                 lav.samp.err.se,lav.coll.err.se,reg.samp.err.se,reg.coll.err.se)

save(all.objects,file=paste0("/Users/isobelridler/Desktop/RosalindHPC/outputGenMIMIC/tablesGenMIMICrep100_3item", nsnps.v,".RData"))

}
