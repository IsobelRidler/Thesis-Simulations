## simulated genomic MIMIC model and data with regsem ----------------------------------------
############### based on Jacobucci et al. code on OSF, adapted to suit my model scenario and simulation ###############

rm(list=ls())

arrayID <- Sys.getenv("SLURM_ARRAY_TASK_ID")
whichRep <- as.numeric(arrayID)
cat("Starting run whichRep = ",whichRep,"\n")
if(is.na(whichRep)) whichRep <- 1

# load R packages
library(lavaan)
library(regsem)
library(tidyverse)

# grid for indexing array jobs
val <- c(0, 0.1, 0.2) #3 collinearities
samp <- c(1000, 4000, 7000, 10000) #4 sample sizes
nsnps <- c(30, 60, 90, 120) #no. of SNPs included
grid <- expand.grid(val = val,samp = samp, nsnps = nsnps)
tot <- length(val)*length(samp)*length(nsnps)
grid <- grid %>% 
  mutate(repl = 1:tot)

# index using whichRep
val1 <- grid[whichRep,"val"]
samp1 <- grid[whichRep,"samp"]
nsnps1 <- grid[whichRep, "nsnps"]
repl <- grid[whichRep,"repl"]

# number of times each cell is repeated (Monte Carlo), 5 for test run, eventually found ~200 sufficient.. use 1000 in final as smaller MCSE
nsReps <- 200

do.one <- function(srep,val1,samp1,nsnps1,repl,nsReps){
  ################################################################################
  ### Data Generation
  ################################################################################
  
  trueRep <- nsReps*(repl - 1) + srep

### sim pop model spec ###

# item loadings
w1 <- 0.6
w2 <- 0.4
w3 <- 0.3
# SNP beta weights
w4 <- 0.05
w5 <- 0.01
w6 <- 0.005
w7 <- 0.001
w8 <- 0

# no. of snps with certain weight
subsnp <- nsnps1/6

sim.list = list()

sim.mod.MIMIC <- paste0(paste0("f1 ~ ", paste0(c(rep(w8, subsnp), rep(w8, subsnp), rep(w7, subsnp), rep(w6, subsnp), rep(w5, subsnp), rep(w4, subsnp)), "*", "x", 1:nsnps1, collapse = " + "), collapse = ""), "\n", paste0("f1 =~ ", paste0(c(rep(w1, 1), rep(w2, 1), rep(w3, 1)), "*", "x", (nsnps1+1):(nsnps1+3), collapse = " + "), collapse = ""), collapse = "")

sim.list[[1]] = sim.mod.MIMIC

xxx = list()
for(i in 1:(nsnps1+3)){
  xxx[i] = paste("x", i, sep = "")
}

xx1 = list()
for(i in 1:nsnps1){
  xx1[i] = paste("x", i, sep = "")
}

jj1 = list()
for(i in 1:(length(xx1))){
  jj1[i] = paste(xx1[i],"~~",1,"*",xx1[i],sep="")
}
sim.list[[2]] = paste(jj1,collapse=";")

kk1 = list()
count=0
for(i in 1:nsnps1){
  for(j in 1:nsnps1){
    if(i != j & j > i){
      count = count+1
      kk1[count] = paste(xx1[i], "~~", val1, "*", xx1[j], sep="")
    }
  }
}
sim.list[[3]] = paste(kk1, collapse=";")
sim.list[[4]] = "f1 ~~ 1*f1"

pop.mod = " "
for(i in 1:length(sim.list)){
  pop.mod = paste(pop.mod, sim.list[[i]], sep="\n")
}

# create data from pop.mod (generates z-score multivariate normal data)
dat <- simulateData(pop.mod,sample.nobs=samp1,seed=trueRep,model.type="sem")

# transform items into ordinal 4-category data (3 splits, equal interval)
# dat$catx1 <- cut(dat$x1, breaks = 4, labels = c(0:3))
## using mutate and across dplyr
# for all items
dat.c <- dat %>% 
  mutate(
    across(
      .cols = c(-(4:(3+nsnps1))),
      .fns  = ~ cut(.x, breaks = 4, labels = c(0:3)),
      .names = "{col}_c"
    )
  )


# transform SNPs into 0,1,2 ordinal categories (0 being q^2, 1 being 2pq, and 2 being p^2 proportions; where p = MAF, q = 1-p and p^2 + 2pq + q^2 = 1) - qnorm(MAF) would give percentage of one allele but not the genotype/haplotype
# choose MAF between 0.01-0.45/0.49 (?)
# dat$snp31 <- cut(dat$x31, breaks = c(-4, qnorm(q^2), -qnorm(p^2), 4), labels = c(0, 1, 2))
# for all snps, but want different random MAFs for each SNP, so need a MAP function (inside or wrapped around?) which randomly takes a number from a vector/list of numbers from 0.01-0.49 - is this possible in purrr or do I need to stick to for loops/lapply?

### TIDYVERSE/PURRR WAYS ####

# most consistent (with how we categorised x1-30) way using function and mutate/across
# p = MAF
snp_cut <- function (x){
  p <- runif(1, 0.01, 0.5)
  q <- 1-p
  cut(x, breaks = c(-6, qnorm(q^2), -qnorm(p^2), 6), labels = c(0, 1, 2))
}

set.seed(1234)
dat.c <- dat.c %>% 
  mutate(
    across(
      .cols = (4:(3+nsnps1)),
      .fns  = ~ snp_cut(.x),
      .names = "{col}_c"
    )
  )


# orrrr, make numeric, then
#scale
dat.n <- dat.c %>%
  mutate_if(is.factor, as.numeric)
dat.n <- as_tibble(scale(dat.n))

### run model spec ###

xxx = list()
for(i in 1:(nsnps1+3)){
  xxx[i] = paste("x", i, "_c", sep = "")
}

xx1 = list()
for(i in 1:nsnps1){
  xx1[i] = paste("x", i, "_c", sep = "")
}

xx2 = list()
for(i in (nsnps1+1):(nsnps1+3)){
  xx2[i] = paste("x", i, "_c", sep = "")
}

cc1 = list()
for(i in 1:nsnps1){
  cc1[i] = paste("c", i, sep = "")
}

cc2 = list()
for(i in (nsnps1+1):(nsnps1+3)){
  cc2[i] = paste("c", i, sep = "")
}

ll3 = list()
for(i in 1:nsnps1){
  ll3[i] = paste(cc1[[i]],"*",xx1[[i]],sep="")
}

ll4 = list()
for(i in (nsnps1+1):(nsnps1+3)){
  ll4[i] = paste(cc2[[i]],"*",xx2[[i]],sep="")
}

ll4 <- ll4 %>% discard(is.null)

load.list3 = paste(ll3, collapse = "+")
load.list4 = paste(ll4, collapse = "+")

run.list = list()
run.list[[1]] = paste0("f1 ~ ", load.list3, collapse = "")
run.list[[2]] = paste0("f1 =~ ", load.list4, collapse = "")
run.list[[3]] = "f1 ~~ 1*f1"

run.mod = " "
for(k in 1:length(run.list)){
  run.mod = paste(run.mod,run.list[[k]],sep="\n")
}

## run model using lavaan
# do.fit: If FALSE, the model is not fit, and the current starting values of the model parameters are preserved
# fixed.x: If TRUE, the exogenous ‘x’ covariates are considered fixed variables and the means, variances and covariances of these variables are fixed to their sample values. If FALSE, they are considered random, and the means, variances and covariances are free parameters. If "default", the value is set depending on the mimic option.
# std.lv: If TRUE, the metric of each latent variable is determined by fixing their variances to 1.0. If FALSE, the metric of each latent variable is determined by fixing the factor loading of the first indicator to 1.0.
# start = simple: all parameter values are set to zero, except the factor loadings (set to one), the variances of latent variables (set to 0.05), and the residual variances of observed variables (set to half the observed variance)
lav.out1 <- try(lavaan::sem(run.mod,dat.s1,fixed.x=TRUE, std.lv=TRUE),silent=FALSE)
## summary(lav.out)

if(inherits(lav.out,"try-error")){
  lav.pvalue <- NA
  lav.coef <- NA
  lasso.sim <- NA
  lav.out <- lavaan::sem(run.mod,dat.n,fixed.x=TRUE,std.lv=TRUE,do.fit=FALSE)
}else{
  lav.pvalue = parameterestimates(lav.out)[c(1:(nsnps1+3), (nsnps1+5):(nsnps1+7)), "pvalue"]
  lav.coef = coef(lav.out)
  
  lasso.sim <- try(cv_regsem(lav.out,n.lambda=20,jump=.01,type="lasso",
                             pars_pen=cc1,optMethod="rsolnp",
                             fit.ret=c("BIC","rmsea","AIC","CAIC","EBIC.5","EBIC.25"),
                             fit.ret2="train",warm.start=FALSE),silent=TRUE)
  if(inherits(lasso.sim,"try-error")){
    lasso.sim <- NA
  }
}



out <- list(lasso.sim = lasso.sim,
            lav.pvalue = lav.pvalue,
            lav.coef = lav.coef,
            val1 = val1,
            samp1 = samp1,
            trueRep = trueRep)
return(out)
}


library(parallel)

detectBatchCPUs <- function() { 
  ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) 
  if (is.na(ncores)) { 
    ncores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) 
  } 
  if (is.na(ncores)) { 
    return(2) # default
  } 
  return(ncores) 
}

ncpus <- detectBatchCPUs()

out <- mclapply(1:nsReps,do.one,val1 = val1,samp1 = samp1,repl = repl,nsReps = nsReps, mc.cores = ncpus)

setwd("/scratch/users/k1620435/output")
saveRDS(out,file = paste0("chain03_SNP_MIMICsim_",whichRep,".rds"))

