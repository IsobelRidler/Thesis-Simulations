############### lasso Try 5, based on Jacobucci et al. code on OSF, adapted to suit my model scenario and simulation ###############

rm(list=ls())

arrayID <- Sys.getenv("SLURM_ARRAY_TASK_ID")
whichRep <- as.numeric(arrayID)
cat("Starting run whichRep = ",whichRep,"\n")
if(is.na(whichRep)) whichRep <- 1

library(lavaan)
library(regsem)
library(tidyverse)

# grid for indexing array jobs
val <- c(.20,.30,.50,.70,.85) #5 collinearities
samp <- c(100, 150, 200, 500, 1000, 2000) #6 sample sizes
grid <- expand.grid(val = val,samp = samp)
grid <- grid %>% 
  mutate(repl = 1:30)


# index using whichRep
val1 <- grid[whichRep,"val"]
samp1 <- grid[whichRep,"samp"]
repl <- grid[whichRep,"repl"]

# number of times each cell is repeated (Monte Carlo), 5 for test run, eventually found ~200 sufficient.. use 1000 in final as smaller MCSE
nsReps <- 1000

# function to apply over array

do.one <- function(srep,val1,samp1,repl,nsReps){
  ################################################################################
  ### Data Generation
  ################################################################################
  
  trueRep <- nsReps*(repl - 1) + srep 
  
  # simulated pop model specification
  
  load.reg.model <- paste0(paste0("f1 =~ ", paste0(paste0("0.3*x", 1:4, collapse = " + "), " + ", paste0("0.7*x", 5:8, collapse = " + "), " + ", paste0("0.5*x", 9:15, collapse = " + "), collapse = " + "), collapse = ""), "\n", paste0("f2 =~ ", paste0(paste0("0.3*x", 16:19, collapse = " + "), " + ", paste0("0.7*x", 20:23, collapse = " + "), " + ", paste0("0.5*x", 24:30, collapse = " + "), collapse = " + "), collapse = " + "), "\n", "y ~ 0.4*f1 + 0.6*f2 + 0*x1 + 0.2*x2 + 0.4*x3 + 0.6*x4 + 0*x5 + 0.2*x6 + 0.4*x7 + 0.6*x8 + 0*x16 + 0.2*x17 + 0.4*x18 + 0.6*x19 + 0*x20 + 0.2*x21 + 0.4*x22 + 0.6*x23", "\n", "f1 ~~ 1*f1", "\n", "f2 ~~ 1*f2", "\n", "f1 ~~ 0.3*f2", collapse = "")
  
  
  sim.list = list()
  
  
  sim.list[[1]] = load.reg.model
  
  xxx = list()
  for(i in 1:30){
    xxx[i] = paste("x", i, sep = "")
  }
  
  
  xx1 = list()
  for(i in 1:15){
    xx1[i] = paste("x", i, sep = "")
  }
  
  xx2 = list()
  for(i in 16:30){
    xx2[i] = paste("x", i, sep = "")
  }
  
  # set collinearities between items for each factor 
  kk1 = list()
  count=0
  for(i in 1:15){
    for(j in 1:15){
      if(i != j & j > i){
        count = count+1
        kk1[count] = paste(xx1[i],"~~",val1,"*",xx1[j],sep="")
      }
    }
  }
  
  kk2 = list()
  count=0
  for(i in 16:30){
    for(j in 16:30){
      if(i != j & j > i){
        count = count+1
        kk2[count] = paste(xx2[i],"~~",val1,"*",xx2[j],sep="")
      }
    }
  }
  
  sim.list[[2]] = paste(kk1,collapse=";")
  sim.list[[3]] = paste(kk2,collapse=";")
  
  pop.mod = " "
  for(i in 1:length(sim.list)){
    pop.mod = paste(pop.mod,sim.list[[i]],sep="\n")
  }
  
  # run model specification
  
  ll1 = list()
  for(i in 1:15){
    ll1[i] = paste("start(0.5)","*",xx1[[i]],sep="")
  }
  ll2 = list()
  for(i in 16:30){
    ll2[i] = paste("start(0.5)","*",xx2[[i]],sep="")
  }
  
  ll2 <- ll2 %>% discard(is.null)
  
  load.list1 = paste(ll1, collapse = "+")
  load.list2 = paste(ll2, collapse = "+")
  
  cc1 = list()
  for(i in 1:15){
    cc1[i] = paste("c", i, sep = "")
  }
  
  cc2 = list()
  for(i in 16:30){
    cc2[i] = paste("c", i, sep = "")
  }
  
  
  ll3 = list()
  for(i in 1:15){
    ll3[i] = paste(cc1[[i]],"*",xx1[[i]],sep="")
  }
  ll4 = list()
  for(i in 16:30){
    ll4[i] = paste(cc2[[i]],"*",xx2[[i]],sep="")
  }
  
  ll4 <- ll4 %>% discard(is.null)
  
  load.list3 = paste(ll3, collapse = "+")
  load.list4 = paste(ll4, collapse = "+")
  
  
  
  run.list = list()
  run.list[[1]] = paste0("f1 =~ ", load.list1, " + ", load.list3, collapse = "")
  run.list[[2]] = paste0("f2 =~ ", load.list2, " + ", load.list4, collapse = "")
  run.list[[3]] = "y ~ start(0.5)*f1 + start(0.5)*f2 + e1*f1 + e2*f2 + d1*x1 + d2*x2 + d3*x3 + d4*x4 + d5*x5 + d6*x6 + d7*x7 + d8*x8 + d9*x16 + d10*x17 + d11*x18 + d12*x19 + d13*x20 + d14*x21 + d15*x22 + d16*x23"
  run.list[[4]] = "f1 ~~ 1*f1"
  run.list[[5]] = "f2 ~~ 1*f2"
  run.list[[6]] = "f1 ~~ e3*f2"
  
  run.mod = " "
  for(k in 1:length(run.list)){
    run.mod = paste(run.mod,run.list[[k]],sep="\n")
  }
  
  # parameters to penalise in regsem
  penparam <- c("d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", "d9", "d10", "d11", "d12", "d13", "d14", "d15", "d16", "e1", "e2")
  
  # create data from pop.mod
  dat <- simulateData(pop.mod,sample.nobs=samp1,seed=trueRep,model.type="lavaan")
  
  ## run model using lavaan
  
  lav.out <- try(lavaan::sem(run.mod,dat,fixed.x=FALSE, std.lv=TRUE),silent=FALSE)
  ## summary(lav.out)
  
  if(inherits(lav.out,"try-error")){
    lav.pvalue <- NA
    lav.coef <- NA
    lasso.sim <- NA
    lav.out <- lavaan::sem(run.mod,dat,fixed.x=FALSE, std.lv = TRUE,do.fit=FALSE)
  }else{
    lav.pvalue = parameterestimates(lav.out)[c(1:48, 51:82), "pvalue"]
    lav.coef = coef(lav.out)
    
    lasso.sim <- try(cv_regsem(lav.out,n.lambda=30,jump=.01,type="lasso",
                               pars_pen=penparam,optMethod="rsolnp",
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
saveRDS(out,file = paste0("chain05e_",whichRep,".rds"))

