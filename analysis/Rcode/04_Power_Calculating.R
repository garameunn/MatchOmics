source("00_library_setup.R")
source("99_utils.R")
source("02_simulation_generator.R") # 
bat.dat <- read.csv("data/batchinfo_full.csv", header = T)

load(file="data/otulist.RData") # otulist (otu (1000/500/200), lbsize (1000))
load(file="data/indiclist.RData") # indiclist
load(file="data/psdifflist.RData") # psdifflist
nsamchar <- c(1000,500,200) # norder=1


#####################################
## Scenario: 1
## Statistic : Power
## Topic : Calculating
#####################################

#### nsamp change
nsamp = 200

targetnum <- which(nsamchar==nsamp)

otutable <- otulist[[targetnum]]
indicator <- indiclist[[targetnum]]
psdiff <- psdifflist[[targetnum]]
lbs = otulist[[4]][1:nsamp]

b = c(0.001,0.002,0.005, 0.01)
k <- c(1,8,61,77)

#### scenario change
scenario="S1"

#### rate, rep, nonlin change
set.seed(1)
rep=500
dimm=dim(otutable)[1]; omics="Metagenomics"
datatable=otutable;rate=0.5;nonlin="comp" # if Scenario is not S1, nonlin = NULL
trsftable=trnsf(dt=datatable, omics=omics, lb.size=lbs);
source("02_simulation_generator.R") 
# h<-x<-i<-1

bat.dat <- read.csv("/home2/nekim/scratch/NOTOfull/batchinfo_full.csv", header = T)

batord <- match(colnames(otutable),bat.dat[,1])
bat.dat1 <- bat.dat[batord,]
rownames(bat.dat1) <- bat.dat1[,1]
cmbatable=ComBat(trsftable,batch=as.character(bat.dat1$Batch),par.prior = T,prior.plots = F) 
limtable <- removeBatchEffect(trsftable, batch=bat.dat1$Batch)

pvalues <- mclapply(1:rep, function(h){
  if(h%%10==1){print(h)}
  sumsig <- (mclapply(1:dimm, function(x){ 
    
    indx =dimm*(h-1)+x 
    f(indx)
    
    
    tablefun <- get("tablemaking_power")
    tablefun1 <- get("tablemaking_power_wo")
    
    
    tble0 <- tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=NULL, scenario=scenario, k=NULL, nonlin=nonlin)  
    tble0_.05 <-  tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.05, scenario=scenario, k=NULL, nonlin=nonlin)
    tble0_.1 <- tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.1, scenario=scenario, k=NULL, nonlin=nonlin)
    tble0_.2 <-  tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.2, scenario=scenario, k=NULL, nonlin=nonlin)
    tble0_.3 <- tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.3, scenario=scenario, k=NULL, nonlin=nonlin)
    
    tble0_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=NULL, scenario=scenario, k=NULL, nonlin=nonlin) 
    tble0_.05_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.05, scenario=scenario, k=NULL, nonlin=nonlin)  
    tble0_.1_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.1, scenario=scenario, k=NULL, nonlin=nonlin)
    tble0_.2_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.2, scenario=scenario, k=NULL, nonlin=nonlin)  
    tble0_.3_wo <- tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.3, scenario=scenario, k=NULL, nonlin=nonlin)
    
    
    # 
    r.sig <- coeff <- r.sig2 <- coeff2 <- c() # withrep nocal (Indep - EX)
    r.sig0.1 <- coeff0.1 <- r.sig20.1 <- coeff20.1 <- c() # withrep cal0.1 
    r.sig0.2 <- coeff0.2 <- r.sig20.2 <- coeff20.2 <- c() # withrep cal0.2 
    r.sig0.3 <- coeff0.3 <- r.sig20.3 <- coeff20.3 <- c() # withrep cal0.3
    r.sig0.05 <- coeff0.05 <- r.sig20.05 <- coeff20.05 <- c() # withrep cal0.05
    
    r.sig_w <- coeff_w <- r.sig2_w <- coeff2_w <- c() # without rep nocal (Indep - EX)
    r.sig0.1_w <- coeff0.1_w <- r.sig20.1_w <- coeff20.1_w <- c() # withrep cal0.1 
    r.sig0.2_w <- coeff0.2_w <- r.sig20.2_w <- coeff20.2_w <- c() # withrep cal0.2 
    r.sig0.3_w <- coeff0.3_w <- r.sig20.3_w <- coeff20.3_w <- c() # withrep cal0.3
    r.sig0.05_w <- coeff0.05_w <- r.sig20.05_w <- coeff20.05_w <- c() # withrep cal0.05
    
    logit.sig <- logit.coef <- logit.sig1 <- logit.coef1 <- combat.sig <- combat.coef <- limm.sig <- limm.coef <- c()
    
    
    for (i in 1:length(b)){
      
      tble <- data.frame(tble0[[1]][[i]])
      tble0.1 <- tble0_.1[[1]][[i]]
      tble0.2 <- tble0_.2[[1]][[i]]
      tble0.3 <- tble0_.3[[1]][[i]]
      tble0.05 <- tble0_.05[[1]][[i]]
      tble_wo <- tble0_wo[[1]][[i]]
      tble0.1_wo <- tble0_.1_wo[[1]][[i]]
      tble0.2_wo <- tble0_.2_wo[[1]][[i]]
      tble0.3_wo <- tble0_.3_wo[[1]][[i]]
      tble0.05_wo <- tble0_.05_wo[[1]][[i]]
      
      tble2 <- data.frame(tble0[[2]][[i]]) # entire data
      tble2$limma <- unlist(limtable[x,])
      tble2$cmbt <- unlist(cmbatable[x,])
      
      
      #######################
      ###### with rep #######
      #######################
      # no caliper
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble, family="binomial", 
                                    corstr = "independence", corr.mat = NULL,id = tble$pair_taxa, weights = tble$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig <- c(r.sig,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff <- c(coeff,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble, family="binomial", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble$pair_taxa, weights = tble$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig2 <- c(r.sig2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff2 <- c(coeff2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      rm(geemfit);rm(geemfit2); 
      
      # 0.2 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2, family="binomial", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.2$pair_taxa, weights = tble0.2$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.2 <- c(r.sig0.2,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.2 <- c(coeff0.2,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2, family="binomial", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.2$pair_taxa, weights = tble0.2$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.2 <- c(r.sig20.2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.2 <- c(coeff20.2,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      
      rm(geemfit);rm(geemfit2); 
      
      # 0.1 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1, family="binomial", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.1$pair_taxa, weights = tble0.1$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.1 <- c(r.sig0.1,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.1 <- c(coeff0.1,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1, family="binomial", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.1$pair_taxa, weights = tble0.1$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.1 <- c(r.sig20.1,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.1 <- c(coeff20.1,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      rm(geemfit);rm(geemfit2);
      
      # 0.3 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3, family="binomial", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.3$pair_taxa, weights = tble0.3$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.3 <- c(r.sig0.3,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.3 <- c(coeff0.3,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3, family="binomial", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.3$pair_taxa, weights = tble0.3$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.3 <- c(r.sig20.3,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.3 <- c(coeff20.3,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      
      
      rm(geemfit);rm(geemfit2);
      
      # 0.05 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05, family="binomial", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.05$pair_taxa, weights = tble0.05$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.05 <- c(r.sig0.05,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.05 <- c(coeff0.05,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05, family="binomial", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.05$pair_taxa, weights = tble0.05$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.05 <- c(r.sig20.05,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.05 <- c(coeff20.05,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      
      
      #########################
      ###### without rep ######
      #########################
      rm(geemfit);rm(geemfit2); 
      # no caliper
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble_wo, family="binomial", 
                                    corstr = "independence", corr.mat = NULL,id = tble_wo$pair_taxa, weights = tble_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig_w <- c(r.sig_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff_w <- c(coeff_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble_wo, family="binomial", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble_wo$pair_taxa, weights = tble_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig2_w <- c(r.sig2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff2_w <- c(coeff2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      rm(geemfit);rm(geemfit2);
      
      # 0.2 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2_wo, family="binomial", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.2_wo$pair_taxa, weights = tble0.2_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.2_w <- c(r.sig0.2_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.2_w <- c(coeff0.2_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.2_wo, family="binomial", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.2_wo$pair_taxa, weights = tble0.2_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.2_w <- c(r.sig20.2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.2_w <- c(coeff20.2_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      
      rm(geemfit);rm(geemfit2);
      
      # 0.1 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1_wo, family="binomial", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.1_wo$pair_taxa, weights = tble0.1_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.1_w <- c(r.sig0.1_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.1_w <- c(coeff0.1_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.1_wo, family="binomial", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.1_wo$pair_taxa, weights = tble0.1_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.1_w <- c(r.sig20.1_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.1_w <- c(coeff20.1_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      # 0.3 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3_wo, family="binomial", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.3_wo$pair_taxa, weights = tble0.3_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.3_w <- c(r.sig0.3_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.3_w <- c(coeff0.3_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.3_wo, family="binomial", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.3_wo$pair_taxa, weights = tble0.3_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.3_w <- c(r.sig20.3_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.3_w <- c(coeff20.3_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      
      
      rm(geemfit);rm(geemfit2);
      
      # 0.05 caliper 
      geemfit <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05_wo, family="binomial", 
                                    corstr = "independence", corr.mat = NULL,id = tble0.05_wo$pair_taxa, weights = tble0.05_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig0.05_w <- c(r.sig0.05_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,summary(geemfit)$p[4],1),1))
      coeff0.05_w <- c(coeff0.05_w,ifelse(sum(class(geemfit)=="try-error")==0, ifelse(geemfit$converged&abs(geemfit$beta[4])<100,geemfit$beta[4],NA),NA))
      
      geemfit2 <-  tryCatch(try(geem(formula = as.formula(paste("case~x1+PS","marker",sep="+")), data = tble0.05_wo, family="binomial", 
                                     corstr = "exchangeable", corr.mat = NULL,id = tble0.05_wo$pair_taxa, weights = tble0.05_wo$mm, sandwich = T, useP = T, maxit = 25, tol = 10^-3),silent = T))
      
      r.sig20.05_w <- c(r.sig20.05_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,summary(geemfit2)$p[4],1),1))
      coeff20.05_w <- c(coeff20.05_w,ifelse(sum(class(geemfit2)=="try-error")==0, ifelse(geemfit2$converged&abs(geemfit2$beta[4])<100,geemfit2$beta[4],NA),NA))
      
      
      rm(geemfit);rm(geemfit2);
      
      
      logit <- glm(paste("case~x1+PS","marker",sep="+"), data=tble2, family = "binomial")
      logit.sig <- c(logit.sig,ifelse(logit$converged,coef(summary(logit))[4,4],1))
      logit.coef <- c(logit.coef,ifelse(logit$converged,coef(summary(logit))[4,1],NA))
      rm(logit)
      
      logit1 <- glm(paste("case~x1","marker",sep="+"), data=tble2, family = "binomial")
      logit.sig1 <- ifelse(logit1$converged,coef(summary(logit1))[3,4],1)
      logit.coef1 <- ifelse(logit1$converged,coef(summary(logit1))[3,1],NA)
      rm(logit1)
      
      
      combat <- glm(paste("case~x1+hetero1+hetero2","cmbt",sep="+"), data=tble2, family = "binomial") 
      combat.sig <- c(combat.sig,ifelse(combat$converged,coef(summary(combat))[5,4],1))
      combat.coef <- c(combat.coef,ifelse(combat$converged,coef(summary(combat))[5,1],NA))
      
      limm <- glm(paste("case~x1+hetero1+hetero2","limma",sep="+"), data=tble2, family = "binomial") 
      limm.sig <- c(limm.sig,ifelse(limm$converged,coef(summary(limm))[5,4],1))
      limm.coef <- c(limm.coef,ifelse(limm$converged,coef(summary(limm))[5,1],NA))

    }
    
    return(data.frame(logit.sig1, logit.sig, combat.sig, limm.sig,
                      r.sig, r.sig2,  r.sig0.3, r.sig20.3, r.sig0.2, r.sig20.2, r.sig0.1, r.sig20.1, r.sig0.05, r.sig20.05,
                      r.sig_w, r.sig2_w, r.sig0.3_w, r.sig20.3_w, r.sig0.2_w, r.sig20.2_w, r.sig0.1_w, r.sig20.1_w, r.sig0.05_w, r.sig20.05_w,  #24
                      logit.coef1, logit.coef, combat.coef, limm.coef,
                      coeff, coeff2, coeff0.3, coeff20.3, coeff0.2, coeff20.2, coeff0.1, coeff20.1, coeff0.05, coeff20.05, 
                      coeff_w, coeff2_w, coeff0.3_w, coeff20.3_w, coeff0.2_w, coeff20.2_w, coeff0.1_w, coeff20.1_w, coeff0.05_w, coeff20.05_w)) #24   
    
    
  },mc.cores=8))
},mc.cores = 6)
head(pvalues[[1]])

save(pvalues, file=paste0("data/Power/NSAMP",nsamp,"_REP",rep,"_RATE",rate,"_",nonlin,"_",scenario,".RData")) 





