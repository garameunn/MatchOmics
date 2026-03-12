source("00_library_setup.R")
source("99_utils.R")
source("02_simulation_generator.R") 
bat.dat <- read.csv("data/batchinfo_full.csv", header = T)

load(file="data/otulist.RData") # otulist (otu (1000/500/200), lbsize (1000))
load(file="data/indiclist.RData") # indiclist
load(file="data/psdifflist.RData") # psdifflist
nsamchar <- c(1000,500,200) # norder=1

EffectSizeFun <- function(omics, trsftable, datatable, indicator, rate, b, x,cl, scenario, k=NULL){
  # browser()
  
  otu0 <- data.frame(t(datatable)[,x])
  # taxa up/down labeling
  taxa_class <- ifelse(otu0 >= median(as.numeric(unlist(otu0))), "up", "down")
  namefull <- paste0(taxa_class,"_",rownames(otu0))
  formatch <- data.frame(taxa_class=as.factor(c(taxa_class)), indicator)
  matchform <- as.formula(paste("taxa_class",paste(colnames(indicator),collapse = "+"),sep = "~"))
  
  # matmat <- matchit(matchform, formatch, method = "nearest", replace = T) 
  matmat <- matchit(matchform, formatch, method = "nearest", replace = T, caliper = cl)
  matr <- matmat$match.matrix
  
  subclass <- c()
  for(i in 1:length(unique(na.omit(matr[,1])))){
    nw <- unique(na.omit(matr[,1]))[i]
    new <- c(names(matr[which(matr[,1]==nw),]),nw)
    newsub <- rep(i,length(new))
    names(newsub) <- new
    subclass <- c(subclass,newsub)
  }
  subclass <- as.factor(subclass)
  
  tmp <- data.frame(taxaclass = rep(c("Var1","Var2"),each=nrow(matmat$match.matrix)), #Var1=up
                    samID = c(rownames(matmat$match.matrix), matmat$match.matrix[,1]))
  # browser()
  
  ## matched table making 
  samdat <- data.frame(ps = matmat$distance[match(names(subclass),names(matmat$distance))], 
                       pair_taxa=subclass, 
                       variable=tmp[match(names(subclass),tmp$samID),"taxaclass"], 
                       value=names(subclass))
  
  # Liability score 
  trsftable_matched <- trsftable[,match(samdat[,4], colnames(trsftable))]
  
  if (scenario=="S1"|scenario=="S2"){
    
    if (scenario=="S1"){
      PPS <- 1/(samdat$ps)
    } else if (scenario=="S2"){
      PPS <- indicator[match(samdat[,4], rownames(indicator)),1]+indicator[match(samdat[,4], rownames(indicator)),2]
    }
    means <- c(1, mean(unlist(trsftable_matched[x,])), mean(PPS), 0) # x1, otu, hetero, epsilon
    sds <- c(1, sd(unlist(trsftable_matched[x,])), sd(PPS), 1) # x1, otu, hetero, epsilon
    
    x1 = rnorm(mean=means[1], sd=sds[1], nrow(samdat))
    epsilon <- rnorm(mean=means[4], sd=sds[4], nrow(samdat))
    v23 <- cov(unlist(trsftable_matched[x,]), PPS) 
    
  } else if (scenario=="S3"){
    
    otu_for_generate0 <- trsftable[k,match(samdat[,4], colnames(trsftable))]
    
    otu_for_generate <- matrix(ncol=ncol(trsftable_matched))
    for (ii in 1:length(k)){
      otu_for_generate <- rbind(otu_for_generate,lm(t(otu_for_generate0)[,ii]~trsftable_matched[x,])$residuals)
    }
    otu_for_generate <- otu_for_generate[-1,]
    
    means <- c(1, mean(unlist(trsftable_matched[x,])), apply(otu_for_generate,1,mean), 0) # x1, otutobetested, residuals, epsilon
    sds <- c(1, sd(unlist(trsftable_matched[x,])), apply(otu_for_generate,1,sd), 1) # x1, otutobetested, residuals, epsilon
    
    x1 = rnorm(mean=means[1], sd=sds[1], nrow(samdat))
    epsilon <- rnorm(mean=means[length(means)], sd=sds[length(sds)], nrow(samdat))
    vresid <- do.call(sum,lapply(data.frame(combn(c(1:length(k)),2)), function(y) cov(otu_for_generate[y[1],], otu_for_generate[y[2],])))
    TT = sum(apply(otu_for_generate,1,var))+2*vresid 
    v23 = cov(unlist(trsftable_matched[x,]), apply(otu_for_generate,2,sum))
    
  }
  
  b1=0.1
  
  listtable <- lapply(1:length(b), function(blen){
    # browser()
    # add effect size
    B=b[[blen]]
    
    if (scenario=="S1"|scenario=="S2"){
      
      if (omics=="Metagenomics"){
        b3=sqrt(1.01/(9*sds[3]^2))    
      } else if (omics=="Proteomics"){
        b3=sqrt(1.01/(99*sds[3]^2)) 
      } else if (omics=="Metabolomics") {
        b3=sqrt(3.03/(7*sds[3]^2)) 
      }
      b2=Re(polyroot(c(-1.01*B+B*(b3*sds[3])^2, -2*B*b3*v23, (1-B)*(sds[2]^2))))[1]
      
      } else if (scenario=="S3"){
      
      if (omics=="Metagenomics"){
        b3=sqrt(1.01/(9*TT)) 
      } else if (omics=="Proteomics"){
        b3=sqrt(1.01/(99*TT)) 
      } else if (omics=="Metabolomics") {
        b3=sqrt(3.03/(7*TT)) 
      }
      
      b2=Re(polyroot(c(-1.01*B+B*TT*(b3^2), -2*B*b3*v23, (1-B)*(sds[2]^2))))[1] 
      
    }
    return(b2)
    
  })
  
  
  listtable1 <- lapply(1:length(b), function(blen){
    ## Entire table making 
    B=b[[blen]]
    if (scenario=="S1"|scenario=="S2"){
      if (scenario=="S1"){
        PPS1 <- 1/matmat$distance 
      } else if (scenario=="S2"){
        PPS1 <- indicator[,1]+indicator[,2]
      }
      means1 <- c(1, mean(unlist(trsftable[x,])), mean(PPS1), 0) # x1, otu, hetero, epsilon
      sds1 <- c(1, sd(unlist(trsftable[x,])), sd(PPS1), 1) # x1, otu, hetero, epsilon
      
      x11 = rnorm(mean=means1[1], sd=sds1[1], length(PPS1))
      epsilon1 <- rnorm(mean=means1[4], sd=sds1[4], length(PPS1))
      v231 <- cov(unlist(trsftable[x,]), PPS1) 
      
      if (omics=="Metagenomics"){
        b3=sqrt(1.01/(9*sds1[3]^2)) 
      } else if (omics=="Proteomics"){
        b3=sqrt(1.01/(99*sds1[3]^2)) 
      } else if (omics=="Metabolomics") {
        b3=sqrt(3.03/(7*sds1[3]^2)) 
      }
      
      b2=Re(polyroot(c(-1.01*B+B*(b3*sds1[3])^2, -2*B*b3*v231, (1-B)*(sds1[2]^2))))[1] 
      
    } else if (scenario=="S3"){
      
      otu_for_generate0 <- trsftable[k,]
      
      otu_for_generate <- matrix(ncol=ncol(trsftable))
      for (ii in 1:length(k)){
        otu_for_generate <- rbind(otu_for_generate,lm(t(otu_for_generate0)[,ii]~trsftable[x,])$residuals)
      }
      otu_for_generate <- otu_for_generate[-1,]
      
      
      ## Entire table making 
      
      means1 <- c(1, mean(unlist(trsftable[x,])), apply(otu_for_generate,1,mean), 0) # x1, otu, hetero, epsilon
      sds1 <- c(1, sd(unlist(trsftable[x,])), apply(otu_for_generate ,1,sd), 1) # x1, otu, hetero, epsilon
      
      x11 = rnorm(mean=means1[1], sd=sds1[1], dim(trsftable)[2])
      epsilon1 <- rnorm(mean=means1[length(means1)], sd=sds1[length(means1)], dim(trsftable)[2])
      vresid1 <- do.call(sum,lapply(data.frame(combn(c(1:length(k)),2)), function(y) cov(otu_for_generate[y[1],], otu_for_generate[y[2],])))
      TT1 = sum(apply(otu_for_generate,1,var))+2*vresid1 # otu for gener = trsftable[1:k,]
      v231 = cov(unlist(trsftable[x,]), apply(otu_for_generate,2,sum))
      
      if (omics=="Metagenomics"){
        b3=sqrt(1.01/(9*TT1))  
      } else if (omics=="Proteomics"){
        b3=sqrt(1.01/(99*TT1))
      } else if (omics=="Metabolomics") {
        b3=sqrt(3.03/(7*TT1))
      }
      
      b2=Re(polyroot(c(-1.01*B+B*TT1*(b3^2), -2*B*b3*v231, (1-B)*(sds1[2]^2))))[1]       
    }
    return(b2)
    
    
  })
  
  listtable1 <- list(listtable, listtable1)
  
  return(listtable1)
  
  
}
EffectSizeFun_wo <- function(omics, trsftable, datatable, indicator, rate, b, x, cl, scenario, k=NULL){
  # browser()
  
  otu0 <- data.frame(t(datatable)[,x])
  # taxa up/down labeling
  taxa_class <- ifelse(otu0 >= median(as.numeric(unlist(otu0))), "up", "down")
  namefull <- paste0(taxa_class,"_",rownames(otu0))
  formatch <- data.frame(taxa_class=as.factor(c(taxa_class)), indicator)
  matchform <- as.formula(paste("taxa_class",paste(colnames(indicator),collapse = "+"),sep = "~"))
  
  matmat <- matchit(matchform, formatch, method = "nearest", replace = F, caliper = cl)
  
  samdat0 <- data.frame(ps = matmat$distance, 
                        pair_taxa=matmat$subclass, 
                        variable=taxa_class[match(names(matmat$subclass),rownames(taxa_class))], 
                        value=names(matmat$subclass))
  samdat <- samdat0[!is.na(samdat0$pair_taxa),]
  
  # Liability score 
  trsftable_matched <- trsftable[,match(samdat[,4], colnames(trsftable))]
  
  if (scenario=="S1"|scenario=="S2"){
    
    if (scenario=="S1"){
      PPS <- 1/(samdat$ps)
    } else if (scenario=="S2"){
      PPS <- indicator[match(samdat[,4], rownames(indicator)),1]+indicator[match(samdat[,4], rownames(indicator)),2]
    }
    means <- c(1, mean(unlist(trsftable_matched[x,])), mean(PPS), 0) # x1, otu, hetero, epsilon
    sds <- c(1, sd(unlist(trsftable_matched[x,])), sd(PPS), 1) # x1, otu, hetero, epsilon
    
    x1 = rnorm(mean=means[1], sd=sds[1], nrow(samdat))
    epsilon <- rnorm(mean=means[4], sd=sds[4], nrow(samdat))
    v23 <- cov(unlist(trsftable_matched[x,]), PPS)
    
  } else if (scenario=="S3"){
    
    otu_for_generate0 <- trsftable[k,match(samdat[,4], colnames(trsftable))]
    
    otu_for_generate <- matrix(ncol=ncol(trsftable_matched))
    for (ii in 1:length(k)){
      otu_for_generate <- rbind(otu_for_generate,lm(t(otu_for_generate0)[,ii]~trsftable_matched[x,])$residuals)
    }
    otu_for_generate <- otu_for_generate[-1,]
    
    means <- c(1, mean(unlist(trsftable_matched[x,])), apply(otu_for_generate,1,mean), 0) # x1, otutobetested, residuals, epsilon
    sds <- c(1, sd(unlist(trsftable_matched[x,])), apply(otu_for_generate,1,sd), 1) # x1, otutobetested, residuals, epsilon
    
    x1 = rnorm(mean=means[1], sd=sds[1], nrow(samdat))
    epsilon <- rnorm(mean=means[length(means)], sd=sds[length(sds)], nrow(samdat))
    vresid <- do.call(sum,lapply(data.frame(combn(c(1:length(k)),2)), function(y) cov(otu_for_generate[y[1],], otu_for_generate[y[2],])))
    TT = sum(apply(otu_for_generate,1,var))+2*vresid # var(sum(residuals))
    v23 = cov(unlist(trsftable_matched[x,]), apply(otu_for_generate,2,sum))
    
  }
  
  b1=0.1
  
  listtable <- lapply(1:length(b), function(blen){
    B=b[[blen]]
    
    if (scenario=="S1"|scenario=="S2"){
      
      if (omics=="Metagenomics"){
        b3=sqrt(1.01/(9*sds[3]^2))
      } else if (omics=="Proteomics"){
        b3=sqrt(1.01/(99*sds[3]^2))
      } else if (omics=="Metabolomics") {
        b3=sqrt(3.03/(7*sds[3]^2))  
      }
      b2=Re(polyroot(c(-1.01*B+B*(b3*sds[3])^2, -2*B*b3*v23, (1-B)*(sds[2]^2))))[1]
      
      } else if (scenario=="S3"){
      
      if (omics=="Metagenomics"){
        b3=sqrt(1.01/(9*TT))
      } else if (omics=="Proteomics"){
        b3=sqrt(1.01/(99*TT))
      } else if (omics=="Metabolomics") {
        b3=sqrt(3.03/(7*TT))
      }
      
      b2=Re(polyroot(c(-1.01*B+B*TT*(b3^2), -2*B*b3*v23, (1-B)*(sds[2]^2))))[1]
      
    }
    return(b2)
    
  })
  
  
  listtable1 <- lapply(1:length(b), function(blen){
    ## Entire table making 
    B=b[[blen]]
    if (scenario=="S1"|scenario=="S2"){
      
      if (scenario=="S1"){
        PPS1 <- 1/matmat$distance # log10 취하면 아예 inflation 이 없어짐 
      } else if (scenario=="S2"){
        PPS1 <- indicator[,1]+indicator[,2]
      }
      means1 <- c(1, mean(unlist(trsftable[x,])), mean(PPS1), 0) # x1, otu, hetero, epsilon
      sds1 <- c(1, sd(unlist(trsftable[x,])), sd(PPS1), 1) # x1, otu, hetero, epsilon
      
      x11 = rnorm(mean=means1[1], sd=sds1[1], length(PPS1))
      epsilon1 <- rnorm(mean=means1[4], sd=sds1[4], length(PPS1))
      v231 <- cov(unlist(trsftable[x,]), PPS1)
      
      if (omics=="Metagenomics"){
        b3=sqrt(1.01/(9*sds1[3]^2))
      } else if (omics=="Proteomics"){
        b3=sqrt(1.01/(99*sds1[3]^2))
      } else if (omics=="Metabolomics") {
        b3=sqrt(3.03/(7*sds1[3]^2))
      }
      
      b2=Re(polyroot(c(-1.01*B+B*(b3*sds1[3])^2, -2*B*b3*v231, (1-B)*(sds1[2]^2))))[1]
      
    } else if (scenario=="S3"){
      
      otu_for_generate0 <- trsftable[k,]
      
      otu_for_generate <- matrix(ncol=ncol(trsftable))
      for (ii in 1:length(k)){
        otu_for_generate <- rbind(otu_for_generate,lm(t(otu_for_generate0)[,ii]~trsftable[x,])$residuals)
      }
      otu_for_generate <- otu_for_generate[-1,]
      
      
      ## Entire table making 
      
      means1 <- c(1, mean(unlist(trsftable[x,])), apply(otu_for_generate,1,mean), 0) # x1, otu, hetero, epsilon
      sds1 <- c(1, sd(unlist(trsftable[x,])), apply(otu_for_generate ,1,sd), 1) # x1, otu, hetero, epsilon
      
      x11 = rnorm(mean=means1[1], sd=sds1[1], dim(trsftable)[2])
      epsilon1 <- rnorm(mean=means1[length(means1)], sd=sds1[length(means1)], dim(trsftable)[2])
      vresid1 <- do.call(sum,lapply(data.frame(combn(c(1:length(k)),2)), function(y) cov(otu_for_generate[y[1],], otu_for_generate[y[2],])))
      TT1 = sum(apply(otu_for_generate,1,var))+2*vresid1 # otu for gener = trsftable[1:k,]
      v231 = cov(unlist(trsftable[x,]), apply(otu_for_generate,2,sum))
      
      if (omics=="Metagenomics"){
        b3=sqrt(1.01/(9*TT1))
      } else if (omics=="Proteomics"){
        b3=sqrt(1.01/(99*TT1))
      } else if (omics=="Metabolomics") {
        b3=sqrt(3.03/(7*TT1))
      }
      
      b2=Re(polyroot(c(-1.01*B+B*TT1*(b3^2), -2*B*b3*v231, (1-B)*(sds1[2]^2))))[1]
       
    }
    return(b2)
    
  })
  
  listtable1 <- list(listtable, listtable1)
  
  return(listtable1)
  
  
}

load(file="/home2/nekim/scratch/matching/Robject/otulist.RData") # otulist (otu (1000/500/200), lbsize (1000))
load(file="/home2/nekim/scratch/matching/Robject/indiclist.RData") # indiclist
load(file="/home2/nekim/scratch/matching/Robject/psdifflist.RData") # psdifflist

k <- c(1,8,61,77)
for (i in nsamchar){
  for (j in c("S1","S2","S3")){
    nsamp=i
    lbs = otulist[[4]][1:nsamp]
    otutable <- otulist[[which(nsamchar==nsamp)]]
    indicator=indiclist[[which(nsamchar==nsamp)]]
    scenario=j
    
    #### rate, rep, nonlin 변경
    set.seed(1)
    rep=1
    dimm=dim(otutable)[1]
    datatable=otutable;rate=0.2;nonlin="inv" # comp
    trsftable=trnsf(dt=datatable, omics="Metagenomics", lb.size=lbs)
    b = c(0.001,0.002,0.005, 0.01)
    omics="Metagenomics"
    
    effectsize <- mclapply(1:rep, function(h){
      if(h%%10==1){print(h)}
      sumsig <- (mclapply(1:dimm, function(x){ #dim(otutable1)[1]
        
        indx =dimm*(h-1)+x 
        f(indx)
        
        if (nonlin == "inv"){
          tablefun <- get("EffectSizeFun")
          tablefun1 <- get("EffectSizeFun_wo")
        } 
        
        tble0 <- unlist(tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=NULL, scenario=scenario, k=k)[[1]])
        tble0_.1 <- unlist(tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.1, scenario=scenario, k=k)[[1]])
        tble0_.2 <-  unlist(tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.2, scenario=scenario, k=k)[[1]])
        tble0_.3 <- unlist(tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.3, scenario=scenario, k=k)[[1]])
        tble0_.05 <-  unlist(tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.05, scenario=scenario, k=k)[[1]])
        
        tble0_wo <- unlist(tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=NULL, scenario=scenario, k=k)[[1]])
        tble0_.1_wo <- unlist(tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.1, scenario=scenario, k=k)[[1]])
        tble0_.2_wo <- unlist(tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.2, scenario=scenario, k=k)[[1]])
        tble0_.3_wo <- unlist(tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.3, scenario=scenario, k=k)[[1]])
        tble0_.05_wo <- unlist(tablefun1(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=0.05, scenario=scenario, k=k)[[1]])
        
        tble0_whole <- unlist(tablefun(omics=omics, trsftable=trsftable, datatable=datatable, indicator=indicator, rate=rate, b=b, x=x,cl=NULL, scenario=scenario, k=k)[[2]]) 
        
        out <- data.frame(tble0_whole,tble0,tble0_.3,tble0_.2,tble0_.1,tble0_.05,
                          tble0_wo,tble0_.3_wo,tble0_.2_wo,tble0_.1_wo,tble0_.05_wo)
        
        return(out)
        
        
      },mc.cores=20))
    },mc.cores = 1) # no need to take 5000 reps. all the same through reps
    assign(paste0("ES_",nsamp,"_",scenario), effectsize[[1]])
    rm(effectsize)
  }
}

ES_*d
ESlist <- list(ES_1000_S1,ES_500_S1,ES_200_S1,ES_1000_S2,ES_500_S2,ES_200_S2,ES_1000_S3,ES_500_S3,ES_200_S3)
b2list <- lapply(1:length(ESlist), function(xx){
  x <- ESlist[[xx]]
  # all taxa
  b2est <- do.call(cbind, nominalb2(x,c(1:length(x))))
  # Top 20 taxa
  if(xx==1|xx==4|xx==7){psdiff<-psdifflist[[1]]} else if (xx==2|xx==5|xx==8){psdiff<-psdifflist[[2]]} else if (xx==3|xx==6|xx==9){psdiff<-psdifflist[[3]]}
  b2est20 <- do.call(cbind, nominalb2(x,order(psdiff, decreasing = T)[1:20]))
  
  colnames(b2est) <- colnames(b2est20) <- c("whole","wrep_nocal","wrep_0.3cal","wrep_0.2cal","wrep_0.1cal","wrep_0.05cal",
                                            "worep_nocal","worep_0.3cal","worep_0.2cal","worep_0.1cal","worep_0.05cal")
  rownames(b2est) <- rownames(b2est20) <- c("var.0001","var.0005","var.001","var.002")# 이건 b variance % 별 값
  
  b2list <- list(b2est, b2est20)
  return(b2list)
})
names(b2list) <- c("S1_1000","S1_500","S1_200","S2_1000","S2_500","S2_200","S3_1000","S3_500","S3_200")
save(b2list, file="data/EffectSizeList.RData") 


