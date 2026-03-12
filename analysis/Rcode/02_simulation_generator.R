source("00_library_setup.R")
source("99_utils.R")

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
f <- function(idx){ 
  set.seed(seeds[idx]) # set to proper seed
}
indc <- function(dt, omics) {
  if(omics=="Metagenomics"){
    shannon <- vegan::diversity(t(dt), index = "shannon")
    simpson <- vegan::diversity(t(dt), index = "simpson")
    d <- data.frame(shannon, simpson)
  } else if (omics%in%c("Proteomics","Metabolomics")){
    pp <- prcomp(t(dt), scale=T)
    d <- data.frame(pp$x[,1:2])
  }
  return(d)
}
trnsf <- function(dt, omics,lb.size){
  if(omics=="Metagenomics"){
    # d <- cpm(dt, log = TRUE, lib.size = apply(dt,2,sum))
    d <- cpm(dt, log = TRUE, lib.size = lb.size) # 전체의 library size
  } else {
    d <- dt
  }
  return(d)
}
getqic1 <- function(model.R, model.Indep){
  browser()
  mu.R <- model.R$fitted.values
  y <- model.R$y
  # Quasi Likelihood for Binomial
  quasi.R <- sum(y*log(mu.R/(1-mu.R)) + log(1-mu.R))
  
  # Trace Term (penalty for model complexity)
  AIinverse <- ginv(model.Indep$naive.variance) #
  
  Vr <- model.R$robust.variance
  
  trace.R <- sum(diag(AIinverse %*% Vr))
  # QIC
  QIC <- (-2)*quasi.R + 2*trace.R 
  return(QIC)
} 
getqic2 <- function(model.R, model.Indep){
  mu.R <- model.R$eta
  y <- model.R$y
  # Quasi Likelihood for Binomial
  quasi.R <- sum(dbinom(y, size = 1, prob = plogis(mu.R), log = TRUE)) # sum(y*log(mu.R/(1-mu.R)) + log(1-mu.R))
  
  # Trace Term (penalty for model complexity)
  AIinverse <- ginv(as.matrix(model.Indep$naiv.var)) #
  
  Vr <- as.matrix(model.R$var)
  
  trace.R <- sum(diag(AIinverse %*% Vr)) 
  # QIC
  QIC <- (-2)*quasi.R + 2*trace.R 
  return(QIC)
}  


tablemaking_t1 <- function(omics, trsftable, datatable, indicator, rate, b, x, cl, scenario, k=NULL, nonlin="inv"){ # datatable 은 행이 marker, 열이 sample 인 테이블이어야 함
  # browser()
  
  otu0 <- data.frame(t(datatable)[,x])
  # taxa up/down labeling
  taxa_class <- ifelse(otu0 >= median(as.numeric(unlist(otu0))), "up", "down")
  namefull <- paste0(taxa_class,"_",rownames(otu0))
  formatch <- data.frame(taxa_class=as.factor(c(taxa_class)), indicator)
  matchform <- as.formula(paste("taxa_class",paste(colnames(indicator),collapse = "+"),sep = "~"))
  
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
      if (nonlin=="inv"){
        
        PPS <- 1/(samdat$ps)
      } else if (nonlin=="comp"){
        if (omics=="Metagenomics"){
          
          PPS <- log(exp(samdat$ps)+log(samdat$ps)*((samdat$ps)^-2-6)/-2) # microbiome
        } else {
          
          PPS <- log(exp(samdat$ps)+log(samdat$ps)*((samdat$ps)^-2-8)/-2)
        }
      }
      
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
  b2=0 
  
  
  
  listtable <- lapply(1:length(b), function(blen){
    
    # add effect size
    B=b[[blen]]
    
    if (scenario=="S1"|scenario=="S2"){
      if (omics=="Metagenomics"){
        b3=sqrt(1.01/(9*sds[3]^2)) # batch effect = 10% coefficient when beta1=0.1   
      } else if (omics=="Proteomics"){
        b3=sqrt(1.01/(99*sds[3]^2)) # batch effect = 1% coefficient when beta1=0.1   
      } else if (omics=="Metabolomics") {
        b3=sqrt(3.03/(7*sds[3]^2)) # batch effect = 30% coefficient when beta1=0.1     
      }
      L=b1*x1 + b2*unlist(trsftable_matched[x,]) + b3*PPS + epsilon
      
      thres <- qnorm(1-rate,
                     mean=b1*means[1]+b2*means[2]+b3*means[3],
                     sd=sqrt((b1*sds[1])^2+(b2*sds[2])^2+(b3*sds[3])^2+2*b2*b3*v23+1)) 
    } else if (scenario=="S3"){
      
      if (omics=="Metagenomics"){
        b3=sqrt(1.01/(9*TT)) # batch effect = 10% coefficient when beta1=0.1      
      } else if (omics=="Proteomics"){
        b3=sqrt(1.01/(99*TT)) # batch effect = 1% coefficient when beta1=0.1   
      } else if (omics=="Metabolomics") {
        b3=sqrt(3.03/(7*TT)) # batch effect = 30% coefficient when beta1=0.1     
      }
      L=b1*x1 + b2*unlist(trsftable_matched[x,]) + b3*apply(otu_for_generate,2,sum) + epsilon
      
      thres <- qnorm(1-rate,
                     mean=b1*means[1]+b2*means[2]+b3*sum(apply(otu_for_generate,1,mean)),
                     sd=sqrt((b1*sds[1])^2+(b2*sds[2])^2+(b3^2)*TT+2*b2*b3*v23+1)) 
      
    }

    case = ifelse(L>=thres,"1","0")
    
    samdat1 <- data.frame(ID=samdat[,4],                                  
                          case=as.numeric(as.character(case[match(samdat[,4],names(case))])), 
                          pair_taxa=samdat$pair_taxa,
                          indOTU = ifelse(samdat[,3]=="Var1","large","small"),
                          x1 = x1,
                          PS = samdat$ps) 
    samdat1 <- samdat1[!duplicated(samdat1$ID),]
    rownames(samdat1) <- samdat1$ID
    # browser()
    
    # final data
    trsftable1 <- trsftable[,match(rownames(samdat1),colnames(trsftable))] 
    matchdata <- data.frame(samdat1, 
                            hetero1 = indicator[match(rownames(samdat1), rownames(indicator)),1],
                            hetero2 = indicator[match(rownames(samdat1), rownames(indicator)),2],
                            marker=unlist(trsftable1[x,]))
    matchdata1 <- matchdata[order(matchdata$pair_taxa),]
    matchdata1[,"pair_taxa"] <- as.factor(matchdata1[,"pair_taxa"])
    
    tbpair <- table(matchdata1$pair_taxa)
    matchdata1$m <- rep(1, nrow(matchdata1))
    matchdata1$m[cumsum(tbpair)] <- tbpair-1
    matchdata1$mm <- matchdata1$m/sum(matchdata1$m)
    matchdata1$mmm <- matmat$weights[match(rownames(matchdata1),names(matmat$weights))] 

    return(matchdata1)
  })
  
  rm(means);rm(sds);rm(x1);rm(epsilon);rm(v23)
  
  
  
  ## Entire table making 
  if (scenario=="S1"|scenario=="S2"){
    if (scenario=="S1"){
      if (nonlin=="inv"){
        
        PPS1 <- 1/matmat$distance
      } else if (nonlin=="comp"){
        
        if (omics=="Metagenomics"){
          PPS1 <- log(exp(matmat$distance)+log(matmat$distance)*((matmat$distance)^-2-6)/-2)  # microbiome
        } else {
          PPS1 <- log(exp(matmat$distance)+log(matmat$distance)*((matmat$distance)^-2-8)/-2) 
        }
      }
      
    } else if (scenario=="S2"){
      PPS1 <- indicator[,1]+indicator[,2]
    }
    means1 <- c(1, mean(unlist(trsftable[x,])), mean(PPS1), 0) # x1, otu, hetero, epsilon
    sds1 <- c(1, sd(unlist(trsftable[x,])), sd(PPS1), 1) # x1, otu, hetero, epsilon
    
    x11 = rnorm(mean=means1[1], sd=sds1[1], length(PPS1))
    epsilon1 <- rnorm(mean=means1[4], sd=sds1[4], length(PPS1))
    v231 <- cov(unlist(trsftable[x,]), PPS1) 
    
    if (omics=="Metagenomics"){
      b3=sqrt(1.01/(9*sds1[3]^2)) # batch 의 효과가 10% 가 되는 coefficient when beta1=0.1   
    } else if (omics=="Proteomics"){
      b3=sqrt(1.01/(99*sds1[3]^2)) # batch 의 효과가 20% 가 되는 coefficient when beta1=0.1   
    } else if (omics=="Metabolomics") {
      b3=sqrt(3.03/(7*sds1[3]^2)) # batch 의 효과가 30% 가 되는 coefficient when beta1=0.1   
    }
    
    L1=b1*x11 + b2*unlist(trsftable[x,]) + b3*PPS1 + epsilon1
    
    thres1 <- qnorm(1-rate,
                    mean=b1*means1[1]+b2*means1[2]+b3*means1[3],
                    sd=sqrt((b1*sds1[1])^2+(b2*sds1[2])^2+(b3*sds1[3])^2+2*b2*b3*v231+1)) 
    
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
    
    L1=b1*x11 + b2*unlist(trsftable[x,]) + b3*apply(otu_for_generate,2,sum) + epsilon1
    
    thres1 <- qnorm(1-rate,
                    mean=b1*means1[1]+b2*means1[2]+b3*sum(apply(otu_for_generate,1,mean)),
                    sd=sqrt((b1*sds1[1])^2+(b2*sds1[2])^2+(b3^2)*TT1+2*b2*b3*v231+1)) # y=1 비율 20%, sum of correlated Normal dist.
    
  }
  
  case1 = ifelse(L1>=thres1,"1","0")
  
  entiredata <- data.frame(ID=names(matmat$distance), 
                           case=as.numeric(as.character(case1)),
                           x1=x11,
                           PS=matmat$distance,
                           hetero1 = indicator[match(names(matmat$distance), rownames(indicator)),1],
                           hetero2 = indicator[match(names(matmat$distance), rownames(indicator)),2],
                           marker=unlist(trsftable[x,]))
  
  meanrep <- mean(table(matmat$match.matrix)[which(table(matmat$match.matrix)>1)])
  
  listtable1 <- append(listtable, list(entiredata, meanrep))
  
  return(listtable1)
  
  
}

tablemaking_t1_wo <- function(omics, trsftable, datatable, indicator, rate, b, x, cl, scenario, k=NULL, nonlin="inv"){
  # browser()
  
  otu0 <- data.frame(t(datatable)[,x])
  # taxa up/down labeling
  taxa_class <- ifelse(otu0 >= median(as.numeric(unlist(otu0))), "up", "down")
  # taxa_class <- ifelse(otu0 < summary(as.numeric(unlist(otu0)))[5], "up", "down")
  namefull <- paste0(taxa_class,"_",rownames(otu0))
  formatch <- data.frame(taxa_class=as.factor(c(taxa_class)), indicator)
  matchform <- as.formula(paste("taxa_class",paste(colnames(indicator),collapse = "+"),sep = "~"))
  
  # matmat <- matchit(matchform, formatch, method = "nearest", replace = T) 
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
      if (nonlin=="inv"){
        
        PPS <- 1/(samdat$ps)
      } else if (nonlin=="comp"){
        if (omics=="Metagenomics"){
          
          PPS <- log(exp(samdat$ps)+log(samdat$ps)*((samdat$ps)^-2-6)/-2) # microbiome
        } else {
          
          PPS <- log(exp(samdat$ps)+log(samdat$ps)*((samdat$ps)^-2-8)/-2)
        }
      }
      
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
  b2=0 # 마커와는 관계 없이 Y 생성해야 TYPE1 error
  
  
  
  listtable <- lapply(1:length(b), function(blen){
    
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
      
      L=b1*x1 + b2*unlist(trsftable_matched[x,]) + b3*PPS + epsilon
      
      thres <- qnorm(1-rate,
                     mean=b1*means[1]+b2*means[2]+b3*means[3],
                     sd=sqrt((b1*sds[1])^2+(b2*sds[2])^2+(b3*sds[3])^2+2*b2*b3*v23+1))
    } else if (scenario=="S3"){
      
      if (omics=="Metagenomics"){
        b3=sqrt(1.01/(9*TT)) # batch 의 효과가 10% 가 되는 coefficient when beta1=0.1   
      } else if (omics=="Proteomics"){
        b3=sqrt(1.01/(99*TT)) # batch 의 효과가 1% 가 되는 coefficient when beta1=0.1   
      } else if (omics=="Metabolomics") {
        b3=sqrt(3.03/(7*TT)) # batch 의 효과가 30% 가 되는 coefficient when beta1=0.1     
      }
      
      L=b1*x1 + b2*unlist(trsftable_matched[x,]) + b3*apply(otu_for_generate,2,sum) + epsilon
      
      thres <- qnorm(1-rate,
                     mean=b1*means[1]+b2*means[2]+b3*sum(apply(otu_for_generate,1,mean)),
                     sd=sqrt((b1*sds[1])^2+(b2*sds[2])^2+(b3^2)*TT+2*b2*b3*v23+1))
      
    }
    
    case = ifelse(L>=thres,"1","0")
    
    samdat1 <- data.frame(ID=samdat[,4],                                  
                          case=as.numeric(as.character(case[match(samdat[,4],names(case))])), 
                          pair_taxa=samdat$pair_taxa,
                          indOTU = ifelse(samdat[,3]=="up","large","small"),
                          x1 = x1,
                          PS = samdat$ps) 
    samdat1 <- samdat1[!duplicated(samdat1$ID),]
    rownames(samdat1) <- samdat1$ID
    
    
    # final data
    trsftable1 <- trsftable[,match(rownames(samdat1),colnames(trsftable))]
    matchdata <- data.frame(samdat1, 
                            hetero1 = indicator[match(rownames(samdat1), rownames(indicator)),1],
                            hetero2 = indicator[match(rownames(samdat1), rownames(indicator)),2],
                            marker=unlist(trsftable1[x,]))
    matchdata1 <- matchdata[order(matchdata$pair_taxa),]
    matchdata1[,"pair_taxa"] <- as.factor(matchdata1[,"pair_taxa"])
    
    tbpair <- table(matchdata1$pair_taxa)
    matchdata1$m <- rep(1, nrow(matchdata1))
    matchdata1$m[cumsum(tbpair)] <- tbpair-1
    matchdata1$mm <- matchdata1$m/sum(matchdata1$m)
    matchdata1$mmm <- matmat$weights[match(rownames(matchdata1),names(matmat$weights))]
    
    return(matchdata1)
    
    
  })
  
  
  
  
  
  ## Entire table making 
  
  if (scenario=="S1"|scenario=="S2"){
    
    if (scenario=="S1"){
      if (nonlin=="inv"){
        
        PPS1 <- 1/matmat$distance
      } else if (nonlin=="comp"){
        
        if (omics=="Metagenomics"){
          PPS1 <- log(exp(matmat$distance)+log(matmat$distance)*((matmat$distance)^-2-6)/-2)  # microbiome
        } else {
          PPS1 <- log(exp(matmat$distance)+log(matmat$distance)*((matmat$distance)^-2-8)/-2) 
        }
      }
      
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
    
    L1=b1*x11 + b2*unlist(trsftable[x,]) + b3*PPS1 + epsilon1
    
    thres1 <- qnorm(1-rate,
                    mean=b1*means1[1]+b2*means1[2]+b3*means1[3],
                    sd=sqrt((b1*sds1[1])^2+(b2*sds1[2])^2+(b3*sds1[3])^2+2*b2*b3*v231+1))
    
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
      b3=sqrt(1.01/(9*TT1)) # batch 의 효과가 10% 가 되는 coefficient when beta1=0.1   
    } else if (omics=="Proteomics"){
      b3=sqrt(1.01/(99*TT1)) # batch 의 효과가 20% 가 되는 coefficient when beta1=0.1   
    } else if (omics=="Metabolomics") {
      b3=sqrt(3.03/(7*TT1)) # batch 의 효과가 30% 가 되는 coefficient when beta1=0.1   
    }
    
    L1=b1*x11 + b2*unlist(trsftable[x,]) + b3*apply(otu_for_generate,2,sum) + epsilon1
    
    thres1 <- qnorm(1-rate,
                    mean=b1*means1[1]+b2*means1[2]+b3*sum(apply(otu_for_generate,1,mean)),
                    sd=sqrt((b1*sds1[1])^2+(b2*sds1[2])^2+(b3^2)*TT1+2*b2*b3*v231+1))
    
  }
  
  case1 = ifelse(L1>=thres1,"1","0")
  # browser()
  entiredata <- data.frame(ID=names(matmat$distance), 
                           case=as.numeric(as.character(case1)),
                           x1=x11,
                           PS=matmat$distance,
                           hetero1 = indicator[match(names(matmat$distance), rownames(indicator)),1],
                           hetero2 = indicator[match(names(matmat$distance), rownames(indicator)),2],
                           marker=unlist(trsftable[x,]))
  
  listtable1 <- append(listtable, list(entiredata))
  
  return(listtable1)
  
  
}




tablemaking_power <- function(omics, trsftable, datatable, indicator, rate, b, x,cl, scenario, k=NULL, nonlin="inv"){
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
      if (nonlin=="inv"){
        
        PPS <- 1/(samdat$ps)
      } else if (nonlin=="comp"){
        if (omics=="Metagenomics"){
          
          PPS <- log(exp(samdat$ps)+log(samdat$ps)*((samdat$ps)^-2-6)/-2) # microbiome
        } else {
          
          PPS <- log(exp(samdat$ps)+log(samdat$ps)*((samdat$ps)^-2-8)/-2)
        }
      }
      
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
    TT = sum(apply(otu_for_generate,1,var))+2*vresid # var(sum(residuals)) = var(apply(dd,2,sum)) 로 편하게 나타내도 됨..
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
      
      L=b1*x1 + b2*unlist(trsftable_matched[x,]) + b3*PPS + epsilon
      
      thres <- qnorm(1-rate,
                     mean=b1*means[1]+b2*means[2]+b3*means[3],
                     sd=sqrt((b1*sds[1])^2+(b2*sds[2])^2+(b3*sds[3])^2+2*b2*b3*v23+1))
    } else if (scenario=="S3"){
      
      if (omics=="Metagenomics"){
        b3=sqrt(1.01/(9*TT)) # batch 의 효과가 10% 가 되는 coefficient when beta1=0.1   
      } else if (omics=="Proteomics"){
        b3=sqrt(1.01/(99*TT)) # batch 의 효과가 1% 가 되는 coefficient when beta1=0.1   
      } else if (omics=="Metabolomics") {
        b3=sqrt(3.03/(7*TT)) # batch 의 효과가 30% 가 되는 coefficient when beta1=0.1     
      }
      
      b2=Re(polyroot(c(-1.01*B+B*TT*(b3^2), -2*B*b3*v23, (1-B)*(sds[2]^2))))[1]
      L=b1*x1 + b2*unlist(trsftable_matched[x,]) + b3*apply(otu_for_generate,2,sum) + epsilon
      
      thres <- qnorm(1-rate,
                     mean=b1*means[1]+b2*means[2]+b3*sum(apply(otu_for_generate,1,mean)),
                     sd=sqrt((b1*sds[1])^2+(b2*sds[2])^2+(b3^2)*TT+2*b2*b3*v23+1))
      
    }
    
    case = ifelse(L>=thres,"1","0")
    
    samdat1 <- data.frame(ID=samdat[,4],                                  
                          case=as.numeric(as.character(case[match(samdat[,4],names(case))])),
                          pair_taxa=samdat$pair_taxa,
                          indOTU = ifelse(samdat[,3]=="Var1","large","small"),
                          x1 = x1,
                          PS = samdat$ps) 
    samdat1 <- samdat1[!duplicated(samdat1$ID),]
    rownames(samdat1) <- samdat1$ID
    
    
    # final data
    trsftable1 <- trsftable[,match(rownames(samdat1),colnames(trsftable))]
    matchdata <- data.frame(samdat1, 
                            hetero1 = indicator[match(rownames(samdat1), rownames(indicator)),1],
                            hetero2 = indicator[match(rownames(samdat1), rownames(indicator)),2],
                            marker=unlist(trsftable1[x,]))
    matchdata1 <- matchdata[order(matchdata$pair_taxa),]
    matchdata1[,"pair_taxa"] <- as.factor(matchdata1[,"pair_taxa"])
    
    tbpair <- table(matchdata1$pair_taxa)
    matchdata1$m <- rep(1, nrow(matchdata1))
    matchdata1$m[cumsum(tbpair)] <- tbpair-1
    matchdata1$mm <- matchdata1$m/sum(matchdata1$m)
    matchdata1$mmm <- matmat$weights[match(rownames(matchdata1),names(matmat$weights))]
    
    return(matchdata1)
    
    
  })
  
  rm(means);rm(sds);rm(x1);rm(epsilon);rm(v23)
  
  listtable1 <- lapply(1:length(b), function(blen){
    ## Entire table making 
    B=b[[blen]]
    
    if (scenario=="S1"|scenario=="S2"){
      if (scenario=="S1"){
        if (nonlin=="inv"){
          
          PPS1 <- 1/matmat$distance
        } else if (nonlin=="comp"){
          
          if (omics=="Metagenomics"){
            PPS1 <- log(exp(matmat$distance)+log(matmat$distance)*((matmat$distance)^-2-6)/-2)  # microbiome
          } else {
            PPS1 <- log(exp(matmat$distance)+log(matmat$distance)*((matmat$distance)^-2-8)/-2) 
          }
        }
        
      } else if (scenario=="S2"){
        PPS1 <- indicator[,1]+indicator[,2]
      }
      means1 <- c(1, mean(unlist(trsftable[x,])), mean(PPS1), 0) # x1, otu, hetero, epsilon
      sds1 <- c(1, sd(unlist(trsftable[x,])), sd(PPS1), 1) # x1, otu, hetero, epsilon
      
      x11 = rnorm(mean=means1[1], sd=sds1[1], length(PPS1))
      epsilon1 <- rnorm(mean=means1[4], sd=sds1[4], length(PPS1))
      v231 <- cov(unlist(trsftable[x,]), PPS1)
      
      if (omics=="Metagenomics"){
        b3=sqrt(1.01/(9*sds1[3]^2)) # batch 의 효과가 10% 가 되는 coefficient when beta1=0.1   
      } else if (omics=="Proteomics"){
        b3=sqrt(1.01/(99*sds1[3]^2)) # batch 의 효과가 20% 가 되는 coefficient when beta1=0.1   
      } else if (omics=="Metabolomics") {
        b3=sqrt(3.03/(7*sds1[3]^2)) # batch 의 효과가 30% 가 되는 coefficient when beta1=0.1   
      }
      
      b2=Re(polyroot(c(-1.01*B+B*(b3*sds1[3])^2, -2*B*b3*v231, (1-B)*(sds1[2]^2))))[1] 
      
      L1=b1*x11 + b2*unlist(trsftable[x,]) + b3*PPS1 + epsilon1
      
      thres1 <- qnorm(1-rate,
                      mean=b1*means1[1]+b2*means1[2]+b3*means1[3],
                      sd=sqrt((b1*sds1[1])^2+(b2*sds1[2])^2+(b3*sds1[3])^2+2*b2*b3*v231+1))
      
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
      L1=b1*x11 + b2*unlist(trsftable[x,]) + b3*apply(otu_for_generate,2,sum) + epsilon1
      
      thres1 <- qnorm(1-rate,
                      mean=b1*means1[1]+b2*means1[2]+b3*sum(apply(otu_for_generate,1,mean)),
                      sd=sqrt((b1*sds1[1])^2+(b2*sds1[2])^2+(b3^2)*TT1+2*b2*b3*v231+1)) 
      
    }
    
    case1 = ifelse(L1>=thres1,"1","0")
    # browser()
    entiredata <- data.frame(ID=names(matmat$distance), 
                             case=as.numeric(as.character(case1)),
                             x1=x11,
                             PS=matmat$distance,
                             hetero1 = indicator[match(names(matmat$distance), rownames(indicator)),1],
                             hetero2 = indicator[match(names(matmat$distance), rownames(indicator)),2],
                             marker=unlist(trsftable[x,]))
    
    return(entiredata)
    
    
  })
  # browser()
  
  
  
  listtable1 <- list(listtable, listtable1)
  
  return(listtable1)
  
  
}

tablemaking_power_wo <- function(omics, trsftable, datatable, indicator, rate, b, x, cl, scenario, k=NULL, nonlin="inv"){
  
  otu0 <- data.frame(t(datatable)[,x])
  # taxa up/down labeling
  taxa_class <- ifelse(otu0 >= median(as.numeric(unlist(otu0))), "up", "down")
  namefull <- paste0(taxa_class,"_",rownames(otu0))
  formatch <- data.frame(taxa_class=as.factor(c(taxa_class)), indicator)
  matchform <- as.formula(paste("taxa_class",paste(colnames(indicator),collapse = "+"),sep = "~"))
  
  # matmat <- matchit(matchform, formatch, method = "nearest", replace = T) 
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
      if (nonlin=="inv"){
        
        PPS <- 1/(samdat$ps)
      } else if (nonlin=="comp"){
        if (omics=="Metagenomics"){
          
          PPS <- log(exp(samdat$ps)+log(samdat$ps)*((samdat$ps)^-2-6)/-2) # microbiome
        } else {
          
          PPS <- log(exp(samdat$ps)+log(samdat$ps)*((samdat$ps)^-2-8)/-2)
        }
      }
      
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
      L=b1*x1 + b2*unlist(trsftable_matched[x,]) + b3*PPS + epsilon
      
      thres <- qnorm(1-rate,
                     mean=b1*means[1]+b2*means[2]+b3*means[3],
                     sd=sqrt((b1*sds[1])^2+(b2*sds[2])^2+(b3*sds[3])^2+2*b2*b3*v23+1))
    } else if (scenario=="S3"){
      
      if (omics=="Metagenomics"){
        b3=sqrt(1.01/(9*TT))    
      } else if (omics=="Proteomics"){
        b3=sqrt(1.01/(99*TT)) 
      } else if (omics=="Metabolomics") {
        b3=sqrt(3.03/(7*TT))
      }
      
      b2=Re(polyroot(c(-1.01*B+B*TT*(b3^2), -2*B*b3*v23, (1-B)*(sds[2]^2))))[1] 
      L=b1*x1 + b2*unlist(trsftable_matched[x,]) + b3*apply(otu_for_generate,2,sum) + epsilon
      
      thres <- qnorm(1-rate,
                     mean=b1*means[1]+b2*means[2]+b3*sum(apply(otu_for_generate,1,mean)),
                     sd=sqrt((b1*sds[1])^2+(b2*sds[2])^2+(b3^2)*TT+2*b2*b3*v23+1)) 
      
    }
    
    case = ifelse(L>=thres,"1","0")
    
    samdat1 <- data.frame(ID=samdat[,4],                                  
                          case=as.numeric(as.character(case[match(samdat[,4],names(case))])), 
                          pair_taxa=samdat$pair_taxa,
                          indOTU = ifelse(samdat[,3]=="up","large","small"),
                          x1 = x1,
                          PS = samdat$ps) 
    samdat1 <- samdat1[!duplicated(samdat1$ID),]
    rownames(samdat1) <- samdat1$ID
    
    
    # final data
    trsftable1 <- trsftable[,match(rownames(samdat1),colnames(trsftable))] 
    matchdata <- data.frame(samdat1, 
                            hetero1 = indicator[match(rownames(samdat1), rownames(indicator)),1],
                            hetero2 = indicator[match(rownames(samdat1), rownames(indicator)),2],
                            marker=unlist(trsftable1[x,]))
    matchdata1 <- matchdata[order(matchdata$pair_taxa),]
    matchdata1[,"pair_taxa"] <- as.factor(matchdata1[,"pair_taxa"])
    
    tbpair <- table(matchdata1$pair_taxa)
    matchdata1$m <- rep(1, nrow(matchdata1))
    matchdata1$m[cumsum(tbpair)] <- tbpair-1
    matchdata1$mm <- matchdata1$m/sum(matchdata1$m)
    matchdata1$mmm <- matmat$weights[match(rownames(matchdata1),names(matmat$weights))]
    return(matchdata1)
    
    
  })
  
  
  listtable1 <- lapply(1:length(b), function(blen){
    ## Entire table making 
    B=b[[blen]]
    
    if (scenario=="S1"|scenario=="S2"){
      
      if (scenario=="S1"){
        if (nonlin=="inv"){
          
          PPS1 <- 1/matmat$distance
        } else if (nonlin=="comp"){
          
          if (omics=="Metagenomics"){
            PPS1 <- log(exp(matmat$distance)+log(matmat$distance)*((matmat$distance)^-2-6)/-2)  # microbiome
          } else {
            PPS1 <- log(exp(matmat$distance)+log(matmat$distance)*((matmat$distance)^-2-8)/-2) 
          }
        }
        
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
      
      L1=b1*x11 + b2*unlist(trsftable[x,]) + b3*PPS1 + epsilon1
      
      thres1 <- qnorm(1-rate,
                      mean=b1*means1[1]+b2*means1[2]+b3*means1[3],
                      sd=sqrt((b1*sds1[1])^2+(b2*sds1[2])^2+(b3*sds1[3])^2+2*b2*b3*v231+1))
      
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
        b3=sqrt(1.01/(9*TT1)) # batch 의 효과가 10% 가 되는 coefficient when beta1=0.1   
      } else if (omics=="Proteomics"){
        b3=sqrt(1.01/(99*TT1)) # batch 의 효과가 20% 가 되는 coefficient when beta1=0.1   
      } else if (omics=="Metabolomics") {
        b3=sqrt(3.03/(7*TT1)) # batch 의 효과가 30% 가 되는 coefficient when beta1=0.1   
      }
      
      b2=Re(polyroot(c(-1.01*B+B*TT1*(b3^2), -2*B*b3*v231, (1-B)*(sds1[2]^2))))[1] 
      L1=b1*x11 + b2*unlist(trsftable[x,]) + b3*apply(otu_for_generate,2,sum) + epsilon1
      
      thres1 <- qnorm(1-rate,
                      mean=b1*means1[1]+b2*means1[2]+b3*sum(apply(otu_for_generate,1,mean)),
                      sd=sqrt((b1*sds1[1])^2+(b2*sds1[2])^2+(b3^2)*TT1+2*b2*b3*v231+1))
      
    }
    
    case1 = ifelse(L1>=thres1,"1","0")
    # browser()
    entiredata <- data.frame(ID=names(matmat$distance), 
                             case=as.numeric(as.character(case1)),
                             x1=x11,
                             PS=matmat$distance,
                             hetero1 = indicator[match(names(matmat$distance), rownames(indicator)),1],
                             hetero2 = indicator[match(names(matmat$distance), rownames(indicator)),2],
                             marker=unlist(trsftable[x,]))
    
    return(entiredata)
    
    
  })
  
  
  
  
  listtable1 <- list(listtable, listtable1)
  
  return(listtable1)
  
  
}



