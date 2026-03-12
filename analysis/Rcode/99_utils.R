# Source Code For Dealing with Simulation Result

options(digits = 3)        
cnames <- c("Unadjusted","Adjusted","Combat","Limma",
            "wrep_nocal_ind","wrep_nocal_ex","wrep_cal0.3_ind","wrep_cal0.3_ex","wrep_cal0.2_ind","wrep_cal0.2_ex","wrep_cal0.1_ind","wrep_cal0.1_ex","wrep_cal0.05_ind","wrep_cal0.05_ex",
            "worep_nocal_ind","worep_nocal_ex","worep_cal0.3_ind","worep_cal0.3_ex","worep_cal0.2_ind","worep_cal0.2_ex","worep_cal0.1_ind","worep_cal0.1_ex","worep_cal0.05_ind","worep_cal0.05_ex")
filtOTU <- function(x){
  cr1 <- which(apply(x, 1, function(x) sum(x==0))/ncol(x) > 0.5)
  cr2 <- which(apply(x, 1, sum)/sum(apply(x, 1, sum)) < 0.001) # 0.01 에서 바꿔봄 
  if (length(c(cr1,cr2))!=0) {return(x[-unique(c(cr1,cr2)),])} else {return(x)}
}
type1 <- function(x, lev){
  t001 <- sum(ifelse(is.na(unlist(x)),1,unlist(x))<0.01)/length(unlist(x))
  t005 <- sum(ifelse(is.na(unlist(x)),1,unlist(x))<0.05)/length(unlist(x))
  t01 <- sum(ifelse(is.na(unlist(x)),1,unlist(x))<0.1)/length(unlist(x))
  t02 <- sum(ifelse(is.na(unlist(x)),1,unlist(x))<0.2)/length(unlist(x))
  return(round(c(t001,t005,t01,t02),4))
}
type1_1 <- function(x){
  t0001 <- sum(ifelse(is.na(x),1,x)<0.001)/length(x)
  t0005 <- sum(ifelse(is.na(x),1,x)<0.005)/length(x)
  t001 <- sum(ifelse(is.na(x),1,x)<0.01)/length(x)
  t005 <- sum(ifelse(is.na(x),1,x)<0.05)/length(x)
  t01 <- sum(ifelse(is.na(x),1,x)<0.1)/length(x)
  return(round(c(t0001,t0005,t001,t005,t01),5))
}
CIfun <- function(pvals, level){
  # browser()
  testp <- prop.test(x=sum(pvals<level, na.rm = T),n=length(pvals),correct=T)
  return(testp$conf.int)
}
powerfold <- function(p, nmodel){
  # browser()
  p0001 <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[1,]))), byrow=T, ncol=nmodel)
  p0005 <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[2,]))), byrow=T, ncol=nmodel)
  p001 <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[3,]))), byrow=T, ncol=nmodel)
  p002 <- matrix(unlist(lapply(p, function(x) lapply(x, function(y) y[4,]))), byrow=T, ncol=nmodel)
  return(list(p0001,p0005,p001,p002))
}
powerfold1 <- function(x){
  t005 <- sum(ifelse(is.na(x),1,x)<0.05)/length(x)
}
powerbias <- function(pw, nmodel){
  # browser()
  pvalues_p <- lapply(pw, function(x) lapply(x, function(xx) xx[,1:nmodel])) # 1:24=sig, 25:48=beta
  pvalues1 <- powerfold(pvalues_p, nmodel)
  powers <- do.call(rbind,lapply(pvalues1, function(xx) apply(xx, 2, powerfold1)))
  conf.int <- do.call(rbind,lapply(pvalues1, function(y) apply(y, 2, FUN=CIfun, level=0.05)))
  powers1 <- rbind(powers, conf.int)
  
  betas_p <- lapply(pw, function(x) lapply(x, function(xx) xx[,((nmodel+1):(2*nmodel))])) 
  betas1 <- powerfold(betas_p, nmodel)
  betamean <- matrix(unlist(lapply(betas1, function(x) { # 그냥 평균
    d <- apply(x,2,mean, na.rm=T)
    return(d)
  })),ncol=nmodel, byrow=T)
  
  rownames(powers1)[5:12] <- paste0("CI_",c(rep(c("var.001","var.002","var.005","var.01"),each=2)))
  rownames(powers1)[1:4] <- paste0("Power_", c("var.001","var.002","var.005","var.01")) # row=marker variance % 가 0.0001~0.002 일 때 power
  rownames(betamean) <- paste0("Beta_", c("var.001","var.002","var.005","var.01"))
  colnames(powers1) <- colnames(betamean) <- cnames
  
  fin <- list(powers1, betamean)
  
  return(fin)
}
relbiasfun <- function(dt, ref, map){
    dt1 <- dt[,-ncol(dt)]
    dtt <- c()
    for (i in 1:ncol(dt1)){
      # browser()
      cols <- (dt1[,i]-ref[,map[i]])/ref[,map[i]]
      # cols <- (exp(dt1[,i])-exp(ref[,map[i]]))/exp(ref[,map[i]])
      dtt <- cbind(dtt, cols)
    }
    rownames(dtt) <- paste("RelBias",rownames(dtt),sep="_")
    colnames(dtt) <- colnames(dt1)
    dtt <- data.frame(dtt)
    dtt$name <- dt$name
    return(dtt)
  }
nominalb2 <- function(blist, diffbig){
  # browser()
  
  bblist <- blist[diffbig]
  
  whole <- apply(do.call(rbind,lapply(bblist, function(y) y[,1])), 2, mean)
  wrep_nocal <- apply(do.call(rbind,lapply(bblist, function(y) y[,2])), 2, mean)
  wrep_0.3cal <- apply(do.call(rbind,lapply(bblist, function(y) y[,3])), 2, mean)
  wrep_0.2cal <- apply(do.call(rbind,lapply(bblist, function(y) y[,4])), 2, mean)
  wrep_0.1cal <- apply(do.call(rbind,lapply(bblist, function(y) y[,5])), 2, mean)
  wrep_0.05cal <- apply(do.call(rbind,lapply(bblist, function(y) y[,6])), 2, mean)
  worep_nocal <- apply(do.call(rbind,lapply(bblist, function(y) y[,7])), 2, mean)
  worep_0.3cal <- apply(do.call(rbind,lapply(bblist, function(y) y[,8])), 2, mean)
  worep_0.2cal <- apply(do.call(rbind,lapply(bblist, function(y) y[,9])), 2, mean)
  worep_0.1cal <- apply(do.call(rbind,lapply(bblist, function(y) y[,10])), 2, mean)
  worep_0.05cal <- apply(do.call(rbind,lapply(bblist, function(y) y[,11])), 2, mean)
  return(list(whole,wrep_nocal,wrep_0.3cal,wrep_0.2cal,wrep_0.1cal,wrep_0.05cal,worep_nocal,worep_0.3cal,worep_0.2cal,worep_0.1cal,worep_0.05cal))
}

inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
filtextreme <- function(otus, k){
  
  # taxa up/down labeling
  taxa_class <- ifelse(otus >= median(as.numeric(unlist(otus))), "up", "down")
  formatch <- data.frame(taxa_class=as.factor(c(taxa_class)), indicator)
  matchform <- as.formula(paste("taxa_class",paste(colnames(indicator),collapse = "+"),sep = "~"))
  
  ps <- matchit(matchform, formatch, method = "nearest")$distance
  
  
  otusup <- mean(ps[which(otus>=median(as.numeric(unlist(otus))))])
  otusdown <- mean(ps[which(otus<median(as.numeric(unlist(otus))))])
  
  return(otusup-otusdown>0.1)
}
pickp <- function(x,otus) {
  xx <- x[otus]
  return(xx)
}

psdiffun <- function(otus, k){ # 매칭되는 otu 크기에 따른 두 군이 각 otu 별로 얼마나 다른 ps 값을 가지고 있는지 계산
  
  # taxa up/down labeling
  taxa_class <- ifelse(otus >= median(as.numeric(unlist(otus))), "up", "down")
  formatch <- data.frame(taxa_class=as.factor(c(taxa_class)), indicator)
  matchform <- as.formula(paste("taxa_class",paste(colnames(indicator),collapse = "+"),sep = "~"))
  
  ps <- matchit(matchform, formatch, method = "nearest")$distance
  
  
  otusup <- mean(ps[which(otus>=median(as.numeric(unlist(otus))))])
  otusdown <- mean(ps[which(otus<median(as.numeric(unlist(otus))))])
  # browser()
  return(otusup-otusdown)
}
# power 결과 정리 코드
nsamchar <- c(1000,500,200) # norder=1
mappingchar <- c(rep(1,4),rep(2:11,each=2))
tablesfun <- function(norder, scene, bin=TRUE, prevalence=0.2, append=TRUE, nonlin){
  nsamp=nsamchar[norder]
  if (bin){
    load(paste0("/data4/nekim/matching/Robject/Power/NSAMP",nsamp,"_REP500_RATE",prevalence,"_",nonlin,"_",scene,".RData"))
  } else {
    load(paste0("/data4/nekim/matching/Robject/Power_conti/NSAMP",nsamp,"_REP500_",nonlin,"_",scene,".RData"))
  }
  pvalues_test <- pvalues
  rm(pvalues)
  # browser()
  ### 1. Power & Bias
  # 전체
  totres <- data.frame(do.call(rbind,powerbias(pvalues_test, nmodel=24)))
  totres$name <- "Overall"
  # Top 20
  psdiff <- psdifflist[[norder]]
  pvalues_pp<-lapply(pvalues_test, pickp, otus=order(psdiff, decreasing = T)[1:20])
  t20res <- data.frame(do.call(rbind,powerbias(pvalues_pp, nmodel=24)))
  t20res$name <- "Top20"
  
  
  ### 2. Beta mean
  combi <- paste(scene, nsamp, sep = "_")
  # 전체
  resbeta <- totres[13:16,]
  refbias <- b2list[[combi]][[1]]
  relb <- relbiasfun(dt=resbeta, ref=refbias, map=mappingchar)  # https://www.biorxiv.org/content/10.1101/385740v1.full.pdf logistic 에서 liability 로 변환하는 공식 여럿이고 이 논문에서는 (6)번 공식 보면 됨 (그 외: https://link.springer.com/article/10.1007/s10519-021-10042-2)
  # Top 20
  resbeta20 <- t20res[13:16,]
  refbias20 <- b2list[[combi]][[2]]
  relb20 <- relbiasfun(dt=resbeta20, ref=refbias20, map=mappingchar)
  
  ### 3. Save to global env.
  tables <- data.frame(rbind(totres, t20res, relb, relb20))

  if (bin){
    globname <- paste("Bin", scene, nsamp, prevalence, nonlin, sep= "_")
  } else {
    globname <- paste("Con", scene, nsamp, sep = "_")
  }
  tables$name2 <- globname
  print(globname)
  assign(globname,tables)

  if (file.exists("/data4/nekim/matching/Routput_240818/Plot/Power/powerbias.csv")){
    write.table(tables, file="/data4/nekim/matching/Routput_240818/Plot/Power/powerbias.csv", append=append, col.names = F, quote=F, sep=",")
  } else {
    write.table(tables, file="/data4/nekim/matching/Routput_240818/Plot/Power/powerbias.csv", append=FALSE, col.names = T, quote=F, sep=",")
  }


  # browser()
  ### 4. plotting - Power
  # 전체
  pd <- totres[1:4,c(2,3,5,6,11,12,15,16,21,22)]
  rownames(pd) <- paste0("Level ",c("1","2","3","4"))
  colnames(pd) <- c(paste0("Model ",c(1:2,rep(2+(1:((ncol(pd)-2)/2)),each=2))))
  colnames(pd)[c(4,6,8,10)] <- paste0(colnames(pd)[c(4,6,8,10)],"-1")
  pdm <- melt(as.matrix(pd))
  colnames(pdm) <- c("X1", "X2","value")
  if (bin){
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Power/Bin_NSAMP",nsamp,"_REP500_RATE",prevalence,"_",nonlin,"_",scene,".png"), width=300,height=450)
  } else {
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Power/Con_NSAMP",nsamp,"_REP500_",nonlin,"_",scene,".png"), width=300,height=450)
  }
  p<-ggplot(pdm, aes(x=X1, y=value, group=X2)) +
    geom_line(aes(color=X2, linetype=X2))+
    geom_point(aes(color=X2))+ theme_bw() + ylim(c(0.05,0.6)) +
    xlab("Effect size")+ylab("Power")+
    labs(color = "Model", linetype="Model")
  print(p)
  dev.off()
  # Top 20
  pd1 <- t20res[1:4,c(2,3,5,6,11,12,15,16,21,22)]
  rownames(pd1) <- paste0("Level ",c("1","2","3","4"))
  colnames(pd1) <- paste0("Model ",c(1:2,rep(2+(1:((ncol(pd1)-2)/2)),each=2)))
  colnames(pd1)[c(4,6,8,10)] <- paste0(colnames(pd1)[c(4,6,8,10)],"-1")
  pdm1 <- melt(as.matrix(pd1))
  colnames(pdm1) <- c("X1", "X2","value")
  if (bin){
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Power/Bin_NSAMP",nsamp,"_REP500_RATE",prevalence,"_",nonlin,"_",scene,"_Top20.png"), width=300,height=450)
  } else {
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Power/Con_NSAMP",nsamp,"_REP500_inv_",scene,"_Top20.png"), width=300,height=450)
  }
  p<-ggplot(pdm1, aes(x=X1, y=value, group=X2)) +
    geom_line(aes(color=X2, linetype=X2))+
    geom_point(aes(color=X2))+ theme_bw() + ylim(c(0.05,0.6)) +
    xlab("Effect size")+ylab("Power")+
    labs(color = "Model", linetype="Model")
  print(p)
  dev.off()


  ### 5. plotting - Relative Bias
  # 전체
  rd <- relb[1:4,c(2,3,5,6,11,12,15,16,21,22)]
  rownames(rd) <- paste0("Level ",c("1","2","3","4"))
  colnames(rd) <- paste0("Model ",c(1:2,rep(2+(1:((ncol(rd)-2)/2)),each=2)))
  colnames(rd)[c(4,6,8,10)] <- paste0(colnames(rd)[c(4,6,8,10)],"-1")
  rdm <- melt(as.matrix(rd))
  colnames(rdm) <- c("X1", "X2","value")
  if (bin) {
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Bias/Bin_NSAMP",nsamp,"_REP500_RATE",prevalence,"_",nonlin,"_",scene,".png"), width=300,height=450)
  } else {
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Bias/Con_NSAMP",nsamp,"_REP500_inv_",scene,".png"), width=300,height=450)
  }
  if (nsamp==1000){
    pp<-ggplot(rdm, aes(x=X1, y=value, group=X2)) +
      geom_line(aes(color=X2, linetype=X2))+
      geom_point(aes(color=X2))+ theme_bw() +ylim(c(-1,1.5)) +
      xlab("Nominal Effect size")+ylab("Relative Bias")+
      labs(color = "Model", linetype="Model")
  } else {
    pp<-ggplot(rdm, aes(x=X1, y=value, group=X2)) +
      geom_line(aes(color=X2, linetype=X2))+
      geom_point(aes(color=X2))+ theme_bw() +ylim(c(-2,2)) +
      xlab("Nominal Effect size")+ylab("Relative Bias")+
      labs(color = "Model", linetype="Model")
  }
  print(pp)
  dev.off()
  
  # Top20
  rd1 <- relb20[1:4,c(2,3,5,6,11,12,15,16,21,22)]
  rownames(rd1) <- paste0("Level ",c("1","2","3","4"))
  colnames(rd1) <- paste0("Model ",c(1:2,rep(2+(1:((ncol(rd1)-2)/2)),each=2)))
  colnames(rd1)[c(4,6,8,10)] <- paste0(colnames(rd1)[c(4,6,8,10)],"-1")
  rdm1 <- melt(as.matrix(rd1))
  colnames(rdm1) <- c("X1", "X2","value")
  if (bin) {
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Bias/Bin_NSAMP",nsamp,"_REP500_RATE",prevalence,"_",nonlin,"_",scene,"_Top20.png"), width=300,height=450)
  } else {
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Bias/Con_NSAMP",nsamp,"_REP500_inv_",scene,"_Top20.png"), width=300,height=450)
  }
  if (nsamp==1000){
    pp<-ggplot(rdm1, aes(x=X1, y=value, group=X2)) +
      geom_line(aes(color=X2, linetype=X2))+
      geom_point(aes(color=X2))+ theme_bw() +ylim(c(-1,1.5)) +
      xlab("Nominal Effect size")+ylab("Relative Bias")+
      labs(color = "Model", linetype="Model")
  } else {
    pp<-ggplot(rdm1, aes(x=X1, y=value, group=X2)) +
      geom_line(aes(color=X2, linetype=X2))+
      geom_point(aes(color=X2))+ theme_bw() +ylim(c(-2,2)) +
      xlab("Nominal Effect size")+ylab("Relative Bias")+
      labs(color = "Model", linetype="Model")
  }
  print(pp)
  dev.off()
}
tablesfun_T1 <- function(norder, scene, bin=TRUE, prevalence=0.2, nonlin="inv", append=TRUE){
  nsamp=nsamchar[norder]
  if (bin){
    load(paste0("/data4/nekim/matching/Robject/T1/NSAMP",nsamp,"_REP5000_RATE",prevalence,"_",nonlin,"_",scene,".RData"))
  } else {
    load(paste0("/data4/nekim/matching/Robject/T1_conti/NSAMP",nsamp,"_REP5000_",nonlin,"_",scene,".RData"))
  }
  pvalues_test <- pvalues_comb
  rm(pvalues_comb)
  # browser()
  ### 1. Type 1 error & Bias
  # 전체
  t1err <- apply(pvalues_test[,1:23], 2, type1_1)
  t1err$name <- "Overall"
  # Top 20
  psdiff <- psdifflist[[norder]]
  tot = nrow(pvalues_test); dimm=tot/5000;
  pvalues_test_pp <- lapply(1:dimm, function(i){
    
    if (i<dimm){
      otui <- totminus1[which(c(1:tot)%%dimm==i),] # ith otu 
    } else {
      otui <- totminus1[which(c(1:tot)%%dimm==0),]
    }
    
    return(otui)
  })
  pvalues_test_otu <- lapply(pvalues_test_pp, function(X) apply(X, 2, type1_1))
  tour_otu <- do.call(rbind, pvalues_test_otu[order(psdiff, decreasing = T)[1:20]]) # top 20 만 합친 pvalues_test 가 되는 것. 
  
  t1err_otu <- apply(tour_otu[,1:23], 2, type1_1)
  t1err_otu$name <- "Top20"

  
  ### 2. Bias
  combi <- paste(scene, nsamp, sep = "_")
  # 전체
  resbeta <- totres[13:16,]
  refbias <- b2list[[combi]][[1]]
  relb <- relbiasfun(dt=resbeta, ref=refbias, map=mappingchar)  # https://www.biorxiv.org/content/10.1101/385740v1.full.pdf logistic 에서 liability 로 변환하는 공식 여럿이고 이 논문에서는 (6)번 공식 보면 됨 (그 외: https://link.springer.com/article/10.1007/s10519-021-10042-2)
  # Top 20
  resbeta20 <- t20res[13:16,]
  refbias20 <- b2list[[combi]][[2]]
  relb20 <- relbiasfun(dt=resbeta20, ref=refbias20, map=mappingchar)
  
  ### 3. Save to global env.
  tables <- data.frame(rbind(totres, t20res, relb, relb20))
  
  if (bin){
    globname <- paste("Bin", scene, nsamp, prevalence, nonlin, sep= "_")
  } else {
    globname <- paste("Con", scene, nsamp, sep = "_")
  }
  tables$name2 <- globname
  print(globname)
  assign(globname,tables)
  
  # if (file.exists("/data4/nekim/matching/Robject/Plot/Power/powerbias.csv")){
  #   write.table(tables, file="/data4/nekim/matching/Robject/Plot/Power/powerbias.csv", append=append, col.names = F, quote=F, sep=",")
  # } else {
  #   write.table(tables, file="/data4/nekim/matching/Robject/Plot/Power/powerbias.csv", append=FALSE, col.names = T, quote=F, sep=",")
  # }
  
  
  # browser()
  ### 4. plotting - Power
  # 전체
  pd <- totres[1:4,c(2,3,5,6,11,12,15,16,21,22)]
  rownames(pd) <- paste0("Level ",c("1","2","3","4"))
  colnames(pd) <- c(paste0("Model ",c(1:2,rep(2+(1:((ncol(pd)-2)/2)),each=2))))
  colnames(pd)[c(4,6,8,10)] <- paste0(colnames(pd)[c(4,6,8,10)],"-1")
  pdm <- melt(as.matrix(pd))
  colnames(pdm) <- c("X1", "X2","value")
  if (bin){
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Power/Bin_NSAMP",nsamp,"_REP500_RATE",prevalence,"_inv_",scene,".png"), width=300,height=450, res=300)
  } else {
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Power/Con_NSAMP",nsamp,"_REP500_inv_",scene,".png"), width=300,height=450, res=300)
  }
  p<-ggplot(pdm, aes(x=X1, y=value, group=X2)) +
    geom_line(aes(color=X2, linetype=X2))+
    geom_point(aes(color=X2))+ theme_bw() + ylim(c(0.05,0.6)) +
    xlab("Effect size")+ylab("Power")+
    labs(color = "Model", linetype="Model")
  print(p)
  dev.off()
  # Top 20
  pd1 <- t20res[1:4,c(2,3,5,6,11,12,15,16,21,22)]
  rownames(pd1) <- paste0("Level ",c("1","2","3","4"))
  colnames(pd1) <- paste0("Model ",c(1:2,rep(2+(1:((ncol(pd1)-2)/2)),each=2)))
  colnames(pd1)[c(4,6,8,10)] <- paste0(colnames(pd1)[c(4,6,8,10)],"-1")
  pdm1 <- melt(as.matrix(pd1))
  colnames(pdm1) <- c("X1", "X2","value")
  if (bin){
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Power/Bin_NSAMP",nsamp,"_REP500_RATE",prevalence,"_inv_",scene,"_Top20.png"), width=300,height=450, res=300)
  } else {
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Power/Con_NSAMP",nsamp,"_REP500_inv_",scene,"_Top20.png"), width=300,height=450, res=300)
  }
  p<-ggplot(pdm1, aes(x=X1, y=value, group=X2)) +
    geom_line(aes(color=X2, linetype=X2))+
    geom_point(aes(color=X2))+ theme_bw() + ylim(c(0.05,0.6)) +
    xlab("Effect size")+ylab("Power")+
    labs(color = "Model", linetype="Model")
  print(p)
  dev.off()
  
  
  ### 5. plotting - Relative Bias
  # 전체
  rd <- relb[1:4,c(2,3,5,6,11,12,15,16,21,22)]
  rownames(rd) <- paste0("Level ",c("1","2","3","4"))
  colnames(rd) <- paste0("Model ",c(1:2,rep(2+(1:((ncol(rd)-2)/2)),each=2)))
  colnames(rd)[c(4,6,8,10)] <- paste0(colnames(rd)[c(4,6,8,10)],"-1")
  rdm <- melt(as.matrix(rd))
  colnames(rdm) <- c("X1", "X2","value")
  if (bin) {
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Bias/Bin_NSAMP",nsamp,"_REP500_RATE",prevalence,"_inv_",scene,".png"), width=300,height=450, res=300)
  } else {
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Bias/Con_NSAMP",nsamp,"_REP500_inv_",scene,".png"), width=300,height=450, res=300)
  }
  pp<-ggplot(rdm, aes(x=X1, y=value, group=X2)) +
    geom_line(aes(color=X2, linetype=X2))+
    geom_point(aes(color=X2))+ theme_bw() +ylim(c(-1.5,1.5)) +
    xlab("Nominal Effect size")+ylab("Relative Bias")+
    labs(color = "Model", linetype="Model")
  print(pp)
  dev.off()
  
  # Top20
  rd1 <- relb20[1:4,c(2,3,5,6,11,12,15,16,21,22)]
  rownames(rd1) <- paste0("Level ",c("1","2","3","4"))
  colnames(rd1) <- paste0("Model ",c(1:2,rep(2+(1:((ncol(rd1)-2)/2)),each=2)))
  colnames(rd1)[c(4,6,8,10)] <- paste0(colnames(rd1)[c(4,6,8,10)],"-1")
  rdm1 <- melt(as.matrix(rd1))
  colnames(rdm1) <- c("X1", "X2","value")
  if (bin) {
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Bias/Bin_NSAMP",nsamp,"_REP500_RATE",prevalence,"_inv_",scene,"_Top20.png"), width=300,height=450, res=300)
  } else {
    png(file=paste0("/data4/nekim/matching/Routput_240818/Plot/Bias/Con_NSAMP",nsamp,"_REP500_inv_",scene,"_Top20.png"), width=300,height=450, res=300)
  }
  pp<-ggplot(rdm1, aes(x=X1, y=value, group=X2)) +
    geom_line(aes(color=X2, linetype=X2))+
    geom_point(aes(color=X2))+ theme_bw() +ylim(c(-1.5,1.5)) +
    xlab("Nominal Effect size")+ylab("Relative Bias")+
    labs(color = "Model", linetype="Model")
  print(pp)
  dev.off()
}
