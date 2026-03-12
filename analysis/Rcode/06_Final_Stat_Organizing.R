## Organizing the final empirical type 1 error rate and statistical power

#####################################
## Scenario: 1,2,3
## Statistic : T1
## Topic : Organizing
#####################################
source("99_utils.R")
load(file="data/EffectSizeList.RData") # b2list 

CIfun1<-function(pvals, level){ 
  # browser()
  if(sum(is.na(pvals))>0) pvals[which(is.na(pvals))] <- 1
  xx <-round((sum(pvals<level)/length(pvals))*5000,0)
  testp <- prop.test(x=xx,n=5000,correct=T)
  return(testp$conf.int)
}

nsample <- c(200, 500, 1000)
scen <- c("S1", "S2","S3") 
rate <- c("0.2", "0.3", "0.5")
rep <- 5000

result_list <- list()

for (nsamp in nsample) {
  for (scenario in scen) {
    if (scenario == "S1") {
      nonlin_values <- c("inv", "comp")
    } else {
      nonlin_values <- NULL
    }
    
    for (rate_value in rate) {
      for (nonlin in if (is.null(nonlin_values)) "" else nonlin_values) {
        file_name <- paste0("NSAMP", nsamp, "_REP", rep, "_RATE", rate_value, "_",
                            if (nonlin != "") nonlin, "_", scenario, ".RData")
        print(file_name)
        load(file = paste0("/data4/nekim/matching/Robject/T1/", file_name))
        
        t1err <- apply(pvalues_comb[,1:23], 2, type1_1)
        rownames(t1err) <- c("p.001","p.005","p.01","p.05","p.1")
        
        conf.int <- do.call(rbind, lapply(c(0.001, 0.005, 0.01, 0.05, 0.1), 
                                          function(y) apply(pvalues_comb[,1:23], 2, FUN = CIfun1, level = y)))
        rownames(conf.int) <- paste0("CI_", c(rep(c("p.001", "p.005", "p.01", "p.05", "p.1"), each = 2)))
        
        t1err1 <- rbind(t1err, conf.int)
        colnames(t1err1) <- c("LogitZ", "LogitZX", "Combat",
                              paste0("Cal", rep(c("No", 0.3, 0.2, 0.1, 0.05), each = 2), c("Ind", "Exp")),
                              paste0("Cal", rep(c("No", 0.3, 0.2, 0.1, 0.05), each = 2), c("Ind", "Exp"), "Wrep"))
        
        result_list[[file_name]] <- t1err1
      }
    }
  }
}

result_list
combined_table <- do.call(rbind, lapply(names(result_list), function(table_name) {
  data <- result_list[[table_name]]
  data <- as.data.frame(data[-c(1,6,7),]) 
  table_name1 <- gsub(pattern = ".RData","",table_name)
  data1 <- data.frame(table_name1, data)
  return(data1)
}))

head(combined_table)
combined_table <- data.frame(title=rep(rownames(result_list[[1]])[-c(1,6,7)],length(result_list)),
                             combined_table)

write.xlsx(combined_table, file = "data/combined_T1err.xlsx", row.names = FALSE)




#####################################
#####################################
## Scenario: 1,2,3
## Statistic : T1_nrow, mean rep
## Topic : Organizing
#####################################

nsample <- c(200, 500, 1000)
scen <- c("S1", "S2","S3") 
rate <- c("0.2", "0.3", "0.5")
rep <- 5000

mean_nsample_list <- list()
mean_rep_list <- list()

for (nsamp in nsample) {
  for (scenario in scen) {
    if (scenario == "S1") {
      nonlin_values <- c("inv", "comp")
    } else {
      nonlin_values <- NULL
    }
    
    for (rate_value in rate) {
      for (nonlin in if (is.null(nonlin_values)) "" else nonlin_values) {
        file_name <- paste0("NSAMP", nsamp, "_REP", rep, "_RATE", rate_value, "_",
                            if (nonlin != "") nonlin, "_", scenario, ".RData")
        print(file_name)
        load(file = paste0("data/T1/", file_name))
        
        mean_nsample_values <- colMeans(pvalues_comb[, 47:56], na.rm = TRUE)
        names(mean_nsample_values) <- c("CalNo", "Cal0.3", "Cal0.2", "Cal0.1", "Cal0.05", 
                                        "CalNoRep", "Cal0.3Rep", "Cal0.2Rep", "Cal0.1Rep", "Cal0.05Rep")
        
        mean_nsample_list[[file_name]] <- mean_nsample_values
        
        mean_rep_values <- colMeans(pvalues_comb[, 57:61], na.rm = TRUE)
        names(mean_rep_values) <- c("CalNoRep", "Cal0.3Rep", "Cal0.2Rep", "Cal0.1Rep", "Cal0.05Rep")
        
        mean_rep_list[[file_name]] <- mean_rep_values
      }
    }
  }
}

# mean_nsample table bind
mean_nsample <- do.call(rbind, lapply(names(mean_nsample_list), function(table_name) {
  data <- mean_nsample_list[[table_name]]
  table_name1 <- gsub(pattern = ".RData", "", table_name)
  data1 <- as.data.frame(matrix(data, ncol = 10, byrow = T)) # 값을 행으로 변환
  colnames(data1) <- names(data)
  data2 <- data.frame(table_name1, data1)
  return(data2)
}))

# mean_rep table bind
mean_rep <- do.call(rbind, lapply(names(mean_rep_list), function(table_name) {
  data <- mean_rep_list[[table_name]]
  table_name1 <- gsub(pattern = ".RData", "", table_name)
  data1 <- as.data.frame(matrix(data, ncol = 5, byrow = T)) # 값을 행으로 변환
  colnames(data1) <- names(data)
  data2 <- data.frame(table_name1, data1)
  return(data2)
}))

# print table
head(mean_nsample)
head(mean_rep)
write.xlsx(mean_nsample, file = "data/combined_nsample.xlsx", row.names = FALSE)
write.xlsx(mean_rep, file = "data/combined_rep.xlsx", row.names = FALSE)






#####################################
#####################################
## Scenario: 1,2,3
## Statistic : Power
## Topic : Organizing
#####################################
source("99_utils.R")
load(file="data/EffectSizeList.RData") # b2list

Binary=TRUE
for (i in c(1,2,3)){
  for (j in c("S1","S2","S3")){
    if (j=="S1"){
      for (nonlin in c("comp","inv")){
        for (Rate in c(0.2,0.3,0.5)){
          tablesfun(norder=i, scene=j, bin=Binary, prevalence=Rate, nonlin=nonlin)
        }
      }
    } else {
      for (Rate in c(0.2,0.3,0.5)){
        tablesfun(norder=i, scene=j, bin=Binary, prevalence=Rate, nonlin="inv")
      }
    }

  }
}


