#COX PH models: recurrence, progression and relapse
#libraries
library(dplyr)
library(glue)
library(stringr)
library(parallel)
library(R.utils)
library(survival)

#argument for parallelization
argParam <- commandArgs(trailingOnly=TRUE)
numchr<- argParam[1] #chromosomes' numbers

#datasets
#Epicuro: dataset with clinical information
#SNPs_info: information about SNPs (MAF and imputation quality)

###### RECURRENCE #######
cox_snp_recurrence <- function(i){
  try({ 
  GenoMod  <- coxph(Surv(TR1,indR1) ~  
                      as.numeric(unlist(geno1.epicuro[,i])) + #SNP
                      age  +
                      gender +
                      TG_recurrence+
                      tumor_size+
                      number_tumor+
                      treatment+
                      area,
                     data=Epicuro)

  sum<- summary(GenoMod)
  output<-tibble('var'=colnames(geno1.epicuro)[i],
                 'coef'=sum$coefficients[1],
                 'p_value'= sum$coefficients[1,5],
                 'CI_L'= sum$conf.int[1,3], #lower bound IC for HR (95%)
                 'CI_U'= sum$conf.int[1,4], #upper bound IC for HR (95%)
                 'se_coef'=sum$coefficients[1,3],
                 'MAF'=SNPs_info$MAF[i-1],
                 'info'=SNPs_info$info[i-1]) %>%
  mutate(coef=as.numeric(coef),
         CI_L=as.numeric(CI_L),
         CI_U=as.numeric(CI_U),
         MAF=as.numeric(MAF),
         info=as.numeric(info))

  return(output)

  }, silent=TRUE)
}


###### PROGRESSION #######
cox_snp_progression <- function(i){
  try({ 
  GenoMod  <- coxph(Surv(TPFS,indPFS) ~  
                      as.numeric(unlist(geno1.epicuro[,i])) + #SNP
                      age  +
                      gender +
                      TG_progression+
                      number_recurrences+
                      number_tumor+
                      treatment+
                      area,
                     data=Epicuro)

  sum<- summary(GenoMod)
  output<-tibble('var'=colnames(geno1.epicuro)[i],
                 'coef'=sum$coefficients[1],
                 'p_value'= sum$coefficients[1,5],
                 'CI_L'= sum$conf.int[1,3], #lower bound IC for HR (95%)
                 'CI_U'= sum$conf.int[1,4], #upper bound IC for HR (95%)
                 'se_coef'=sum$coefficients[1,3],
                 'MAF'=SNPs_info$MAF[i-1],
                 'info'=SNPs_info$info[i-1]) %>%
  mutate(coef=as.numeric(coef),
         CI_L=as.numeric(CI_L),
         CI_U=as.numeric(CI_U),
         MAF=as.numeric(MAF),
         info=as.numeric(info))

  return(output)

  }, silent=TRUE)
}

###### RELAPSE #######
cox_snp_relapse <- function(i){
  try({ 
  GenoMod  <- coxph(Surv(TEFS,indEFS) ~  
                      as.numeric(unlist(geno1.epicuro[,i])) + #SNP
                      age  +
                      gender +
                      TG_recurrence+
                      tumor_size+
                      number_tumor+
                      treatment+
                      area,
                     data=Epicuro)

  sum<- summary(GenoMod)
  output<-tibble('var'=colnames(geno1.epicuro)[i],
                 'coef'=sum$coefficients[1],
                 'p_value'= sum$coefficients[1,5],
                 'CI_L'= sum$conf.int[1,3], #lower bound IC for HR (95%)
                 'CI_U'= sum$conf.int[1,4], #upper bound IC for HR (95%)
                 'se_coef'=sum$coefficients[1,3],
                 'MAF'=SNPs_info$MAF[i-1],
                 'info'=SNPs_info$info[i-1]) %>%
  mutate(coef=as.numeric(coef),
         CI_L=as.numeric(CI_L),
         CI_U=as.numeric(CI_U),
         MAF=as.numeric(MAF),
         info=as.numeric(info))

  return(output)

  }, silent=TRUE)
}


result_recurrence<-mclapply(2:ncol(geno1.epicuro), #from 2: first column is id
         cox_snp_recurrence,
         mc.cores=16) 
result_recurrence<-do.call(rbind,result_recurrence) %>%
  as_tibble()

result_progression<-mclapply(2:ncol(geno1.epicuro), #from 2: first column is id
         cox_snp_progression,
         mc.cores=16) 
result_progression<-do.call(rbind,result_progression) %>%
  as_tibble()

result_relapse<-mclapply(2:ncol(geno1.epicuro), #from 2: first column is id
         cox_snp_relapse,
         mc.cores=16) 
result_relapse<-do.call(rbind,result_relapse) %>%
  as_tibble()




