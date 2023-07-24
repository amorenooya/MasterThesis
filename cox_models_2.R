library(dplyr)
library(glue)
library(stringr)
library(parallel)
library(R.utils)
library(survival)


argParam <- commandArgs(trailingOnly=TRUE)
numchr<- argParam[1] #chromosomes' numbers


setwd('/home/amorenoo/GWASP/NMIBC/EPICURO')

load(glue('geno1.epicuro_{numchr}.Rdata')) #imputed SNPS: tidy data
load(glue('geno_epicuro_info_{numchr}.Rdata')) #info and MAF of each SNP
load('epicuro_genomic_CNMIBC_n.RData') #clinical data
load('outputEPICURO.prcomp.RData') #PCA components

setwd('/home/amorenoo/GWASP/NMIBC/EPICURO/output_models')

Epicuro_Genomic_CNMIBC<-as_tibble(Epicuro_Genomic_CNMIBC)

geno1.epicuro<- geno1.epicuro %>% filter(id %in% Epicuro_Genomic_CNMIBC$EPICURO.ID) 
#modelo con variables geneticas

#recurrence

cox_snp_recurrence <- function(i){

  try({ #que continue aunque haya error: ok?

  GenoMod  <- coxph(Surv(TR1,indR1) ~  
                      as.numeric(unlist(geno1.epicuro[,i])) + #SNP: unlist to avoid error
                      #very important to convert it to numeric/factor! when transposing to obtain geno.epicuro
                      #all SNPs are converted to character
                      #de momento numerica y sin redondear
                      age  +
                      gender +
                      TG_final_recurrence+
                      tumor_size+
                      number_tumor+
                      treatment_recurrence+
                      area,
                     data=Epicuro_Genomic_CNMIBC)

  sum<- summary(GenoMod)

  output<-tibble('var'=colnames(geno1.epicuro)[i],
                 'SNP_new'=geno_epicuro_info$SNP_new[i-1],
                 'coef'=sum$coefficients[1],
                 'p_value'= sum$coefficients[1,5], #1 because it is the first variable
                 'CI_L'= sum$conf.int[1,3], #lower bound IC for HR (95%)
                 'CI_U'= sum$conf.int[1,4], #upper bound IC for HR (95%)
                 'se_coef'=sum$coefficients[1,3],
                 'MAF'=geno_epicuro_info$MAF[i-1],
                 'info'=geno_epicuro_info$info[i-1]) %>%
  mutate(coef=as.numeric(coef),
         CI_L=as.numeric(CI_L),
         CI_U=as.numeric(CI_U),
         MAF=as.numeric(MAF),
         info=as.numeric(info))

  return(output)

  }, silent=TRUE)
}


result<-mclapply(2:ncol(geno1.epicuro), #from 2: first column is id
         cox_snp_recurrence,
         mc.cores=16) #list: for each SNP we have a value (tibble) in the list

result<-do.call(rbind,result) %>%
  as_tibble()

rm(geno1.epicuro_rec)

save(result, file=glue('cox_recurrence_{numchr}.Rdata')) #result: chr21, chr22, 2:100 rows

#for example: to filter SNPs in chr21:
# result %>% filter(str_detect(var, 'chr21'))

#progression

cox_snp_progression <- function(i){
  
  try({ #que continue aunque haya error: ok?
    
    GenoMod  <- coxph(Surv(TPFS,indPFS) ~  
                        as.numeric(unlist(geno1.epicuro[,i])) + #SNP: unlist to avoid error
                        #very important to convert it to numeric/factor! when transposing to obtain geno.epicuro
                        #all SNPs are converted to character
                        #de momento numerica y sin redondear
                        age  + 
                        gender + 
                        TG_final_progression+
                        treatment_progression+
                        number_tumor+
                        nRec_new+
                        area,
                      data=Epicuro_Genomic_CNMIBC))
    
    sum<- summary(GenoMod)
    
    output<-tibble('var'=colnames(geno1.epicuro)[i],
                   'SNP_new'=geno_epicuro_info$SNP_new[i-1],
                   'coef'=sum$coefficients[1],
                   'p_value'= sum$coefficients[1,5], #1 because it is the first variable
                   'CI_L'= sum$conf.int[1,3], #lower bound IC for HR (95%)
                   'CI_U'= sum$conf.int[1,4], #upper bound IC for HR (95%)
                   'se_coef'=sum$coefficients[1,3],
                   'MAF'=geno_epicuro_info$MAF[i-1],
                   'info'=geno_epicuro_info$info[i-1]) %>%
      mutate(coef=as.numeric(coef),
             CI_L=as.numeric(CI_L),
             CI_U=as.numeric(CI_U),
             MAF=as.numeric(MAF),
             info=as.numeric(info))
    
    return(output)
    
  }, silent=TRUE)

}

result<-mclapply(2:ncol(geno1.epicuro), #from 2: first column is id
                 cox_snp_progression,
                 mc.cores=16) #list: for each SNP we have a value (tibble) in the list

result<-do.call(rbind,result) %>%
  as_tibble()

save(result, file=glue('cox_progression_{numchr}.Rdata'))
