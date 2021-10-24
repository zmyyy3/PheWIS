setwd("/Users/zmyyy/Desktop/rotation_2/PheWAS/debug")
#library("PheWAS")
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)
#read files
age <- read.table("age.txt", header = T)
sex <- read.table("sex.txt", header = T)
phe <- read.table("phe.txt", header = T)
gen <- read.table("gen.txt", header = T)

#prepare data
exposers = age
exp=exposers
covariates=sex
cov=covariates
adjustments=list(NA)

#
id=intersect(intersect(names(phe),names(gen)),names(exp))
phenotypes=names(phe)
phenotypes=phenotypes[!(phenotypes %in% id)]
genotypes=names(gen)
genotypes=genotypes[!(genotypes %in% id)]
exposers=names(exp)
exposers=exposers[!(exposers %in% id)]


#
if(length(phenotypes)<1 || length(genotypes)<1 ||  length(exposers)<1)
{stop("Either phenotypes, genotypes or exposers contained no non-shared columns, yielding no variables for analysis after the merge.")}

data=merge((merge(phe,gen,by=id)),exp, by=id)
if(!is.null(names(covariates)[-1])) {
  covariates=names(covariates)
  if(sum(id %in% covariates)!=length(id)){stop(paste("The shared ID column(s) do not all exist in covariates: ",id))}
  covariates=covariates[!(covariates %in% id)]
  data=merge(data,cov,by=id)
}
if(!is.null(names(adjustments)[-1])) {
  adjustments=names(adjustments)
  if(sum(id %in% adjustments)!=length(id)){stop(paste("The shared ID column(s) do not all exist in adjustments: ",id))}
  adjustments=as.list(c(NA,adjustments[!(adjustments %in% id)]))
  data=merge(data,adjustment,by=id)}

full_list=data.frame(t(expand.grid(phenotypes,genotypes,exposers,adjustments,stringsAsFactors=F)),stringsAsFactors=F)
result=lapply(full_list,FUN=phe_as)
              ,data, covariates)

phe.gen=full_list[,1]
phe_o=phe.gen[[1]]
phe=phe_o
gen=phe.gen[[2]]
gens=gen
exp=phe.gen[[3]]  ## MZ
ex=exp  ## MZ
adjustment=phe.gen[[4]]  ## MZ

#Subset the data
d=data[,na.omit(unlist(c(gen,phe,exp,covariates,adjustment)))]  ##MZ
#Turn adjustment into a string, if not NA
if(!is.na(adjustment[1])) {adjustment=paste(adjustment,collapse=",")}else {adjustment=NA_character_} #Make sure it is a character NA for aggregation
#Alter phe_o if necessary for the regression formula
if(suppressWarnings(!is.na(as.numeric(phe_o)))) {
  phe=paste0("pheno_",phe_o)
  names(d)[2]=phe
}
#Exclude the exclusions for the target phenotype
d=d[!is.na(d[[phe]]),]
n_no_snp=sum(is.na(d[[gen]]))
#Exclude rows with missing data
d=na.omit(d)
n_total=nrow(d)
n_cases=NA_integer_
n_controls=NA_integer_
allele_freq=NA_real_
HWE_pval=NA_real_
or=NA_real_
se=NA_real_
p=NA_real_
beta=NA_real_
type=NA_character_
note=""
model=NA
if (length(covariates)==0) {model_name=sprintf("No model: %s ~ %s",phe,paste0(names(d)[c(1,3)],collapse = " * "))
}else{ if (length(covariates)==1){
  model_name=sprintf("No model: %s ~ %s",phe,paste(paste0(names(d)[c(1,3)],collapse = " * "),"+",names(d)[c(-1,-2,-3)]))
}else{ model_name=sprintf("No model: %s ~ %s",phe,paste(paste0(names(d)[c(1,3)],collapse = " * "),"+",paste0(names(d)[4:c(2+length(names(cov)))],collapse = "+")))
}
}
if(n_total<min.records) {
  note=paste(note,"[Error: <", min.records, " complete records]")
} else if(length(unique(na.omit(d[[phe]])))<=1 | length(unique(na.omit(d[[gen]]))) <=1) {
  note=paste(note,"[Error: non-varying phenotype or genotype]")
} else {
  if(additive.genotypes) {
    if(class(d[[gen]]) %in% c("numeric","integer")){
      allele_freq=sum(d[[gen]])/(2*n_total)
    }
    if(class(d[[gen]]) %in% c("numeric","integer") & sum(!(na.omit(d[[gen]]) %in% 0:2))==0) {
      P=allele_freq
      Q=1-allele_freq
      AA=sum(d[[gen]]==2)
      xAA=P^2*n_total
      Aa=sum(d[[gen]]==1)
      xAa=2*P*Q*n_total
      aa=sum(d[[gen]]==0)
      xaa=Q^2*n_total
      HWE_pval=pchisq((AA-xAA)^2/(xAA)+(Aa-xAa)^2/(xAa)+(aa-xaa)^2/(xaa),1)
    } else {note=paste(note,"[Warning: Genotype is not coded 0,1,2, but additive.genotypes was TRUE.]")}
  } 
}






if(n_cases<min.records|n_controls<min.records){output=data.frame(phenotype=phe_o,snp=gens,exposers=exp,
                                                                 covariates=covariates,  type=type,
                                                                 n_total=n_total, n_cases=n_cases, n_controls=n_controls,
                                                                 beta_phenotype_genotype, SE_phenotype_genotype,
                                                                 OR_phenotype_genotype,
                                                                 p_phenotype_genotype, 
                                                                 beta_phenotype_exposer, SE_phenotype_exposer,
                                                                 OR_phenotype_exposer,
                                                                 p_phenotype_exposer,
                                                                 beta_phenotype_gen_exp, SE_phenotype_gen_exp,
                                                                 OR_phenotype_gen_exp,
                                                                 p_phenotype_gen_exp,
                                                                 HWE_p=HWE_pval,allele_freq=allele_freq,n_no_snp=n_no_snp, 
                                                                 note=note, stringsAsFactors=F)
} else { output=data.frame(phenotype=phe_o,snp=gens,exposers=exp,
                           covariates=covariates,  type=type,
                           n_total=n_total, n_cases=n_cases, n_controls=n_controls,
                           beta_phenotype_genotype=beta_phenotype_genotype, SE_phenotype_genotype=se_phenotype_genotype,
                           OR_phenotype_genotype=or_phenotype_genotype,
                           p_phenotype_genotype=p_phenotype_genotype, 
                           beta_phenotype_exposer=beta_phenotype_exposer, SE_phenotype_exposer=se_phenotype_exposer,
                           OR_phenotype_exposer=or_phenotype_exposer,
                           p_phenotype_exposer=p_phenotype_exposer,
                           beta_phenotype_gen_exp=beta_phenotype_gen_exp, SE_phenotype_gen_exp=se_phenotype_gen_exp,
                           OR_phenotype_gen_exp=or_phenotype_gen_exp,
                           p_phenotype_gen_exp=p_phenotype_gen_exp,
                           HWE_p=HWE_pval,allele_freq=allele_freq,n_no_snp=n_no_snp, 
                           note=note, stringsAsFactors=F) }  ## MZ
 