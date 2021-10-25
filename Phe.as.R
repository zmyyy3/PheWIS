phe_as <-
  function(phe.gen, additive.genotypes=T,min.records=20,return.models=F,confint.level=NA, my.data, my.covariates) {
    if(!missing(my.data)) data=my.data
    if(!missing(my.covariates)) covariates=my.covariates
    #Retrieve the targets for this loop
    phe_o=phe.gen[[1]]
    phe=phe_o
    gen=phe.gen[[2]]
    gens=gen
    exp=phe.gen[[3]]  ## MZ
    ex=exp  ## MZ
    #covariates=phe.gen[[4]]
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
    beta_phenotype_genotype=NA_real_  ##MZ
    se_phenotype_genotype=NA_real_
    or_phenotype_genotype=NA_real_
    p_phenotype_genotype=NA_real_
    beta_phenotype_exposure=NA_real_
    se_phenotype_exposure=NA_real_
    or_phenotype_exposure=NA_real_
    p_phenotype_exposure=NA_real_
    beta_phenotype_gen_exp=NA_real_
    se_phenotype_gen_exp=NA_real_
    or_phenotype_gen_exp=NA_real_
    p_phenotype_gen_exp=NA_real_  ##MZ
    type=NA_character_
    note=""
    model=NA
    if (length(covariates)==0) {model_name=sprintf("No model: %s ~ %s",phe,paste0(names(d)[c(1,3)],collapse = " * "))
     }else{ if (length(covariates)==1){
      model_name=sprintf("No model: %s ~ %s",phe,paste(paste0(names(d)[c(1,3)],collapse = " * "),"+",names(d)[c(-1,-2,-3)]))
    }else{ model_name=sprintf("No model: %s ~ %s",phe,paste(paste0(names(d)[c(1,3)],collapse = " * "),"+",paste0(names(d)[4:c(2+length(names(cov)))],collapse = "+")))
    }  ## MZ
    }
    if(n_total< min.records) {
      note=paste(note,"[Error: <",min.records," complete records]")
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
    
      #Check if genotype was available
      #Check if phenotype is logical (boolean)
      if(class(d[[phe]]) %in% c("logical")) {
        type = "logistic"
        #Create the logistic model
        n_cases=sum(d[[phe]])
        n_controls=n_total-n_cases
        if(n_cases<min.records|n_controls<min.records) {note=paste(note,"[Error: <",min.records," cases or controls]")}else{ 
          if (length(covariates)==0){model=glm(as.formula(sprintf("%s ~ %s",phe,paste0(names(d)[c(1,3)],collapse = " * "))), data=d, family=binomial)
          } else{ if (length(covariates)==1){
            model=glm(as.formula(sprintf("%s ~ %s",phe,paste(paste0(names(d)[c(1,3)],collapse = " * "),"+",names(d)[c(-1,-2,-3)]))), data=d, family=binomial)
          }else{ model=glm(as.formula(sprintf("%s ~ %s",phe,paste(paste0(names(d)[c(1,3)],collapse = " * "),"+",paste0(names(d)[4:c(2+length(names(cov)))],collapse = "+")))), data=d, family=binomial)
          }  ## MZ
          }  ## MZ
          ## MZ
          modsum= summary(model)
          model_name=paste0(as.character(terms(model))[c(2,1,3)],collapse=" ")
          #If the models did not converge, report NA values instead.
          if(model$converged) {
            #Find the rows with results that gets merged across all loops
            gen_list=sort(unique(c(grep(gen,row.names(modsum$coef)),grep(exp,row.names(modsum$coef))))) ## MZ
            gens=row.names(modsum$coef)[gen_list][1]   ## MZ
            or_phenotype_genotype=exp(modsum$coef[gen_list,1][1])   ## MZ
            beta_phenotype_genotype=modsum$coef[gen_list,1][1]   ## MZ
            se_phenotype_genotype=modsum$coef[gen_list,2][1]  ## MZ
            p_phenotype_genotype=modsum$coef[gen_list,4][1]   ## MZ
            or_phenotype_exposure=exp(modsum$coef[gen_list,1][2])  ## MZ
            beta_phenotype_exposure=modsum$coef[gen_list,1][2]  ## MZ
            se_phenotype_exposure=modsum$coef[gen_list,2][2]  ## MZ
            p_phenotype_exposure=modsum$coef[gen_list,4][2]  ## MZ
            or_phenotype_gen_exp=exp(modsum$coef[gen_list,1][3]) ## MZ
            beta_phenotype_gen_exp=modsum$coef[gen_list,1][3]  ## MZ
            se_phenotype_gen_exp=modsum$coef[gen_list,2][3]  ## MZ
            p_phenotype_gen_exp=modsum$coef[gen_list,4][3] ## MZ
          } else {
            note=paste(note,"[Error: The model did not converge]")
          }
        }
      } else {
        type = "linear"
        if(n_total<min.records) {
          note=paste(note,"[Error: <",min.records," records with phenotype and genotype]")
        } else {
          model = glm(as.formula(paste(phe," ~ .", sep="", collapse="")), data=d)
          modsum= summary(model)
          model_name=paste0(as.character(terms(model))[c(2,1,3)],collapse=" ")
          
          #If the models did not converge, report NA values instead.
          if(model$converged) {
            #Find the rows with results that gets merged across all loops
            gen_list=grep(gen,row.names(modsum$coef))
            gens=row.names(modsum$coef)[gen_list]
            beta=modsum$coef[gen_list,1]
            se=modsum$coef[gen_list,2]
            p=modsum$coef[gen_list,4]
          } else {
            note=paste(note,"[Error: The model did not converge]")
          }
        }
      }
    }  
    
   output=data.frame(phenotype=phe_o,snp=gens,exposures=exp,
                  covariates=covariates,  type=type,
                  n_total=n_total, n_cases=n_cases, n_controls=n_controls,
                  beta_phenotype_genotype=beta_phenotype_genotype, 
                  SE_phenotype_genotype=se_phenotype_genotype,
                  OR_phenotype_genotype=or_phenotype_genotype,
                  p_phenotype_genotype=p_phenotype_genotype, 
                  beta_phenotype_exposure=beta_phenotype_exposure, 
                  SE_phenotype_exposure=se_phenotype_exposure,
                  OR_phenotype_exposure=or_phenotype_exposure,
                  p_phenotype_exposure=p_phenotype_exposure,
                  beta_phenotype_gen_exp=beta_phenotype_gen_exp, 
                  SE_phenotype_gen_exp=se_phenotype_gen_exp,
                  OR_phenotype_gen_exp=or_phenotype_gen_exp,
                  p_phenotype_gen_exp=p_phenotype_gen_exp,
                  HWE_p=HWE_pval,allele_freq=allele_freq,n_no_snp=n_no_snp, 
                  note=note, stringsAsFactors=F)  ## MZ
    
    #Add confidenested.
   if(!is.na(confint.level)) {
     if(!is.na(model)[1]){
       suppressMessages(conf<-confint(model,c("(Intercept)",gens),level=confint.level))
       lower=conf[-1,1]
       upper=conf[-1,2]
       if(type=="logistic") {
         lower=exp(lower)
         upper=exp(upper)
       }
     } else {
       lower=NA_real_
       upper=NA_real_
     }
     output$lower=lower
     output$upper=upper
     
     output=output[,c("phenotype","snp","adjustment","beta","SE",
                      "lower","upper","OR","p","type",
                      "n_total","n_cases","n_controls",
                      "HWE_p","allele_freq","n_no_snp","note")]
   }
    
    #If the complete models were requested, add them as well.
    if(return.models) {
      attributes(output)$model=model
      attributes(output)$model_name=model_name
    }
    attributes(output)$successful.phenotype=ifelse(is.na(p_phenotype_genotype),NA,phe_o)
    attributes(output)$successful.genotype=ifelse(is.na(p_phenotype_genotype),NA,gen)
    #Return this to the loop to be merged.
    output
}


