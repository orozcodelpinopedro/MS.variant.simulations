#Pedro Orozco del Pino
#resample from a vcf with replacement 
#input: raw file, snp file, prevalence, effect size, ncases
#output: beta estimates
#arg1<-"DATA/CHR_221EUR.raw"
#arg2<-"DATA/snp221.frq"
#arg3<-0.1
#arg4<-0.0954
#arg5<-30000
#arg6<-"221"
#arg7<-"X_sim_DEBUG"
#arg8<-1
#setwd("/Users/pedro/Box Sync/Zöllner Lab/Project 1/")
# # 
args = commandArgs(trailingOnly=TRUE)
#args<-c(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
file1<-args[1] #the name of the raw file from population 1
code<-args[6]
SCRATCH<-"/net/wonderland/home/porozco/GreatLakes/scratch/"
fileAFR<-paste(SCRATCH,"out_f/CHR_",code,"AFR.raw",sep="")
fileAMR<-paste(SCRATCH,"out_f/CHR_",code,"AMR.raw",sep="")
fileEAS<-paste(SCRATCH,"out_f/CHR_",code,"EAS.raw",sep="")
fileSAS<-paste(SCRATCH,"out_f/CHR_",code,"SAS.raw",sep="")


file2<-args[2] #the name of the file with the selected snps
beta<-as.numeric(args[4]) #The effect size in log(ODS) scale
library(data.table)
library(tidyverse)
snps<-as.data.frame(fread(file2))
col<-paste(snps$SNP,"_",sep="")
k<-dim(snps)[1]
n<-as.numeric(args[5])#sample size for cases and controls
#######################################################################################################
#######################################################################################################
############################################################
##..................------------------....................##
############################################################
############ Sampling cases with replacement
get.cases2<-function(dat,G,beta,beta0,P.G,P.A,n.cases){
  P.AG<-exp(G*beta+beta0)/(1+exp(G*beta+beta0)) #probability of beign affected given Genotype
  ############ 
  P.GA<-P.AG*P.G/P.A # probability of genotype given status of disease is case
  P.GA<-P.GA/sum(unique(P.GA))
  N<-merge(data.frame(dat$IID,G,P.GA),data.frame(table(G)),by.x=names(G),by.y="G",sort=F)[,-1]
  
  ind<-data.frame(FID=sample(N[,1],size=n.cases,prob=N[,2]/N[,3],replace=T))
  CASES<-merge(dat,ind)
  CASES$PHENOTYPE=1
  return(CASES)
}
############################################################
##..................------------------....................##
############################################################
############ Sampling controls with replacement
get.controls2<-function(dat,G,beta,beta0,P.G,P.A,n.controls){
  P.AG<-1/(1+exp(G*beta+beta0)) #probability of beign affected given Genotype
  ############ 
  P.GA<-P.AG*P.G/(1-P.A) # probability of genotype given status of disease is case
  P.GA<-P.GA/sum(unique(P.GA))
  N<-merge(data.frame(dat$IID,G,P.GA),data.frame(table(G)),by.x=names(G),by.y="G",sort=F)[,-1]
  ind<-data.frame(FID=sample(N[,1],size=n.controls,prob=N[,2]/N[,3],replace=T))
  CONTROLS<-merge(dat,ind)
  CONTROLS$PHENOTYPE=0
  return(CONTROLS)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################

dat1<-as.data.frame(fread(file1))
NAM<-names(dat1)[-c(1:6)]
dat1<-dat1[,c(1:6,(6+grep("rs",NAM)))]
NAM<-names(dat1)[-c(1:6)]
names(dat1)[-c(1:6)]<-substr(NAM,1,regexpr("_",NAM))
G1<-as.data.frame(dat1[,names(dat1)%in%col])
names(G1)<-snps$SNP  
P.A1<-as.numeric(args[3]) # probability of being affected, i.e. prevalence of the disease in Population 1
allele.fr1<-colSums(G1)/(dim(G1)[1]*2)
#To be revisited when k>1
P.G1<-unlist(lapply(G1,function(x){dbinom(x,2,allele.fr1)}))  #P.G probability of genotype
f.b<-function(x){sum((P.A1-exp(x+beta*G1)/(1+exp(x+beta*G1)*P.G1))^2)} #obtain intercept 
beta01<-optimize(f.b,c(-10,10))$minimum
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#Simulate data from population 2
# print("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")
# print("-*-*-*-*-*  AFRICAN   POPULATION      -*-*-*-*-*-*-*-*-*-*")
dat2<-as.data.frame(fread(fileAFR))
NAM<-names(dat2)[-c(1:6)]
dat2<-dat2[,c(1:6,(6+grep("rs",NAM)))]
NAM<-names(dat2)[-c(1:6)]
names(dat2)[-c(1:6)]<-substr(NAM,1,regexpr("_",NAM))
G2<-as.data.frame(dat2[,names(dat2)%in%col])
names(G2)<-snps$SNP  
P.A2<-as.numeric(args[3]) # probability of being affected, i.e. prevalence of the disease
# P.A2<-as.numeric(arg3) # probability of being affected, i.e. prevalence of the disease
allele.fr2<-colSums(G2)/(dim(G2)[1]*2)
#To be revisited when k>1
P.G2<-unlist(lapply(G2,function(x){dbinom(x,2,allele.fr2)}))  #P.G probability of genotype
f.b<-function(x){sum((P.A2-exp(x+beta*G2)/(1+exp(x+beta*G2)*P.G2))^2)} #obtain intercept 
beta02<-optimize(f.b,c(-10,10))$minimum



###############################################################################################################################################################################################################################################################################
## print("-*-*-*-*-*  AMERICAN   POPULATION      -*-*-*-*-*-*-*-*-*-*")
dat3<-as.data.frame(fread(fileAMR))
NAM<-names(dat3)[-c(1:6)]
dat3<-dat3[,c(1:6,(6+grep("rs",NAM)))]
NAM<-names(dat3)[-c(1:6)]
names(dat3)[-c(1:6)]<-substr(NAM,1,regexpr("_",NAM))
G3<-as.data.frame(dat3[,names(dat3)%in%col])
names(G3)<-snps$SNP  
P.A3<-as.numeric(args[3]) # probability of being affected, i.e. prevalence of the disease
# P.A2<-as.numeric(arg3) # probability of being affected, i.e. prevalence of the disease
allele.fr3<-colSums(G3)/(dim(G3)[1]*2)
#To be revisited when k>1
P.G3<-unlist(lapply(G3,function(x){dbinom(x,2,allele.fr3)}))  #P.G probability of genotype
f.b<-function(x){sum((P.A3-exp(x+beta*G3)/(1+exp(x+beta*G3)*P.G3))^2)} #obtain intercept 
beta03<-optimize(f.b,c(-10,10))$minimum



###############################################################################################################################################################################################################################################################################
### print("-*-*-*-*-*  EAST ASIAN   POPULATION      -*-*-*-*-*-*-*-*-*-*")
dat4<-as.data.frame(fread(fileEAS))
NAM<-names(dat4)[-c(1:6)]
dat4<-dat4[,c(1:6,(6+grep("rs",NAM)))]
NAM<-names(dat4)[-c(1:6)]
names(dat4)[-c(1:6)]<-substr(NAM,1,regexpr("_",NAM))
G4<-as.data.frame(dat4[,names(dat4)%in%col])
names(G4)<-snps$SNP  
P.A4<-as.numeric(args[3]) # probability of being affected, i.e. prevalence of the disease
# P.A2<-as.numeric(arg3) # probability of being affected, i.e. prevalence of the disease
allele.fr4<-colSums(G4)/(dim(G4)[1]*2)
#To be revisited when k>1
P.G4<-unlist(lapply(G4,function(x){dbinom(x,2,allele.fr4)}))  #P.G probability of genotype
f.b<-function(x){sum((P.A4-exp(x+beta*G4)/(1+exp(x+beta*G4)*P.G4))^2)} #obtain intercept 
beta04<-optimize(f.b,c(-10,10))$minimum



###############################################################################################################################################################################################################################################################################
#### print("-*-*-*-*-*  SOUTH ASIAN   POPULATION      -*-*-*-*-*-*-*-*-*-*")
dat5<-as.data.frame(fread(fileSAS))
NAM<-names(dat5)[-c(1:6)]
dat5<-dat5[,c(1:6,(6+grep("rs",NAM)))]
NAM<-names(dat5)[-c(1:6)]
names(dat5)[-c(1:6)]<-substr(NAM,1,regexpr("_",NAM))
G5<-as.data.frame(dat5[,names(dat5)%in%col])
names(G5)<-snps$SNP  
P.A5<-as.numeric(args[3]) # probability of being affected, i.e. prevalence of the disease
# P.A2<-as.numeric(arg3) # probability of being affected, i.e. prevalence of the disease
allele.fr5<-colSums(G5)/(dim(G5)[1]*2)
#To be revisited when k>1
P.G5<-unlist(lapply(G5,function(x){dbinom(x,2,allele.fr5)}))  #P.G probability of genotype
f.b<-function(x){sum((P.A5-exp(x+beta*G5)/(1+exp(x+beta*G5)*P.G5))^2)} #obtain intercept 
beta05<-optimize(f.b,c(-10,10))$minimum



###############################################################################################################################################################################################################################################################################
###############################################################################################################################################################################################################################################################################
###############################################################################################################################################################################################################################################################################
###############################################################################################################################################################################################################################################################################
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                                                               REPETITIONS
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


B<-as.numeric(args[8])# Number of repetitions
k<-1
betac.1<-c(NULL);beta.1<-c(NULL)
betac.2<-c(NULL);beta.2<-c(NULL)
for(h in 1:B){
Time0<-Sys.time()	
        ############################################################
        ##..................------------------....................##
        ############################################################
        ############ Population1
	CASES<-get.cases2(dat1,G1,beta,beta01,P.G1,P.A1,as.numeric(n))#Number of cases
	
	
	CONTROLS<-get.controls2(dat1,G1,beta,beta01,P.G1,P.A1,as.numeric(n))#Number of cases
	
	
	ALL<-rbind(CASES,CONTROLS)
	j<-1
	p.tres<-5e-8#p-value threshold genome wide
	z.tres<-abs(qnorm(p.tres))
	p.cases<-colSums(CASES[,-c(1:6)],na.rm = T)/(2*n)
	p.controls<-colSums(CONTROLS[,-c(1:6)],na.rm = T)/(2*n)
	unval.cases<-which(p.cases==0|p.cases==1)
	unval.control<-which(p.controls==0|p.controls==1)
	unval<-unique(c(unval.cases,unval.control))
	p.cases<-p.cases[-unval]
	p.controls<-p.controls[-unval]
	sd.cases<-p.cases*(1-p.cases)
	sd.controls<-p.controls*(1-p.controls)
	Z.scores<-abs(sqrt(2*n)*(p.cases-p.controls)/sqrt(sd.controls+sd.cases))
	p.val<-1-pnorm(abs(Z.scores))
	t.snps<-names(ALL)[-c(1:6)]
	t.snps<-t.snps[-unval]
	
	#DON'T INCLUDE Casual SNP
	significant<-which(Z.scores>z.tres&t.snps!=col)
	K<-length(significant)
	#Set of tag snps excluding the causal SNP
	tag.snp<-t.snps[significant][which.max(Z.scores[significant])]
	#Is the causal SNP the most significant?
	# Msig<-which(p.val<p.tres)#all significant snps, including causal snp
	allele.fr1<-mean(ALL[,col],na.rm=T)/2
	if(K>0){
	  print("or causal or tag or both")
	  flag=ifelse(Z.scores[col]>Z.scores[tag.snp],"causal",
	                 ifelse(Z.scores[col]<Z.scores[tag.snp],"tag",
	                        "causal/tag"))}
	else
	  if(Z.scores[col]>z.tres){
	    print("causal/only")
	    flag="causal/only"
	  }else{
	    print("none")
	    flag="none"
	  }
	AF.cases.t<-mean(CASES[,tag.snp],na.rm=T)/2
        AF.controls.t<-mean(CONTROLS[,tag.snp],na.rm=T)/2
        AF.cases.c<-mean(CASES[,col],na.rm=T)/2
        AF.controls.c<-mean(CONTROLS[,col],na.rm=T)/2	
	if(K>0){
	
		###########################################################
		##..................-----------------....................##
		###########################################################
		############ Population1
		allele.fr1.t<-mean(ALL[,tag.snp],na.rm=T)/2
		cor.1<-cor(ALL[,col],ALL[,tag.snp],use = "complete.obs")^2
		betac.1[k]<-glm(ALL$PHENOTYPE~ALL[,col],family=binomial)$coefficients[2]
		beta.1[k]<-glm(ALL$PHENOTYPE~ALL[,tag.snp],family=binomial)$coefficients[2]
		############################################################
		##..................------------------....................##
		############################################################
		############ AFRICAN Population
		
        	CASES<-get.cases2(dat2,G2,beta,beta02,P.G2,P.A2,as.numeric(n))
        	# CASES<-get.cases(dat2,G2,beta,beta02,P.G2,P.A2,as.numeric(n))
	  
  
    
        	CONTROLS<-get.controls2(dat2,G2,beta,beta02,P.G2,P.A2,as.numeric(n))
        	# CONTROLS<-get.controls(dat2,G2,beta,beta02,P.G2,P.A2,as.numeric(n))
    
  
    
		ALL<-rbind(CASES,CONTROLS)

		allele.fr2.t<-mean(ALL[,tag.snp],na.rm=T)/2
		allele.fr2<-mean(ALL[,col],na.rm=T)/2
       	 	if(sum(names(ALL)==tag.snp)>0){
			mc<-glm(ALL$PHENOTYPE~ALL[,col],family=binomial)
			mt<-glm(ALL$PHENOTYPE~ALL[,tag.snp],family=binomial)
  		      	betac.2[k]<-mc$coefficients[2]
  		        beta.2[k]<-mt$coefficients[2]
			r2.c<-1 - mc$deviance / mc$null.deviance
			r2.t<-1 - mt$deviance / mt$null.deviance
			file.out<-as.character(args[7])#file to write out
			cor.2<-cor(ALL[,col],ALL[,tag.snp],use = "complete.obs")^2
			print(paste("al escribir caso 1",snps$SNP))
			write.table(as.data.frame(cbind(snps$SNP,beta,betac.1[k],betac.2[k],allele.fr1,allele.fr2,tag.snp,allele.fr1.t,allele.fr2.t,
				beta.1[k],beta.2[k],flag,"CT",cor.1,cor.2,"AFR",K,Z.scores[col],Z.scores[tag.snp],AF.cases.c,AF.controls.c,AF.cases.t,AF.controls.t,r2.c,r2.t)),file=args[7],col.names=F,row.names=F,quote=F,append=TRUE)
  		        
       	 	}
			############################################################
		##..................------------------....................##
		############################################################
		############ AMERICAN Population
		
        	CASES<-get.cases2(dat3,G3,beta,beta03,P.G3,P.A3,as.numeric(n))
        	# CASES<-get.cases(dat2,G2,beta,beta02,P.G2,P.A2,as.numeric(n))
	  
  
    
        	CONTROLS<-get.controls2(dat3,G3,beta,beta03,P.G3,P.A3,as.numeric(n))
        	# CONTROLS<-get.controls(dat2,G2,beta,beta02,P.G2,P.A2,as.numeric(n))
    
  
    
		ALL<-rbind(CASES,CONTROLS)
	
		allele.fr2.t<-mean(ALL[,tag.snp],na.rm=T)/2
		allele.fr2<-mean(ALL[,col],na.rm=T)/2
       	 	if(sum(names(ALL)==tag.snp)>0){
			mc<-glm(ALL$PHENOTYPE~ALL[,col],family=binomial)
                        mt<-glm(ALL$PHENOTYPE~ALL[,tag.snp],family=binomial)
                        betac.2[k]<-mc$coefficients[2]
                        beta.2[k]<-mt$coefficients[2]
                        r2.c<-1 - mc$deviance / mc$null.deviance
                        r2.t<-1 - mt$deviance / mt$null.deviance
			file.out<-as.character(args[7])#file to write out
			cor.2<-cor(ALL[,col],ALL[,tag.snp],use = "complete.obs")^2
			write.table(as.data.frame(cbind(snps$SNP,beta,betac.1[k],betac.2[k],allele.fr1,allele.fr2,tag.snp,allele.fr1.t,allele.fr2.t,
				beta.1[k],beta.2[k],flag,"CT",cor.1,cor.2,"AMR",K,Z.scores[col],Z.scores[tag.snp],AF.cases.c,AF.controls.c,AF.cases.t,AF.controls.t,r2.c,r2.t)),file=args[7],col.names=F,row.names=F,quote=F,append=TRUE)
  		        
       	 	}
		##..................------------------....................##
		############################################################
		############ EAST ASIAN Population
		
        	CASES<-get.cases2(dat4,G4,beta,beta04,P.G4,P.A4,as.numeric(n))
        	# CASES<-get.cases(dat2,G2,beta,beta02,P.G2,P.A2,as.numeric(n))
	  
  
    
        	CONTROLS<-get.controls2(dat4,G4,beta,beta04,P.G4,P.A4,as.numeric(n))
        	# CONTROLS<-get.controls(dat2,G2,beta,beta02,P.G2,P.A2,as.numeric(n))
    
  
    
		ALL<-rbind(CASES,CONTROLS)
	
		allele.fr2.t<-mean(ALL[,tag.snp],na.rm=T)/2
		allele.fr2<-mean(ALL[,col],na.rm=T)/2
       	 	if(sum(names(ALL)==tag.snp)>0){
			mc<-glm(ALL$PHENOTYPE~ALL[,col],family=binomial)
                        mt<-glm(ALL$PHENOTYPE~ALL[,tag.snp],family=binomial)
                        betac.2[k]<-mc$coefficients[2]
                        beta.2[k]<-mt$coefficients[2]
                        r2.c<-1 - mc$deviance / mc$null.deviance
                        r2.t<-1 - mt$deviance / mt$null.deviance
			file.out<-as.character(args[7])#file to write out
			cor.2<-cor(ALL[,col],ALL[,tag.snp],use = "complete.obs")^2
			write.table(as.data.frame(cbind(snps$SNP,beta,betac.1[k],betac.2[k],allele.fr1,allele.fr2,tag.snp,allele.fr1.t,allele.fr2.t,
				beta.1[k],beta.2[k],flag,"CT",cor.1,cor.2,"EAS",K,Z.scores[col],Z.scores[tag.snp],AF.cases.c,AF.controls.c,AF.cases.t,AF.controls.t,r2.c,r2.t)),file=args[7],col.names=F,row.names=F,quote=F,append=TRUE)
  		        
       	 	}
		##..................------------------....................##
		############################################################
		############ SOUTH ASIAN Population
		
        	CASES<-get.cases2(dat5,G5,beta,beta05,P.G5,P.A5,as.numeric(n))
        	# CASES<-get.cases(dat2,G2,beta,beta02,P.G2,P.A2,as.numeric(n))
	  
  
    
        	CONTROLS<-get.controls2(dat5,G5,beta,beta05,P.G5,P.A5,as.numeric(n))
        	# CONTROLS<-get.controls(dat2,G2,beta,beta02,P.G2,P.A2,as.numeric(n))
    
  
    
		ALL<-rbind(CASES,CONTROLS)
	
		allele.fr2.t<-mean(ALL[,tag.snp],na.rm=T)/2
		allele.fr2<-mean(ALL[,col],na.rm=T)/2
       	 	if(sum(names(ALL)==tag.snp)>0){
			mc<-glm(ALL$PHENOTYPE~ALL[,col],family=binomial)
                        mt<-glm(ALL$PHENOTYPE~ALL[,tag.snp],family=binomial)
                        betac.2[k]<-mc$coefficients[2]
                        beta.2[k]<-mt$coefficients[2]
                        r2.c<-1 - mc$deviance / mc$null.deviance
                        r2.t<-1 - mt$deviance / mt$null.deviance
			file.out<-as.character(args[7])#file to write out
			cor.2<-cor(ALL[,col],ALL[,tag.snp],use = "complete.obs")^2
			write.table(as.data.frame(cbind(snps$SNP,beta,betac.1[k],betac.2[k],allele.fr1,allele.fr2,tag.snp,allele.fr1.t,allele.fr2.t,
				beta.1[k],beta.2[k],flag,"CT",cor.1,cor.2,"SAS",K,Z.scores[col],Z.scores[tag.snp],AF.cases.c,AF.controls.c,AF.cases.t,AF.controls.t,r2.c,r2.t)),file=args[7],col.names=F,row.names=F,quote=F,append=TRUE)
  		        
       	 	}
		
		##..................------------------....................##
		############################################################
		############ EUROPEAN Population
		
        	CASES<-get.cases2(dat1,G1,beta,beta01,P.G1,P.A1,as.numeric(n))
        	# CASES<-get.cases(dat2,G2,beta,beta02,P.G2,P.A2,as.numeric(n))
	  
  
    
        	CONTROLS<-get.controls2(dat1,G1,beta,beta01,P.G1,P.A1,as.numeric(n))
        	# CONTROLS<-get.controls(dat2,G2,beta,beta02,P.G2,P.A2,as.numeric(n))
    
  
    
		ALL<-rbind(CASES,CONTROLS)
	
		allele.fr2.t<-mean(ALL[,tag.snp],na.rm=T)/2
		allele.fr2<-mean(ALL[,col],na.rm=T)/2
       	 	if(sum(names(ALL)==tag.snp)>0){
			mc<-glm(ALL$PHENOTYPE~ALL[,col],family=binomial)
                        mt<-glm(ALL$PHENOTYPE~ALL[,tag.snp],family=binomial)
                        betac.2[k]<-mc$coefficients[2]
                        beta.2[k]<-mt$coefficients[2]
                        r2.c<-1 - mc$deviance / mc$null.deviance
                        r2.t<-1 - mt$deviance / mt$null.deviance
			file.out<-as.character(args[7])#file to write out
			cor.2<-cor(ALL[,col],ALL[,tag.snp],use = "complete.obs")^2
			write.table(as.data.frame(cbind(snps$SNP,beta,betac.1[k],betac.2[k],allele.fr1,allele.fr2,tag.snp,allele.fr1.t,allele.fr2.t,
				beta.1[k],beta.2[k],flag,"CT",cor.1,cor.2,"EUR",K,Z.scores[col],Z.scores[tag.snp],AF.cases.c,AF.controls.c,AF.cases.t,AF.controls.t,r2.c,r2.t)),file=args[7],col.names=F,row.names=F,quote=F,append=TRUE)

       	 	}
		
	


  
    
	}else{
	  if(Z.scores[col]>z.tres){	
	  
	    betac.1[k]<-glm(ALL$PHENOTYPE~ALL[,col],family=binomial)$coefficients[2]

	    CASES<-get.cases2(dat2,G2,beta,beta02,P.G2,P.A2,as.numeric(n))
	    CONTROLS<-get.controls2(dat2,G2,beta,beta02,P.G2,P.A2,as.numeric(n))
		  ALL<-rbind(CASES,CONTROLS)
		  allele.fr2<-mean(ALL[,col],na.rm=T)/2
		  mc<-glm(ALL$PHENOTYPE~ALL[,col],family=binomial)
                  betac.2[k]<-mc$coefficients[2]
                  r2.c<-1 - mc$deviance / mc$null.deviance
                  r2.t<-NA
		  file.out<-as.character(args[7])#file to write out
		  print(paste("al escribir caso 2",snps$SNP))
		  write.table(as.data.frame(cbind(snps$SNP,beta,betac.1[k],betac.2[k],allele.fr1,
		                                  allele.fr2,NA,NA,NA,NA,NA,flag,"C",1,1,"AFR",K,Z.scores[col],NA,AF.cases.c,AF.controls.c,AF.cases.t,AF.controls.t,r2.c,r2.t)),
		              file=args[7],col.names=F,row.names=F,quote=F,append=TRUE)
		  
	    CASES<-get.cases2(dat3,G3,beta,beta03,P.G3,P.A3,as.numeric(n))
	    CONTROLS<-get.controls2(dat3,G3,beta,beta03,P.G3,P.A3,as.numeric(n))
		  ALL<-rbind(CASES,CONTROLS)
		  allele.fr2<-mean(ALL[,col],na.rm=T)/2
		  mc<-glm(ALL$PHENOTYPE~ALL[,col],family=binomial)
                  betac.2[k]<-mc$coefficients[2]
                  r2.c<-1 - mc$deviance / mc$null.deviance
                  r2.t<-NA
		  file.out<-as.character(args[7])#file to write out
		  write.table(as.data.frame(cbind(snps$SNP,beta,betac.1[k],betac.2[k],allele.fr1,
		                                  allele.fr2,NA,NA,NA,NA,NA,flag,"C",1,1,"AMR",K,Z.scores[col],NA,AF.cases.c,AF.controls.c,AF.cases.t,AF.controls.t,r2.c,r2.t)),
		              file=args[7],col.names=F,row.names=F,quote=F,append=TRUE)

	    CASES<-get.cases2(dat4,G4,beta,beta04,P.G4,P.A4,as.numeric(n))
	    CONTROLS<-get.controls2(dat4,G4,beta,beta04,P.G4,P.A4,as.numeric(n))
		  ALL<-rbind(CASES,CONTROLS)
		  allele.fr2<-mean(ALL[,col],na.rm=T)/2
		  mc<-glm(ALL$PHENOTYPE~ALL[,col],family=binomial)
                  betac.2[k]<-mc$coefficients[2]
                  r2.c<-1 - mc$deviance / mc$null.deviance
                  r2.t<-NA
		  file.out<-as.character(args[7])#file to write out
		  write.table(as.data.frame(cbind(snps$SNP,beta,betac.1[k],betac.2[k],allele.fr1,
		                                  allele.fr2,NA,NA,NA,NA,NA,flag,"C",1,1,"EAS",K,Z.scores[col],NA,AF.cases.c,AF.controls.c,AF.cases.t,AF.controls.t,r2.c,r2.t)),
		              file=args[7],col.names=F,row.names=F,quote=F,append=TRUE)
		  
	    CASES<-get.cases2(dat5,G5,beta,beta05,P.G5,P.A5,as.numeric(n))
	    CONTROLS<-get.controls2(dat5,G5,beta,beta05,P.G5,P.A5,as.numeric(n))
		  ALL<-rbind(CASES,CONTROLS)
		  allele.fr2<-mean(ALL[,col],na.rm=T)/2
		  mc<-glm(ALL$PHENOTYPE~ALL[,col],family=binomial)
                  betac.2[k]<-mc$coefficients[2]
                  r2.c<-1 - mc$deviance / mc$null.deviance
                  r2.t<-NA
		  file.out<-as.character(args[7])#file to write out
		  write.table(as.data.frame(cbind(snps$SNP,beta,betac.1[k],betac.2[k],allele.fr1,
		                                  allele.fr2,NA,NA,NA,NA,NA,flag,"C",1,1,"SAS",K,Z.scores[col],NA,AF.cases.c,AF.controls.c,AF.cases.t,AF.controls.t,r2.c,r2.t)),
		              file=args[7],col.names=F,row.names=F,quote=F,append=TRUE)
		  
	    CASES<-get.cases2(dat1,G1,beta,beta01,P.G1,P.A1,as.numeric(n))
	    CONTROLS<-get.controls2(dat1,G1,beta,beta01,P.G1,P.A1,as.numeric(n))
		  ALL<-rbind(CASES,CONTROLS)
		  allele.fr2<-mean(ALL[,col],na.rm=T)/2
		  mc<-glm(ALL$PHENOTYPE~ALL[,col],family=binomial)
                  betac.2[k]<-mc$coefficients[2]
                  r2.c<-1 - mc$deviance / mc$null.deviance
                  r2.t<-NA
		  file.out<-as.character(args[7])#file to write out
		  write.table(as.data.frame(cbind(snps$SNP,beta,betac.1[k],betac.2[k],allele.fr1,
		                                  allele.fr2,NA,NA,NA,NA,NA,flag,"C",1,1,"EUR",K,Z.scores[col],NA,AF.cases.c,AF.controls.c,AF.cases.t,AF.controls.t,r2.c,r2.t)),
		              file=args[7],col.names=F,row.names=F,quote=F,append=TRUE)

	  }
	
	}
	
	
  Tb<-Sys.time()-Time0 
}
