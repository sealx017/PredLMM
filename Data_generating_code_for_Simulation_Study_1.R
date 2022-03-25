## Data generating procedure for Simulation Study 1 ##

### The technical details about the simulation strategy are provided in the main text.
### Below, we simulate 100 datasets having 8000 individuals and 13000 genes. The true heritability value is 0.2.

library(data.table)
library(reshape2)

for(step in c(1:100)){

set.seed(step)
  
T=80 # the following seven parameters control the Balding-Nichols model
t0=20
t1=50
t2=60
t3=70
t4=80
popsize=15000

nmarkers=13000 # the number of SNPs
npeople=2000 # there are 4 groups of individuals each of size npeople



freq1A=rep(0,nmarkers)
freq2A=rep(0,nmarkers)
freq3A=rep(0,nmarkers)
freq4A=rep(0,nmarkers)
marker1A=matrix(0,nrow=npeople,ncol=nmarkers)
marker2A=matrix(0,nrow=npeople,ncol=nmarkers)
marker3A=matrix(0,nrow=npeople,ncol=nmarkers)
marker4A=matrix(0,nrow=npeople,ncol=nmarkers)
Pop1=matrix(0,nrow=npeople,ncol=nmarkers)
Pop2=matrix(0,nrow=npeople,ncol=nmarkers)
Pop3=matrix(0,nrow=npeople,ncol=nmarkers)
Pop4=matrix(0,nrow=npeople,ncol=nmarkers)
p0s<-c()
p1s<-c()
p2s<-c()
p3s<-c()
p4s<-c()
p5s<-c()
p6s<-c()
for(i in 1:nmarkers){ # allele frequencies of SNPs for four different groups of individuals
  p0=rnorm(1,0.5,.0004)
  p1=rnorm(1,p0, sqrt(p0*(1-p0)*T/(2*popsize)))
  p3=rnorm(1,p1,sqrt(p1*(1-p1)*t1/(2*popsize)))
  p4=rnorm(1,p1,sqrt(p1*(1-p1)*t1/(2*popsize)))
  p2=rnorm(1,p0, sqrt(p0*(1-p0)*T/(2*popsize)))
  p5=rnorm(1,p2,sqrt(p2*(1-p2)*t2/(2*popsize)))
  p6=rnorm(1,p2,sqrt(p2*(1-p2)*t2/(2*popsize)))
  p0s<-c(p0s,p0)
  p1s<-c(p1s,p1)
  p2s<-c(p2s,p2)
  p3s<-c(p3s,p3)
  p4s<-c(p4s,p4)
  p5s<-c(p5s,p5)
  p6s<-c(p6s,p6)
  case1<-c(p3,p4,p5,p6)
  freq=case1
  
  freq1A[i]=freq[1]
  freq2A[i]=freq[2]
  freq3A[i]=freq[3]
  freq4A[i]=freq[4]
}
if(all(c(freq1A,freq2A,freq3A,freq4A)>0 & c(freq1A,freq2A,freq3A,freq4A)<1)){ # simulating raw and standardized SNP data
  for(i in 1:length(freq1A)){marker1A[,i]=rbinom(npeople,2,freq1A[i]);Pop1[,i]= (marker1A[,i] -2*freq1A[i])/(sqrt(2*freq1A[i]*(1-freq1A[i])))}
  for(i in 1:length(freq2A)){marker2A[,i]=rbinom(npeople,2,freq2A[i]);Pop2[,i]= (marker2A[,i] -2*freq2A[i])/(sqrt(2*freq2A[i]*(1-freq2A[i])))}
  for(i in 1:length(freq3A)){marker3A[,i]=rbinom(npeople,2,freq3A[i]);Pop3[,i]= (marker3A[,i] -2*freq3A[i])/(sqrt(2*freq3A[i]*(1-freq3A[i])))}
  for(i in 1:length(freq4A)){marker4A[,i]=rbinom(npeople,2,freq4A[i]);Pop4[,i]= (marker4A[,i] -2*freq4A[i])/(sqrt(2*freq4A[i]*(1-freq4A[i])))}
  
  
  ncausalsnps = 300 # the number of SNPs which have causal effect
  samplecausalsnps=sample(1:nmarkers,ncausalsnps)
  heritability = 0.2 # the true value of the heritability
  
  beta=rnorm(length(samplecausalsnps),0,sqrt(heritability/(length(samplecausalsnps)))) # fixed effect-size based on the true heritability
  
  # simulating phenotype data based on the causal SNPs
  CausalGeno1=Pop1[,samplecausalsnps]
  CausalFreq1=freq1A[samplecausalsnps] 
  phenotype1=CausalGeno1%*%beta + rnorm(npeople,0,sqrt(1-heritability))
  CausalGeno2=Pop2[,samplecausalsnps]
  CausalFreq2=freq2A[samplecausalsnps]
  phenotype2=CausalGeno2%*%beta + rnorm(npeople,0,sqrt(1-heritability))
  CausalGeno3=Pop3[,samplecausalsnps]
  CausalFreq3=freq3A[samplecausalsnps]
  phenotype3=CausalGeno3%*%beta + rnorm(npeople,0,sqrt(1-heritability))
  CausalGeno4=Pop4[,samplecausalsnps]
  CausalFreq4=freq4A[samplecausalsnps]
  phenotype4=CausalGeno4%*%beta + rnorm(npeople,0,sqrt(1-heritability))
  
  
  phenotype=c(phenotype1,phenotype2,phenotype3,phenotype4)
  genotype=rbind(marker1A,marker2A,marker3A,marker4A)
  require(genetics)
  popCom = genotype
  l = 1
  j = 1
  newdat = NULL
  assigner<-function(j){
    while(j <= ncol(popCom)){
      vec =  as.character(unlist(strsplit(as.character(as.genotype.allele.count(popCom[,j], alleles=c("C","T")) ),"[/]")))
      newdat <<-cbind(newdat,matrix(vec,nrow=nrow(popCom),ncol=2,byrow=T))
      print(j)
      j = j+1
      l = l+2
    }}
  assigner(1)
  Genotype_data = newdat
  
  # create ped and map file for plink
  npeople1=npeople*4
  IndivID=c(paste(1:npeople1,sep=""))
  FamID=c(paste(1:npeople1,sep=""))
  FatherID=rep(0,npeople1)
  MotherID=rep(0,npeople1)
  Sex=c(rep(1,npeople1/2),rep(2,npeople1/2))
  plink.ped=data.frame(FamID,IndivID, FatherID,MotherID, Sex,phenotype, Genotype_data)
  plink.fam=data.frame(FamID,IndivID, FatherID,MotherID, Sex,phenotype)
  
  snpnames=c(paste("snp",1:ncol(genotype),sep=""))
  chrom=rep(1,ncol(genotype))
  bpdis=sort.default(sample(1:902606,ncol(genotype)))
  plink.map=data.frame(chrom,snpnames,bpdis)
  
  tgeno = t(as.matrix(genotype))
  
  # the folder-names where ped, binary, phenotype and GRM files will be saved
  
  pedfiles_location = "/Data/pedfiles/"
  binaryfiles_location = "/Data/binaryfiles/"
  grmfiles_location = "/Data/GRM/"
  phenfiles_location = "/Data/phen_data/"
  filename = "Heritability_0_2_8k_indiv_13ksnps" # the first part of the name of every file, followed by the dataset's number and extension
  
  write.table(plink.fam[,c(1,2,6)],paste0(pedfiles_location,filename,step,'.pheno'),quote=F,col.names=F,row.names=F) # writing plink files
  write.table(plink.fam,paste0(pedfiles_location,filename,step,'.fam'),quote=F,col.names = F,row.names = F)
  write.table(plink.map,paste0(pedfiles_location,filename,step,'.map'),quote=F,col.names = F,row.names = F)
  write.table(plink.ped,paste0(pedfiles_location,filename,step,'.ped'),quote=F,col.names = F,row.names = F)
  
  phen = data.frame(FamID,IndivID,phenotype)
  colnames(phen) = c("FID","IID","phen")
  write.table(phen,paste0(phenfiles_location,filename,step,'.csv'),quote=F,row.names = F) # writing phenotype file
  
  # next plink and GCTA softwares are used to respectively construct the binary and the GRM files
  
  plink_location = "/plink/plink_linux_x86_64/plink" # location of the plink software
  gcta_location = "/gcta_1.91.4beta/gcta64" # location of the GCTA software
  system(paste(plink_location," -file ",pedfiles_location,filename,step," --out ",binaryfiles_location,filename,step," --make-bed",sep=""))
  system(paste(gcta_location," --bfile ",binaryfiles_location,filename,step," --make-grm --out ",grmfiles_location,filename,step," --thread-num 10",sep=""))

}
}