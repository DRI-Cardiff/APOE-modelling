#Valentina Escott-Price
#Dementia Research Institute (DRI) at Cardiff, UK
#02.12.2020


library(SDMTools)
library(ROCR)
library(pROC)
library(ResourceSelection)
library(bdpv)


Nc <- 10000  # number of cases
Nn <- 10000  # number of controls
Nc_o<-Nc
status <- c(rep(1, Nc), rep(0, Nn))  # case/non-case status
age<- c(rnorm(Nc,mean=79,sd=8), rnorm(Nn,mean=54.3, sd=14.5))
gwas<-read.table(file="GWAS_SNPs.txt", header=T, sep=" ")

# model APOE
APOE <- array(dim=Nc+Nn)
genAPOE <- matrix(nrow=Nc+Nn, ncol=2)
lorAPOE<-log(gwas$OR[37])

#matched population controls
#e4
f0<-0.14;f1<-0.351;o1<-f1*(1-f0)/f0/(1-f1); o1
genAPOE[1:Nc, 1] <- rbinom(Nc, 1, f1)+rbinom(Nc, 1, f1);
genAPOE[1:Nc, 1] <- rbinom(Nc, 1, f1)+rbinom(Nc, 1, f1);
genAPOE[(Nc+1):(Nc+Nn), 1] <- rbinom(Nn, 1, f0)+rbinom(Nn, 1, f0)
#e2
f0<-0.084;f1<-0.039;o2<-f1*(1-f0)/f0/(1-f1);o2
genAPOE[1:Nc, 2] <- rbinom(Nc, 1, f1)+rbinom(Nc, 1, f1);
genAPOE[(Nc+1):(Nc+Nn), 2] <- rbinom(Nn, 1, f0)+rbinom(Nn, 1, f0)
#e33
a<-which((genAPOE[1:Nc, 1]+genAPOE[1:Nc, 2])==0); length(a)/Nc
a<-which((genAPOE[(Nc+1):(Nc+Nn),1]+genAPOE[(Nc+1):(Nc+Nn),2])==0); length(a)/Nn

#young controls (54-)
#e4
a<-which(genAPOE[(Nc+1):(Nc+Nn), 1]==2); nc1<-round(length(a)*0.91) #91% of e4e4 are cases at 68 years old
a<-which(genAPOE[(Nc+1):(Nc+Nn), 1]==1); nc2<-round(length(a)*0.47) #47% of het e4 are cases at 76 years old
a<-which(genAPOE[(Nc+1):(Nc+Nn), 1]==0); nc3<-round(length(a)*0.2) #20% of non-e4 are cases at 84 years old
Nc<-Nc+nc1+nc2+nc3; Nc
Nn<-Nn-nc1-nc2-nc3; Nn
status1 <- c(rep(1, Nc), rep(0, Nn))  # "TRUE" case/control status
Nc_o-Nn
  
# model APOE
APOE <- array(dim=Nc+Nn)
genAPOE <- matrix(nrow=Nc+Nn, ncol=2)

#SIMULATE ACTUAL cases and controls
#e4
f1<-0.35
genAPOE[1:Nc_o, 1] <- rbinom(Nc_o, 1, f1)+rbinom(Nc_o, 1, f1);
f1<-0.355
genAPOE[(1+Nc_o):Nc, 1] <- rbinom((Nc-Nc_o), 1, f1)+rbinom((Nc-Nc_o), 1, f1);
f0<-0.12
genAPOE[(Nc+1):(Nc+Nn), 1] <- rbinom(Nn, 1, f0)+rbinom(Nn, 1, f0)
for (i in (1:(Nc+Nn))) APOE[i] <- lorAPOE*genAPOE[i,1] # polygenic risk score

gwas<-gwas[-37,] #exclude APOE SNP
M<-nrow(gwas)
OdR<-gwas$OR
f0<-gwas$AltAlleleFreq
f1<-OdR*f0/(1-f0*(1-OdR))
lor <- log(f1*(1-f0)/f0/(1-f1))  # log odds ratio

# create random genotypes based on allele frequencies
gen <- matrix(nrow=Nc+Nn, ncol=M)
for (m in (1:M)) 
{
  gen[1:Nc, m] <- rbinom(Nc, 1, f1[m])+rbinom(Nc, 1, f1[m]);
  gen[(Nc+1):(Nc+Nn), m] <- rbinom(Nn, 1, f0[m])+rbinom(Nn, 1, f0[m])
}

GRS <- array(dim=Nc+Nn)
for (i in (1:(Nc+Nn))) GRS[i] <- sum(lor*gen[i,])  # polygenic risk score


gwas<-read.table(file="GWAS_SNPs.txt", header=T, sep=" ")
M<-nrow(gwas)
OdR<-gwas$OR
f0<-gwas$AltAlleleFreq
f1<-OdR*f0/(1-f0*(1-OdR))
lor <- log(f1*(1-f0)/f0/(1-f1))  # log odds ratio

gen1 <- matrix(nrow=Nc+Nn, ncol=M)
gen1[,1:36]<-gen[,1:36]
gen1[,37]<-genAPOE[,1]
gen1[,38:40]<-gen[,37:39]

GRSwithAPOE <- array(dim=Nc+Nn)
for (i in (1:(Nc+Nn))) GRSwithAPOE[i] <- sum(lor*gen1[i,])  # polygenic risk score

PRS<-GRS
PRSwithAPOE<-GRSwithAPOE
#Simulation of additional PRS
Nmar<-  c(500,  1000, 2500,  2000,  2000, 2000)
Bmean_o  <- c(0.005, 0.0025, 0.002, 0.001, 0.0005, 0.00025)
Bmean    <- c(0.0025,0.001,0.001, 0.0005, 0.00025, 0.0001)
Bsdev <-c(0.02,  0.015, 0.01,   0.01,   0.01,   0.01)
RAN   <-c(0.99,   0.95,   0.9,  0.85, 0.8, 0.75)
mafL <-0.01
mafU <- 0.45
sum(Nmar)

for (SC in 1:length(Nmar))
{
  KK<-floor(Nmar[SC]*RAN[SC])
  KKminorAllele<-floor(KK*0.7)

  lor_o<-rnorm(KKminorAllele, Bmean_o[SC], Bsdev[SC])
  lor<-rnorm(KKminorAllele, Bmean[SC], Bsdev[SC])
  f0<- runif(KKminorAllele, mafL, mafU)
  f1_o<-f0*exp(lor_o)/(1+f0*(exp(lor_o)-1))
  f1<-f0*exp(lor)/(1+f0*(exp(lor)-1))
  gen <- matrix(nrow=Nc+Nn, ncol=KKminorAllele)
  for (m in (1:KKminorAllele)) {
    gen[1:Nc_o, m] <- rbinom(Nc_o, 1, f1_o[m])+rbinom(Nc_o, 1, f1_o[m]);
    gen[(1+Nc_o):Nc, m] <- rbinom((Nc-Nc_o), 1, f1[m])+rbinom((Nc-Nc_o), 1, f1[m]);
    gen[(Nc+1):(Nc+Nn), m] <- rbinom(Nn, 1, f0[m])+rbinom(Nn, 1, f0[m])
  }
  for (i in (1:(Nc+Nn))) PRS[i] <- PRS[i]+sum(lor_o*gen[i,])
  for (i in (1:(Nc+Nn))) PRSwithAPOE[i] <- PRSwithAPOE[i]+sum(lor*gen[i,])
  
  KKmajorAllele<-KK-KKminorAllele
  lor<-rnorm(KKmajorAllele, Bmean[SC], Bsdev[SC])
  lor_o<-rnorm(KKmajorAllele, Bmean_o[SC], Bsdev[SC])
  f0<- runif(KKmajorAllele, (1-mafU), (1-mafL))
  f1_o<-f0*exp(lor_o)/(1+f0*(exp(lor_o)-1))
  f1<-f0*exp(lor)/(1+f0*(exp(lor)-1))
  gen <- matrix(nrow=Nc+Nn, ncol=KKmajorAllele)
  for (m in (1:KKmajorAllele)) {
    gen[1:Nc_o, m] <- rbinom(Nc_o, 1, f1_o[m])+rbinom(Nc_o, 1, f1_o[m]);
    gen[(1+Nc_o):Nc, m] <- rbinom((Nc-Nc_o), 1, f1[m])+rbinom((Nc-Nc_o), 1, f1[m]);
    gen[(Nc+1):(Nc+Nn), m] <- rbinom(Nn, 1, f0[m])+rbinom(Nn, 1, f0[m])
  }
  for (i in (1:(Nc+Nn))) PRS[i] <- PRS[i]+sum(lor_o*gen[i,])  # polygenic risk score
  for (i in (1:(Nc+Nn))) PRSwithAPOE[i] <- PRSwithAPOE[i]+sum(lor*gen[i,])  # polygenic risk score
  
  #RANDOM
  KK<-Nmar[SC]-KK
  lor<-rnorm(KK, Bmean[SC], Bsdev[SC])
  f1<-runif(KK, mafL, mafU)
  f0<-f1
  gen <- matrix(nrow=Nc+Nn, ncol=KK)
  for (m in (1:KK)) {
    gen[1:Nc, m] <- rbinom(Nc, 1, f1[m])+rbinom(Nc, 1, f1[m]);
    gen[(Nc+1):(Nc+Nn), m] <- rbinom(Nn, 1, f0[m])+rbinom(Nn, 1, f0[m])
  }
  for (i in (1:(Nc+Nn))) PRS[i] <- PRS[i]+sum(lor*gen[i,])
  for (i in (1:(Nc+Nn))) PRSwithAPOE[i] <- PRSwithAPOE[i]+sum(lor*gen[i,])
}#SCENARIO

res<-matrix(0.5,ncol=2,nrow=6)
resR2<-matrix(0,ncol=2,nrow=6)
N<-Nc+Nn

model<-glm(status1~APOE, family=binomial);res[1,1]<-auc(model$y, model$fitted.values);
l0<-model$null.deviance; df0<-model$df.null; l1 <- model$deviance;df1 <- model$df.residual
resR2[1,1]<-(1 - exp((l1 - l0)/N))/(1 - exp( - l0/N))	    #Nagelkerke

model<-glm(status1~GRS, family=binomial); res[2,1]<-auc(model$y, model$fitted.values);
l0<-model$null.deviance; df0<-model$df.null; l1 <- model$deviance;df1 <- model$df.residual
resR2[2,1]<-(1 - exp((l1 - l0)/N))/(1 - exp( - l0/N))

model<-glm(status1~GRSwithAPOE, family=binomial); res[3,1]<-auc(model$y, model$fitted.values);
l0<-model$null.deviance; df0<-model$df.null; l1 <- model$deviance;df1 <- model$df.residual
resR2[3,1]<-(1 - exp((l1 - l0)/N))/(1 - exp( - l0/N))

model<-glm(status1~PRS, family=binomial); res[4,1]<-auc(model$y, model$fitted.values);
l0<-model$null.deviance; df0<-model$df.null; l1 <- model$deviance;df1 <- model$df.residual
resR2[4,1]<-(1 - exp((l1 - l0)/N))/(1 - exp( - l0/N))

model<-glm(status1~PRSwithAPOE, family=binomial); res[5,1]<-auc(model$y, model$fitted.values);
l0<-model$null.deviance; df0<-model$df.null; l1 <- model$deviance;df1 <- model$df.residual
resR2[5,1]<-(1 - exp((l1 - l0)/N))/(1 - exp( - l0/N))

model<-glm(status1~APOE+PRS, family=binomial); res[6,1]<-auc(model$y, model$fitted.values);
l0<-model$null.deviance; df0<-model$df.null; l1 <- model$deviance;df1 <- model$df.residual
resR2[6,1]<-(1 - exp((l1 - l0)/N))/(1 - exp( - l0/N))



model<-glm(status~APOE, family=binomial);res[1,2]<-auc(model$y, model$fitted.values);
l0<-model$null.deviance; df0<-model$df.null; l1 <- model$deviance;df1 <- model$df.residual
resR2[1,2]<-(1 - exp((l1 - l0)/N))/(1 - exp( - l0/N))

model<-glm(status~GRS, family=binomial); res[2,2]<-auc(model$y, model$fitted.values);
l0<-model$null.deviance; df0<-model$df.null; l1 <- model$deviance;df1 <- model$df.residual
resR2[2,2]<-(1 - exp((l1 - l0)/N))/(1 - exp( - l0/N))

model<-glm(status~GRSwithAPOE, family=binomial); res[3,2]<-auc(model$y, model$fitted.values);
l0<-model$null.deviance; df0<-model$df.null; l1 <- model$deviance;df1 <- model$df.residual
resR2[3,2]<-(1 - exp((l1 - l0)/N))/(1 - exp( - l0/N))

model<-glm(status~PRS, family=binomial); res[4,2]<-auc(model$y, model$fitted.values);
l0<-model$null.deviance; df0<-model$df.null; l1 <- model$deviance;df1 <- model$df.residual
resR2[4,2]<-(1 - exp((l1 - l0)/N))/(1 - exp( - l0/N))

model<-glm(status~PRSwithAPOE, family=binomial); res[5,2]<-auc(model$y, model$fitted.values);
l0<-model$null.deviance; df0<-model$df.null; l1 <- model$deviance;df1 <- model$df.residual
resR2[5,2]<-(1 - exp((l1 - l0)/N))/(1 - exp( - l0/N))

model<-glm(status~APOE+PRS, family=binomial); res[6,2]<-auc(model$y, model$fitted.values);
l0<-model$null.deviance; df0<-model$df.null; l1 <- model$deviance;df1 <- model$df.residual
resR2[6,2]<-(1 - exp((l1 - l0)/N))/(1 - exp( - l0/N))
res

barplot(t(resR2[-2,]), beside = TRUE, ylim=c(0,0.33), main = "Variance Explained (R2)",
        xlab = "", ylab="R2", col = c("grey", "black"),
        names.arg=c("APOE", "ORS", "PRS.no.APOE", "PRS", "APOE+PRS.no.APOE"),
        cex.names=0.8)
abline(h=max(resR2), lty=2)
legend("topleft",bty="n", c("True case/control","True case, Population control"),
       fill = c("grey", "black"))

barplot(t(res[-2,]), beside = TRUE, ylim=c(0,1), main = "Prediction accuracy",
        xlab = "", ylab="AUC", col = c("grey", "black"),
        names.arg=c("APOE", "ORS", "PRS.no.APOE", "PRS", "APOE+PRS.no.APOE"),
        cex.names=0.8)
abline(h=max(res), lty=2)
legend("topleft",bty="n", c("True case/control","True case, Population control"),
       fill = c("grey", "black"))

