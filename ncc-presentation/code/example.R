library(survival)
# all the R functions needed
source("ncc-presentation/code/FUN.R"); 

# two data needed
data0 <- read.csv("ncc-presentation/code/data0.csv",header=T)
NCC.dat <- read.csv("ncc-presentation/code/NCC.dat.csv",header=T)

# arguments needed to specify in execute.fun
m0 <- c(3,1) # the number of controls selected in Phase II and III subcohort
a0 <- c(0,1) # the constraints for the matching variable using in the nested case control study

## the name of the covariates you will use in different models
cov.list = list(Z=list(nm=c("cov1"),ncc="S0"),
                B=list(nm=c("cov1","cov2"),ncc="S1"),
                G=list(nm=c("cov1","cov3"),ncc="S2"),
                "Full"=list(nm=c("cov1","cov2","cov3"), ncc="S2")) 
t0 = 5 # the pre-determined prediction time
u0 = 0.1 # the cut-off value for false positive rate to calculate true postive rate, postive predictive values, negative predictive values

## accomodate marker-dependent censoring which requires the covarites that the distribution of censoring time might depend on
yes.DepCen=T
cov.DepCen=list(Cov=c("cov1"),Stage = "S0") 

## specify which model you want to use: GLM or Cox PH model. If you want to use Cox PH model, you can specify: model = list("type"="Cox","family"=NULL)
model = list("type"="GLM","family"=binomial(link=cloglog)) 

result = execute.fun(data0=data0,
                     NCC.dat=NCC.dat,
                     cov.list=cov.list,
                     model=model,
                     u0=u0,
                     t0=t0,
                     execute.rtn="ALL",
                     yes.DepCen,
                     cov.DepCen,
                     control=list(yes.cv=T,yes.IncV=T))
