execute.fun <- function(data0,
                        NCC.dat,
                        cov.list,model,
                        u0,
                        t0,
                        execute.rtn="ALL",
                        yes.DepCen=F,
                        cov.DepCen=list(Cov=NULL,Stage=NULL),
                        control=list(yes.cv=F,yes.IncV=F)){
  
  n.S = length(cov.list) # number of the stages
  model.type = model$type
  model.family = model$family
  yes.cv = control$yes.cv
  yes.IncV = control$yes.IncV
  
  id = as.character(data0$id); nn = length(id) # id
  xi = as.numeric(data0$time) # event time which might be censored
  di = as.numeric(data0$delta) # indicator of censoring: = 0, censored
  cov.mat = data0[,-match(c("id","time","delta"),colnames(data0)),drop=F]
  
  col.nm = colnames(NCC.dat)
  ## indicate whether sampled to NCC cohort as a case
  vi.case = as.numeric(NCC.dat$case) 
  ## indicate whether sampled to NCC as a case or a control 
  vi.S = as.matrix(NCC.dat[,which(substr(col.nm,1,1)=="S"),drop=F]) 
  
  if (sum(substr(col.nm,1,1)=="M")>0)
      Mi = as.matrix(NCC.dat[,which(substr(col.nm,1,1)=="M"),drop=F]) else
          Mi = NULL # Matching variables
  
  ## === calculate the NCC weighting === ##
  WGT.NCC.out = WGT.NCC.FUN(data=cbind(id,xi,di,vi.case,vi.S),
                            M.dat=Mi,
                            a.match=a0,
                            m.match=m0)
  WGT.NCC.out$S0 = list("wgt"=rep(1,nn),
                        "I0kj"=NULL,
                        "m"=Inf)
  WGT.NCC.out$S1$wgt = WGT.NCC.out$S1$wgt
  WGT.NCC.out$S2$wgt = WGT.NCC.out$S2$wgt
  
  ## === calculate the censoring weighting === ##
  if (yes.DepCen) {
    Z.DepCen = cov.mat[,match(cov.DepCen$Cov,colnames(cov.mat))] 
    wgt = WGT.NCC.out[[match(cov.DepCen$Stage,names(WGT.NCC.out))]]$wgt
  } else {
    Z.DepCen=NULL
    wgt = NULL
  }
  WGT.CEN.out = WGT.CEN.FUN(data=cbind(xi,di,wgt,Z.DepCen),
                            t0=t0,yes.exp=T,
                            yes.DepCen=yes.DepCen)
  wgt.cen = WGT.CEN.out$wgt
  LamCexp = WGT.CEN.out$LamC.exp
  
  output.acc = output.beta = output.cv = U.acc = as.list(1:n.S)
  names(output.acc) = names(output.beta)= names(output.cv)= names(U.acc) = names(cov.list)
  for (s in 1:n.S){
    nm.s = names(cov.list)[[s]]
    cov.s = cov.list[[s]]
    var.nm.s = cov.s$nm
    ncc.nm.s = cov.s$ncc
    wgt.ncc.out = WGT.NCC.out[[match(ncc.nm.s,names(WGT.NCC.out))]]
    wgt.ncc.out$tj = WGT.NCC.out$tj
    wgt.ncc.out$pi.tj = WGT.NCC.out$pi.tj
    wgt.ncc = wgt.ncc.out$wgt
    vi = 1*(wgt.ncc!=0); 
    cat("Model",nm.s,sum(vi),"\n")
    
    ## === markers === ##
    var = cov.mat[,match(var.nm.s,colnames(cov.mat)),drop=F]
    mydata.ncc = cbind(id,vi.case,wgt.ncc,wgt.cen,xi,di,var)[vi==1,]
    var.ncc = mydata.ncc[,-c(1:6),drop=F]
    if (model.type=="GLM") var.ncc = cbind(1,var.ncc)
    if (yes.cv==T){
      cv.out = ROC.DIPW.CV.632(data.ncc=mydata.ncc,u0=u0,t0=t0,model=model)
      acc.cv = apply(cv.out$ACC.u0,1,mean,trim=0.1,na.rm=T)
      roc.cv = cbind(cv.out$roc[,1],apply(cv.out$roc[,-1],1,mean,trim=0.1,na.rm=T))
      output.cv[[s]] = list("acc"=acc.cv,"roc"=roc.cv)
    } else output.cv[[s]]= NA
    
    betafit = try(NCC.Fit.FUN(data.ncc=mydata.ncc,
                              model=model.type,
                              model.family=model.family,
                              t0=t0,
                              rtn="ALL"),
                  silent=T)
    if (any(substr(betafit,1,5)!="Error")){
      betahat = betafit[,1]; 
      score = as.matrix(var.ncc)%*%betahat
      mydata.ncc = cbind(mydata.ncc[,c(1:6)], score, var.ncc)
      out = try(ROC.DIPW.FUN(data.ncc=mydata.ncc,
                             u0=u0, t0=t0,
                             rtn=execute.rtn,
                             var.control=list(wgt.ncc.dat=wgt.ncc.out,
                                              betafit=betafit,
                                              data0=cbind(di,vi),
                                              model=model.type,
                                              family=model.family,
                                              LamCexp=LamCexp)),
                silent=T)
      if (any(substr(out,1,5)!="Error")) yes.error=F else yes.error =T
    } else {yes.error=T}
    
    if (yes.error==F){				
      if (execute.rtn=="ALL"){ 
        beta.out = cbind(betahat,
                         sqrt(out$var$VAR.Beta[,1]),
                         sqrt(out$var$VAR.Beta[,2]))
        colnames(beta.out) = c("est","se-naive","se-adj")
        acc.out = cbind(out$est$ACC.u0,
                        sqrt(out$var$VAR.ACC[,1]),
                        sqrt(out$var$VAR.ACC[,2]))
        colnames(acc.out) = c("est","se-naive","se-adj")
        rownames(acc.out) = c("AUC",paste("TPR",u0,sep="-"),
                              paste("NPV",u0,sep="-"),
                              paste("PPV",u0,sep="-"))
        output.beta[[s]] = beta.out
        output.acc[[s]] = acc.out
        U.acc[[s]] = out$var$Ui.ACC
      } else{
        output.beta[[s]]=betahat
        output.acc[[s]]=out$ACC.u0
        U.acc[[s]]=NA
      } 
    } else {
      output.beta[[s]] = output.acc[[s]] = U.acc[[s]] = NA
    } # end if (yes.error==F)	
    
  } # end for s
  
  if (yes.IncV==T){
    output.IncV = as.list(1:((n.S-1)*n.S/2))
    nm.IncV = NULL
    index = 0
    for (s1 in 1:(n.S-1)){
      nm.s1 = names(cov.list)[s1]; U.s1 = U.acc[[s1]]
      if (yes.cv==T) acc.s1 = output.cv[[s1]]$acc else acc.s1 = output.acc[[s1]][,1]; 
      for (s2 in (s1+1):n.S){
        nm.s2 = names(cov.list)[s2];U.s2 = U.acc[[s2]]
        if (yes.cv==T) acc.s2 = output.cv[[s2]]$acc else acc.s2 = output.acc[[s2]][,1]
        cat("Model",paste(nm.s2,"-vs-",nm.s1,sep=""),"\n")
        nm.IncV = c(nm.IncV,paste(nm.s2,"-vs-",nm.s1,sep=""))
        index = index + 1
        
        IncV.est = acc.s2-acc.s1
        if (execute.rtn=="ALL"){					
          if (any(U.s1!="Error") & any(U.s2!="Error")) {
            var.out = NULL
            for (kk in 1:length(U.s1)){
              var.est = VarAdj.FUN(datA=U.s2[[kk]],datB=U.s1[[kk]])
              var.out = rbind(var.out,var.est)
            } # end for kk
            IncV.out = cbind(IncV.est,sqrt(var.out)); colnames(IncV.out) = c("est","se-naive","se-adj"); rownames(IncV.out) = c("AUC",paste("TPR",u0,sep="-"),paste("NPV",u0,sep="-"),paste("PPV",u0,sep="-"))
            output.IncV[[index]]=IncV.out	
          } else {output.IncV[[index]]=IncV.est}	
        } else {output.IncV[[index]]=IncV.est}
      } # end for s2
    } # end for s1
    names(output.IncV) = nm.IncV
    output=list("beta"=output.beta,"acc"=output.acc,"cv"=output.cv,"IncV"=output.IncV)	
  } else {
    output= list("beta"=output.beta,"acc"=output.acc,"cv"=output.cv)
  }
  return(output)	
}

sum.I <- function(yy,FUN,Yi,Vi=NULL)
{
  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
  # for each distinct ordered failure time t[j], number of Xi < t[j]
  pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')    
  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
  if (!is.null(Vi)) {
    ## if FUN contains '=', tmpind is the order of decending
    if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
    ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    return(rbind(0,Vi)[pos+1,])
  } else return(pos)
}

VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

BACC.fun <- function(xk,ck,type="FPR",u0){
  N = length(xk)
  scl = seq(min(ck),max(ck),length.out=1000)
  St0.Fcl = sum.I(scl,">=",ck,1*(xk==0))/N
  Ft0.Fcl = sum.I(scl,">=",ck,1*(xk==1))/N
  Fcl = sum.I(scl,">=",ck)/N
  
  St0 = max(St0.Fcl); Ft0 = 1-St0## St0 = P(T> t0); Ft0 = P(T<=t0)
  FPR.cl= (St0-St0.Fcl)/St0     ## P(Y> cl|T> t0)
  TPR.cl= (Ft0-Ft0.Fcl)/Ft0     ## P(Y> cl|T<=t0)
  NPV.cl= St0.Fcl/Fcl           ## P(T> t0|Y<=cl)
  PPV.cl= (Ft0-Ft0.Fcl)/(1-Fcl); if (sum(Fcl==1)>0) PPV.cl[Fcl==1] = NA  ## P(T<=t0|Y> cl)
  AUC = sum(TPR.cl*(FPR.cl-c(FPR.cl[-1],0)))	
  acc.cl = data.frame("cutoff"=scl,"FPR"=FPR.cl,"TPR"=TPR.cl,"NPV"=NPV.cl, "PPV"=PPV.cl)
  acc.ul = acc.cl; ind0 = match(type,names(acc.cl)); 
  ul = acc.ul[,ind0]; acc.ul = acc.ul[order(ul),]; ul = sort(ul); 
  if(ind0==1){
    indu = sum.I(u0,">=",ul)
  } else {
    indu = sum.I(u0,">",ul)
  } #%%%%
  c.u0 = acc.ul[indu,1]; 
  return(c("AUC"=AUC, unlist(acc.ul[indu,-c(1,ind0)])))
}

WGT.NCC.FUN <- function(data,M.dat=NULL,a.match=NULL,m.match,pi1=1,wp=1){
  
  ## ====================================================== ##
  ## Calculate the NCC sampling weights
  ## pi1 = prob. of the subject sampled as a case
  ## could incoporate multiple-stage sampling of controls
  ## ======================================================= ##
  N = nrow(data); n.S = length(m.match)
  data = data[,-1,drop=F]; mode(data) = "numeric"; xi = data[,1]; di = data[,2]; vi.case = data[,3]; vi.S = data[,-c(1:3),drop=F]; vi.S1 = vi.S[,1]
  ind.case = (1:N)[vi.case==1] ## all cases
  ti1 = xi[ind.case]; n1 = length(ti1)
  Ii1j = VTM(xi,n1)>=ti1 ## I.ij = I(Tj>=Ti)
  if (!is.null(M.dat)){
    p.M = ncol(M.dat)
    M1 = M.dat[ind.case,,drop=F]
    for (ll in 1:p.M){Ii1j = Ii1j*(abs(VTM(M.dat[,ll],n1)-M1[,ll])<=a.match[ll])} ## I.ij = I(Tj>=Ti,|Mj-Mi|<=a.match)
  }
  n.Ri1 = apply(Ii1j,1,sum); pi.i1 = n.Ri1/N
  # browser()
  ## stage 1
  vi.S1 = vi.S[,1]; m.S1 = m.match[1]
  I0k.i1.S1 = t(1*Ii1j[,vi.S1==1 & vi.case==0 & di==0])
  wk.S1 = rep(0,N); wk.S1[vi.case==1] = 1/pi1; 
  junk = t(I0k.i1.S1)*log(1-m.S1/pmax(n.Ri1*wp-1,m.S1)); junk[is.na(junk)]= 0; 
  wk.S1[vi.S1==1 & vi.case==0 & di==0]=1/(1-exp(apply(junk,2,sum,na.rm=T)))
  wk.S1.out = list("wgt"=wk.S1,"I0kj"=I0k.i1.S1,"m"=m.S1) # I0lj: n0 x n1
  wk.out = list("tj"=ti1,"pi.tj"=pi.i1,"S1"=wk.S1.out)
  if (n.S>1){ # stage 2
    vi.S2 = vi.S[,2]; m.S2 = m.match[2]
    I0k.i1.S2 = t(1*Ii1j[,vi.S2==1 & vi.case==0 & di==0])
    wk.S2 = rep(0,N); wk.S2[vi.case==1]=1/pi1
    junk = t(I0k.i1.S2)*log(1-m.S2/pmax(n.Ri1*wp-1,m.S2)); junk[is.na(junk)]=0
    wk.S2[vi.S2==1 & vi.case==0 & di==0] = 1/(1-exp(apply(junk,2,sum,na.rm=T)))
    wk.S2.out = list("wgt"=wk.S2,"I0kj"=I0k.i1.S2,"m"=m.S2)
    wk.out = list("tj"=ti1,"pi.tj"=pi.i1,"S1"=wk.S1.out,"S2"=wk.S2.out)
  }	 
  return(wk.out)
}


## This is changed by Michelle on July 17th, 2014
##    allow the censoring depending on markers
WGT.CEN.FUN <- function(data,t0,yes.exp=F,yes.DepCen=F){
  ## ================================================================ ##
  ## ============== KM Estimator of Censoring Survival ============== ##
  ## ================================================================ ##
  Ghat.FUN <- function(tt, data.fit, wgt=NULL, yes.DepCen){
    if (yes.DepCen){
      cox.fit = coxph(Surv(time,delta)~.,data=data.fit,weights=wgt)
      new.data = data.frame(matrix(0,nrow=1,ncol=(ncol(data.fit)-2))); colnames(new.data)=colnames(data.fit)[-c(1,2)]
      surv.fit = survfit(cox.fit,newdata=new.data,se.fit=F)
    } else{
      cox.fit = coxph(Surv(time,delta)~1,data=data.fit) 
      surv.fit = survfit(cox.fit,se.fit=F,type="fl")
    }     
    surv.ti = surv.fit$time; surv.fit.surv = surv.fit$surv
    surv.til = sort(surv.fit$time); surv.fit.surv = surv.fit.surv[order(surv.ti)]
    tmpind = sum.I(tt,">=",surv.til) + 1
    c(1,surv.fit.surv)[tmpind]
  }
  
  Ti <- data[,1]; Di <- data[,2]; N = length(Ti)
  
  if (yes.DepCen) {
    wgt <- data[,3]
    Zi = as.matrix(data[,-c(1:3),drop=F])
    data.fit <- data.frame(time=Ti,delta=1-Di,Zi) 
    beta = coxph(Surv(time,delta)~.,data=data.fit,weights=wgt)$coef
  }  else {
    wgt <- NULL
    data.fit <- data.frame(time=Ti,delta=1-Di)
  }
  
  S0.Ti = Ghat.FUN(Ti,data.fit,wgt,yes.DepCen); S0.t0 = Ghat.FUN(t0,data.fit,wgt,yes.DepCen)
  Wi <- rep(0,length(Ti))
  if (yes.DepCen){
    Ghat.Ti = S0.Ti^(exp(Zi%*%beta)); Ghat.t0 = S0.t0^(exp(Zi%*%beta))
    Wi[Ti <= t0] <- Di[Ti<=t0]/Ghat.Ti[Ti<=t0]  
    Wi[Ti >  t0] <- 1/Ghat.t0[Ti > t0]
  } else{
    Ghat.Ti = S0.Ti; Ghat.t0 = S0.t0
    Wi[Ti <= t0] <- Di[Ti<=t0]/S0.Ti[Ti<=t0]
    Wi[Ti >  t0] <- 1/S0.t0
  }
  Wi[is.na(Wi)==T] <- 0
  
  if (yes.exp){
    tj = sort(Ti[Di==1&Ti<=t0]); tj = c(t0,tj); pi.xi = sum.I(Ti,"<=",Ti)/N
    I.xi.le.tj = 1*(Ti <= VTM(tj,N)); I.Di0 = 1*(Di==0)
    pi.tj = sum.I(tj,"<=",Ti)/N
    
    Til = sort(Ti); Ghat.Til = Ghat.Ti[order(Ti)]
    dGhat.Til = c(-Ghat.Til[-1],1) - (-Ghat.Til)
    junk = rep(0,N)
    junk[order(Ti)] = dGhat.Til/Ghat.Til
    
    term1 = matrix(0,nrow=N,ncol=length(tj))
    term1[(I.xi.le.tj*I.Di0)==1] = matrix(1/pi.xi,nrow=N,ncol=length(tj))[(I.xi.le.tj*I.Di0)==1]  
    term2 = sum.I(Ti,">=",Ti,I.xi.le.tj*junk/pi.xi)
    LamC.exp = term1 - term2
    return(list("wgt"=Wi,"Ghat.Xi"=Ghat.Ti,"LamC.exp"=as.matrix(rbind(tj,LamC.exp)),ncol=length(tj))) 
  } else return(Wi)   
}

NCC.Fit.FUN <- function(data.ncc,model="COX",model.family=NULL,t0=NULL,rtn="EST")
{
  ## data.ncc: ID, Vi.case, wgtk, wgtk.C, time, delta, covariate 
  data.ncc = data.ncc[,-c(1,2)]
  data.ncc = as.matrix(data.ncc); mode(data.ncc)="numeric"
  wgtk = data.ncc[,1]; wgtk.cen = data.ncc[,2]; 
  xk = data.ncc[,3]; dk = data.ncc[,4]; yk = data.ncc[,-(1:4),drop=F]; 
  if(model=="COX"){
    betahat = coxph(Surv(xk,dk)~yk,weights=wgtk,robust=T,subset=wgtk>0); betahat = summary(betahat)$coef[,c(1,4),drop=F]
  } else {	
    ## GLM fit with intercept 
    betahat = glm(1*(xk <= t0) ~ yk, weight=wgtk*wgtk.cen,family=model.family); betahat = summary(betahat)$coef[,c(1,2),drop=F]
  }
  if(rtn=="EST"){ betahat = betahat[,1]}
  betahat 
}

PI.k.FUN <- function(tt,ebyi,xi,yi,wgti,k0=0)
{
  out = ebyi*wgti; yi=as.matrix(yi); py = ncol(yi); 
  if(k0==1){out=out*yi}
  if(k0==2){out=c(out)*yi[,rep(1:py,py)]*yi[,rep(1:py,rep(py,py))]}
  as.matrix(sum.I(tt,"<=",xi,out))
}

Cox.Beta.EXP <- function(data.ncc)
{
  ## data.ncc: ID, V1,wgtk, wgtk.C, time, delta, score, covariate
  data.ncc = data.ncc[,-c(1,2)]
  data.ncc = na.omit(as.matrix(data.ncc)); wgtk = data.ncc[,1]; xk = data.ncc[,3]; dk = data.ncc[,4]; 
  ebyk = exp(data.ncc[,5]); yk = data.ncc[,-(1:5),drop=F]; py = ncol(yk)
  tmpind = dk==1; tj = xk[tmpind]; 
  pi0.tj   = c(PI.k.FUN(tj,ebyk,xk,yk,wgtk,k0=0))/sum(wgtk) 
  pi1.tj   =   PI.k.FUN(tj,ebyk,xk,yk,wgtk,k0=1)/sum(wgtk)    ## tj x py matrix ##
  Ahat = PI.k.FUN(tj,ebyk,xk,yk,wgtk,k0=2)/sum(wgtk)*pi0.tj - pi1.tj[,rep(1:py,py)]*pi1.tj[,rep(1:py,rep(py,py))]
  Ahat = matrix(apply(Ahat/pi0.tj^2*wgtk[tmpind],2,sum),ncol=py)/sum(wgtk)
  term1 = matrix(0,ncol=py,nrow=length(dk))
  term1[tmpind,] = yk[tmpind,] - pi1.tj/pi0.tj
  term2 = (sum.I(xk,">=",tj,wgtk[tmpind]/pi0.tj)*ebyk*yk - sum.I(xk,">=",tj,wgtk[tmpind]*pi1.tj/pi0.tj^2)*ebyk)/sum(wgtk)
  betaexp = (term1 - term2)%*%solve(Ahat)  ## betahat - beta = \sum_k wgtk * betaexp_k
  betaexp
}

GLM.LINK <- function(glm.family){
  link0 = NULL
  if (glm.family$family=='binomial'){
    switch(glm.family$link,
           "logit"={
             link0$g <- function(eta){exp(eta)/(1+exp(eta))}
             link0$dmu <- function(eta){exp(eta)/(1+exp(eta))^2}
           },
           "cloglog"={
             link0$g <- function(eta){1-exp(-exp(eta))}
             link0$dmu <- function(eta){exp(-exp(eta))*exp(eta)}
           },
           "probit"={
             link0$g <- function(eta){pnorm(eta)}
             link0$dmu <- function(eta){dnorm(eta)}
           }
    )
  } # end if (glm.family$family=='binomial')
  return(link0)
}

GLM.Beta.EXP <- function(data.ncc,family,t0)
{
  link0 = GLM.LINK(family)
  ## data.ncc: ID, V1, wgtk, wgtk.C, time, delta, score, covariate
  data.ncc = data.ncc[,-c(1,2)]
  data.ncc = as.matrix(data.ncc)
  data.ncc = na.omit(as.matrix(data.ncc)); N = nrow(data.ncc)
  wgtk = data.ncc[,1]; wgt.cen = data.ncc[,2]; xk = data.ncc[,3]; dk = data.ncc[,4]; sk = data.ncc[,5]; 
  yk = data.ncc[,-c(1:5),drop=F]; py = ncol(yk)
  muk = link0$g(sk); dmuk = link0$dmu(sk)
  Ahat = apply(wgtk*wgt.cen*dmuk*yk[,rep(1:py,py)]*yk[,rep(1:py,rep(py,py))],2,sum,na.rm=T)/sum(wgtk*wgt.cen)
  Ahat = matrix(Ahat,ncol=py)
  # browser()
  betaexp = ((1*(xk<=t0)-muk)*yk)%*%solve(Ahat) ## betahat - beta = \sum_k wgtk * wgt.cen * betaexp_k
  betaexp
}


Kern.FUN <- function(zz,zi,bw,kern0="gauss") 
{ 
  ## returns an (n x nz) matrix ##
  out = (VTM(zz,length(zi))- zi)/bw
  switch(kern0,
         "epan"= 0.75*(1-out^2)*(abs(out)<=1)/bw,
         "gauss"= dnorm(out)/bw
  )
}

dACC.FUN <- function(uu0, uu=SE.yy, A.uu = Sp.yy, bw=NULL)
{
  data = cbind(uu,A.uu); data=na.omit(data); data=data[apply(abs(data),1,sum)<Inf,]
  uu=data[,1]; A.uu=data[,2]; n.uu = length(uu)
  A.uu = A.uu[order(uu)]; uu = sort(uu)
  if(is.null(bw)){bw=1.06*min(sd(uu),IQR(uu)/1.34)*n.uu^(-0.5)}
  Ki.u0 = Kern.FUN(uu0, uu[-1], bw) ## n.uu x n.u0
  c(t(A.uu[-1]-A.uu[-n.uu])%*%Ki.u0)
}


wRk.Lam0exp.FUN <- function(data.ncc,wRk,LamCexp)
{  
  ## Expansion of \sum_k (\what_Ck - w_Ck) * Rk ##
  ## data.ncc : wgt.ncc, wgt.cen, xk, dk
  data.ncc = as.matrix(data.ncc); wgtk = data.ncc[,1]; wgtk.cen = data.ncc[,2]; 
  xk = data.ncc[,3]; dk = data.ncc[,4]; nc = length(dk)
  t0=LamCexp[1,1]; tj = LamCexp[1,-1]; LamCexp.tj = LamCexp[-1,,drop=F]; nj = length(tj)
  LamCexp.t0 = LamCexp.tj[,1]; xi.R.t0 = sum.I(t0,"<",xk,wgtk.cen*wRk*wgtk)/sum(wgtk.cen*wgtk) #1 x nR %%%%%
  if (nj>0) {
    LamCexp.tj = LamCexp.tj[,-1,drop=F] ## N x nj
    xi.R.tj = matrix(sum.I(tj,">=",xk,wgtk.cen*wRk*dk*wgtk),nrow=nj,ncol=ncol(wRk))/sum(wgtk.cen*wgtk) #nj x nR
    LamCexp.tj%*%(xi.R.tj-rbind(0,xi.R.tj[-nj,,drop=F])) + LamCexp.t0%*%matrix(xi.R.t0,nrow=1)# N x nR
  } else { LamCexp.t0%*%matrix(xi.R.t0,nrow=1) }
}

dACC.Beta.FUN <- function(data.ncc,type,u0,t0,betahat,eps)
{
  # data.ncc: ID, V1, wgtk, wgtk.c, Xi, Di, Score, Yi
  pp = length(betahat); d.acc.u = NULL
  if(length(eps)==1){eps = eps/apply(data.ncc[,-(1:7),drop=F],2,sd,na.rm=T)}
  for(ip in 1:pp) {
    betapnew = betahat + eps[ip]*diag(pp)[,ip]
    betamnew = betahat - eps[ip]*diag(pp)[,ip]
    scorepnew = as.matrix(data.ncc[,-(1:7),drop=F])%*%betapnew; datapnew = data.ncc; datapnew[,7] = scorepnew
    scoremnew = as.matrix(data.ncc[,-(1:7),drop=F])%*%betamnew; datamnew = data.ncc; datamnew[,7] = scoremnew
    accpnew = ROC.DIPW.FUN(data.ncc=datapnew,type=type,u0=u0,t0=t0,rtn="EST",sigmoid=F)
    accmnew = ROC.DIPW.FUN(data.ncc=datamnew,type=type,u0=u0,t0=t0,rtn="EST",sigmoid=F)
    d.acc.u = cbind(d.acc.u, (accpnew$ACC.u0-accmnew$ACC.u0)/(2*eps[ip]))
  } # end for ip
  return(d.acc.u)
}




VarAdj.FUN <- function(datA,datB=NULL,adj=NULL){ 
  
  adj.index <- function(v,adj=adj){
    out = rep(0,length(v))
    index = which(is.na(v)==F & v>0)
    junk = v[index]
    junk2 = log(junk)-median(log(junk))
    if (is.null(adj)) out[index]=rep(1,length(junk)) else out[index]=1*(junk2<=adj*sd(junk2))
    out
  }	
  ## always use datB with larger number of matching number if m0A != m0B ##
  if(!is.null(datB)){if(datB$m0<datA$m0){tmpdat=datB;datB=datA;datA=tmpdat}}
  pi.tjA=datA$pi.tj; pi.dA=datA$pi.d; I0kjA=datA$I0kj; m0A=datA$m0; m1A=datA$m1; cexp.nccA=datA$cexp.ncc; datA=datA$dat;
  id.nccA = cexp.nccA[,1]; nomiss.nccA=cexp.nccA[,2]; d.nccA=cexp.nccA[,3]; cexp.nccA=cexp.nccA[,-(1:3),drop=F]
  datA0 = matrix(0,ncol=ncol(datA),nrow=length(nomiss.nccA)); datA0[nomiss.nccA==1,] = datA
  wkA=datA0[,1];  wckA=datA0[,2]; wcrkA=datA0[,-(1:2),drop=F];   	
  if (!is.null(datB)){
    pi.tjB=datB$pi.tj; I0kjB=datB$I0kj; m0B=datB$m0; cexp.nccB=datB$cexp.ncc; datB=datB$dat
    id.nccB = cexp.nccB[,1]; nomiss.nccB=cexp.nccB[,2]; d.nccB=cexp.nccB[,3]; cexp.nccB = cexp.nccB[,-c(1:3),drop=F]
    datB0 = matrix(0,ncol=ncol(datB),nrow=length(nomiss.nccB)); datB0[nomiss.nccB==1,] = datB
    wkB=datB0[,1]; wckB=datB0[,2]; wcrkB=datB0[,-c(1,2),drop=F]
    AinB = match(id.nccA,id.nccB)
  } # end if (!is.null(datB))
  
  L = ncol(wcrkA); var.out = NULL
  for (l in 1:L){
    # browser()
    wcrkA0 = wcrkA[,l]; cexp.nccA0 = cexp.nccA[,l]
    v.naive0 = wkA*((wkA-1)*wcrkA0^2+(wcrkA0+cexp.nccA0)^2)
    naive.ind = adj.index(v.naive0,adj=adj)
    var.naive0 = sum(v.naive0[naive.ind==1],na.rm=T)/sum(wkA[naive.ind==1],na.rm=T)
    
    if (m0A < Inf){
      eta1.tjzjA = sum((wcrkA0*wkA*(wkA-1))[naive.ind==1 & d.nccA==1],na.rm=T)/sum(wkA[naive.ind==1],na.rm=T)
      var1.adjA = (m1A*eta1.tjzjA^2/pi.dA^2)/sum(wkA[naive.ind==1],na.rm=T)
      eta0.tjzjA = t(I0kjA[naive.ind[d.nccA==0]==1,,drop=F])%*%((wcrkA0*wkA*(wkA-1))[d.nccA==0 & naive.ind==1])/sum(wkA[naive.ind==1],na.rm=T)
      var0.adjA = sum(m0A*eta0.tjzjA^2/pi.tjA^2,na.rm=T)/sum(wkA[naive.ind==1],na.rm=T)
      var.adjA0 = var1.adjA+var0.adjA
      out = c(var.naive0, var.naive0 - var.adjA0)/sum(wkA,na.rm=T)
    } else {
      out = c(var.naive0,var.naive0)/sum(wkA,na.rm=T)
    } 
    
    if (!is.null(datB)){
      wcrkB0 = wcrkB[,l]; cexp.nccB0 = cexp.nccB[,l]
      if (m0B > m0A){
        v.naive0 = wkA*(wkB[AinB]-1)*(wcrkB0[AinB]-wcrkA0)^2 + wkA*(wkA-wkB[AinB])*wcrkA0^2 + wkA*(wcrkA0-wcrkB0[AinB]+cexp.nccA0-cexp.nccB0[AinB])^2; 
        naive.ind = adj.index(v.naive0,adj=adj)				
        var.naive0 = sum(v.naive0[naive.ind==1],na.rm=T)/sum(wkA[naive.ind==1],na.rm=T)
        
        if (m0B<Inf){
          eta1.tjzjA = sum((wcrkA0*wkA*(wkA-1))[naive.ind==1 & d.nccA==1],na.rm=T)/sum(wkA[naive.ind==1],na.rm=T)
          var1.adjA = (m1A*eta1.tjzjA^2/pi.dA^2)/sum(wkA[naive.ind==1],na.rm=T)
          eta0.tjzjA = t(I0kjA[naive.ind[d.nccA==0]==1,,drop=F])%*%((wcrkA0*wkA*(wkA-1))[d.nccA==0 & naive.ind==1])/sum(wkA[naive.ind==1],na.rm=T) 
          naive.indB = rep(1,length(wcrkB0)); naive.indB[AinB]=naive.ind
          eta0.tjzjB = t(I0kjB[naive.indB[d.nccB==0]==1,,drop=F])%*%((wcrkB0*wkB*(wkB-1))[d.nccB==0 & naive.indB==1])/sum(wkB[naive.indB==1],na.rm=T)
          var0.adj = sum(m0A*(eta0.tjzjB-eta0.tjzjA)^2/pi.tjB^2,na.rm=T)/sum(wkA[naive.ind==1],na.rm=T) + sum((m0B-m0A)*eta0.tjzjB^2/pi.tjB,na.rm=T)/sum(wkB[naive.indB==1],na.rm=T)
          var.adj0 = var1.adjA + var0.adj
        } else {var.adj0 = var.adjA0} 
      } else {
        ## m0B = m0A
        
        v.naive0 = wkA*(wkA-1)*(wcrkA0-wcrkB0)^2 + wkA*(wcrkA0-wcrkB0+cexp.nccA0-cexp.nccB0)^2
        naive.ind = adj.index(v.naive0,adj=adj)
        var.naive0 = sum(v.naive0[naive.ind==1],na.rm=T)/sum(wkA[naive.ind==1],na.rm=T)
        
        eta1.tjzjA = sum((wcrkA0*wkA*(wkA-1))[naive.ind==1 & d.nccA==1],na.rm=T)/sum(wkA[naive.ind==1],na.rm=T)
        var1.adjA = (m1A*eta1.tjzjA^2/pi.dA^2)/sum(wkA[naive.ind==1],na.rm=T)
        naive.indB = rep(1,length(wcrkB0)); naive.indB[AinB]=naive.ind
        eta0.tjzjB = t(I0kjB[naive.indB[d.nccB==0]==1,,drop=F])%*%((wcrkB0*wkB*(wkB-1))[d.nccB==0 & naive.indB==1])/sum(wkB[naive.indB==1],na.rm=T)
        var0.adj = sum(m0B*(eta0.tjzjB+eta0.tjzjA)^2/pi.tjA^2,na.rm=T)/sum(wkB[naive.indB==1],na.rm=T); 
        var.adj0 = var1.adjA + var0.adj
        var.adj0 = var.adj0*(var.naive0>var.adj0)
      }
      out = c(var.naive0,var.naive0 - var.adj0)/sum(wkB,na.rm=T)	
    } # end if (!is.null(datB))
    var.out = rbind(var.out,out)
  } # end for l
  var.out  
}

ROC.DIPW.FUN <- function(data.ncc,type="FPR",u0,t0,rtn="EST",sigmoid=F,var.control=list(wgt.ncc.dat=NULL,betafit=NULL,data0=NULL,model="COX",family=NULL,LamCexp=NULL)) 
{  	
  ## ===================================================================##
  ## 
  ## ===================================================================##
  # browser()
  ## == data.ncc.f: ID, Vi.case, wgt.ncc, wgtk.cen, Xi, Di, score, Yi == ##
  data.ncc.f = data.ncc; id.ncc = data.ncc[,1]; vi.case.ncc = data.ncc[,2]
  ## == data.ncc: wgt.ncc, wgt.cen, Xi, Di, score, Yi == ##	
  data.ncc = data.ncc[,-c(1,2)]; nomiss.ck = !is.na(data.ncc[,5]);
  data.ncc = as.matrix(data.ncc); data.ncc0 = data.ncc; ## data.ncc0: with missing score
  ## == data.ncc: without missing score == ##
  data.ncc = na.omit(data.ncc) 
  wgtk = data.ncc[,1]; wgtk.c = data.ncc[,2]; xk = data.ncc[,3]; dk = data.ncc[,4]; ck = data.ncc[,5]  
  scl = sort(ck); nv = length(ck)
  
  if (!sigmoid){
    ## empirical estimate
    St0.Fcl = sum.I(scl,">=",ck,wgtk*wgtk.c*(xk >= t0))/sum(wgtk.c*wgtk)
    Ft0.Fcl = sum.I(scl,">=",ck,wgtk*wgtk.c*(xk <  t0))/sum(wgtk.c*wgtk)
    Fcl = sum.I(scl,">=",ck,wgtk*wgtk.c)/sum(wgtk.c*wgtk) ## I(score <= c)
  } else {
    scl = seq(min(scl),max(scl),length=1000)
    bw=1.06*min(sd(ck),IQR(ck)/1.34)*nv^(-0.4)
    CK.mat = pnorm((VTM(scl,nv)-ck)/bw)
    St0.Fcl = c(t(wgtk*wgtk.c*(xk >= t0))%*%CK.mat)/sum(wgtk.c*wgtk)
    Ft0.Fcl = c(t(wgtk*wgtk.c*(xk <  t0))%*%CK.mat)/sum(wgtk.c*wgtk)
    Fcl = c(t(wgtk*wgtk.c)%*%CK.mat)/sum(wgtk.c*wgtk)
  }
  St0 = max(St0.Fcl); Ft0 = 1-St0## St0 = P(T> t0); Ft0 = P(T<=t0)
  FPR.cl= (St0-St0.Fcl)/St0     ## P(Y> cl|T> t0)
  TPR.cl= (Ft0-Ft0.Fcl)/Ft0     ## P(Y> cl|T<=t0)
  NPV.cl= St0.Fcl/Fcl           ## P(T> t0|Y<=cl)
  PPV.cl= (Ft0-Ft0.Fcl)/(1-Fcl); if (sum(Fcl==1)>0) PPV.cl[Fcl==1] = 1  ## P(T<=t0|Y> cl)
  AUC = sum(TPR.cl*(FPR.cl-c(FPR.cl[-1],0)))
  nm.acc = c("FPR","TPR","NPV","PPV")
  acc.cl = data.frame("cutoff"=scl,"FPR"=FPR.cl,"TPR"=TPR.cl,"NPV"=NPV.cl, "PPV"=PPV.cl)    
  # browser()
  
  ## transform from c by TPR/FPR
  acc.ul = acc.cl; ind0 = match(type,names(acc.cl)); 
  ul = acc.ul[,ind0]; acc.ul = acc.ul[order(ul),]; ul = sort(ul); # sorted by ul
  indu = sum.I(u0,">=",ul)
  c.u0 = acc.ul[indu,1]; 
  summary.est = c("AUC"=AUC, unlist(acc.ul[indu,-c(1,ind0)]))
  
  uu = seq(0.001,0.999,by=0.001)
  FPR.scl = sort(FPR.cl); TPR.scl = TPR.cl[order(FPR.cl)]
  ind.u = sum.I(uu,">=",FPR.scl)
  roc.uu = TPR.scl[ind.u]
  
  est.out = list("ACC.u0"=summary.est,"ALL"=acc.cl,"roc"=cbind(uu,roc.uu),"St0"=St0)
  
  if(rtn=="EST"){ # Point estimate only
    return(est.out)
  } else { # Variance estimate
    wgt.ncc.dat = var.control$wgt.ncc.dat; betafit = var.control$betafit; data0 = var.control$data0; model = var.control$model; LamCexp = var.control$LamCexp; family=var.control$family
    nu0 = length(u0); nck = length(ck)
    acc.u0 = acc.ul[indu,]
    tmpind = sum.I(c.u0,">=",scl); F.u0 = Fcl[tmpind]
    U.ACC.u0 = as.list(1:4); names(U.ACC.u0)=c("FPR","TPR","PPV","NPV")
    I.ck.u0 = 1*(ck >= VTM(c.u0,nck)); ## I(score > c)
    U.ACC.u0$FPR = (xk > t0)*(I.ck.u0 - VTM(acc.u0$FPR,nck))/St0 ## exp for FPRhat(c)-FPR(c)
    U.ACC.u0$TPR = (xk <= t0)*(I.ck.u0 - VTM(acc.u0$TPR,nck))/(1-St0) ## exp for TPRhat(c)-TPR(c)
    U.ACC.u0$PPV = I.ck.u0*(1*(xk<=t0) - VTM(acc.u0$PPV,nck))/(1-VTM(F.u0,nck)) ## exp for PPVhat(c)-PPV(c)
    U.ACC.u0$NPV = (1-I.ck.u0)*(1*(xk>t0) - VTM(acc.u0$NPV,nck))/VTM(F.u0,nck)## exp for NPVhat(c)-NPV(c) 
    U.AUC = (xk<=t0)/(1-St0)*(1-FPR.cl[rank(ck)]-AUC)+(xk>t0)/St0*(TPR.cl[rank(ck)]-AUC)    		
    U.ACC = as.list(1:4)
    names(U.ACC) = c("AUC",nm.acc[-(ind0-1)])
    var.ACC = NULL    		
    Di = data0[,1]; Vi = data0[,2]    
    
    data.fun = function(Uexp,Yes.Cexp=T){
      exp.cen.k = wRk.Lam0exp.FUN(data.ncc,Uexp,LamCexp)[Vi==1,,drop=F]
      # browser()
      exp.cen.k = cbind(id.ncc,nomiss.ck,data.ncc0[,4],exp.cen.k*Yes.Cexp)# indicator of missing score, Dk, (if Yes.Cexp=F, 0; if Yes.Cexp=T, exp.cen.k)
      tmpdata=cbind(wgtk,wgtk.c^Yes.Cexp,Uexp*wgtk.c^Yes.Cexp) # wgt.ncc, (if Yes.Cexp=F, 1; if Yes.Cexp=T, wgt.cen); (if Yes.Cexp=F, Uexp; if Yes.Cexp=T, Uexp*wgt.cen)
      list("data"=tmpdata,"pi.tj"=wgt.ncc.dat$pi.tj,"pi.d"= mean(Di), "I0kj"=wgt.ncc.dat$I0kj,"m0"=wgt.ncc.dat$m,"m1"=sum(vi.case.ncc),"cexp.ncc"=exp.cen.k)
    }
    ## ==== variance for beta ==== ##
    ## ==== expansion for betahat ==== ##
    if (!is.null(betafit)){   				
      switch(model,
             "COX"={
               beta.exp = Cox.Beta.EXP(data.ncc.f)
               U.beta = data.fun(beta.exp,Yes.Cexp=F)
             },
             "GLM"={
               beta.exp = GLM.Beta.EXP(data.ncc.f,family=family,t0=t0)
               U.beta = data.fun(beta.exp,Yes.Cexp=T)
             }
      )# end for switch
      var.beta = VarAdj.FUN(datA=U.beta,datB=NULL)	
      tmpeps = betafit[,2]
      d.acc = dACC.Beta.FUN(data.ncc=data.ncc.f,type=type,u0=u0,t0=t0,betahat=betafit[,1],eps=tmpeps)
      acc.beta.exp = beta.exp%*%t(d.acc)
      
      
    } else { # if beta is true
      U.beta = NULL; acc.beta.exp = 0
    } # end for lese
    ## ==== variance for ACC ==== ##
    for (kk in 1:length(U.ACC)){
      tmpnm = names(U.ACC)[kk]
      if (tmpnm=="AUC"){
        U.ACC[[kk]]= as.matrix(U.AUC,ncol=1); U.acc.beta = acc.beta.exp[,1,drop=F]
      } else {
        dACC.hat = dACC.FUN(u0,uu=ul,A.uu=acc.ul[,match(tmpnm,names(acc.ul))])
        U.ACC[[kk]] = as.matrix(U.ACC.u0[[match(tmpnm,names(U.ACC.u0))]] - VTM(dACC.hat,nck)*U.ACC.u0[[match(type,names(U.ACC.u0))]],ncol=nu0)
        junk = acc.beta.exp[,-1]
        U.acc.beta = junk[,(1:nu0)+(kk-2)*nu0,drop=F]   
      } # end for else 
      
      switch(model,
             "COX"={
               U.ACC[[kk]] = data.fun(U.ACC[[kk]])
               U.ACC[[kk]]$data[,-c(1,2)] = U.ACC[[kk]]$data[,-c(1,2)] + U.acc.beta
             },
             "GLM"={
               U.ACC[[kk]]=U.ACC[[kk]]+U.acc.beta
               U.ACC[[kk]]=data.fun(U.ACC[[kk]])
             }
      )	# end for switch
      tmpvar = VarAdj.FUN(datA=U.ACC[[kk]],datB=NULL)
      var.ACC = rbind(var.ACC,tmpvar)
    } # end for kk
    var.out = list("Ui.ACC"=U.ACC,"Ui.Beta"=U.beta,"VAR.ACC"=var.ACC,"VAR.Beta"=var.beta)
    list("est"=est.out,"var"=var.out)
  } # end for else	
}
ROC.DIPW.CV.632 <-function(data.ncc,type="FPR",u0,t0,n.rep=200,model) ## cross-validation
{ 
  model.type = model$type; model.family=model$family
  ## data.ncc: ID, V1, wgtk, wgtk.C, time, delta, covariate
  thetahat.all = try(NCC.Fit.FUN(data.ncc=data.ncc,model=model.type,model.family=model.family,t0=t0,rtn="EST"),silent=T)	
  if(substring(thetahat.all[1],1,5)!="Error"){
    if (model.type=="GLM") score = as.matrix(cbind(1,data.ncc[,-c(1:6),drop=F]))%*%thetahat.all else score = as.matrix(data.ncc[,-c(1:6),drop=F])%*%thetahat.all
    tmpdata = cbind(data.ncc[,1:6],score)
    rocout.all=ROC.DIPW.FUN(data.ncc=tmpdata,type=type,u0=u0,t0=t0,rtn="EST")
    acc.sum.all = rocout.all$ACC.u0;
    St0.all = rocout.all$St0
    roc.uu.all = rocout.all$roc[,2]
    uu = rocout.all$roc[,1]
  }
  acc.sum.cv = St0.cv = roc.uu.cv = NULL; set.seed(702)
  for(rep in 1:n.rep){
    set.seed(1000+rep)
    n0 = nrow(data.ncc)
    ind.t = sort(sample(1:n0,n0,replace=T))
    ind.v=setdiff(1:n0,unique(ind.t))
    thetahat.tr = try(NCC.Fit.FUN(data.ncc[ind.t,],model=model.type,model.family=model.family,t0=t0,rtn="EST"),silent=T)
    if(substring(thetahat.tr[1],1,5)!="Error"){
      if (model.type=="GLM") score = as.matrix(cbind(1,data.ncc[ind.v,-c(1:6),drop=F]))%*%thetahat.tr else score = as.matrix(data.ncc[ind.v,-c(1:6),drop=F])%*%thetahat.tr
      tmpdata = cbind(data.ncc[ind.v,1:6],score)
      rocout=ROC.DIPW.FUN(data.ncc=tmpdata,type=type,u0=u0,t0=t0,rtn="EST");  
      acc.sum = rocout$ACC.u0; St0 = rocout$St0; roc.uu = rocout$roc[,2]
      
      acc.sum.cv = cbind(acc.sum.cv,0.632*acc.sum+0.368*acc.sum.all)
      St0.cv = c(St0.cv, 0.632*St0+0.368*St0.all)
      roc.uu.cv = cbind(roc.uu.cv,0.632*roc.uu+0.368*roc.uu.all)
    } # end for if
  }  # end for rep
  return(list("ACC.u0"=acc.sum.cv,"St0"=St0.cv,"roc"=cbind(uu,roc.uu.cv)))
}
