#######################################
# AD DISTANCES OBTAINED BY SIMULATION
#######################################
library(parallel)

#######################################
# This program writes many files
# You must indicate in which directory you want these files
# changing the following path:
setwd(".") 

# Initial parameters
must=0
sigmast=1
ite = 1000000
size = 300
mosaic = 49
parallel=TRUE

# Random generation function
rsepd2=function(n,must=0,sigmast=1,alpha=0.5,p=2){
  kp=1/(2*gamma(1+1/p)*p^(1/p))
  u=runif(n)
  w=rgamma(n,shape=1/p,scale=1)
  y=ifelse(u>alpha,1-alpha,-alpha)*(w^(1/p))/gamma(1+1/p)
  y=y/kp
  must+sigmast*y
}
# Distribution function
psepd2=function(x, must=0, sigmast=1, alpha=0.5, p=2){
  x=(x-must)/sigmast
  alpaux=ifelse(x<=0,alpha,1-alpha)
  alpha+sign(x)*alpaux*ifelse(alpaux>0,pgamma(((abs(x)/(2*alpaux))^p)/p,shape=1/p),0)
}

if (parallel){
  mc=detectCores(logical=TRUE)-1
  cl=makeCluster(mc)
}

cell<-function(par){
  alpha=par[1]
  p=par[2]
  cat("\nAlpha=",alpha," p=",p," time=", as.character(Sys.time()))
  
  calcStat=function(iter){
      samp = matrix(rsepd2(size*iter,must,sigmast,alpha,p),nrow=iter)  
      if (!parallel){
        apply(samp,1,simul,alpha=alpha,p=p)
      } else {
        junk=clusterExport(cl, c("alpha","p"),envir = environment())
        parApply(cl,samp,1,simul,alpha=alpha,p=p)
      }
    }
  
  filename=paste0("Sam_Size_",size,"_p_",p,"_alpha_",alpha,"_Ite_",ite,".csv")
  if (!file.exists(filename)){
      if (size<200) {
        result=calcStat(ite)
      } else{
        result=c(apply(matrix(rep(ite/100,100)),1,calcStat))
      }
      write.csv2(result,filename, row.names=FALSE)
  }
}

simul<-function(x,alpha,p){
  equ=function(par,sam){
    dif=par-sort(sam)
    xmen=sum((dif[dif>0])^(p-1))
    xmas=sum((-dif[dif<0])^(p-1))
    (1-alpha)^p*xmen-alpha^p*xmas
  }
  # AD distance function
  AD<-function (x, distr.fun, mu, sigma, alpha, p,...) 
  {
    x <- sort(distr.fun(x, mu, sigma, alpha, p))
    tmp <- x * (1 - rev(x))
    tmp <- (2 * seq(x) - 1) * log(tmp+.Machine$double.eps)
    tmp <- -mean(tmp) - length(x)
  }
  if (p==1){
    must2=quantile(x, probs = alpha, type =4)                           
  }else{
    must2=uniroot(equ,c(min(x),max(x)),tol=1e-16,sam=x)$root
  }
  dif=must2-sort(x)
  xmen=sum((dif[dif>0])^p)/size
  xmas=sum((-dif[dif<0])^p)/size
  sigst2=((xmas/(1-alpha)^p+xmen/alpha^p))^(1/p)/2
  AD(x, psepd2, mu=must2, sigma=sigst2, alpha=alpha, p=p)
}

j=(log10(50)/log10(2))^(1/(mosaic/2-0.5))
lp=round(c(1,2^(1/j^((mosaic/2-1.5):1)),2,2^(j^((1:(mosaic/2-0.5))))),4)
la = round(seq(from=0, to=0.5, by=1/48),4)
la[1] = 0.01

if (parallel) junk=clusterExport(cl, c("simul","size","psepd2"))

apply(expand.grid(la,lp),1,cell)

if (parallel) stopCluster(cl)
end = Sys.time()

