################################################################
#    Supplementary material that accompanies the article:      #
# "Concurrent Goodness-of-fit Test for Unimodal Distributions" #
#    Published in:                                             #
################################################################

# The working directory have to be the same where the file CriticalValues.csv is

# You can change the values of the following parameters:
sample="RANDOM" # sample: DATA=from file, RANDOM=random generation
   # If sample="DATA" The sample has to be in an ASCCI file with all the values in a column without header and named sample.dat)
   # If sample!="RANDOM" The data are randomly generated
size    = 30   # Size for random samples
loc_p   = 0    # Location parameter for random samples
scale_p = 1    # Scale parameter for random samples
alpha   = 0.5  # Alpha parameter for random samples
p       = 2    # p parameter for random samples
set.seed(123)  # You can change the seed (or eliminate it) for random samples


#########################################################
# The following functions are needed

# Random generation function
rsepd2=function(n,must,sigmast,alpha,p){
  kp=1/(2*gamma(1+1/p)*p^(1/p))
  u=runif(n)
  w=rgamma(n,shape=1/p,scale=1)
  y=ifelse(u>alpha,1-alpha,-alpha)*(w^(1/p))/gamma(1+1/p)
  y=y/kp
  must+sigmast*y
}
# Distribution function
psepd2=function(x, must, sigmast, alpha, p){
  x=(x-must)/sigmast
  alpaux=ifelse(x<=0,alpha,1-alpha)
  alpha+sign(x)*alpaux*ifelse(alpaux>0,pgamma(((abs(x)/(2*alpaux))^p)/p,shape=1/p),0)
}
#
##################################################################
# Loading and preparing the critical values for each cell and size

y = read.csv2("CriticalValues.csv")

tab<-lapply(split(y[,-c(1,2)],list(y$q,y$n)),function(el){
  da<-xtabs(VC~alpha+p,el)
  rbind(da,da[24:1,])})


###################################################
samp=switch(sample,
  #DATA   = list(values=unlist(read.table("sample.dat",header=F))),
  DATA   = as.numeric(read.table("sample.dat",header=F)[,1]),
  RANDOM = rsepd2(size,loc_p, scale_p, alpha, p)
)

size=length(samp)

###################################################
#Building the table with the reference quantiles for each cell

s <- c(5,6,7,8,9,10,11,12,13,14,15,20,25,30,35,40,45,50,55,60,65,70,75,80,90,100,150,200,250,300,400,500,1000)

if (size<1000){
  pos=findInterval(size,s)
  ref90<-(tab[[paste0("90.",s[pos])]]*(s[pos+1]-size)+tab[[paste0("90.",s[pos+1])]]*(size-s[pos]))/(s[pos+1]-s[pos])
  ref95<-(tab[[paste0("95.",s[pos])]]*(s[pos+1]-size)+tab[[paste0("95.",s[pos+1])]]*(size-s[pos]))/(s[pos+1]-s[pos])
  ref99<-(tab[[paste0("99.",s[pos])]]*(s[pos+1]-size)+tab[[paste0("99.",s[pos+1])]]*(size-s[pos]))/(s[pos+1]-s[pos])
} else {
  ref90<-tab[["90.1000"]]
  ref95<-tab[["95.1000"]]
  ref99<-tab[["99.1000"]]
}
ref<-cbind(ref90,ref95,ref99)
dim(ref)<-c(49,49,3)
###################################################

###################################################
#Calculation of the distance for each cell

distan<-function(par,samp){
  alpha=par[1]
  p=par[2]
  equ=function(par,sam,alpha=alpha,p=p){
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
    must2=quantile(samp, probs = alpha, type =4)                           
  }else{
    must2=uniroot(equ,c(min(samp),max(samp)),tol=1e-16,sam=samp,alpha=alpha,p=p)$root
  }
  dif=must2-sort(samp)
  xmen=sum((dif[dif>0])^p)/size
  xmas=sum((-dif[dif<0])^p)/size
  sigst2=((xmas/(1-alpha)^p+xmen/alpha^p))^(1/p)/2
  AD(samp, psepd2, mu=must2, sigma=sigst2, alpha=alpha, p=p)
}

j=(log10(50)/log10(2))^(1/(49/2-0.5))
lp=round(c(1,2^(1/j^((49/2-1.5):1)),2,2^(j^((1:(49/2-0.5))))),4)
la = round(seq(from=0, to=1, by=1/48),4)
la[1] = 0.01
la[49] = 0.99

dist=matrix(apply(expand.grid(la,lp),1,distan,samp=samp),nrow=49)

##################################################
# Plotting the Mosaic

plot (NULL, xlim=c(0.5,49+0.5), ylim=c(0.5,49+0.5), yaxs="i", xaxs="i", xaxt='n', yaxt='n', xlab="",ylab="")
for (i in 1:49) {
  for (j in 1:49) {
    rect(i-0.5,49-j+0.5,i+0.5,49-j+1.5, lty=3, border = "grey", col=gray(findInterval(dist[i,j],ref[i,j,])/3))
  }
}

# If the sample has been generated randomly with the parameters of a mosaic box, a red box is painted in that box
if (sample == "RANDOM") {
  cx = findInterval(alpha,la)
  cy = 49-findInterval(p,lp)+1
  if (alpha %in% la & p %in% lp) rect(cx-0.5, cy-0.5, cx+0.5, cy+0.5, lty=1, border = "red" )
}

