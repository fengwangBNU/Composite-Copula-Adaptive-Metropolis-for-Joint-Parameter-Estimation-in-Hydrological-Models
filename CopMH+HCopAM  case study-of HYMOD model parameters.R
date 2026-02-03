## ============================================================
# The code below is used to sample and compare the posterior distributions
# of HYMOD model parameters (five parameters in total) for the Yuanjiaan catchment
# using three MCMC algorithms: AM, CopMH, and HCopAM.
# The implementation illustrates how CopMH and HCopAM can be applied in practice,
# and can be readily adapted for use in other hydrological models or related studies.
# ============================================================



# ============================================================
# Preparation
# ============================================================

library(rtop)

library(ggplot2)
excess<-function(x_loss,cmax,bexp,Pval,PETval)
{
  output<-list()
  
  xn_prev = x_loss;
  if(((bexp+1)*(xn_prev)/cmax)>=1){
    ct_prev = cmax
  } else{
    ct_prev = cmax*(1-(1-((bexp+1)*(xn_prev)/cmax))^(1/(bexp+1)))
  };
  #ct_prev = cmax*(1-(1-((bexp+1)*(xn_prev)/cmax))^(1/(bexp+1)));
  UT1 = max((Pval-cmax+ct_prev),0.0);
  Pval = Pval-UT1;
  dummy = min(((ct_prev+Pval)/cmax),1);
  xn = (cmax/(bexp+1))*(1-(1-dummy)^(bexp+1));
  UT2 = max(Pval-(xn-xn_prev),0);
  evap = min(xn,PETval);
  xn = xn-evap;
  output$UT1<-UT1
  output$UT2<-UT2
  output$x_loss<-xn
  return(output)
  
}

linres<-function(x_slow,inflow,Rs)
{
  x_slow = (1-Rs)*x_slow + (1-Rs)*inflow;
  outflow = (Rs/(1-Rs))*x_slow;
  results <- list()
  results$outflow <-outflow
  results$x_slow <-x_slow
  return(results)
}

hymod_daily<-function(x,Extra,x_loss, x_slow,x_quick,t) 
{
  PET = Extra$E; 
  Precip = Extra$pre; 
  F = Extra$F;
  # Define the parameters
  cmax = x[1]; 
  bexp = x[2]; 
  alpha = x[3]; 
  Rs = x[4]; 
  Rq = x[5];
  outflow = c()
  
  Pval = Precip[t]
  PETval = PET[t]
  
  output<-list()
  
  output=excess(x_loss,cmax,bexp,Pval,PETval);
  
  UT1<-output$UT1
  UT2<-output$UT2
  x_loss<-output$x_loss
  
  # Partition UT1 and UT2 into quick and slow flow component
  UQ = alpha*UT2 + UT1
  US = (1-alpha)*UT2
  
  # Route slow flow component with single linear reservoir
  inflow = US; 
  
  temp<-list()
  temp <- linres(x_slow,inflow,Rs);
  x_slow<-temp$x_slow
  outflow<-temp$outflow
  #[x_slow,outflow] = linres(x_slow,inflow,outflow,Rs); 
  QS = outflow;
  
  inflow = UQ
  
  k = 1;
  while (k < 4){
    tempQ<-list()
    tempQ<-linres(x_quick[k],inflow,Rq)
    x_quick[k] <- tempQ$x_slow
    outflow<-tempQ$outflow
    inflow<-outflow	   
    k = k+1;
    
  }
  Qsim = (QS + outflow)*F;
  
  SimResults<-list()
  SimResults$Qsim <- Qsim
  SimResults$x_loss <- x_loss
  SimResults$x_slow <-x_slow
  SimResults$x_quick<-x_quick
  
  return(SimResults)
  
}


gape <-function(Qobs,Qsim){
  e <- c()
  e <- Qobs-Qsim
  return(e)}



Likelihood <- function(e,theta){
  mu <- theta[1]
  sigma <- theta[2]
  N <- length(e)
  logL <- -0.5*N*log(2*pi) - N*log(sigma) - sum(0.5*(e - mu)^2/sigma^2)
  return (logL)}


NSE <-function(Qobs,Qsim){
  E <- c()
  E <- 1 - (sum((Qobs-Qsim)^2)/sum((Qobs - mean(Qobs))^2))
  return(E)
}


#############################################################
path = 'E:/Case/'
setwd(path)


Filename <- paste("Yuanjiaan.txt",sep="")
rawdata <- read.table(Filename, header = FALSE) # pre ET, Qobs
data_all<-as.matrix(rawdata) #as.matrix()ת???ɾ???

T1 = 1981; 
T2 = 1982; 
T3 = 1983;
T4 = 1984;
row_index = 1;
temp<-matrix(0,1,6)
bound<-matrix( , ,6)
for (j in 1:nrow(data_all)){
  if((data_all[j,1]==T1)|(data_all[j,1]==T2)|(data_all[j,1]==T3)|(data_all[j,1]==T4))
  {
    for (k in 1:6)
    {temp[1,k] = data_all[j,k]}
    
    bound<-rbind(bound,temp)
    row_index = row_index+1
  }
}
bound = bound[2:nrow(bound),]


FileName <- paste("Yuanjiaan1981-1984",".txt",sep="")	
write.table(bound,file=FileName,row.names = FALSE, col.names = FALSE)



ExtraF = 1661 * (1000 * 1000 ) / (1000 * 60 * 60 * 24)
MaxT = 365+365+365+366
Qobs <- matrix(bound[,4],MaxT,1)
P <- matrix(bound[,6],MaxT,1)
E<-matrix(bound[,5],MaxT,1)

Extra <-list()
Extra$E <- E
Extra$pre<-P
Extra$F<-ExtraF
#??ֵ????

#####################################################################################
N=2000#  N denotes the length of the MCMC chain and can be adjusted according to practical considerations, such as model complexity
x_loss <-c()
x_slow<-c()
x_loss = 0.05;
x_slow = 2.3503/(0.05*ExtraF)
Tx_quick = matrix(0,3,1);
Q_prediction <-c()
x=matrix(0,N,5)
e=matrix(0,N,(MaxT-50))
Like_e=matrix(0,N,1)
nash=matrix(0,N,1)
ratio=matrix(0,N,1)
##
Par.minn = c(100,0.01,0.10,0.01,0.3)#ParRange.minn
Par.maxn = c(700,15,0.99,0.20,0.7)#ParRange.maxn

#******************************************************
#********************First 1000 sample **********************************
t=1;
for(i in 1:5){x[1,i]=runif(1,Par.minn[i],Par.maxn[i])}
x_loss <-c()
x_slow<-c()
x_loss = 0.05;
x_slow = 2.3503/(0.05*ExtraF)
Tx_quick = matrix(0,3,1);
Q_prediction <-c() 
while(t<=MaxT+1){
  dailySim <-list()
  dailySim <-hymod_daily(x[1,],Extra,x_loss,x_slow,Tx_quick[1:3,1],t)
  Q_prediction[t]<-dailySim$Qsim
  x_loss<-dailySim$x_loss
  x_slow<-dailySim$x_slow
  Tx_quick[1:3,1]<-matrix(dailySim$x_quick,3,1)
  t = t+1
}

Q_prediction_1=Q_prediction[1:MaxT]

e[1,] = gape(Qobs[51:MaxT],Q_prediction[51:MaxT]);
Like_e[1,]= Likelihood(e[1,],theta=c(0,sd(e[1,])))
ratio[1,]=1

nash[1,]=NSE(Qobs[51:MaxT],Q_prediction[51:MaxT])

y=c()

for (k in 1:1000) {
  for (i in 1:5) {  y[i]=rnorm(1,x[k,i],0.05*abs(x[k,i]))}

  for (i in 1:5) {  if(y[i]<=Par.minn[i]|y[i]>=Par.maxn[i]){y[i]=x[k,i]}}
  
  x_loss <-c()
  x_slow<-c()
  x_loss = 0.05;
  x_slow = 2.3503/(0.05*ExtraF)
  Tx_quick = matrix(0,3,1);
  Q_prediction <-c() 
  t=1;
  while(t<=MaxT+1)
  {
    dailySim <-list()
    dailySim <-hymod_daily(y,Extra,x_loss,x_slow,Tx_quick[1:3,1],t)
    Q_prediction[t]<-dailySim$Qsim
    x_loss<-dailySim$x_loss
    x_slow<-dailySim$x_slow
    Tx_quick[1:3,1]<-matrix(dailySim$x_quick,3,1)
    t = t+1
  }
  
  e[k+1,] = gape(Qobs[51:MaxT],Q_prediction[51:MaxT]);
  Like_e[k+1,] = Likelihood(e[k+1,] ,theta=c(0,sd(e[k+1,])))
  ratio[k+1]=exp(Like_e[k+1,]-Like_e[k,])
  if (runif(1)<ratio[k+1,])
  {x[k+1,]=y
  }else {x[k+1,]=x[k,]}
  
  nash[k+1,]=NSE(Qobs[51:MaxT],Q_prediction[51:MaxT])
  print(k) 
  }


#*********************************AM*************************************
#**********************************************************************
library(MASS) 
library(mvtnorm) 
library(copula) 
library(CDVineCopulaConditional)

timestart<-Sys.time()

y1=y2=c()
for (k in 1000:N) {
  d=5 #d is the number of HYMOD parameters
  if(k%%500 == 0){
    DFN=x[(k-999):k,]
    sigma1=2.38^2/d*(cov(DFN)+10^(-6)*diag(d))#
  }
  
  mean1=c(x[k,1],x[k,2],x[k,3],x[k,4],x[k,5])
  y=rmvnorm(n=1, mean1, sigma1)
 
  for (i in 1:5) {  if(y[i]<=Par.minn[i]|y[i]>=Par.maxn[i]){y[i]=x[k,i]}}
  
  x_loss <-c()
  x_slow<-c()
  x_loss = 0.05;
  x_slow = 2.3503/(0.05*ExtraF)
  Tx_quick = matrix(0,3,1);
  Q_prediction <-c() 
  for (t in 1:MaxT) {
    dailySim <-list()
    dailySim <-hymod_daily(y,Extra,x_loss,x_slow,Tx_quick[1:3,1],t)
    Q_prediction[t]<-dailySim$Qsim
    x_loss<-dailySim$x_loss
    x_slow<-dailySim$x_slow
    Tx_quick[1:3,1]<-matrix(dailySim$x_quick,3,1)
  }
  
  e[k+1,] = gape(Qobs[51:MaxT],Q_prediction[51:MaxT]);
  Like_e[k+1,] = Likelihood(e[k+1,] ,theta=c(0,sd(e[k+1,])))
  ratio[k+1]=exp(Like_e[k+1,]-Like_e[k,])
  if (runif(1)<ratio[k+1,])
  {x[k+1,]=y
  }else {x[k+1,]=x[k,]}
  nash[k+1,]=NSE(Qobs[51:MaxT],Q_prediction[51:MaxT])
  print(k)
}

timesend<-Sys.time()
timesend-timestart 


x_AM=x
nash_AM=nash
library("coda")
#*********************************CopMH*************************************
#**********************************************************************
    
timestart<-Sys.time()
library(CDVineCopulaConditional)
for (k in 1000:N) {
  d=5 #d is the number of HYMOD parameters 
  if (k %% 500 == 0) {
    DFN = x[(k-999):k, ]      
    norm_x = matrix(0, 1000, d) 
    for (t in 1:d) { 
      norm_x[, t] = pnorm(DFN[, t], mean(DFN[, t]), sd(DFN[, t]))
    }

  RVM <- CDVineCondFit(norm_x,Nx = 5,familyset=c(1,3,4,5,6),treecrit = "BIC",type = "CVine",selectioncrit = "AIC")
     }  
  
  
  y_copula=CDVineCondSim(RVM,runif(1,0,1))

  for (t in 1:5) {y[t]=qnorm(y_copula[t],mean(DFN[,t]),sd(DFN[,t])) } 
  timesend<-Sys.time()
  timestart - timesend
  for (i in 1:5) {  if(y[i]<=Par.minn[i]|y[i]>=Par.maxn[i]){y[i]=x[k,i]}}
  x_loss <-c()
  x_slow<-c()
  x_loss = 0.05;
  x_slow = 2.3503/(0.05*ExtraF)
  Tx_quick = matrix(0,3,1);
  Q_prediction <-c() 
  t=1;
  while(t<=MaxT+1)
  {
    dailySim <-list()
    dailySim <-hymod_daily(y,Extra,x_loss,x_slow,Tx_quick[1:3,1],t)
    Q_prediction[t]<-dailySim$Qsim
    x_loss<-dailySim$x_loss
    x_slow<-dailySim$x_slow
    Tx_quick[1:3,1]<-matrix(dailySim$x_quick,3,1)
    t = t+1
  }
  
  e[k+1,] = gape(Qobs[51:MaxT],Q_prediction[51:MaxT]);
  Like_e[k+1,] = Likelihood(e[k+1,] ,theta=c(0,sd(e[k+1,])))
  ratio[k+1]=exp(Like_e[k+1,]-Like_e[k,])
  if (runif(1)<ratio[k+1,])
  {x[k+1,]=y
  }else {x[k+1,]=x[k,]}
  nash[k+1,]=NSE(Qobs[51:MaxT],Q_prediction[51:MaxT])
  print(k) 
}

timesend<-Sys.time()
timestart - timesend

x_CopMH=x
nash_CopMH=nash

##***************************************HCopAM***************************************
##**********************************************************************************************
##HCopAM   illustrated by the AM–CopMH hybrid scheme with a 5:5 sampling proportion
    
timestart<-Sys.time()
y1=y2=c()

(k in 1000:N) {
  d=5
  i              #d is the number of HYMOD parametersf(k%%100 5= 0){
    DFN=x[(k-999):k,]##da sigma1=2.38^2/d*(cov(DFN)+10^(-6)*diag(d))#
   norm_x=matrix(0,1000,5)#??? for (t in 1:5) { norm_x[,t]=pnorm(DFN[,t],mean(DFN[,t]),sd(DFN[,t]))}
    fam1RVM <- CDVineCondFit(norm_x,Nx = 5,familyset=c(1,3,4,5,6),treecrit = "BIC",type = "CVine",selectioncrit = "AIC")
  
  mean1=c(x[k,1],x[k,2],x[k,3],x[k,4],x[k,5])
  y1=rmvnorm(n=1, mean1, sigma1)
  y_r y_copula=CDVineCondSim(RVM,runif(1,0,1))
  
  (t in 1:5) {y2[t]=qnorm(y_copula[t],mean(DFN[,t]),sd(DFN[,t])) } 
  y=0.5*y1+0.5*y2
  #?�or (i in 1:5) {  if(y[i]<=Par.minn[i]|y[i]>=Par.maxn[i]){y[i]=x[k,i]}}
  
  x_loss <-c()
  x_slow<-c()
  x_loss = 0.05;
  x_slow = 2.3503/(0.05*ExtraF)
  Tx_quick = matrix(0,3,1);
  Q_prediction <-c() 
  t=1;
  while(t<=MaxT+1)
  {
    dailySim <-list()
    dailySim <-hymod_daily(y,Extra,x_loss,x_slow,Tx_quick[1:3,1],t)
    Q_prediction[t]<-dailySim$Qsim
    x_loss<-dailySim$x_loss
    x_slow<-dailySim$x_slow
    Tx_quick[1:3,1]<-matrix(dailySim$x_quick,3,1)
    t = t+1
  }
  
  e[k+1,] = gape(Qobs[51:MaxT],Q_prediction[51:MaxT]);
  Like_e[k+1,] = Likelihood(e[k+1,] ,theta=c(0,sd(e[k+1,])))
  ratio[k+1]=exp(Like_e[k+1,]-Like_e[k,])
  if (runif(1)<ratio[k+1,])
  {x[k+1,]=y
  }else {x[k+1,]=x[k,]}
  nash[k+1,]=NSE(Qobs[51:MaxT],Q_prediction[51:MaxT])
  #prit(k) 
}

#???CoHCoppuash_AMCoHCoppuh
##########################################
#wrie.table(nash_AM,"nash_AM_SAC",row.names = FALSE)
write.table(nash_Cop,"naMHsh_Copula_MH",row.names = FALSE)
#wrie.table(nash_AMCoHCoppush_AMCoHCoppu",row.names = FALSE)

#wrie.table(x_AM,"x_AMHymod_parameter__SAC",row.names = FALSE)
write.table(x_Cop,"x_MHCoHymod_parameter_pula_MH",row.names = FALSE)
#wrie.table(x_AMCoHCopAMAMHymod_parameter_CoHCopAM",row.names = FALSE)


#***#
####################Parameter visualization, using Cmax as an example##########################
par(mfrow=c(3,2))
plot(1000:N,x_AM[1000:N,1],ylab="Cmax",xlab="",ylim=c(Par.minn[1],Par.maxn[1]),main="AM",col=1,type="l",cex.lab = 1.5,cex.axis = 1.3,bty = 'l')
plot(1000:N,x_CopMH[1000:N,1],ylab="Cmax",xlab="",ylim=c(Par.minn[1],Par.maxn[1]),main="CopMH",col=5,type="l",cex.lab = 1.5,cex.axis = 1.3,bty = 'l')
plot(1000:N,x_HCopAM[1000:N,1],ylab="Cmax",xlab="",ylim=c(Par.minn[1],Par.maxn[1]),main="HCopAM",col=5,type="l",cex.lab = 1.5,cex.axis = 1.3,bty = 'l')

#####********************************************************************************************#########
####*********************model performance visualization，, using NSE as an example******************####
par(mfrow=c(2,3))
plot(1000:N,nash_AM[1000:N,1],ylab="NSE",ylim=c(0,1),xlab="",main="AM",col=1,type="p",cex.lab = 1.5,cex.axis = 1.3,bty = 'l')
plot(1000:N,nash_CopMH[1000:N,1],ylab="NSE",ylim=c(0,1),xlab="",main="CopHM",col=2,type="p",cex.lab = 1.5,cex.axis = 1.3,bty = 'l')
plot(1000:N,nash_HCopAM[1000:N,1],ylab="NSE",ylim=c(0,1),xlab="",main="HCopAM",col=3,type="p",cex.lab = 1.5,cex.axis = 1.3,bty = 'l')

#########********************************************************************************########
par(mfrow=c(3,5))
autocorr.plot(x_AM[1000:N,],auto.layout = TRUE,cex.lab = 1.7,cex.axis = 1.8,bty = 'l')
autocorr.plot(x_CopMH[1000:N,],auto.layout = TRUE,cex.lab = 1.7,cex.axis = 1.8,bty = 'l')
autocorr.plot(x_HCopAM[1000:N,], auto.layout = TRUE,cex.lab = 1.7,cex.axis = 1.8,bty = 'l')
#****************************************************