###############################################################################
#    Metapopulation regulation acts at multiple spatial scales:               #
#     insights from a century of seabird colony census data                   #
#                                                                             #
#  Jana W.E. Jeglinski, Sarah Wanless, Stuart Murray, Robert T. Barrett,      #
#  Arnthor Gardarsson, Mike P. Harris, Jochen Dierschke, Hallvard Strøm,      #
#        Svein-Håkon Lorentsen and Jason Matthiopoulos                        #
#                                                                             #
#               JAGS code for gannet metapopulation model                     #       
###############################################################################

##### required R libraries  ####----

library(runjags)
library(coda)

#### Data requirements  ####----

# Data is available on Enlighten under  http://dx.doi.org/10.5525/gla.researchdata.1286

# Pd: matrix of colony census data for each colony (columns) and year (rows) of the time series, for forecasting add empty 
# rows 
# harvest: matrix with harvest information or colonies and years of time series, expressed as 0 = no harvest, 1 = harvest
# reg: a vector of regions
# cc: a vector of terrestrial carrying capacity estimates based on expert judgement
# nYrs: number of years of time series
# N: number of colonies


#### Model statement ####----

# Note:
# The model contains several switches, i.e. two subsequent code lines that each represent alternative scenarios
# of metapopulation regulation
# Un-comment the respective line of code to run a particular scenario
#
# We investigate 5 different regulatory scenarios with this model:
# scenario 1 - closed populations with local density dependence
# scenario 2: metapopulation with local density dependence and equipartitioning of immigrants
# scenario 3: metapopulation with local density dependence and conspecific attraction in immigration term
# scenario 4: metapopulation with local or regional density dependence and equipartitioning of immmigrants
# scenario 5: metapopulation with local or regional density dependence and conspecific attraction in immigration term
# the model code below is un-commented to run the most complex scenario (model 5)


popMod <- "model{

### MAIN LOOP ###

for(n in 1:N)
{
  for(t in 5:(nYrs-1))
  {
  
    ### PROCESS MODEL ###
    
    logit(r[n,t])<-100-max(a1[reg[n]],a2[n])*P[n,t]                     # Switch 1: Recruitment probability dampened by  regional or local density dependance
    #logit(r[n,t])<-100-(a2[n]*P[n,t])                                  # Switch 1: Recruitment probability dampened by local density dependence
    
    
    logit(b[n,t])<-eps[t]                                               # Fecundity with effect of annual environmental perturbations
    y[n,t]<-sj*b[n,t-4]*P[n,t-4]*(aH^harvest[n,t-4])                    # state: Number of young at colony n and year t, as function of harvest (proportion), fecundity (rate) and survival (proportion)
    

    #w[n,t]<-min(1,P[n,t])                                              # Switch 2: equal redistribution to all existing colonies
    w[n,t]<-P[n,t]                                                      # Switch 2 : redistribution according to size of receiving colony/ conspecific attraction
    W[n,t]<-w[n,t]/wtot[t]                                              # redistribution function weighted by w[n,t]
    
    R[n,t]<-r[n,t]*((1-aI)*y[n,t]+aI*ytot[t]*W[n,t])                    # state: number of recruits at colony n in year t
    
    lambda[n,t]<-sa*P[n,t]+R[n,t]                                       # growth rate based on adult and juvenile survival, recruitment and immigration

    #P[n,t+1]~dpois(lambda[n,t])                                        # Switch 3: poisson form of growth model
    P[n,t+1]~dpois(lambda[n,t] * h[n,t])                                # Switch 3: negative binomial form of growth model
    h[n,t] ~ dgamma(theta,theta)                                        # Switch 3: required for negative binomial
    
    
    ### OBSERVATION MODEL ###
    Pd[n,t+1]~dnorm(P[n,t+1],1/(0.05*P[n,t+1]+1)^2)                     # Accounting for 10 % observation error in census data
    
  }
}

for (t in 5:(nYrs-1))
{
  ytot[t]<-sum(y[,t])                                       # Total number of young in metapop in year t
  wtot[t]<-sum(w[,t])                                       # Total number of extant colonies in year t (needed for distribution of emigrants)
 
}




for(n in 1:N) {R[n,4]<-1}                      # equal recruitment for all colonies at time step 4 to generate last year R as starting value

## Colony size vectors for monitoring

 P1<-P[1,1:nYrs]
 P2<-P[2,1:nYrs]
 P3<-P[3,1:nYrs]
 P4<-P[4,1:nYrs]
 P5<-P[5,1:nYrs]
 P6<-P[6,1:nYrs]
 P7<-P[7,1:nYrs]
 P8<-P[8,1:nYrs]
 P9<-P[9,1:nYrs]
 P10<-P[10,1:nYrs]
 P11<-P[11,1:nYrs]
 P12<-P[12,1:nYrs]
 P13<-P[13,1:nYrs]
 P14<-P[14,1:nYrs]
 P15<-P[15,1:nYrs]
 P16<-P[16,1:nYrs]
 P17<-P[17,1:nYrs]
 P18<-P[18,1:nYrs]
 P19<-P[19,1:nYrs]
 P20<-P[20,1:nYrs]
 P21<-P[21,1:nYrs]
 P22<-P[22,1:nYrs]
 P23<-P[23,1:nYrs]
 P24<-P[24,1:nYrs]
 P25<-P[25,1:nYrs]
 P26<-P[26,1:nYrs]
 P27<-P[27,1:nYrs]
 P28<-P[28,1:nYrs]
 P29<-P[29,1:nYrs]
 P30<-P[30,1:nYrs]
 P31<-P[31,1:nYrs]
 P32<-P[32,1:nYrs]
 P33<-P[33,1:nYrs]
 P34<-P[34,1:nYrs]
 P35<-P[35,1:nYrs]
 P36<-P[36,1:nYrs]
 P37<-P[37,1:nYrs]
 P38<-P[38,1:nYrs]
 P39<-P[39,1:nYrs]
 P40<-P[40,1:nYrs]
 P41<-P[41,1:nYrs]
 P42<-P[42,1:nYrs]
 P43<-P[43,1:nYrs]
 P44<-P[44,1:nYrs]
 P45<-P[45,1:nYrs]
 P46<-P[46,1:nYrs]
 P47<-P[47,1:nYrs]
 P48<-P[48,1:nYrs]
 P49<-P[49,1:nYrs]
 P50<-P[50,1:nYrs]
 P51<-P[51,1:nYrs]
 P52<-P[52,1:nYrs]
 P53<-P[53,1:nYrs]


### PRIORS ###

## adult survival
sa0~dbeta(12,12)                 # min and max values based on Wanless et al. (2006)
samax<-0.918+(2*0.023)
samin<-0.918-(2*0.023)
sa<-samin+sa0*(samax-samin)

## immature survival
sj0~dbeta(12,12)                 # min and max values based on Wanless et al. (2006)
sjmax<-0.279+(2*0.05)
sjmin<-0.279-(2*0.05)
sj<-sjmin+sj0*(sjmax-sjmin)

## immigration
aI0~dbeta(12,12)                 # min and max values based on BTO ring recovery data
aImax<-0.52                      
aImin<-0.32
aI<-aImin+aI0*(aImax-aImin)

#aI <- 0                           # set aI to zero when running model for closed populations

## harvest
aH~dbeta(1,1)      

## negative binomial growth
 theta <- max(10000-theta0,1000)     # shrinkage tendency to poisson
 theta0 ~ dexp(1/200) 

## regional density dependance
# for scenarios without regional density dependence comment loop out (required for switch 1)

for(nr in 1:nreg) 
{
  
  a1[nr]~dgamma(6.25,125)                        # a1 is regional density dependent parameter
  K1[nr]<-1/a1[nr]*(100-log(rstar/(1-rstar)))    # K1 is regional carrying capacity
  
}


## local density dependance

a0<--0.58                          # mean fecundity (half for each pair) based on JNCC breeding success data 
logit(bstar)<-a0                   # Baseline fecundity
rstar<-(1-sa)/sj*bstar             # Constant related to demographic replacement for steady state

for(n in 1:N) 
{
  K0[n]~dbeta(9.2,13.8)           
  Kmax[n]<-1.2*cc[n]
  Kmin[n]<-0.8*cc[n]
  K[n]<-Kmin[n]+K0[n]*(Kmax[n]-Kmin[n])      # K is local carrying capacity
  
  a2[n]<-1/K[n]*(100-log(rstar/(1-rstar)))   # a2 is local density dependent parameter
  
  }

### end of PRIORS ###


# stochastic annual fecundity fluctuations in each year 

  for(t in 1:(nYrs))
  {
    eps[t]~dnorm(a0,15.4)           #integrates mean fecundity and standard deviation based on long term data from 13 colonies
  }


### INITIALIZATION ###

for(n in 1:N){
  for (t in 1:4){
    y[n,t]~dunif(0,5000)
    logit(b[n,t])<-eps[t] 
    w[n,t]<-min(1,P[n,t])         #using simplest form of immigrant distribution for initialisation
    }
}

for (n in 1:N){

  for (t in 1:5){
    P[n,t]~dunif(0,5000)
    Pd[n,t]~dnorm(P[n,t],1/(0.05*P[n,t]+1)^2)
    
    }
}




#data# nYrs,N,Pd,harvest,reg,cc

#monitor#  K,K1,a1,a2,aH,sa,sj,aI

#monitor#  P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,P23,P24,P25,P26,P27

#monitor# P28,P29,P30,P31,P32,P33,P34,P35,P36,P37,P38,P39,P40,P41,P42,P43,P44,P45,P46,P47,P48,P49,P50,P51,P52,P53
  
  }"
  


