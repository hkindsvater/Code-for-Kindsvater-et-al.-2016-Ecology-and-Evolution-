 #this is the code used to generate the plot in Box 1
 library(plyr)
library(fields)

#Model parameters:
Npop=200 #initial pop size (not important)
Tmax=25 #max years
Tfishing=1 #when we start fishing
Trecovery=Tmax #stop fishing 
  
  
   #Life history parameters      	
	
 	  	sizeatbirth=1
    Amax=4 #max age

 	amat=2 #age at maturation - this varies depending on the matration function chosen
	Linf=15 #max size
	fishingsize = .75*Linf #size-selective fishery threshold size 
       Lmat=1 #length at maturity
       k=0.1 # Growth function param
 	  kappa=0.2 #natural mortality coef #params from Lorenzen et al. 2005
	 eps=0.18 #natural mortality exp
	  
	  c=1   #egg prod coef 
	  d=3 #egg &sperm prod exp
	  g = 7000000#sperm prod coef
	  v = 50 #mass-at-size coef
	  wexp = .3  #mass-at-size exp
	  h=1 #sperm fertilization effectiveness

#params that can be scalar or vector:
 natmort=0.2 #instantaneous natural mortality rate 
 alpha = 0.0001#recruitment function param
 beta = 0.000000001#recruitment function param - strenght of density dependence

       q=0.8 #maturation ogive steepness (1 = very steep)
  
    Fvec = seq(0, 1.2, by = 0.3) #Fishing mortality vector
    
   
      
    #set up matrices
N=matrix(nrow=Amax, ncol=Tmax, data=0)
 
W=matrix(nrow=Amax, ncol=1, data=0)
L=matrix(nrow=Amax, ncol=1, data=0)
E=matrix(nrow=Amax, ncol=Tmax, data=0)
S=matrix(nrow=1, ncol=Tmax, data=0)
P=matrix(nrow=1, ncol=Tmax)
B=matrix(nrow=1, ncol=Tmax)
Y=matrix(nrow=1, ncol=Tmax) 
 
 #set up vectors to hold age-specific variables 
Pmale=rep(0, Amax) #probability of sex change
fert=rep(0, Tmax) #fert probability
pmat=rep(0, Amax) #maturation probability
eggs=rep(0, Amax) #fecundity
mu=rep(0, Amax) #natural mortality
select=rep(0, Amax) #fishery selectivity

N[1, 1]=Npop #initial population size
L[1]=sizeatbirth #larval size at recruitment
 sums=matrix(nrow=Tmax, ncol=length(Fvec))

for(Fish in 1:length(Fvec)) {
  
  mu_f=Fvec[Fish]
  
#First: Define the age-specific growth, mortality, and maturation function
  
  for (a in 1:(Amax-1)) {
  	
  	
#GROWTH 
		 
	L[a+1]= Linf*(1-exp(-k)) + (L[a])*exp(-k) #length at age
     W[a]= v*L[a]^wexp #weight at age
      #natural mortality
  
  mu[a]=0.05 #for simplicity, fixed at 0.05
         
     pmat[a]=ifelse( a < amat, 0, 1)     
     
     
 #FISHING SELECTIVITY    
     #JUST AFTER AGE AT MATURATION
     
    select[a] = ifelse( a < amat+1, 0, 1)  #knife edge gear selectivity right after maturation  
     
      

 #FECUNDITY
      # eggs[a]= c*W[a]^2 #fecundity is a function of individual mass-at-age 
       eggs[a]=100000

   } #end 1st age loop

plot(1:Amax, eggs, type="l")


#Now for these age specific rates, simulate pop dynamics through time       
 for(t in 1:(Tmax-1)) {
	 
   
E[,t]=N[, t]*pmat*eggs  #assuming spawning occurs between 1 t and the next
 
P[t]= sum(E[,t])  #assuming fertilization is 100%
  
 N[1,t+1]= alpha*P[t]/(1+beta*P[t]) #this is the N_0 class that is born and enters in the next time step.. 
 
     for (age in 1:(Amax-1)) {
      
   Fishing= if (Tfishing < t & Trecovery > t) {
   				Fishing = select[age]*mu_f } else {
   				Fishing = 0
                     }
    
 	N[age+1,t+1] = N[age, t]*exp(-mu[age]-Fishing)  #assuming mortality happens after spawning

 	B[t]=sum(N[,t]*W)
 	 
 			
    } #end second age loop
        } #end t loop
      sums[,Fish]=apply(N, 2, sum) #calculate popsize (numbers) at each time
        } #end Fish loop
 
 
 #plot steady state for each value of fishing:
 quartz() 
 matplot(1:Tmax, sums/100, type="l", lwd=2, xlab="Time", ylab="", ylim=c(0, 4500), lty = 1, col = c(1, "gray30", "gray50", "gray60", "gray80"), las=1)
     