##this code is used to calculte the various metrics of reproductive value discussed in Box 2 and Appendix 3, including relative fitness (plotted in figure 4). The values for Yellowfin Tuna are included here as an example. 

#rm(list=ls(all=TRUE))
library(reshape)
  library(plyr)
library(fields)

Tmax = 500 
Npop = 1
   #Life history parameters      	
 
 #  NATURAL MORTALITY - AS A FUNCTION OF AGE (Data from Hampton 2000)
	mu=c(6.1, 1.5, 0.68, .44, .69, 1.5, 1.5, 1.5, 1.5)
	Amax=length(mu) #max age
 	amat=3.5 #age at maturation - this varies depending on the matration function chosen
	Linf=180 #max size
	k=0.4 # Growth function param
    c= 10  #egg prod coef 1 mod; 10 very high, .004 shark     
  q=2 #maturation ogive steepness (1 = very steep) #DOESN'T VARY
  recruitsize=4 #DOESN'T VARY 

#FECUNDITY ASSUMED TO BE PROPORTIONAL TO MASS: 
	  
	  d=3 #egg &sperm prod exp
	  v = 50 #mass-at-size coef
	  wexp = .3  #mass-at-size exp
	  h=1 #sperm fertilization effectiveness
	  
	        #recruitment
 alpha = 0.002 #recruitment function param
 beta = 0.0000000001#recruitment function param - strenght of density dependence
 
   #MATRICES
  N=matrix(nrow=Amax, ncol=Tmax, data=0)
E=matrix(nrow=Amax, ncol=Tmax, data=0)
S=matrix(nrow=1, ncol=Tmax, data=0)
P=matrix(nrow=1, ncol=Tmax)

  #set up vectors to hold age-specific variables 

pmat=rep(0, Amax) #maturation probability
eggs=rep(0, Amax) #fecundity
W=matrix(nrow=Amax, ncol=1, data=1)
L=matrix(nrow=Amax, ncol=1, data=1)


N[1, 1]=Npop #initial population size
L[1]=recruitsize #larval size at recruitment
 
 #First: Define the age-specific growth, mortality, and maturation function
  
  for (a in 1:(Amax-1)) {
  	
#GROWTH 
		 
	L[a+1]= Linf*(1-exp(-k)) + (L[a])*exp(-k) #length at age
     W[a]= v*L[a]^wexp #weight at age
     
  
###  AGE-DEP MATURATION:
       pmat[a]= 1/(1+exp(-q*(a-amat))) 
     #pmat[a]=ifelse( a < amat, 0, 1)     #knife edge maturation
     

 #FECUNDITY
       eggs[a]= c*W[a]^2 #fecundity is a function of individual mass-at-age 
 
   } #end 1st age loop

   #NOW FOLLOW COHORT THRU TIME: 
       
 for(t in 1:(Tmax-1)) {
	 
   
E[,t]=N[, t]*pmat*eggs  #assuming spawning occurs between 1 t and the next
 
P[t]= sum(E[,t])  #assuming fertilization is 100%
  
 N[1,t+1]= alpha*P[t]/(1+beta*P[t]) #this is the N_0 class that is born and enters in the next time step.
 
     for (age in 1:(Amax-1)) {
      
   	N[age+1,t+1] = N[age, t]*exp(-mu[age])
 	#assuming mortality happens after spawning
		
     } #end second age loop
    
   
       } #end t loop
 
 #now calculate different metrics of RV: 
 
 surv_a=cumprod(exp(-mu))  
 RV_fisher = sum(eggs*pmat/surv_a)
  plot(1:(Amax), RV_fisher  )
 
   
   
EP_a=(E[-Amax,t]*pmat[-Amax]) /sum(E[-Amax,t]*pmat[-Amax])   
   RV_ham=rep(0, Amax)
   RV_hamperage=rep(0, Amax)
   RV_fisher=rep(0, Amax)
   
   for (i in Amax : 2) {
   	RV_fisher[i-1] = sum(eggs[i]*pmat[i], eggs[i-1]*pmat[i-1])/surv_a[i] 
   	  
   	   RV_ham[i-1]=RV_ham[i]+EP_a[i-1]
   	
   	     }
   	   
   	   CF_yellowfin=E[-Amax, t]*pmat[-Amax]/sum(E[-Amax,t]*pmat[-Amax])   
  
  	     
 #plotting  	     
   	quartz() 
   	par(mfrow=c(3, 1))    
   	 plot(1:(Amax-1),RV_fisher[-Amax], type="l", ylab="V(a)", xlab="Age", las=1, main = "Yellowfin tuna")    
	 plot(1:(Amax-1), RV_ham[-Amax], type="l", ylab="W(a)", xlab="Age", las=1)  
	 plot(1:(Amax-1), CF_yellowfin, type="l", ylab = "w(a)", xlab="Age", las=1)   
 
 quartz()
 par(mfrow=c(3,2))
   plot(1:Amax, L, type="l", main="von Bert growth", ylab="Size", xlab="Age")
    legend("bottomright", legend=c( paste("Linf = ", Linf), paste("von Bert k = ", k)))

   plot(1:(Amax-1), eggs[-Amax]*pmat[-Amax], type="l", main="Age-specific Fecundity", ylab="Number of female progeny", xlab= "Age")
    legend("bottomright", legend=c(paste("Amat = ", amat)))

  plot(1:(Amax-1), pmat[-Amax], type="l", main="Maturity", ylab="Probability Mature", xlab= "Age")

   plot(1:(Amax-1), exp(-mu[-Amax]), type="l", main= "Age-specific survival", ylab = "Annual survival", xlab="Age")
   
   plot(1:Tmax, colSums(N), type="l", main = "Population Dynamics", ylab = "Abundance", xlab = "Time")
   
   plot(1:(Amax-1),  CF_yellowfin, type="l", ylab = "Current Fitness", main = "Relative Fitness at each age", xlab="Age")
 