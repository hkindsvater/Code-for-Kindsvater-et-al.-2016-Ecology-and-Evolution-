
#This code can be modified to make Figure 1. 

#Simulating a population with logistic growth
Tmax <- 80
  #K <- 10000
    N0 <- 1

 
 a = .001 #this is the proportion by which r declines as a function of n - IE the biology it captures is the effect of crowding on births/deaths. 
 
 
 #For upper line on Fig 1 inset = a= 0.001 and fec = 5 and mu = 4.8; for all other lines N2-N4 a=0.002)
 
 
 #if a changes, the "K" changes a lot, IE doubling 'a' halves the population size. 

# varying r (population growth rate ): for the baseline case (N2) fec = 5 and mu = 4.8; for N3 (increased births) fec = 5.1; for N4 (increased deaths) mu = 4.9
 
 N <- matrix(nrow = Tmax, ncol = length(a), data = N0)# set up matrix to be filled
B <- matrix(nrow = Tmax, ncol = length(a), data = 0)# set up matrix to be filled
fec=5.0 
mu=4.8  
r=fec-mu
i=1
 	
	 for (t in 1:(Tmax - 1))  {
	 	 
		 r=fec-mu 
		K=r[i]/a
		# N[t+1,i] <- N[t]+r*N[t,i]*(1-(N[t,i]/K)) 
		
		 N[t+1,i] = K*N[t,i] / (N[t,i] + (K-N[t,i])*exp(-r[i]*t)  )
		 
		} #end t loop

N1=N[, i]

Nvec=1:max(N[,i])
 dN2= r*Nvec - ((r/K)*Nvec^2) 
  
##FIGURE 1 PANELS  
  
 plot(1:max(N1), dN1, type="l", ylab ="", xlab = "Population Size (N)", lwd = 3, ylim=c(0, max(dN3)+3), las=1)  
  lines(dN2, lwd=2)
  lines(dN3, lwd=1)
  lines(dN4)

quartz()
plot(1:50, N1[1:50], type="l", lwd=3, xlab="", ylab="", ylim=c(0, max(N1)+20), xlim=c(0, Tmax-25),las =1)
   lines(N2[1:50], lwd=2)
  lines(N3[1:50])
  lines(N4[1:50])
 
 
  