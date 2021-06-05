# FARRGAP
# Foster Care Model
# Carolyn LaPrete
# 12/08/2020

####################################################################################
# Rate changes (during pandemic) - changed rates for 1.5 years
####################################################################################

# Load a package that solves differential equations
require(deSolve)

# If you don't have this package on your laptop, try this!
# install.packages("deSolve")

# This defines the differential equations
rhs <- function(t,x,parms){
# t is time, x is the vector of state variables, parms is an annoying
# vector of parameter that I ignore
  P <- x["P"]  #Pull P out of the x vector
  Rh <- x["Rh"]  #Pull Rh out of the x vector
  Rs <- x["Rs"]  #Pull Rs out of the x vector
  G <- x["G"]  #Pull G out of the x vector
  F <- x["F"]  #Pull F out of the x vector
  A <- x["A"]  #Pull A out of the x vector
  dP <- P*(-sigmah-sigmas-gamma)+Rh*(rhoh+alpha)+Rs*(rhos+alpha)+G*(epsilon+alpha)+F*alpha+A*alpha
  dRh <- P*sigmah+Rh*(-rhoh-alpha)
  dRs <- P*(sigmas)+Rs*(-rhos-alpha-miur)
  dG <- P*gamma+G*(-phi-alpha-epsilon)+F*w
  dF <- G*phi+F*(-w-miuf-alpha)
  dA <- Rs*miur+F*miuf-A*alpha 
  return(list(c(dP,dRh,dRs,dG,dF,dA)))
}

# Define parameters
Npop <- 7.3e+07           #population of kids in the US under 18
Nf <- 424000              #population of kids in foster care
Nef <- 251000/365         #population of kids entering foster care per day = same for hidden
Nxf <- 248000/365         #population of kids exiting foster care per day
alpha <- 1/(18*365)       #rate of anything->Ao per day

sigmah <- Nef*1.3/Npop    #rate of P->Rh
rhoh <- 0.4*Nxf/Nf        #rate of Rh->P

ptors <-0.17              #percentage in kinship care, used for sigmas and initial conditions
sigmas <- ptors*Nef/Npop  #rate of P->Rs (times rate of leaving P - entering system)
rhos <- 0.27*Nxf/Nf       #rate of Rs->P  (times rate of leaving system)

ptog <- 0.53              #percentage in foster/group, used for gamma and initial conditions
gamma <- ptog*Nef/Npop    #rate of P->G (times rate of leaving P - entering system)
epsilon <- 0.44*Nxf/Nf    #rate of G->P  (times rate of leaving system)

miur <- 0.32*Nxf/Nf       #rate of Rs->A  (times rate of leaving system)
miuf <- 0.32*Nxf/Nf       #rate of F->A  (times rate of leaving system)
phi <- 0.51               #rate of G->F
w <- 0.61                 #rate of F->G


# Define times to solve
tmax <- 1.5*365
times <- seq(from=0,to=tmax,by=0.1)


#Define initial condition and give it the right name
Rhinit <- Nef*365
Rsinit <- 0.32*Nf         #originally ptors*Nf but ptors changed - use old one, goes with data pre pandemic
Ginit <- 0.68*(1-0.68)*Nf #originally ptog*(1-phi)*Nf
Finit <- 0.68*0.68*Nf     #originally ptog*phi*Nf
Ainit <- 0
Pinit <- Npop-Nf
init <- c(Pinit,Rhinit,Rsinit,Ginit,Finit,Ainit)
names(init) <- c("P","Rh","Rs","G","F","A")

# Solve the equation numerically
FARRGAPout <- as.data.frame(ode(y = init,times = times, func = rhs))

print(tail(FARRGAPout,1))


# Graph (during pandemic)
par(mfrow=c(1,1),mar=c(5,6,4,2))
plot(F ~ time,FARRGAPout,type="l",lty=1,lwd=2, ylim=c(2000,560000),
	ylab="Number of Kids",xlab="Time (days)",
     cex.axis=1,cex.lab=1, col="dark green")
lines(Rs ~ time,FARRGAPout,lty=1,lwd=2,col="blue",
     cex.axis=1,cex.lab=1)
lines(Rh ~ time,FARRGAPout,lty=1,lwd=2,col="light blue",
     cex.axis=1,cex.lab=1)
lines(G ~ time,FARRGAPout,lty=1,lwd=2,col="red",
     cex.axis=1,cex.lab=1)
lines(A ~ time,FARRGAPout,lty=1,lwd=2,col="purple",
     cex.axis=1,cex.lab=1)

legend("topleft",c("Foster Care","Kinship Care (in-system)",
	"Kinship Care (hidden)","Group Homes (and floaters)", "Adopted"),cex=0.7,
      col=c("dark green","blue","light blue","red", "purple"),
      lty=c(1,1,1,1,1),lwd=c(2,2,2,2,2))
title(main="The Foster System (During 1.5 Year Pandemic)")

