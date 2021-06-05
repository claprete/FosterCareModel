# FARRGAP
# Foster Care Model
# Carolyn LaPrete
# 12/08/2020

####################################################################################
# Changed rates for 1.5 yr (pandemic) -> dump into system (post-pandemic)
####################################################################################

# Load a package that solves differential equations
require(deSolve)

####################################################################################
# Start with graphing no pandemic
####################################################################################

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

sigmah <- Nef/Npop        #rate of P->Rh
rhoh <- 0.4*Nxf/Nf        #rate of Rh->P

ptors <-0.32              #percentage in kinship care, used for sigmas and initial conditions
sigmas <- ptors*Nef/Npop  #rate of P->Rs (times rate of leaving P - entering system)
rhos <- 0.47*Nxf/Nf       #rate of Rs->P  (times rate of leaving system)

ptog <- 0.68              #percentage in foster/group, used for gamma and initial conditions
gamma <- ptog*Nef/Npop    #rate of P->G (times rate of leaving P - entering system)
epsilon <- 0.92*Nxf/Nf    #rate of G->P  (times rate of leaving system)

miur <- 0.32*Nxf/Nf       #rate of Rs->A  (times rate of leaving system)
miuf <- 0.32*Nxf/Nf       #rate of F->A  (times rate of leaving system)
phi <- 0.68               #rate of G->F
w <- 0.51                 #rate of F->G (0.51 for eq in 1 year, 0.68 for eq in many)


# Define times to solve
tmax <- 10*365
times <- seq(from=0,to=tmax,by=0.1)


#Define initial condition and give it the right name
Rhinit <- 250000
Rsinit <- ptors*Nf
Ginit <- ptog*(1-phi)*Nf
Finit <- ptog*phi*Nf
Ainit <- 0
Pinit <- Npop-Nf
init <- c(Pinit,Rhinit,Rsinit,Ginit,Finit,Ainit)
names(init) <- c("P","Rh","Rs","G","F","A")

# Solve the equation numerically
FARRGAPout <- as.data.frame(ode(y = init,times = times, func = rhs))

print(tail(FARRGAPout,1))

# Start graph with base
par(mfrow=c(1,1),mar=c(5,6,4,2))
plot(F ~ time,FARRGAPout,type="l",lty=1,lwd=2, ylim=c(75000,300000),
  ylab="Number of Kids",xlab="Time (days)",
     cex.axis=1,cex.lab=1, col="light green")
lines(Rs ~ time,FARRGAPout,lty=1,lwd=2,col="light blue",
     cex.axis=1,cex.lab=1)
lines(G ~ time,FARRGAPout,lty=1,lwd=2,col="pink",
     cex.axis=1,cex.lab=1)

####################################################################################
# Start of new code - re run for pandemic changed rates and post-pandemic surge
####################################################################################

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
  dP <- P*(-sigmah(t)-sigmas(t)-gamma(t))+Rh*(rhoh+alpha)+Rs*(rhos(t)+alpha)+G*(epsilon(t)+alpha)+F*alpha+A*alpha
  dRh <- P*sigmah(t)+Rh*(-rhoh-alpha)
  dRs <- P*(sigmas(t))+Rs*(-rhos(t)-alpha-miur)
  dG <- P*gamma(t)+G*(-phi(t)-alpha-epsilon(t))+F*w(t)
  dF <- G*phi(t)+F*(-w(t)-miuf-alpha)
  dA <- Rs*miur+F*miuf-A*alpha 
  return(list(c(dP,dRh,dRs,dG,dF,dA)))   #For some reason, this has to be a "list"
}

# Define parameters
Npop <- 7.3e+07           #population of kids in the US under 18
Nf <- 424000              #population of kids in foster care
Nef <- function(t) {(t<=365)*(250000/365)+(t>365 & t<7*365)*((-25000*sin((2*pi*t)/(6*365))+250000)/365)+(t>=7*365)*(250000/365)}
                          #population of kids entering foster care per day = same for hidden
Nxf <- 248000/365         #population of kids exiting foster care per day
alpha <- 1/(18*365)       #rate of anything->Ao per day

sigmah <- function(t) {(t<=365)*(1)+(t>365 & t<2.5*365)*(1.3)+(t>=2.5*365 & t<3.5*365)*(0.15*sin((pi/365)*t)+1.15)+(t>=3.5*365)*(1)}*Nef(t)/Npop
                          #rate of P->Rh
rhoh <- 0.4*Nxf/Nf        #rate of Rh->P

ptors <- function(t) {(t<=365)*(0.32)+(t>365 & t<2.5*365)*(0.17)+(t>=2.5*365 & t<3.5*365)*(0.075*-sin((pi/365)*t)+0.245)+(t>=3.5*365)*(0.32)}
                          #percentage in kinship care, used for sigmas and initial conditions
sigmas <- function(t) ptors(t)*Nef(t)/Npop
                          #rate of P->Rs (times rate of leaving P - entering system)
rhos <- function(t) {(t<=365)*(0.47)+(t>365 & t<2.5*365)*(0.27)+(t>=2.5*365 & t<3.5*365)*(0.1*-sin((pi/365)*t)+0.37)+(t>=3.5*365)*(0.47)}*Nxf/Nf
                          #rate of Rs->P  (times rate of leaving system)

ptog <- function(t) {(t<=365)*(0.68)+(t>365 & t<2.5*365)*(0.53)+(t>=2.5*365 & t<3.5*365)*(0.075*-sin((pi/365)*t)+0.605)+(t>=3.5*365)*(0.68)}
                          #percentage in foster/group, used for gamma and initial conditions
gamma <- function(t) ptog(t)*Nef(t)/Npop
                          #rate of P->G (times rate of leaving P - entering system)
epsilon <- function(t) {(t<=365)*(0.92)+(t>365 & t<2.5*365)*(0.44)+(t>=2.5*365 & t<3.5*365)*(0.24*-sin((pi/365)*t)+0.68)+(t>=3.5*365)*(0.92)}*Nxf/Nf 
                          #rate of G->P  (times rate of leaving system)

miur <- 0.32*Nxf/Nf       #rate of Rs->A  (times rate of leaving system)
miuf <- 0.32*Nxf/Nf       #rate of F->A  (times rate of leaving system)
phi <- function(t) {(t<=365)*(0.68)+(t>365 & t<2.5*365)*(0.51)+(t>=2.5*365 & t<3.5*365)*(0.075*-sin((pi/365)*t)+0.585)+(t>=3.5*365)*(0.68)}
                          #rate of G->F
w <- function(t) {(t<=365)*(0.51)+(t>365 & t<2.5*365)*(0.61)+(t>=2.5*365 & t<3.5*365)*(0.05*sin((pi/365)*t)+0.56)+(t>=3.5*365)*(0.51)}
                          #rate of F->G


#Define initial condition and give it the right name
Rhinit <- 250000
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


# Continue graphing with surge
lines(F ~ time,FARRGAPout,lty=1,lwd=2,col="dark green",
     cex.axis=1,cex.lab=1)
lines(Rs ~ time,FARRGAPout,lty=1,lwd=2,col="blue",
     cex.axis=1,cex.lab=1)
lines(G ~ time,FARRGAPout,lty=1,lwd=2,col="red",
     cex.axis=1,cex.lab=1)

legend("topright",c("(Pre) Foster Care","(Pre) Kinship Care (in-system)","(Pre) Group Homes (and floaters)",
  "Foster Care","Kinship Care (in-system)","Group Homes (and floaters)"),cex=0.7,
      col=c("light green","light blue","pink","dark green","blue","red"),
      lty=c(1,1,1,1,1,1),lwd=c(2,2,2,2,2,2))
title(main="The Foster System\n(Pandemic Rates to Post-Pandemic Surge)")

