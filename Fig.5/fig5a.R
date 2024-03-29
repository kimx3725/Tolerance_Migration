# ODE cycles - equations 
  
# load libraries
#install.packages(pacman)
pacman::p_load(deSolve,reshape2,akima,lattice,ggplot2,dplyr,tidyr,gridExtra,reshape2,envalysis)


# First 6 months - Enivornment A together (both M and R)

# First 6 months (Season 1) Envir A 
open.sir.model1 <- function (t, x, params) {
  
  # first extract the state variables
  X1 <- x[1] # Resident Susceptible Individual
  X2 <- x[2] # Migrant Susceptible Individual
  Y1 <- x[3] # Resident Infected Individual
  Y2 <- x[4] # Migrant Infected Individual
  
  # Parameters
  beta <- params["beta"] #Transmission rate
  mu <- params["mu"] # Death rate
  gamma <- params["gamma"] # Recovery rate
  alpha <- params["alpha"] # Parasite Induced Death rate
  
  # system of differential equations for SI model
  dX1dt <- -beta*x[1]*(x[3]+x[4])+(gamma*x[3])-(mu*x[1]) # dfq of x1
  dX2dt <- -beta*x[2]*(x[3]+x[4])+(gamma*x[4])-(mu*x[2]) # dfq of x2
  dY1dt <- beta*x[1]*(x[3]+x[4])-(alpha+gamma+mu)*x[3] # dfq of y1
  dY2dt <- beta*x[2]*(x[3]+x[4])-(alpha+gamma+mu)*x[4] # dfq of y2
  list(c(dX1dt,dX2dt,dY1dt,dY2dt))
}

# Specify parameter rates 
times1 <- seq(from=0,to=6) # Time period 1 yr


# Second 6 months - Enivornment A together (Only R)

# First 6 months (Season 1) Envir A 
open.sir.model2 <- function (t, x, params) {
  
  # first extract the state variables
  X1 <- x[1] # Resident Susceptible Individual
  Y1 <- x[2] # Resident Infected Individual
  
  # Parameters
  beta <- params["beta"] #Transmission rate
  mu <- params["mu"] # Death rate
  gamma <- params["gamma"] # Recovery rate
  alpha <- params["alpha"] # Parasite Induced Death rate
  
  # system of differential equations for SI model
  dX1dt <- -beta*x[1]*(x[2])+(gamma*x[2])-(mu*x[1]) # dfq of x1
  dY1dt <- beta*x[1]*(x[2])-(alpha+gamma+mu)*x[2] # dfq of y1
  list(c(dX1dt,dY1dt))
}

# Specify parameter rates 
times2 <- seq(from=0,to=6) # Time period 1 yr


# Second 6 months - Environment B (Only M)

# First 6 months (Season 1) Envir B 
open.sir.model3 <- function (t, x, params) {
  
  # first extract the state variables
  X2 <- x[1] # Migrant Susceptible Individual
  Y2 <- x[2] # Migrant Infected Individual
  
  # Parameters
  beta <- params["beta"] #Transmission rate
  mu <- params["mu"] # Death rate
  gamma <- params["gamma"] # Recovery rate
  alpha <- params["alpha"] # Parasite Induced Death rate 
  
  # system of differential equations for SI model
  dX2dt <- -beta*x[1]*(x[2])+(gamma*x[2])-(mu*x[1]) # dfq of x2
  dY2dt <- beta*x[1]*(x[2])-(alpha+gamma+mu)*x[2] # dfq of y2
  list(c(dX2dt,dY2dt))
}

# Specify parameter rates 
times3 <- seq(from=0,to=6) # Time period 6 months 


# Specify the parameter values 

# initial parameter inputs
#beta <- 0.01 # tramsmission rate
gammaA <- 0.3 # recovery rate in environment A
gammaB <- 0.7 # recovery rate in environment B 
mu <- 0.01 # mortality rate
#alpha <- 0.2 # parasite induced death rate
phi=3 # maximum per capita fecundity rate for each individual
z=0.0001# density dependent fecundity coefficient 
cm=0.1 # costs of migration
#ct=0.5 # costs of tolerance



# 100 years loop 
cTbetaoutput <- vector(mode="list", length=121) # create 20 list dummy vectors for both migrants and residents
alphaoutput <- vector(mode="list", length=121) # create 20 list dummy vectors for both migrants and residents
betaoutput <- vector(mode="list", length=121) # create 20 list dummy vectors for both migrants and residents
betavals <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) # beta variations
alphavals <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) # alpha variations
for (i in 1:length(alphavals)){
  
  for (g in 1:length(betavals)){
    
    X10 <- 150 # initial susceptible population (Resident) 
    X20 <- 150 # initial susceptible population (Migrant)
    Y10 <- 100 # initial infected population (Resident)
    Y20 <- 100 # initial infected population (Migrant)
    for (t in 1:100){
      # Environment A (M+R) together - first 6 months
      params1 <- c(beta=betavals[g],gamma=gammaA,mu=mu,alpha=alphavals[i]) # parameters
      xstart1 <- c(X1=X10,X2=X20,Y1=Y10,Y2=Y20) # initial conditions
      cTbeta1 <- as.data.frame(ode(func = open.sir.model1,y=xstart1,times=times1,parms=params1,method="ode23")) 
      #method = "ode45",atol=1e-14,rtol=1e-14
      
      # Environment A (only R) - Second 6 months
      params2 <- c(beta=betavals[g],gamma=gammaA,mu=mu,alpha=alphavals[i]) # parameters
      xstart2 <- c(X1=cTbeta1[7,2],Y1=cTbeta1[7,4]) # initial conditions
      cTbeta2 <- as.data.frame(ode(func = open.sir.model2,y=xstart2,times=times2,parms=params2,method="ode23"))
      
      # Environment B for migrants - Second 6 months 
      params3 <- c(beta=betavals[g],gamma=gammaB,mu=mu,alpha=alphavals[i]) 
      xstart3 <- c(X2=cTbeta1[7,3],Y2=cTbeta1[7,5]) # Initial survived population from env A
      cTbeta3 <- as.data.frame(ode(func = open.sir.model3,y=xstart3,times=times3,parms=params3,method="ode23"))
      
      # the annual survivors 
      
      X1hat <- cTbeta2[7,2] # resdient susceptible
      X2hat <- cTbeta3[7,2] # migrant susceptible
      Y1hat <- cTbeta2[7,3] # resident infected
      Y2hat <- cTbeta3[7,3] # migrant susceptible
      
      cT <- 1-alphavals[i]
      
      bRmax <- phi*(X1hat+Y1hat)*(1-cT) # maximum fecundity rate of a resident individual 
      bMmax <- phi*(X2hat+Y2hat)*(1-cm)*(1-cT) # maximum fecundity rate of a migrant individual 
      
      b1 <- bRmax*(exp(-z*(bRmax+bMmax))) # actual fecundity rate of a resident individual
      b2 <- bMmax*(exp(-z*(bRmax+bMmax))) # actual fecundity rate of a migrant individual 
      
      # after fecundity - 2nd yr initial populations 
      X10 <- cTbeta2[7,2] + b1
      X20 <- cTbeta3[7,2] + b2
      Y10 <- cTbeta2[7,3]
      Y20 <- cTbeta3[7,3]
    }
    cTbetaresult <- data.frame(X10,X20,Y10,Y20)
    
    cTbetaoutput[[(g-1)*11+i]]<-as.data.frame(cTbetaresult)
    betaoutput[[(g-1)*11+i]]<-as.data.frame(betavals[g])
    alphaoutput[[(g-1)*11+i]]<-as.data.frame(alphavals[i])}} # CT loop



AcTbeta <- data.frame(cTbetaoutput[[1]]$X10,cTbetaoutput[[2]]$X10,cTbetaoutput[[3]]$X10,cTbetaoutput[[4]]$X10,cTbetaoutput[[5]]$X10,cTbetaoutput[[6]]$X10,cTbetaoutput[[7]]$X10,cTbetaoutput[[8]]$X10,cTbetaoutput[[9]]$X10,cTbetaoutput[[10]]$X10,cTbetaoutput[[11]]$X10,cTbetaoutput[[12]]$X10,cTbetaoutput[[13]]$X10,cTbetaoutput[[14]]$X10,cTbetaoutput[[15]]$X10,cTbetaoutput[[16]]$X10,cTbetaoutput[[17]]$X10,cTbetaoutput[[18]]$X10,cTbetaoutput[[19]]$X10,cTbetaoutput[[20]]$X10,cTbetaoutput[[21]]$X10,cTbetaoutput[[22]]$X10,cTbetaoutput[[23]]$X10,cTbetaoutput[[24]]$X10,cTbetaoutput[[25]]$X10,cTbetaoutput[[26]]$X10,cTbetaoutput[[27]]$X10,cTbetaoutput[[28]]$X10,cTbetaoutput[[29]]$X10,cTbetaoutput[[30]]$X10,cTbetaoutput[[31]]$X10,cTbetaoutput[[32]]$X10,cTbetaoutput[[33]]$X10,cTbetaoutput[[34]]$X10,cTbetaoutput[[35]]$X10,cTbetaoutput[[36]]$X10,cTbetaoutput[[37]]$X10,cTbetaoutput[[38]]$X10,cTbetaoutput[[39]]$X10,cTbetaoutput[[40]]$X10,cTbetaoutput[[41]]$X10,cTbetaoutput[[42]]$X10,cTbetaoutput[[43]]$X10,cTbetaoutput[[44]]$X10,cTbetaoutput[[45]]$X10,cTbetaoutput[[46]]$X10,cTbetaoutput[[47]]$X10,cTbetaoutput[[48]]$X10,cTbetaoutput[[49]]$X10,cTbetaoutput[[50]]$X10,cTbetaoutput[[51]]$X10,cTbetaoutput[[52]]$X10,cTbetaoutput[[53]]$X10,cTbetaoutput[[54]]$X10,cTbetaoutput[[55]]$X10,cTbetaoutput[[56]]$X10,cTbetaoutput[[57]]$X10,cTbetaoutput[[58]]$X10,cTbetaoutput[[59]]$X10,cTbetaoutput[[60]]$X10,cTbetaoutput[[61]]$X10,cTbetaoutput[[62]]$X10,cTbetaoutput[[63]]$X10,cTbetaoutput[[64]]$X10,cTbetaoutput[[65]]$X10,cTbetaoutput[[66]]$X10,cTbetaoutput[[67]]$X10,cTbetaoutput[[68]]$X10,cTbetaoutput[[69]]$X10,cTbetaoutput[[70]]$X10,cTbetaoutput[[71]]$X10,cTbetaoutput[[72]]$X10,cTbetaoutput[[73]]$X10,cTbetaoutput[[74]]$X10,cTbetaoutput[[75]]$X10,cTbetaoutput[[76]]$X10,cTbetaoutput[[77]]$X10,cTbetaoutput[[78]]$X10,cTbetaoutput[[79]]$X10,cTbetaoutput[[80]]$X10,cTbetaoutput[[81]]$X10,cTbetaoutput[[82]]$X10,cTbetaoutput[[83]]$X10,cTbetaoutput[[84]]$X10,cTbetaoutput[[85]]$X10,cTbetaoutput[[86]]$X10,cTbetaoutput[[87]]$X10,cTbetaoutput[[88]]$X10,cTbetaoutput[[89]]$X10,cTbetaoutput[[90]]$X10,cTbetaoutput[[91]]$X10,cTbetaoutput[[92]]$X10,cTbetaoutput[[93]]$X10,cTbetaoutput[[94]]$X10,cTbetaoutput[[95]]$X10,cTbetaoutput[[96]]$X10,cTbetaoutput[[97]]$X10,cTbetaoutput[[98]]$X10,cTbetaoutput[[99]]$X10,cTbetaoutput[[100]]$X10,cTbetaoutput[[101]]$X10,cTbetaoutput[[102]]$X10,cTbetaoutput[[103]]$X10,cTbetaoutput[[104]]$X10,cTbetaoutput[[105]]$X10,cTbetaoutput[[106]]$X10,cTbetaoutput[[107]]$X10,cTbetaoutput[[108]]$X10,cTbetaoutput[[109]]$X10,cTbetaoutput[[110]]$X10,cTbetaoutput[[111]]$X10,cTbetaoutput[[112]]$X10,cTbetaoutput[[113]]$X10,cTbetaoutput[[114]]$X10,cTbetaoutput[[115]]$X10,cTbetaoutput[[116]]$X10,cTbetaoutput[[117]]$X10,cTbetaoutput[[118]]$X10,cTbetaoutput[[119]]$X10,cTbetaoutput[[120]]$X10,cTbetaoutput[[121]]$X10)

AcTbeta <- melt(AcTbeta, variable.name = "population", value.name = "X1")

AcTbeta$Y1 <- c(cTbetaoutput[[1]]$Y10,cTbetaoutput[[2]]$Y10,cTbetaoutput[[3]]$Y10,cTbetaoutput[[4]]$Y10,cTbetaoutput[[5]]$Y10,cTbetaoutput[[6]]$Y10,cTbetaoutput[[7]]$Y10,cTbetaoutput[[8]]$Y10,cTbetaoutput[[9]]$Y10,cTbetaoutput[[10]]$Y10,cTbetaoutput[[11]]$Y10,cTbetaoutput[[12]]$Y10,cTbetaoutput[[13]]$Y10,cTbetaoutput[[14]]$Y10,cTbetaoutput[[15]]$Y10,cTbetaoutput[[16]]$Y10,cTbetaoutput[[17]]$Y10,cTbetaoutput[[18]]$Y10,cTbetaoutput[[19]]$Y10,cTbetaoutput[[20]]$Y10,cTbetaoutput[[21]]$Y10,cTbetaoutput[[22]]$Y10,cTbetaoutput[[23]]$Y10,cTbetaoutput[[24]]$Y10,cTbetaoutput[[25]]$Y10,cTbetaoutput[[26]]$Y10,cTbetaoutput[[27]]$Y10,cTbetaoutput[[28]]$Y10,cTbetaoutput[[29]]$Y10,cTbetaoutput[[30]]$Y10,cTbetaoutput[[31]]$Y10,cTbetaoutput[[32]]$Y10,cTbetaoutput[[33]]$Y10,cTbetaoutput[[34]]$Y10,cTbetaoutput[[35]]$Y10,cTbetaoutput[[36]]$Y10,cTbetaoutput[[37]]$Y10,cTbetaoutput[[38]]$Y10,cTbetaoutput[[39]]$Y10,cTbetaoutput[[40]]$Y10,cTbetaoutput[[41]]$Y10,cTbetaoutput[[42]]$Y10,cTbetaoutput[[43]]$Y10,cTbetaoutput[[44]]$Y10,cTbetaoutput[[45]]$Y10,cTbetaoutput[[46]]$Y10,cTbetaoutput[[47]]$Y10,cTbetaoutput[[48]]$Y10,cTbetaoutput[[49]]$Y10,cTbetaoutput[[50]]$Y10,cTbetaoutput[[51]]$Y10,cTbetaoutput[[52]]$Y10,cTbetaoutput[[53]]$Y10,cTbetaoutput[[54]]$Y10,cTbetaoutput[[55]]$Y10,cTbetaoutput[[56]]$Y10,cTbetaoutput[[57]]$Y10,cTbetaoutput[[58]]$Y10,cTbetaoutput[[59]]$Y10,cTbetaoutput[[60]]$Y10,cTbetaoutput[[61]]$Y10,cTbetaoutput[[62]]$Y10,cTbetaoutput[[63]]$Y10,cTbetaoutput[[64]]$Y10,cTbetaoutput[[65]]$Y10,cTbetaoutput[[66]]$Y10,cTbetaoutput[[67]]$Y10,cTbetaoutput[[68]]$Y10,cTbetaoutput[[69]]$Y10,cTbetaoutput[[70]]$Y10,cTbetaoutput[[71]]$Y10,cTbetaoutput[[72]]$Y10,cTbetaoutput[[73]]$Y10,cTbetaoutput[[74]]$Y10,cTbetaoutput[[75]]$Y10,cTbetaoutput[[76]]$Y10,cTbetaoutput[[77]]$Y10,cTbetaoutput[[78]]$Y10,cTbetaoutput[[79]]$Y10,cTbetaoutput[[80]]$Y10,cTbetaoutput[[81]]$Y10,cTbetaoutput[[82]]$Y10,cTbetaoutput[[83]]$Y10,cTbetaoutput[[84]]$Y10,cTbetaoutput[[85]]$Y10,cTbetaoutput[[86]]$Y10,cTbetaoutput[[87]]$Y10,cTbetaoutput[[88]]$Y10,cTbetaoutput[[89]]$Y10,cTbetaoutput[[90]]$Y10,cTbetaoutput[[91]]$Y10,cTbetaoutput[[92]]$Y10,cTbetaoutput[[93]]$Y10,cTbetaoutput[[94]]$Y10,cTbetaoutput[[95]]$Y10,cTbetaoutput[[96]]$Y10,cTbetaoutput[[97]]$Y10,cTbetaoutput[[98]]$Y10,cTbetaoutput[[99]]$Y10,cTbetaoutput[[100]]$Y10,cTbetaoutput[[101]]$Y10,cTbetaoutput[[102]]$Y10,cTbetaoutput[[103]]$Y10,cTbetaoutput[[104]]$Y10,cTbetaoutput[[105]]$Y10,cTbetaoutput[[106]]$Y10,cTbetaoutput[[107]]$Y10,cTbetaoutput[[108]]$Y10,cTbetaoutput[[109]]$Y10,cTbetaoutput[[110]]$Y10,cTbetaoutput[[111]]$Y10,cTbetaoutput[[112]]$Y10,cTbetaoutput[[113]]$Y10,cTbetaoutput[[114]]$Y10,cTbetaoutput[[115]]$Y10,cTbetaoutput[[116]]$Y10,cTbetaoutput[[117]]$Y10,cTbetaoutput[[118]]$Y10,cTbetaoutput[[119]]$Y10,cTbetaoutput[[120]]$Y10,cTbetaoutput[[121]]$Y10)

AcTbeta$X2 <- c(cTbetaoutput[[1]]$X20,cTbetaoutput[[2]]$X20,cTbetaoutput[[3]]$X20,cTbetaoutput[[4]]$X20,cTbetaoutput[[5]]$X20,cTbetaoutput[[6]]$X20,cTbetaoutput[[7]]$X20,cTbetaoutput[[8]]$X20,cTbetaoutput[[9]]$X20,cTbetaoutput[[10]]$X20,cTbetaoutput[[11]]$X20,cTbetaoutput[[12]]$X20,cTbetaoutput[[13]]$X20,cTbetaoutput[[14]]$X20,cTbetaoutput[[15]]$X20,cTbetaoutput[[16]]$X20,cTbetaoutput[[17]]$X20,cTbetaoutput[[18]]$X20,cTbetaoutput[[19]]$X20,cTbetaoutput[[20]]$X20,cTbetaoutput[[21]]$X20,cTbetaoutput[[22]]$X20,cTbetaoutput[[23]]$X20,cTbetaoutput[[24]]$X20,cTbetaoutput[[25]]$X20,cTbetaoutput[[26]]$X20,cTbetaoutput[[27]]$X20,cTbetaoutput[[28]]$X20,cTbetaoutput[[29]]$X20,cTbetaoutput[[30]]$X20,cTbetaoutput[[31]]$X20,cTbetaoutput[[32]]$X20,cTbetaoutput[[33]]$X20,cTbetaoutput[[34]]$X20,cTbetaoutput[[35]]$X20,cTbetaoutput[[36]]$X20,cTbetaoutput[[37]]$X20,cTbetaoutput[[38]]$X20,cTbetaoutput[[39]]$X20,cTbetaoutput[[40]]$X20,cTbetaoutput[[41]]$X20,cTbetaoutput[[42]]$X20,cTbetaoutput[[43]]$X20,cTbetaoutput[[44]]$X20,cTbetaoutput[[45]]$X20,cTbetaoutput[[46]]$X20,cTbetaoutput[[47]]$X20,cTbetaoutput[[48]]$X20,cTbetaoutput[[49]]$X20,cTbetaoutput[[50]]$X20,cTbetaoutput[[51]]$X20,cTbetaoutput[[52]]$X20,cTbetaoutput[[53]]$X20,cTbetaoutput[[54]]$X20,cTbetaoutput[[55]]$X20,cTbetaoutput[[56]]$X20,cTbetaoutput[[57]]$X20,cTbetaoutput[[58]]$X20,cTbetaoutput[[59]]$X20,cTbetaoutput[[60]]$X20,cTbetaoutput[[61]]$X20,cTbetaoutput[[62]]$X20,cTbetaoutput[[63]]$X20,cTbetaoutput[[64]]$X20,cTbetaoutput[[65]]$X20,cTbetaoutput[[66]]$X20,cTbetaoutput[[67]]$X20,cTbetaoutput[[68]]$X20,cTbetaoutput[[69]]$X20,cTbetaoutput[[70]]$X20,cTbetaoutput[[71]]$X20,cTbetaoutput[[72]]$X20,cTbetaoutput[[73]]$X20,cTbetaoutput[[74]]$X20,cTbetaoutput[[75]]$X20,cTbetaoutput[[76]]$X20,cTbetaoutput[[77]]$X20,cTbetaoutput[[78]]$X20,cTbetaoutput[[79]]$X20,cTbetaoutput[[80]]$X20,cTbetaoutput[[81]]$X20,cTbetaoutput[[82]]$X20,cTbetaoutput[[83]]$X20,cTbetaoutput[[84]]$X20,cTbetaoutput[[85]]$X20,cTbetaoutput[[86]]$X20,cTbetaoutput[[87]]$X20,cTbetaoutput[[88]]$X20,cTbetaoutput[[89]]$X20,cTbetaoutput[[90]]$X20,cTbetaoutput[[91]]$X20,cTbetaoutput[[92]]$X20,cTbetaoutput[[93]]$X20,cTbetaoutput[[94]]$X20,cTbetaoutput[[95]]$X20,cTbetaoutput[[96]]$X20,cTbetaoutput[[97]]$X20,cTbetaoutput[[98]]$X20,cTbetaoutput[[99]]$X20,cTbetaoutput[[100]]$X20,cTbetaoutput[[101]]$X20,cTbetaoutput[[102]]$X20,cTbetaoutput[[103]]$X20,cTbetaoutput[[104]]$X20,cTbetaoutput[[105]]$X20,cTbetaoutput[[106]]$X20,cTbetaoutput[[107]]$X20,cTbetaoutput[[108]]$X20,cTbetaoutput[[109]]$X20,cTbetaoutput[[110]]$X20,cTbetaoutput[[111]]$X20,cTbetaoutput[[112]]$X20,cTbetaoutput[[113]]$X20,cTbetaoutput[[114]]$X20,cTbetaoutput[[115]]$X20,cTbetaoutput[[116]]$X20,cTbetaoutput[[117]]$X20,cTbetaoutput[[118]]$X20,cTbetaoutput[[119]]$X20,cTbetaoutput[[120]]$X20,cTbetaoutput[[121]]$X20)

AcTbeta$Y2 <- c(cTbetaoutput[[1]]$Y20,cTbetaoutput[[2]]$Y20,cTbetaoutput[[3]]$Y20,cTbetaoutput[[4]]$Y20,cTbetaoutput[[5]]$Y20,cTbetaoutput[[6]]$Y20,cTbetaoutput[[7]]$Y20,cTbetaoutput[[8]]$Y20,cTbetaoutput[[9]]$Y20,cTbetaoutput[[10]]$Y20,cTbetaoutput[[11]]$Y20,cTbetaoutput[[12]]$Y20,cTbetaoutput[[13]]$Y20,cTbetaoutput[[14]]$Y20,cTbetaoutput[[15]]$Y20,cTbetaoutput[[16]]$Y20,cTbetaoutput[[17]]$Y20,cTbetaoutput[[18]]$Y20,cTbetaoutput[[19]]$Y20,cTbetaoutput[[20]]$Y20,cTbetaoutput[[21]]$Y20,cTbetaoutput[[22]]$Y20,cTbetaoutput[[23]]$Y20,cTbetaoutput[[24]]$Y20,cTbetaoutput[[25]]$Y20,cTbetaoutput[[26]]$Y20,cTbetaoutput[[27]]$Y20,cTbetaoutput[[28]]$Y20,cTbetaoutput[[29]]$Y20,cTbetaoutput[[30]]$Y20,cTbetaoutput[[31]]$Y20,cTbetaoutput[[32]]$Y20,cTbetaoutput[[33]]$Y20,cTbetaoutput[[34]]$Y20,cTbetaoutput[[35]]$Y20,cTbetaoutput[[36]]$Y20,cTbetaoutput[[37]]$Y20,cTbetaoutput[[38]]$Y20,cTbetaoutput[[39]]$Y20,cTbetaoutput[[40]]$Y20,cTbetaoutput[[41]]$Y20,cTbetaoutput[[42]]$Y20,cTbetaoutput[[43]]$Y20,cTbetaoutput[[44]]$Y20,cTbetaoutput[[45]]$Y20,cTbetaoutput[[46]]$Y20,cTbetaoutput[[47]]$Y20,cTbetaoutput[[48]]$Y20,cTbetaoutput[[49]]$Y20,cTbetaoutput[[50]]$Y20,cTbetaoutput[[51]]$Y20,cTbetaoutput[[52]]$Y20,cTbetaoutput[[53]]$Y20,cTbetaoutput[[54]]$Y20,cTbetaoutput[[55]]$Y20,cTbetaoutput[[56]]$Y20,cTbetaoutput[[57]]$Y20,cTbetaoutput[[58]]$Y20,cTbetaoutput[[59]]$Y20,cTbetaoutput[[60]]$Y20,cTbetaoutput[[61]]$Y20,cTbetaoutput[[62]]$Y20,cTbetaoutput[[63]]$Y20,cTbetaoutput[[64]]$Y20,cTbetaoutput[[65]]$Y20,cTbetaoutput[[66]]$Y20,cTbetaoutput[[67]]$Y20,cTbetaoutput[[68]]$Y20,cTbetaoutput[[69]]$Y20,cTbetaoutput[[70]]$Y20,cTbetaoutput[[71]]$Y20,cTbetaoutput[[72]]$Y20,cTbetaoutput[[73]]$Y20,cTbetaoutput[[74]]$Y20,cTbetaoutput[[75]]$Y20,cTbetaoutput[[76]]$Y20,cTbetaoutput[[77]]$Y20,cTbetaoutput[[78]]$Y20,cTbetaoutput[[79]]$Y20,cTbetaoutput[[80]]$Y20,cTbetaoutput[[81]]$Y20,cTbetaoutput[[82]]$Y20,cTbetaoutput[[83]]$Y20,cTbetaoutput[[84]]$Y20,cTbetaoutput[[85]]$Y20,cTbetaoutput[[86]]$Y20,cTbetaoutput[[87]]$Y20,cTbetaoutput[[88]]$Y20,cTbetaoutput[[89]]$Y20,cTbetaoutput[[90]]$Y20,cTbetaoutput[[91]]$Y20,cTbetaoutput[[92]]$Y20,cTbetaoutput[[93]]$Y20,cTbetaoutput[[94]]$Y20,cTbetaoutput[[95]]$Y20,cTbetaoutput[[96]]$Y20,cTbetaoutput[[97]]$Y20,cTbetaoutput[[98]]$Y20,cTbetaoutput[[99]]$Y20,cTbetaoutput[[100]]$Y20,cTbetaoutput[[101]]$Y20,cTbetaoutput[[102]]$Y20,cTbetaoutput[[103]]$Y20,cTbetaoutput[[104]]$Y20,cTbetaoutput[[105]]$Y20,cTbetaoutput[[106]]$Y20,cTbetaoutput[[107]]$Y20,cTbetaoutput[[108]]$Y20,cTbetaoutput[[109]]$Y20,cTbetaoutput[[110]]$Y20,cTbetaoutput[[111]]$Y20,cTbetaoutput[[112]]$Y20,cTbetaoutput[[113]]$Y20,cTbetaoutput[[114]]$Y20,cTbetaoutput[[115]]$Y20,cTbetaoutput[[116]]$Y20,cTbetaoutput[[117]]$Y20,cTbetaoutput[[118]]$Y20,cTbetaoutput[[119]]$Y20,cTbetaoutput[[120]]$Y20,cTbetaoutput[[121]]$Y20)

colnames(AcTbeta)[colnames(AcTbeta)=="value"] <- "X1"
colnames(AcTbeta)[colnames(AcTbeta)=="variable"] <- "gammaA"
AcTbeta$alpha <-alphaoutput
AcTbeta$beta <-betaoutput
AcTbeta$cT <- c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0)
AcTbeta



AcTbeta2 <- as.data.frame(lapply(AcTbeta, unlist))
AcTbeta2 <- subset(AcTbeta2, select = -c(population))
AcTbeta3 <- subset(AcTbeta2, select = -c(alpha))
AcTbeta3



AcTbeta3[which(AcTbeta3$X1<1),1]<-0
AcTbeta3[which(AcTbeta3$Y1<1),2]<-0
AcTbeta3[which(AcTbeta3$X2<1),3]<-0
AcTbeta3[which(AcTbeta3$Y2<1),4]<-0
AcTbeta3



library(dplyr)
df_AcTbeta<- mutate(AcTbeta3, M = (X2+Y2)/(X1+Y1+X2+Y2),R = (X1+Y1)/(X1+Y1+X2+Y2)) # create new columns with Migrant and Resident ratios
df_AcTbeta <- df_AcTbeta %>% select(beta,cT,M,R)
df_AcTbeta



df2_AcTbeta <- df_AcTbeta %>% filter(cT == 0.9) %>% filter(beta == c(0,0.1,0.2,0.3)) %>% mutate(mig = "M", red= "R") %>% pivot_longer(M:R) %>% select(beta,cT, name, value)
df2_AcTbeta






