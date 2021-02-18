# ODE cycles - equations 
# load libraries
#install.packages(pacman)
pacman::p_load(deSolve,reshape2,akima,lattice,dplyr)

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
beta <- 0.01 # tramsmission rate
gammaA <- 0.3 # recovery rate in environment A
gammaB <- 0.5 # recovery rate in environment B 
mu <- 0.01 # mortality rate
#alpha <- 0.2 # parasite induced death rate
phi=3 # maximum per capita fecundity rate for each individual
z=0.0007 # density dependent fecundity coefficient 
#cm=0.3 # costs of migration
ct=0 # costs of tolerance (let ct == 0 in this case to see the infection factor drives partial migration)

# 100 years loop 
alphacmoutput <- vector(mode="list", length=121) # create 20 list dummy vectors for both migrants and residents
alphaoutput <- vector(mode="list", length=121) # create 20 list dummy vectors for both migrants and residents
cMoutput <- vector(mode="list", length=121) # create 20 list dummy vectors for both migrants and residents
alphavals <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) # recovery rate in environment A
cMvals <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) # recovery rate in environment B
for (i in 1:length(alphavals)){
  print(alphavals[i])
  for (g in 1:length(cMvals)){
    print(cMvals[g])
    X10 <- 150 # initial susceptible population (Resident) 
    X20 <- 150 # initial susceptible population (Migrant)
    Y10 <- 100 # initial infected population (Resident)
    Y20 <- 100 # initial infected population (Migrant)
    for (t in 1:100){
      # Environment A (M+R) together - first 6 months
      params1 <- c(beta=beta,gamma=gammaA,mu=mu,alpha=alphavals[i]) # parameters
      xstart1 <- c(X1=X10,X2=X20,Y1=Y10,Y2=Y20) # initial conditions
      alphacM1 <- as.data.frame(ode(func = open.sir.model1,y=xstart1,times=times1,parms=params1,method="ode23"))
      
      # Environment A (only R) - Second 6 months
      params2 <- c(beta=beta,gamma=gammaA,mu=mu,alpha=alphavals[i]) # parameters
      xstart2 <- c(X1=alphacM1[7,2],Y1=alphacM1[7,4]) # initial conditions
      alphacM2 <- as.data.frame(ode(func = open.sir.model2,y=xstart2,times=times2,parms=params2,method="ode23"))
      
      # Environment B for migrants - Second 6 months 
      params3 <- c(beta=beta,gamma=gammaB,mu=mu,alpha=alphavals[i]) 
      xstart3 <- c(X2=alphacM1[7,3],Y2=alphacM1[7,5]) # Initial survived population from env A
      alphacM3 <- as.data.frame(ode(func = open.sir.model3,y=xstart3,times=times3,parms=params3,method="ode23"))
      
      # Calculation of fecundity rates in 2nd year 
      
      # the annual survivors 
      
      X1hat <- alphacM2[7,2] # resdient susceptible
      X2hat <- alphacM3[7,2] # migrant susceptible
      Y1hat <- alphacM2[7,3] # resident infected
      Y2hat <- alphacM3[7,3] # migrant susceptible 
      
      bRmax <- phi*(X1hat+Y1hat)*(1-ct) # maximum fecundity rate of a resident individual 
      bMmax <- phi*(X2hat+Y2hat)*(1-cMvals[g])*(1-ct) # maximum fecundity rate of a migrant individual 
      
      b1 <- bRmax*(exp(-z*(bRmax+bMmax))) # actual fecundity rate of a resident individual
      b2 <- bMmax*(exp(-z*(bRmax+bMmax))) # actual fecundity rate of a migrant individual 
      
      # after fecundity - 2nd yr initial populations 
      X10 <- alphacM2[7,2] + b1
      X20 <- alphacM3[7,2] + b2
      Y10 <- alphacM2[7,3]
      Y20 <- alphacM3[7,3]
    }
    alphacmresult <- data.frame(X10,X20,Y10,Y20)
    print(alphacmresult)
    
    alphacmoutput[[(g-1)*11+i]]<-as.data.frame(alphacmresult)
    cMoutput[[(g-1)*11+i]]<-as.data.frame(cMvals[g])
    alphaoutput[[(g-1)*11+i]]<-as.data.frame(alphavals[i])}} # CT loop

cmalpha <- data.frame(alphacmoutput[[1]]$X10,alphacmoutput[[2]]$X10,alphacmoutput[[3]]$X10,alphacmoutput[[4]]$X10,alphacmoutput[[5]]$X10,alphacmoutput[[6]]$X10,alphacmoutput[[7]]$X10,alphacmoutput[[8]]$X10,alphacmoutput[[9]]$X10,alphacmoutput[[10]]$X10,alphacmoutput[[11]]$X10,alphacmoutput[[12]]$X10,alphacmoutput[[13]]$X10,alphacmoutput[[14]]$X10,alphacmoutput[[15]]$X10,alphacmoutput[[16]]$X10,alphacmoutput[[17]]$X10,alphacmoutput[[18]]$X10,alphacmoutput[[19]]$X10,alphacmoutput[[20]]$X10,alphacmoutput[[21]]$X10,alphacmoutput[[22]]$X10,alphacmoutput[[23]]$X10,alphacmoutput[[24]]$X10,alphacmoutput[[25]]$X10,alphacmoutput[[26]]$X10,alphacmoutput[[27]]$X10,alphacmoutput[[28]]$X10,alphacmoutput[[29]]$X10,alphacmoutput[[30]]$X10,alphacmoutput[[31]]$X10,alphacmoutput[[32]]$X10,alphacmoutput[[33]]$X10,alphacmoutput[[34]]$X10,alphacmoutput[[35]]$X10,alphacmoutput[[36]]$X10,alphacmoutput[[37]]$X10,alphacmoutput[[38]]$X10,alphacmoutput[[39]]$X10,alphacmoutput[[40]]$X10,alphacmoutput[[41]]$X10,alphacmoutput[[42]]$X10,alphacmoutput[[43]]$X10,alphacmoutput[[44]]$X10,alphacmoutput[[45]]$X10,alphacmoutput[[46]]$X10,alphacmoutput[[47]]$X10,alphacmoutput[[48]]$X10,alphacmoutput[[49]]$X10,alphacmoutput[[50]]$X10,alphacmoutput[[51]]$X10,alphacmoutput[[52]]$X10,alphacmoutput[[53]]$X10,alphacmoutput[[54]]$X10,alphacmoutput[[55]]$X10,alphacmoutput[[56]]$X10,alphacmoutput[[57]]$X10,alphacmoutput[[58]]$X10,alphacmoutput[[59]]$X10,alphacmoutput[[60]]$X10,alphacmoutput[[61]]$X10,alphacmoutput[[62]]$X10,alphacmoutput[[63]]$X10,alphacmoutput[[64]]$X10,alphacmoutput[[65]]$X10,alphacmoutput[[66]]$X10,alphacmoutput[[67]]$X10,alphacmoutput[[68]]$X10,alphacmoutput[[69]]$X10,alphacmoutput[[70]]$X10,alphacmoutput[[71]]$X10,alphacmoutput[[72]]$X10,alphacmoutput[[73]]$X10,alphacmoutput[[74]]$X10,alphacmoutput[[75]]$X10,alphacmoutput[[76]]$X10,alphacmoutput[[77]]$X10,alphacmoutput[[78]]$X10,alphacmoutput[[79]]$X10,alphacmoutput[[80]]$X10,alphacmoutput[[81]]$X10,alphacmoutput[[82]]$X10,alphacmoutput[[83]]$X10,alphacmoutput[[84]]$X10,alphacmoutput[[85]]$X10,alphacmoutput[[86]]$X10,alphacmoutput[[87]]$X10,alphacmoutput[[88]]$X10,alphacmoutput[[89]]$X10,alphacmoutput[[90]]$X10,alphacmoutput[[91]]$X10,alphacmoutput[[92]]$X10,alphacmoutput[[93]]$X10,alphacmoutput[[94]]$X10,alphacmoutput[[95]]$X10,alphacmoutput[[96]]$X10,alphacmoutput[[97]]$X10,alphacmoutput[[98]]$X10,alphacmoutput[[99]]$X10,alphacmoutput[[100]]$X10,alphacmoutput[[101]]$X10,alphacmoutput[[102]]$X10,alphacmoutput[[103]]$X10,alphacmoutput[[104]]$X10,alphacmoutput[[105]]$X10,alphacmoutput[[106]]$X10,alphacmoutput[[107]]$X10,alphacmoutput[[108]]$X10,alphacmoutput[[109]]$X10,alphacmoutput[[110]]$X10,alphacmoutput[[111]]$X10,alphacmoutput[[112]]$X10,alphacmoutput[[113]]$X10,alphacmoutput[[114]]$X10,alphacmoutput[[115]]$X10,alphacmoutput[[116]]$X10,alphacmoutput[[117]]$X10,alphacmoutput[[118]]$X10,alphacmoutput[[119]]$X10,alphacmoutput[[120]]$X10,alphacmoutput[[121]]$X10)
cmalpha <- melt(cmalpha, variable.name = "population", value.name = "X1")
cmalpha$Y1 <- c(alphacmoutput[[1]]$Y10,alphacmoutput[[2]]$Y10,alphacmoutput[[3]]$Y10,alphacmoutput[[4]]$Y10,alphacmoutput[[5]]$Y10,alphacmoutput[[6]]$Y10,alphacmoutput[[7]]$Y10,alphacmoutput[[8]]$Y10,alphacmoutput[[9]]$Y10,alphacmoutput[[10]]$Y10,alphacmoutput[[11]]$Y10,alphacmoutput[[12]]$Y10,alphacmoutput[[13]]$Y10,alphacmoutput[[14]]$Y10,alphacmoutput[[15]]$Y10,alphacmoutput[[16]]$Y10,alphacmoutput[[17]]$Y10,alphacmoutput[[18]]$Y10,alphacmoutput[[19]]$Y10,alphacmoutput[[20]]$Y10,alphacmoutput[[21]]$Y10,alphacmoutput[[22]]$Y10,alphacmoutput[[23]]$Y10,alphacmoutput[[24]]$Y10,alphacmoutput[[25]]$Y10,alphacmoutput[[26]]$Y10,alphacmoutput[[27]]$Y10,alphacmoutput[[28]]$Y10,alphacmoutput[[29]]$Y10,alphacmoutput[[30]]$Y10,alphacmoutput[[31]]$Y10,alphacmoutput[[32]]$Y10,alphacmoutput[[33]]$Y10,alphacmoutput[[34]]$Y10,alphacmoutput[[35]]$Y10,alphacmoutput[[36]]$Y10,alphacmoutput[[37]]$Y10,alphacmoutput[[38]]$Y10,alphacmoutput[[39]]$Y10,alphacmoutput[[40]]$Y10,alphacmoutput[[41]]$Y10,alphacmoutput[[42]]$Y10,alphacmoutput[[43]]$Y10,alphacmoutput[[44]]$Y10,alphacmoutput[[45]]$Y10,alphacmoutput[[46]]$Y10,alphacmoutput[[47]]$Y10,alphacmoutput[[48]]$Y10,alphacmoutput[[49]]$Y10,alphacmoutput[[50]]$Y10,alphacmoutput[[51]]$Y10,alphacmoutput[[52]]$Y10,alphacmoutput[[53]]$Y10,alphacmoutput[[54]]$Y10,alphacmoutput[[55]]$Y10,alphacmoutput[[56]]$Y10,alphacmoutput[[57]]$Y10,alphacmoutput[[58]]$Y10,alphacmoutput[[59]]$Y10,alphacmoutput[[60]]$Y10,alphacmoutput[[61]]$Y10,alphacmoutput[[62]]$Y10,alphacmoutput[[63]]$Y10,alphacmoutput[[64]]$Y10,alphacmoutput[[65]]$Y10,alphacmoutput[[66]]$Y10,alphacmoutput[[67]]$Y10,alphacmoutput[[68]]$Y10,alphacmoutput[[69]]$Y10,alphacmoutput[[70]]$Y10,alphacmoutput[[71]]$Y10,alphacmoutput[[72]]$Y10,alphacmoutput[[73]]$Y10,alphacmoutput[[74]]$Y10,alphacmoutput[[75]]$Y10,alphacmoutput[[76]]$Y10,alphacmoutput[[77]]$Y10,alphacmoutput[[78]]$Y10,alphacmoutput[[79]]$Y10,alphacmoutput[[80]]$Y10,alphacmoutput[[81]]$Y10,alphacmoutput[[82]]$Y10,alphacmoutput[[83]]$Y10,alphacmoutput[[84]]$Y10,alphacmoutput[[85]]$Y10,alphacmoutput[[86]]$Y10,alphacmoutput[[87]]$Y10,alphacmoutput[[88]]$Y10,alphacmoutput[[89]]$Y10,alphacmoutput[[90]]$Y10,alphacmoutput[[91]]$Y10,alphacmoutput[[92]]$Y10,alphacmoutput[[93]]$Y10,alphacmoutput[[94]]$Y10,alphacmoutput[[95]]$Y10,alphacmoutput[[96]]$Y10,alphacmoutput[[97]]$Y10,alphacmoutput[[98]]$Y10,alphacmoutput[[99]]$Y10,alphacmoutput[[100]]$Y10,alphacmoutput[[101]]$Y10,alphacmoutput[[102]]$Y10,alphacmoutput[[103]]$Y10,alphacmoutput[[104]]$Y10,alphacmoutput[[105]]$Y10,alphacmoutput[[106]]$Y10,alphacmoutput[[107]]$Y10,alphacmoutput[[108]]$Y10,alphacmoutput[[109]]$Y10,alphacmoutput[[110]]$Y10,alphacmoutput[[111]]$Y10,alphacmoutput[[112]]$Y10,alphacmoutput[[113]]$Y10,alphacmoutput[[114]]$Y10,alphacmoutput[[115]]$Y10,alphacmoutput[[116]]$Y10,alphacmoutput[[117]]$Y10,alphacmoutput[[118]]$Y10,alphacmoutput[[119]]$Y10,alphacmoutput[[120]]$Y10,alphacmoutput[[121]]$Y10)
cmalpha$X2 <- c(alphacmoutput[[1]]$X20,alphacmoutput[[2]]$X20,alphacmoutput[[3]]$X20,alphacmoutput[[4]]$X20,alphacmoutput[[5]]$X20,alphacmoutput[[6]]$X20,alphacmoutput[[7]]$X20,alphacmoutput[[8]]$X20,alphacmoutput[[9]]$X20,alphacmoutput[[10]]$X20,alphacmoutput[[11]]$X20,alphacmoutput[[12]]$X20,alphacmoutput[[13]]$X20,alphacmoutput[[14]]$X20,alphacmoutput[[15]]$X20,alphacmoutput[[16]]$X20,alphacmoutput[[17]]$X20,alphacmoutput[[18]]$X20,alphacmoutput[[19]]$X20,alphacmoutput[[20]]$X20,alphacmoutput[[21]]$X20,alphacmoutput[[22]]$X20,alphacmoutput[[23]]$X20,alphacmoutput[[24]]$X20,alphacmoutput[[25]]$X20,alphacmoutput[[26]]$X20,alphacmoutput[[27]]$X20,alphacmoutput[[28]]$X20,alphacmoutput[[29]]$X20,alphacmoutput[[30]]$X20,alphacmoutput[[31]]$X20,alphacmoutput[[32]]$X20,alphacmoutput[[33]]$X20,alphacmoutput[[34]]$X20,alphacmoutput[[35]]$X20,alphacmoutput[[36]]$X20,alphacmoutput[[37]]$X20,alphacmoutput[[38]]$X20,alphacmoutput[[39]]$X20,alphacmoutput[[40]]$X20,alphacmoutput[[41]]$X20,alphacmoutput[[42]]$X20,alphacmoutput[[43]]$X20,alphacmoutput[[44]]$X20,alphacmoutput[[45]]$X20,alphacmoutput[[46]]$X20,alphacmoutput[[47]]$X20,alphacmoutput[[48]]$X20,alphacmoutput[[49]]$X20,alphacmoutput[[50]]$X20,alphacmoutput[[51]]$X20,alphacmoutput[[52]]$X20,alphacmoutput[[53]]$X20,alphacmoutput[[54]]$X20,alphacmoutput[[55]]$X20,alphacmoutput[[56]]$X20,alphacmoutput[[57]]$X20,alphacmoutput[[58]]$X20,alphacmoutput[[59]]$X20,alphacmoutput[[60]]$X20,alphacmoutput[[61]]$X20,alphacmoutput[[62]]$X20,alphacmoutput[[63]]$X20,alphacmoutput[[64]]$X20,alphacmoutput[[65]]$X20,alphacmoutput[[66]]$X20,alphacmoutput[[67]]$X20,alphacmoutput[[68]]$X20,alphacmoutput[[69]]$X20,alphacmoutput[[70]]$X20,alphacmoutput[[71]]$X20,alphacmoutput[[72]]$X20,alphacmoutput[[73]]$X20,alphacmoutput[[74]]$X20,alphacmoutput[[75]]$X20,alphacmoutput[[76]]$X20,alphacmoutput[[77]]$X20,alphacmoutput[[78]]$X20,alphacmoutput[[79]]$X20,alphacmoutput[[80]]$X20,alphacmoutput[[81]]$X20,alphacmoutput[[82]]$X20,alphacmoutput[[83]]$X20,alphacmoutput[[84]]$X20,alphacmoutput[[85]]$X20,alphacmoutput[[86]]$X20,alphacmoutput[[87]]$X20,alphacmoutput[[88]]$X20,alphacmoutput[[89]]$X20,alphacmoutput[[90]]$X20,alphacmoutput[[91]]$X20,alphacmoutput[[92]]$X20,alphacmoutput[[93]]$X20,alphacmoutput[[94]]$X20,alphacmoutput[[95]]$X20,alphacmoutput[[96]]$X20,alphacmoutput[[97]]$X20,alphacmoutput[[98]]$X20,alphacmoutput[[99]]$X20,alphacmoutput[[100]]$X20,alphacmoutput[[101]]$X20,alphacmoutput[[102]]$X20,alphacmoutput[[103]]$X20,alphacmoutput[[104]]$X20,alphacmoutput[[105]]$X20,alphacmoutput[[106]]$X20,alphacmoutput[[107]]$X20,alphacmoutput[[108]]$X20,alphacmoutput[[109]]$X20,alphacmoutput[[110]]$X20,alphacmoutput[[111]]$X20,alphacmoutput[[112]]$X20,alphacmoutput[[113]]$X20,alphacmoutput[[114]]$X20,alphacmoutput[[115]]$X20,alphacmoutput[[116]]$X20,alphacmoutput[[117]]$X20,alphacmoutput[[118]]$X20,alphacmoutput[[119]]$X20,alphacmoutput[[120]]$X20,alphacmoutput[[121]]$X20)
cmalpha$Y2 <- c(alphacmoutput[[1]]$Y20,alphacmoutput[[2]]$Y20,alphacmoutput[[3]]$Y20,alphacmoutput[[4]]$Y20,alphacmoutput[[5]]$Y20,alphacmoutput[[6]]$Y20,alphacmoutput[[7]]$Y20,alphacmoutput[[8]]$Y20,alphacmoutput[[9]]$Y20,alphacmoutput[[10]]$Y20,alphacmoutput[[11]]$Y20,alphacmoutput[[12]]$Y20,alphacmoutput[[13]]$Y20,alphacmoutput[[14]]$Y20,alphacmoutput[[15]]$Y20,alphacmoutput[[16]]$Y20,alphacmoutput[[17]]$Y20,alphacmoutput[[18]]$Y20,alphacmoutput[[19]]$Y20,alphacmoutput[[20]]$Y20,alphacmoutput[[21]]$Y20,alphacmoutput[[22]]$Y20,alphacmoutput[[23]]$Y20,alphacmoutput[[24]]$Y20,alphacmoutput[[25]]$Y20,alphacmoutput[[26]]$Y20,alphacmoutput[[27]]$Y20,alphacmoutput[[28]]$Y20,alphacmoutput[[29]]$Y20,alphacmoutput[[30]]$Y20,alphacmoutput[[31]]$Y20,alphacmoutput[[32]]$Y20,alphacmoutput[[33]]$Y20,alphacmoutput[[34]]$Y20,alphacmoutput[[35]]$Y20,alphacmoutput[[36]]$Y20,alphacmoutput[[37]]$Y20,alphacmoutput[[38]]$Y20,alphacmoutput[[39]]$Y20,alphacmoutput[[40]]$Y20,alphacmoutput[[41]]$Y20,alphacmoutput[[42]]$Y20,alphacmoutput[[43]]$Y20,alphacmoutput[[44]]$Y20,alphacmoutput[[45]]$Y20,alphacmoutput[[46]]$Y20,alphacmoutput[[47]]$Y20,alphacmoutput[[48]]$Y20,alphacmoutput[[49]]$Y20,alphacmoutput[[50]]$Y20,alphacmoutput[[51]]$Y20,alphacmoutput[[52]]$Y20,alphacmoutput[[53]]$Y20,alphacmoutput[[54]]$Y20,alphacmoutput[[55]]$Y20,alphacmoutput[[56]]$Y20,alphacmoutput[[57]]$Y20,alphacmoutput[[58]]$Y20,alphacmoutput[[59]]$Y20,alphacmoutput[[60]]$Y20,alphacmoutput[[61]]$Y20,alphacmoutput[[62]]$Y20,alphacmoutput[[63]]$Y20,alphacmoutput[[64]]$Y20,alphacmoutput[[65]]$Y20,alphacmoutput[[66]]$Y20,alphacmoutput[[67]]$Y20,alphacmoutput[[68]]$Y20,alphacmoutput[[69]]$Y20,alphacmoutput[[70]]$Y20,alphacmoutput[[71]]$Y20,alphacmoutput[[72]]$Y20,alphacmoutput[[73]]$Y20,alphacmoutput[[74]]$Y20,alphacmoutput[[75]]$Y20,alphacmoutput[[76]]$Y20,alphacmoutput[[77]]$Y20,alphacmoutput[[78]]$Y20,alphacmoutput[[79]]$Y20,alphacmoutput[[80]]$Y20,alphacmoutput[[81]]$Y20,alphacmoutput[[82]]$Y20,alphacmoutput[[83]]$Y20,alphacmoutput[[84]]$Y20,alphacmoutput[[85]]$Y20,alphacmoutput[[86]]$Y20,alphacmoutput[[87]]$Y20,alphacmoutput[[88]]$Y20,alphacmoutput[[89]]$Y20,alphacmoutput[[90]]$Y20,alphacmoutput[[91]]$Y20,alphacmoutput[[92]]$Y20,alphacmoutput[[93]]$Y20,alphacmoutput[[94]]$Y20,alphacmoutput[[95]]$Y20,alphacmoutput[[96]]$Y20,alphacmoutput[[97]]$Y20,alphacmoutput[[98]]$Y20,alphacmoutput[[99]]$Y20,alphacmoutput[[100]]$Y20,alphacmoutput[[101]]$Y20,alphacmoutput[[102]]$Y20,alphacmoutput[[103]]$Y20,alphacmoutput[[104]]$Y20,alphacmoutput[[105]]$Y20,alphacmoutput[[106]]$Y20,alphacmoutput[[107]]$Y20,alphacmoutput[[108]]$Y20,alphacmoutput[[109]]$Y20,alphacmoutput[[110]]$Y20,alphacmoutput[[111]]$Y20,alphacmoutput[[112]]$Y20,alphacmoutput[[113]]$Y20,alphacmoutput[[114]]$Y20,alphacmoutput[[115]]$Y20,alphacmoutput[[116]]$Y20,alphacmoutput[[117]]$Y20,alphacmoutput[[118]]$Y20,alphacmoutput[[119]]$Y20,alphacmoutput[[120]]$Y20,alphacmoutput[[121]]$Y20)
colnames(cmalpha)[colnames(cmalpha)=="value"] <- "X1"
colnames(cmalpha)[colnames(cmalpha)=="variable"] <- "Alpha"
cmalpha$Alpha <- alphaoutput
cmalpha$cM <-cMoutput
cmalpha

df_cmalpha <- as.data.frame(lapply(cmalpha, unlist))
df_cmalpha <- subset(df_cmalpha, select = -c(population))
library(dplyr)
df_cmalpha<- mutate(df_cmalpha, M = (X2+Y2)/(X1+Y1+X2+Y2),R = (X1+Y1)/(X1+Y1+X2+Y2)) # create new columns with Migrant and Resident ratios
df_cmalpha