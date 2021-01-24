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
z=0.01 # density dependent fecundity coefficient 
#cm=0.3 # costs of migration
#ct=0.5 # costs of tolerance



# 100 years loop 
cTcMoutput <- vector(mode="list", length=121) # create 20 list dummy vectors for both migrants and residents
alphaoutput <- vector(mode="list", length=121) # create 20 list dummy vectors for both migrants and residents
cMoutput <- vector(mode="list", length=121) # create 20 list dummy vectors for both migrants and residents
cMvals <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) # cM variations
alphavals <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) # alpha variations
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
      cTcM1 <- as.data.frame(ode(func = open.sir.model1,y=xstart1,times=times1,parms=params1,method="ode23")) 
      #method = "ode45",atol=1e-14,rtol=1e-14
      
      # Environment A (only R) - Second 6 months
      params2 <- c(beta=beta,gamma=gammaA,mu=mu,alpha=alphavals[i]) # parameters
      xstart2 <- c(X1=cTcM1[7,2],Y1=cTcM1[7,4]) # initial conditions
      cTcM2 <- as.data.frame(ode(func = open.sir.model2,y=xstart2,times=times2,parms=params2,method="ode23"))
      
      # Environment B for migrants - Second 6 months 
      params3 <- c(beta=beta,gamma=gammaB,mu=mu,alpha=alphavals[i]) 
      xstart3 <- c(X2=cTcM1[7,3],Y2=cTcM1[7,5]) # Initial survived population from env A
      cTcM3 <- as.data.frame(ode(func = open.sir.model3,y=xstart3,times=times3,parms=params3,method="ode23"))
      
      # Calculation of fecundity rates in 2nd year 
      #  params4 <- c(phi=phi,z=z,cM=cMvals[g])
      #  xstart4 <- c(X1=cTcM2[7,2],X2=cTcM3[7,2],Y1=cTcM2[7,3],Y2=cTcM3[7,3])
      #  b <- as.data.frame(ode(func = open.sir.model4,y=xstart4,times=times4,parms=params4,method="ode23"))
      
      # the annual survivors 
      
      X1hat <- cTcM2[7,2] # resdient susceptible
      X2hat <- cTcM3[7,2] # migrant susceptible
      Y1hat <- cTcM2[7,3] # resident infected
      Y2hat <- cTcM3[7,3] # migrant susceptible
      
      #cT <- if (alphavals[i]==1){0}else if (alphavals[i]==0.9){0.1}else if (alphavals[i]==0.8){0.2}else if (alphavals[i]==0.7){0.3}else if (alphavals[i]==0.6){0.4}else if (alphavals[i]==0.5){0.5}else if (alphavals[i]==0.4){0.6}else if (alphavals[i]==0.3){0.7}else if (alphavals[i]==0.2){0.8}else if (alphavals[i]==0.1){0.9}else if (alphavals[i]==0){1} # cT variations
      
      cT <- 1-alphavals[i]
      
      bRmax <- phi*(X1hat+Y1hat)*(1-cT) # maximum fecundity rate of a resident individual 
      bMmax <- phi*(X2hat+Y2hat)*(1-cMvals[g])*(1-cT) # maximum fecundity rate of a migrant individual 
      
      b1 <- bRmax*(exp(-z*(bRmax+bMmax))) # actual fecundity rate of a resident individual
      b2 <- bMmax*(exp(-z*(bRmax+bMmax))) # actual fecundity rate of a migrant individual 
      
      # after fecundity - 2nd yr initial populations 
      X10 <- cTcM2[7,2] + b1
      X20 <- cTcM3[7,2] + b2
      Y10 <- cTcM2[7,3]
      Y20 <- cTcM3[7,3]
    }
    cTcMresult <- data.frame(X10,X20,Y10,Y20)
    print(cTcMresult)
    
    cTcMoutput[[(g-1)*11+i]]<-as.data.frame(cTcMresult)
    cMoutput[[(g-1)*11+i]]<-as.data.frame(cMvals[g])
    alphaoutput[[(g-1)*11+i]]<-as.data.frame(alphavals[i])}} # CT loop


cTcM <- data.frame(cTcMoutput[[1]]$X10,cTcMoutput[[2]]$X10,cTcMoutput[[3]]$X10,cTcMoutput[[4]]$X10,cTcMoutput[[5]]$X10,cTcMoutput[[6]]$X10,cTcMoutput[[7]]$X10,cTcMoutput[[8]]$X10,cTcMoutput[[9]]$X10,cTcMoutput[[10]]$X10,cTcMoutput[[11]]$X10,cTcMoutput[[12]]$X10,cTcMoutput[[13]]$X10,cTcMoutput[[14]]$X10,cTcMoutput[[15]]$X10,cTcMoutput[[16]]$X10,cTcMoutput[[17]]$X10,cTcMoutput[[18]]$X10,cTcMoutput[[19]]$X10,cTcMoutput[[20]]$X10,cTcMoutput[[21]]$X10,cTcMoutput[[22]]$X10,cTcMoutput[[23]]$X10,cTcMoutput[[24]]$X10,cTcMoutput[[25]]$X10,cTcMoutput[[26]]$X10,cTcMoutput[[27]]$X10,cTcMoutput[[28]]$X10,cTcMoutput[[29]]$X10,cTcMoutput[[30]]$X10,cTcMoutput[[31]]$X10,cTcMoutput[[32]]$X10,cTcMoutput[[33]]$X10,cTcMoutput[[34]]$X10,cTcMoutput[[35]]$X10,cTcMoutput[[36]]$X10,cTcMoutput[[37]]$X10,cTcMoutput[[38]]$X10,cTcMoutput[[39]]$X10,cTcMoutput[[40]]$X10,cTcMoutput[[41]]$X10,cTcMoutput[[42]]$X10,cTcMoutput[[43]]$X10,cTcMoutput[[44]]$X10,cTcMoutput[[45]]$X10,cTcMoutput[[46]]$X10,cTcMoutput[[47]]$X10,cTcMoutput[[48]]$X10,cTcMoutput[[49]]$X10,cTcMoutput[[50]]$X10,cTcMoutput[[51]]$X10,cTcMoutput[[52]]$X10,cTcMoutput[[53]]$X10,cTcMoutput[[54]]$X10,cTcMoutput[[55]]$X10,cTcMoutput[[56]]$X10,cTcMoutput[[57]]$X10,cTcMoutput[[58]]$X10,cTcMoutput[[59]]$X10,cTcMoutput[[60]]$X10,cTcMoutput[[61]]$X10,cTcMoutput[[62]]$X10,cTcMoutput[[63]]$X10,cTcMoutput[[64]]$X10,cTcMoutput[[65]]$X10,cTcMoutput[[66]]$X10,cTcMoutput[[67]]$X10,cTcMoutput[[68]]$X10,cTcMoutput[[69]]$X10,cTcMoutput[[70]]$X10,cTcMoutput[[71]]$X10,cTcMoutput[[72]]$X10,cTcMoutput[[73]]$X10,cTcMoutput[[74]]$X10,cTcMoutput[[75]]$X10,cTcMoutput[[76]]$X10,cTcMoutput[[77]]$X10,cTcMoutput[[78]]$X10,cTcMoutput[[79]]$X10,cTcMoutput[[80]]$X10,cTcMoutput[[81]]$X10,cTcMoutput[[82]]$X10,cTcMoutput[[83]]$X10,cTcMoutput[[84]]$X10,cTcMoutput[[85]]$X10,cTcMoutput[[86]]$X10,cTcMoutput[[87]]$X10,cTcMoutput[[88]]$X10,cTcMoutput[[89]]$X10,cTcMoutput[[90]]$X10,cTcMoutput[[91]]$X10,cTcMoutput[[92]]$X10,cTcMoutput[[93]]$X10,cTcMoutput[[94]]$X10,cTcMoutput[[95]]$X10,cTcMoutput[[96]]$X10,cTcMoutput[[97]]$X10,cTcMoutput[[98]]$X10,cTcMoutput[[99]]$X10,cTcMoutput[[100]]$X10,cTcMoutput[[101]]$X10,cTcMoutput[[102]]$X10,cTcMoutput[[103]]$X10,cTcMoutput[[104]]$X10,cTcMoutput[[105]]$X10,cTcMoutput[[106]]$X10,cTcMoutput[[107]]$X10,cTcMoutput[[108]]$X10,cTcMoutput[[109]]$X10,cTcMoutput[[110]]$X10,cTcMoutput[[111]]$X10,cTcMoutput[[112]]$X10,cTcMoutput[[113]]$X10,cTcMoutput[[114]]$X10,cTcMoutput[[115]]$X10,cTcMoutput[[116]]$X10,cTcMoutput[[117]]$X10,cTcMoutput[[118]]$X10,cTcMoutput[[119]]$X10,cTcMoutput[[120]]$X10,cTcMoutput[[121]]$X10)
cTcM <- melt(cTcM, variable.name = "population", value.name = "X1")
cTcM$Y1 <- c(cTcMoutput[[1]]$Y10,cTcMoutput[[2]]$Y10,cTcMoutput[[3]]$Y10,cTcMoutput[[4]]$Y10,cTcMoutput[[5]]$Y10,cTcMoutput[[6]]$Y10,cTcMoutput[[7]]$Y10,cTcMoutput[[8]]$Y10,cTcMoutput[[9]]$Y10,cTcMoutput[[10]]$Y10,cTcMoutput[[11]]$Y10,cTcMoutput[[12]]$Y10,cTcMoutput[[13]]$Y10,cTcMoutput[[14]]$Y10,cTcMoutput[[15]]$Y10,cTcMoutput[[16]]$Y10,cTcMoutput[[17]]$Y10,cTcMoutput[[18]]$Y10,cTcMoutput[[19]]$Y10,cTcMoutput[[20]]$Y10,cTcMoutput[[21]]$Y10,cTcMoutput[[22]]$Y10,cTcMoutput[[23]]$Y10,cTcMoutput[[24]]$Y10,cTcMoutput[[25]]$Y10,cTcMoutput[[26]]$Y10,cTcMoutput[[27]]$Y10,cTcMoutput[[28]]$Y10,cTcMoutput[[29]]$Y10,cTcMoutput[[30]]$Y10,cTcMoutput[[31]]$Y10,cTcMoutput[[32]]$Y10,cTcMoutput[[33]]$Y10,cTcMoutput[[34]]$Y10,cTcMoutput[[35]]$Y10,cTcMoutput[[36]]$Y10,cTcMoutput[[37]]$Y10,cTcMoutput[[38]]$Y10,cTcMoutput[[39]]$Y10,cTcMoutput[[40]]$Y10,cTcMoutput[[41]]$Y10,cTcMoutput[[42]]$Y10,cTcMoutput[[43]]$Y10,cTcMoutput[[44]]$Y10,cTcMoutput[[45]]$Y10,cTcMoutput[[46]]$Y10,cTcMoutput[[47]]$Y10,cTcMoutput[[48]]$Y10,cTcMoutput[[49]]$Y10,cTcMoutput[[50]]$Y10,cTcMoutput[[51]]$Y10,cTcMoutput[[52]]$Y10,cTcMoutput[[53]]$Y10,cTcMoutput[[54]]$Y10,cTcMoutput[[55]]$Y10,cTcMoutput[[56]]$Y10,cTcMoutput[[57]]$Y10,cTcMoutput[[58]]$Y10,cTcMoutput[[59]]$Y10,cTcMoutput[[60]]$Y10,cTcMoutput[[61]]$Y10,cTcMoutput[[62]]$Y10,cTcMoutput[[63]]$Y10,cTcMoutput[[64]]$Y10,cTcMoutput[[65]]$Y10,cTcMoutput[[66]]$Y10,cTcMoutput[[67]]$Y10,cTcMoutput[[68]]$Y10,cTcMoutput[[69]]$Y10,cTcMoutput[[70]]$Y10,cTcMoutput[[71]]$Y10,cTcMoutput[[72]]$Y10,cTcMoutput[[73]]$Y10,cTcMoutput[[74]]$Y10,cTcMoutput[[75]]$Y10,cTcMoutput[[76]]$Y10,cTcMoutput[[77]]$Y10,cTcMoutput[[78]]$Y10,cTcMoutput[[79]]$Y10,cTcMoutput[[80]]$Y10,cTcMoutput[[81]]$Y10,cTcMoutput[[82]]$Y10,cTcMoutput[[83]]$Y10,cTcMoutput[[84]]$Y10,cTcMoutput[[85]]$Y10,cTcMoutput[[86]]$Y10,cTcMoutput[[87]]$Y10,cTcMoutput[[88]]$Y10,cTcMoutput[[89]]$Y10,cTcMoutput[[90]]$Y10,cTcMoutput[[91]]$Y10,cTcMoutput[[92]]$Y10,cTcMoutput[[93]]$Y10,cTcMoutput[[94]]$Y10,cTcMoutput[[95]]$Y10,cTcMoutput[[96]]$Y10,cTcMoutput[[97]]$Y10,cTcMoutput[[98]]$Y10,cTcMoutput[[99]]$Y10,cTcMoutput[[100]]$Y10,cTcMoutput[[101]]$Y10,cTcMoutput[[102]]$Y10,cTcMoutput[[103]]$Y10,cTcMoutput[[104]]$Y10,cTcMoutput[[105]]$Y10,cTcMoutput[[106]]$Y10,cTcMoutput[[107]]$Y10,cTcMoutput[[108]]$Y10,cTcMoutput[[109]]$Y10,cTcMoutput[[110]]$Y10,cTcMoutput[[111]]$Y10,cTcMoutput[[112]]$Y10,cTcMoutput[[113]]$Y10,cTcMoutput[[114]]$Y10,cTcMoutput[[115]]$Y10,cTcMoutput[[116]]$Y10,cTcMoutput[[117]]$Y10,cTcMoutput[[118]]$Y10,cTcMoutput[[119]]$Y10,cTcMoutput[[120]]$Y10,cTcMoutput[[121]]$Y10)
cTcM$X2 <- c(cTcMoutput[[1]]$X20,cTcMoutput[[2]]$X20,cTcMoutput[[3]]$X20,cTcMoutput[[4]]$X20,cTcMoutput[[5]]$X20,cTcMoutput[[6]]$X20,cTcMoutput[[7]]$X20,cTcMoutput[[8]]$X20,cTcMoutput[[9]]$X20,cTcMoutput[[10]]$X20,cTcMoutput[[11]]$X20,cTcMoutput[[12]]$X20,cTcMoutput[[13]]$X20,cTcMoutput[[14]]$X20,cTcMoutput[[15]]$X20,cTcMoutput[[16]]$X20,cTcMoutput[[17]]$X20,cTcMoutput[[18]]$X20,cTcMoutput[[19]]$X20,cTcMoutput[[20]]$X20,cTcMoutput[[21]]$X20,cTcMoutput[[22]]$X20,cTcMoutput[[23]]$X20,cTcMoutput[[24]]$X20,cTcMoutput[[25]]$X20,cTcMoutput[[26]]$X20,cTcMoutput[[27]]$X20,cTcMoutput[[28]]$X20,cTcMoutput[[29]]$X20,cTcMoutput[[30]]$X20,cTcMoutput[[31]]$X20,cTcMoutput[[32]]$X20,cTcMoutput[[33]]$X20,cTcMoutput[[34]]$X20,cTcMoutput[[35]]$X20,cTcMoutput[[36]]$X20,cTcMoutput[[37]]$X20,cTcMoutput[[38]]$X20,cTcMoutput[[39]]$X20,cTcMoutput[[40]]$X20,cTcMoutput[[41]]$X20,cTcMoutput[[42]]$X20,cTcMoutput[[43]]$X20,cTcMoutput[[44]]$X20,cTcMoutput[[45]]$X20,cTcMoutput[[46]]$X20,cTcMoutput[[47]]$X20,cTcMoutput[[48]]$X20,cTcMoutput[[49]]$X20,cTcMoutput[[50]]$X20,cTcMoutput[[51]]$X20,cTcMoutput[[52]]$X20,cTcMoutput[[53]]$X20,cTcMoutput[[54]]$X20,cTcMoutput[[55]]$X20,cTcMoutput[[56]]$X20,cTcMoutput[[57]]$X20,cTcMoutput[[58]]$X20,cTcMoutput[[59]]$X20,cTcMoutput[[60]]$X20,cTcMoutput[[61]]$X20,cTcMoutput[[62]]$X20,cTcMoutput[[63]]$X20,cTcMoutput[[64]]$X20,cTcMoutput[[65]]$X20,cTcMoutput[[66]]$X20,cTcMoutput[[67]]$X20,cTcMoutput[[68]]$X20,cTcMoutput[[69]]$X20,cTcMoutput[[70]]$X20,cTcMoutput[[71]]$X20,cTcMoutput[[72]]$X20,cTcMoutput[[73]]$X20,cTcMoutput[[74]]$X20,cTcMoutput[[75]]$X20,cTcMoutput[[76]]$X20,cTcMoutput[[77]]$X20,cTcMoutput[[78]]$X20,cTcMoutput[[79]]$X20,cTcMoutput[[80]]$X20,cTcMoutput[[81]]$X20,cTcMoutput[[82]]$X20,cTcMoutput[[83]]$X20,cTcMoutput[[84]]$X20,cTcMoutput[[85]]$X20,cTcMoutput[[86]]$X20,cTcMoutput[[87]]$X20,cTcMoutput[[88]]$X20,cTcMoutput[[89]]$X20,cTcMoutput[[90]]$X20,cTcMoutput[[91]]$X20,cTcMoutput[[92]]$X20,cTcMoutput[[93]]$X20,cTcMoutput[[94]]$X20,cTcMoutput[[95]]$X20,cTcMoutput[[96]]$X20,cTcMoutput[[97]]$X20,cTcMoutput[[98]]$X20,cTcMoutput[[99]]$X20,cTcMoutput[[100]]$X20,cTcMoutput[[101]]$X20,cTcMoutput[[102]]$X20,cTcMoutput[[103]]$X20,cTcMoutput[[104]]$X20,cTcMoutput[[105]]$X20,cTcMoutput[[106]]$X20,cTcMoutput[[107]]$X20,cTcMoutput[[108]]$X20,cTcMoutput[[109]]$X20,cTcMoutput[[110]]$X20,cTcMoutput[[111]]$X20,cTcMoutput[[112]]$X20,cTcMoutput[[113]]$X20,cTcMoutput[[114]]$X20,cTcMoutput[[115]]$X20,cTcMoutput[[116]]$X20,cTcMoutput[[117]]$X20,cTcMoutput[[118]]$X20,cTcMoutput[[119]]$X20,cTcMoutput[[120]]$X20,cTcMoutput[[121]]$X20)
cTcM$Y2 <- c(cTcMoutput[[1]]$Y20,cTcMoutput[[2]]$Y20,cTcMoutput[[3]]$Y20,cTcMoutput[[4]]$Y20,cTcMoutput[[5]]$Y20,cTcMoutput[[6]]$Y20,cTcMoutput[[7]]$Y20,cTcMoutput[[8]]$Y20,cTcMoutput[[9]]$Y20,cTcMoutput[[10]]$Y20,cTcMoutput[[11]]$Y20,cTcMoutput[[12]]$Y20,cTcMoutput[[13]]$Y20,cTcMoutput[[14]]$Y20,cTcMoutput[[15]]$Y20,cTcMoutput[[16]]$Y20,cTcMoutput[[17]]$Y20,cTcMoutput[[18]]$Y20,cTcMoutput[[19]]$Y20,cTcMoutput[[20]]$Y20,cTcMoutput[[21]]$Y20,cTcMoutput[[22]]$Y20,cTcMoutput[[23]]$Y20,cTcMoutput[[24]]$Y20,cTcMoutput[[25]]$Y20,cTcMoutput[[26]]$Y20,cTcMoutput[[27]]$Y20,cTcMoutput[[28]]$Y20,cTcMoutput[[29]]$Y20,cTcMoutput[[30]]$Y20,cTcMoutput[[31]]$Y20,cTcMoutput[[32]]$Y20,cTcMoutput[[33]]$Y20,cTcMoutput[[34]]$Y20,cTcMoutput[[35]]$Y20,cTcMoutput[[36]]$Y20,cTcMoutput[[37]]$Y20,cTcMoutput[[38]]$Y20,cTcMoutput[[39]]$Y20,cTcMoutput[[40]]$Y20,cTcMoutput[[41]]$Y20,cTcMoutput[[42]]$Y20,cTcMoutput[[43]]$Y20,cTcMoutput[[44]]$Y20,cTcMoutput[[45]]$Y20,cTcMoutput[[46]]$Y20,cTcMoutput[[47]]$Y20,cTcMoutput[[48]]$Y20,cTcMoutput[[49]]$Y20,cTcMoutput[[50]]$Y20,cTcMoutput[[51]]$Y20,cTcMoutput[[52]]$Y20,cTcMoutput[[53]]$Y20,cTcMoutput[[54]]$Y20,cTcMoutput[[55]]$Y20,cTcMoutput[[56]]$Y20,cTcMoutput[[57]]$Y20,cTcMoutput[[58]]$Y20,cTcMoutput[[59]]$Y20,cTcMoutput[[60]]$Y20,cTcMoutput[[61]]$Y20,cTcMoutput[[62]]$Y20,cTcMoutput[[63]]$Y20,cTcMoutput[[64]]$Y20,cTcMoutput[[65]]$Y20,cTcMoutput[[66]]$Y20,cTcMoutput[[67]]$Y20,cTcMoutput[[68]]$Y20,cTcMoutput[[69]]$Y20,cTcMoutput[[70]]$Y20,cTcMoutput[[71]]$Y20,cTcMoutput[[72]]$Y20,cTcMoutput[[73]]$Y20,cTcMoutput[[74]]$Y20,cTcMoutput[[75]]$Y20,cTcMoutput[[76]]$Y20,cTcMoutput[[77]]$Y20,cTcMoutput[[78]]$Y20,cTcMoutput[[79]]$Y20,cTcMoutput[[80]]$Y20,cTcMoutput[[81]]$Y20,cTcMoutput[[82]]$Y20,cTcMoutput[[83]]$Y20,cTcMoutput[[84]]$Y20,cTcMoutput[[85]]$Y20,cTcMoutput[[86]]$Y20,cTcMoutput[[87]]$Y20,cTcMoutput[[88]]$Y20,cTcMoutput[[89]]$Y20,cTcMoutput[[90]]$Y20,cTcMoutput[[91]]$Y20,cTcMoutput[[92]]$Y20,cTcMoutput[[93]]$Y20,cTcMoutput[[94]]$Y20,cTcMoutput[[95]]$Y20,cTcMoutput[[96]]$Y20,cTcMoutput[[97]]$Y20,cTcMoutput[[98]]$Y20,cTcMoutput[[99]]$Y20,cTcMoutput[[100]]$Y20,cTcMoutput[[101]]$Y20,cTcMoutput[[102]]$Y20,cTcMoutput[[103]]$Y20,cTcMoutput[[104]]$Y20,cTcMoutput[[105]]$Y20,cTcMoutput[[106]]$Y20,cTcMoutput[[107]]$Y20,cTcMoutput[[108]]$Y20,cTcMoutput[[109]]$Y20,cTcMoutput[[110]]$Y20,cTcMoutput[[111]]$Y20,cTcMoutput[[112]]$Y20,cTcMoutput[[113]]$Y20,cTcMoutput[[114]]$Y20,cTcMoutput[[115]]$Y20,cTcMoutput[[116]]$Y20,cTcMoutput[[117]]$Y20,cTcMoutput[[118]]$Y20,cTcMoutput[[119]]$Y20,cTcMoutput[[120]]$Y20,cTcMoutput[[121]]$Y20)
colnames(cTcM)[colnames(cTcM)=="value"] <- "X1"
colnames(cTcM)[colnames(cTcM)=="variable"] <- "gammaA"
cTcM$cM <-cMoutput
cTcM$alpha <-alphaoutput
cTcM$cT <- c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0)
cTcM

cTcM2 <- as.data.frame(lapply(cTcM, unlist))
cTcM2 <- subset(cTcM2, select = -c(population))
cTcM3 <- subset(cTcM2, select = -c(alpha))

cTcM3[which(cTcM3$X1<1),1]<-0
cTcM3[which(cTcM3$Y1<1),2]<-0
cTcM3[which(cTcM3$X2<1),3]<-0
cTcM3[which(cTcM3$Y2<1),4]<-0

library(dplyr)
cTcM3<- mutate(cTcM3, M = (X2+Y2)/(X1+Y1+X2+Y2),R = (X1+Y1)/(X1+Y1+X2+Y2)) # create new columns with Migrant and Resident ratios
cTcM3