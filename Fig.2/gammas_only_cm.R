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
#gammaA <- 0.3 # recovery rate in environment A
#gammaB <- 0.5 # recovery rate in environment B 
mu <- 0.01 # mortality rate
alpha <- 0.2 # parasite induced death rate
phi=3 # maximum per capita fecundity rate for each individual
z=0.0005 # density dependent fecundity coefficient 
cm=0.3 # costs of migration
ct=0 # costs of tolerance (let ct == 0 in this case to see the infection factor drives partial migration)

# 100 years loop 
gammasoutput <- vector(mode="list", length=121) # create 20 list dummy vectors for both migrants and residents
gammaAoutput <- vector(mode="list", length=121) # create 20 list dummy vectors for both migrants and residents
gammaBoutput <- vector(mode="list", length=121) # create 20 list dummy vectors for both migrants and residents
gammaAvals <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) # recovery rate in environment A
gammaBvals <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) # recovery rate in environment B
for (i in 1:length(gammaAvals)){
  print(gammaAvals[i])
  for (g in 1:length(gammaBvals)){
    print(gammaBvals[g])
    X10 <- 150 # initial susceptible population (Resident) 
    X20 <- 150 # initial susceptible population (Migrant)
    Y10 <- 100 # initial infected population (Resident)
    Y20 <- 100 # initial infected population (Migrant)
    for (t in 1:100){
      # Environment A (M+R) together - first 6 months
      params1 <- c(beta=beta,gamma=gammaAvals[i],mu=mu,alpha=alpha) # parameters
      xstart1 <- c(X1=X10,X2=X20,Y1=Y10,Y2=Y20) # initial conditions
      gamma1 <- as.data.frame(ode(func = open.sir.model1,y=xstart1,times=times1,parms=params1,method="ode45"))
      
      # Environment A (only R) - Second 6 months
      params2 <- c(beta=beta,gamma=gammaAvals[i],mu=mu,alpha=alpha) # parameters
      xstart2 <- c(X1=gamma1[7,2],Y1=gamma1[7,4]) # initial conditions
      gamma2 <- as.data.frame(ode(func = open.sir.model2,y=xstart2,times=times2,parms=params2,method="ode45"))
      
      # Environment B for migrants - Second 6 months 
      params3 <- c(beta=beta,gamma=gammaBvals[g],mu=mu,alpha=alpha) 
      xstart3 <- c(X2=gamma1[7,3],Y2=gamma1[7,5]) # Initial survived population from env A
      gamma3 <- as.data.frame(ode(func = open.sir.model3,y=xstart3,times=times3,parms=params3,method="ode45"))
      
      # Calculation of fecundity rates in 2nd year 
      
      # the annual survivors 
      
      X1hat <- gamma2[7,2] # resdient susceptible
      X2hat <- gamma3[7,2] # migrant susceptible
      Y1hat <- gamma2[7,3] # resident infected
      Y2hat <- gamma3[7,3] # migrant susceptible 
      
      bRmax <- phi*(X1hat+Y1hat)*(1-ct) # maximum fecundity rate of a resident individual 
      bMmax <- phi*(X2hat+Y2hat)*(1-cm)*(1-ct) # maximum fecundity rate of a migrant individual 
      
      b1 <- bRmax*(exp(-z*(bRmax+bMmax))) # actual fecundity rate of a resident individual
      b2 <- bMmax*(exp(-z*(bRmax+bMmax))) # actual fecundity rate of a migrant individual 
      
      # after fecundity - 2nd yr initial populations 
      X10 <- gamma2[7,2] + b1
      X20 <- gamma3[7,2] + b2
      Y10 <- gamma2[7,3]
      Y20 <- gamma3[7,3]
    }
    gammasresult <- data.frame(X10,X20,Y10,Y20)
    print(gammasresult)
    
    gammasoutput[[(g-1)*11+i]]<-as.data.frame(gammasresult)
    gammaBoutput[[(g-1)*11+i]]<-as.data.frame(gammaBvals[g])
    gammaAoutput[[(g-1)*11+i]]<-as.data.frame(gammaAvals[i])}} # CT loop

gammas <- data.frame(gammasoutput[[1]]$X10,gammasoutput[[2]]$X10,gammasoutput[[3]]$X10,gammasoutput[[4]]$X10,gammasoutput[[5]]$X10,gammasoutput[[6]]$X10,gammasoutput[[7]]$X10,gammasoutput[[8]]$X10,gammasoutput[[9]]$X10,gammasoutput[[10]]$X10,gammasoutput[[11]]$X10,gammasoutput[[12]]$X10,gammasoutput[[13]]$X10,gammasoutput[[14]]$X10,gammasoutput[[15]]$X10,gammasoutput[[16]]$X10,gammasoutput[[17]]$X10,gammasoutput[[18]]$X10,gammasoutput[[19]]$X10,gammasoutput[[20]]$X10,gammasoutput[[21]]$X10,gammasoutput[[22]]$X10,gammasoutput[[23]]$X10,gammasoutput[[24]]$X10,gammasoutput[[25]]$X10,gammasoutput[[26]]$X10,gammasoutput[[27]]$X10,gammasoutput[[28]]$X10,gammasoutput[[29]]$X10,gammasoutput[[30]]$X10,gammasoutput[[31]]$X10,gammasoutput[[32]]$X10,gammasoutput[[33]]$X10,gammasoutput[[34]]$X10,gammasoutput[[35]]$X10,gammasoutput[[36]]$X10,gammasoutput[[37]]$X10,gammasoutput[[38]]$X10,gammasoutput[[39]]$X10,gammasoutput[[40]]$X10,gammasoutput[[41]]$X10,gammasoutput[[42]]$X10,gammasoutput[[43]]$X10,gammasoutput[[44]]$X10,gammasoutput[[45]]$X10,gammasoutput[[46]]$X10,gammasoutput[[47]]$X10,gammasoutput[[48]]$X10,gammasoutput[[49]]$X10,gammasoutput[[50]]$X10,gammasoutput[[51]]$X10,gammasoutput[[52]]$X10,gammasoutput[[53]]$X10,gammasoutput[[54]]$X10,gammasoutput[[55]]$X10,gammasoutput[[56]]$X10,gammasoutput[[57]]$X10,gammasoutput[[58]]$X10,gammasoutput[[59]]$X10,gammasoutput[[60]]$X10,gammasoutput[[61]]$X10,gammasoutput[[62]]$X10,gammasoutput[[63]]$X10,gammasoutput[[64]]$X10,gammasoutput[[65]]$X10,gammasoutput[[66]]$X10,gammasoutput[[67]]$X10,gammasoutput[[68]]$X10,gammasoutput[[69]]$X10,gammasoutput[[70]]$X10,gammasoutput[[71]]$X10,gammasoutput[[72]]$X10,gammasoutput[[73]]$X10,gammasoutput[[74]]$X10,gammasoutput[[75]]$X10,gammasoutput[[76]]$X10,gammasoutput[[77]]$X10,gammasoutput[[78]]$X10,gammasoutput[[79]]$X10,gammasoutput[[80]]$X10,gammasoutput[[81]]$X10,gammasoutput[[82]]$X10,gammasoutput[[83]]$X10,gammasoutput[[84]]$X10,gammasoutput[[85]]$X10,gammasoutput[[86]]$X10,gammasoutput[[87]]$X10,gammasoutput[[88]]$X10,gammasoutput[[89]]$X10,gammasoutput[[90]]$X10,gammasoutput[[91]]$X10,gammasoutput[[92]]$X10,gammasoutput[[93]]$X10,gammasoutput[[94]]$X10,gammasoutput[[95]]$X10,gammasoutput[[96]]$X10,gammasoutput[[97]]$X10,gammasoutput[[98]]$X10,gammasoutput[[99]]$X10,gammasoutput[[100]]$X10,gammasoutput[[101]]$X10,gammasoutput[[102]]$X10,gammasoutput[[103]]$X10,gammasoutput[[104]]$X10,gammasoutput[[105]]$X10,gammasoutput[[106]]$X10,gammasoutput[[107]]$X10,gammasoutput[[108]]$X10,gammasoutput[[109]]$X10,gammasoutput[[110]]$X10,gammasoutput[[111]]$X10,gammasoutput[[112]]$X10,gammasoutput[[113]]$X10,gammasoutput[[114]]$X10,gammasoutput[[115]]$X10,gammasoutput[[116]]$X10,gammasoutput[[117]]$X10,gammasoutput[[118]]$X10,gammasoutput[[119]]$X10,gammasoutput[[120]]$X10,gammasoutput[[121]]$X10)
gammas <- melt(gammas, variable.name = "population", value.name = "X1")
gammas$Y1 <- c(gammasoutput[[1]]$Y10,gammasoutput[[2]]$Y10,gammasoutput[[3]]$Y10,gammasoutput[[4]]$Y10,gammasoutput[[5]]$Y10,gammasoutput[[6]]$Y10,gammasoutput[[7]]$Y10,gammasoutput[[8]]$Y10,gammasoutput[[9]]$Y10,gammasoutput[[10]]$Y10,gammasoutput[[11]]$Y10,gammasoutput[[12]]$Y10,gammasoutput[[13]]$Y10,gammasoutput[[14]]$Y10,gammasoutput[[15]]$Y10,gammasoutput[[16]]$Y10,gammasoutput[[17]]$Y10,gammasoutput[[18]]$Y10,gammasoutput[[19]]$Y10,gammasoutput[[20]]$Y10,gammasoutput[[21]]$Y10,gammasoutput[[22]]$Y10,gammasoutput[[23]]$Y10,gammasoutput[[24]]$Y10,gammasoutput[[25]]$Y10,gammasoutput[[26]]$Y10,gammasoutput[[27]]$Y10,gammasoutput[[28]]$Y10,gammasoutput[[29]]$Y10,gammasoutput[[30]]$Y10,gammasoutput[[31]]$Y10,gammasoutput[[32]]$Y10,gammasoutput[[33]]$Y10,gammasoutput[[34]]$Y10,gammasoutput[[35]]$Y10,gammasoutput[[36]]$Y10,gammasoutput[[37]]$Y10,gammasoutput[[38]]$Y10,gammasoutput[[39]]$Y10,gammasoutput[[40]]$Y10,gammasoutput[[41]]$Y10,gammasoutput[[42]]$Y10,gammasoutput[[43]]$Y10,gammasoutput[[44]]$Y10,gammasoutput[[45]]$Y10,gammasoutput[[46]]$Y10,gammasoutput[[47]]$Y10,gammasoutput[[48]]$Y10,gammasoutput[[49]]$Y10,gammasoutput[[50]]$Y10,gammasoutput[[51]]$Y10,gammasoutput[[52]]$Y10,gammasoutput[[53]]$Y10,gammasoutput[[54]]$Y10,gammasoutput[[55]]$Y10,gammasoutput[[56]]$Y10,gammasoutput[[57]]$Y10,gammasoutput[[58]]$Y10,gammasoutput[[59]]$Y10,gammasoutput[[60]]$Y10,gammasoutput[[61]]$Y10,gammasoutput[[62]]$Y10,gammasoutput[[63]]$Y10,gammasoutput[[64]]$Y10,gammasoutput[[65]]$Y10,gammasoutput[[66]]$Y10,gammasoutput[[67]]$Y10,gammasoutput[[68]]$Y10,gammasoutput[[69]]$Y10,gammasoutput[[70]]$Y10,gammasoutput[[71]]$Y10,gammasoutput[[72]]$Y10,gammasoutput[[73]]$Y10,gammasoutput[[74]]$Y10,gammasoutput[[75]]$Y10,gammasoutput[[76]]$Y10,gammasoutput[[77]]$Y10,gammasoutput[[78]]$Y10,gammasoutput[[79]]$Y10,gammasoutput[[80]]$Y10,gammasoutput[[81]]$Y10,gammasoutput[[82]]$Y10,gammasoutput[[83]]$Y10,gammasoutput[[84]]$Y10,gammasoutput[[85]]$Y10,gammasoutput[[86]]$Y10,gammasoutput[[87]]$Y10,gammasoutput[[88]]$Y10,gammasoutput[[89]]$Y10,gammasoutput[[90]]$Y10,gammasoutput[[91]]$Y10,gammasoutput[[92]]$Y10,gammasoutput[[93]]$Y10,gammasoutput[[94]]$Y10,gammasoutput[[95]]$Y10,gammasoutput[[96]]$Y10,gammasoutput[[97]]$Y10,gammasoutput[[98]]$Y10,gammasoutput[[99]]$Y10,gammasoutput[[100]]$Y10,gammasoutput[[101]]$Y10,gammasoutput[[102]]$Y10,gammasoutput[[103]]$Y10,gammasoutput[[104]]$Y10,gammasoutput[[105]]$Y10,gammasoutput[[106]]$Y10,gammasoutput[[107]]$Y10,gammasoutput[[108]]$Y10,gammasoutput[[109]]$Y10,gammasoutput[[110]]$Y10,gammasoutput[[111]]$Y10,gammasoutput[[112]]$Y10,gammasoutput[[113]]$Y10,gammasoutput[[114]]$Y10,gammasoutput[[115]]$Y10,gammasoutput[[116]]$Y10,gammasoutput[[117]]$Y10,gammasoutput[[118]]$Y10,gammasoutput[[119]]$Y10,gammasoutput[[120]]$Y10,gammasoutput[[121]]$Y10)
gammas$X2 <- c(gammasoutput[[1]]$X20,gammasoutput[[2]]$X20,gammasoutput[[3]]$X20,gammasoutput[[4]]$X20,gammasoutput[[5]]$X20,gammasoutput[[6]]$X20,gammasoutput[[7]]$X20,gammasoutput[[8]]$X20,gammasoutput[[9]]$X20,gammasoutput[[10]]$X20,gammasoutput[[11]]$X20,gammasoutput[[12]]$X20,gammasoutput[[13]]$X20,gammasoutput[[14]]$X20,gammasoutput[[15]]$X20,gammasoutput[[16]]$X20,gammasoutput[[17]]$X20,gammasoutput[[18]]$X20,gammasoutput[[19]]$X20,gammasoutput[[20]]$X20,gammasoutput[[21]]$X20,gammasoutput[[22]]$X20,gammasoutput[[23]]$X20,gammasoutput[[24]]$X20,gammasoutput[[25]]$X20,gammasoutput[[26]]$X20,gammasoutput[[27]]$X20,gammasoutput[[28]]$X20,gammasoutput[[29]]$X20,gammasoutput[[30]]$X20,gammasoutput[[31]]$X20,gammasoutput[[32]]$X20,gammasoutput[[33]]$X20,gammasoutput[[34]]$X20,gammasoutput[[35]]$X20,gammasoutput[[36]]$X20,gammasoutput[[37]]$X20,gammasoutput[[38]]$X20,gammasoutput[[39]]$X20,gammasoutput[[40]]$X20,gammasoutput[[41]]$X20,gammasoutput[[42]]$X20,gammasoutput[[43]]$X20,gammasoutput[[44]]$X20,gammasoutput[[45]]$X20,gammasoutput[[46]]$X20,gammasoutput[[47]]$X20,gammasoutput[[48]]$X20,gammasoutput[[49]]$X20,gammasoutput[[50]]$X20,gammasoutput[[51]]$X20,gammasoutput[[52]]$X20,gammasoutput[[53]]$X20,gammasoutput[[54]]$X20,gammasoutput[[55]]$X20,gammasoutput[[56]]$X20,gammasoutput[[57]]$X20,gammasoutput[[58]]$X20,gammasoutput[[59]]$X20,gammasoutput[[60]]$X20,gammasoutput[[61]]$X20,gammasoutput[[62]]$X20,gammasoutput[[63]]$X20,gammasoutput[[64]]$X20,gammasoutput[[65]]$X20,gammasoutput[[66]]$X20,gammasoutput[[67]]$X20,gammasoutput[[68]]$X20,gammasoutput[[69]]$X20,gammasoutput[[70]]$X20,gammasoutput[[71]]$X20,gammasoutput[[72]]$X20,gammasoutput[[73]]$X20,gammasoutput[[74]]$X20,gammasoutput[[75]]$X20,gammasoutput[[76]]$X20,gammasoutput[[77]]$X20,gammasoutput[[78]]$X20,gammasoutput[[79]]$X20,gammasoutput[[80]]$X20,gammasoutput[[81]]$X20,gammasoutput[[82]]$X20,gammasoutput[[83]]$X20,gammasoutput[[84]]$X20,gammasoutput[[85]]$X20,gammasoutput[[86]]$X20,gammasoutput[[87]]$X20,gammasoutput[[88]]$X20,gammasoutput[[89]]$X20,gammasoutput[[90]]$X20,gammasoutput[[91]]$X20,gammasoutput[[92]]$X20,gammasoutput[[93]]$X20,gammasoutput[[94]]$X20,gammasoutput[[95]]$X20,gammasoutput[[96]]$X20,gammasoutput[[97]]$X20,gammasoutput[[98]]$X20,gammasoutput[[99]]$X20,gammasoutput[[100]]$X20,gammasoutput[[101]]$X20,gammasoutput[[102]]$X20,gammasoutput[[103]]$X20,gammasoutput[[104]]$X20,gammasoutput[[105]]$X20,gammasoutput[[106]]$X20,gammasoutput[[107]]$X20,gammasoutput[[108]]$X20,gammasoutput[[109]]$X20,gammasoutput[[110]]$X20,gammasoutput[[111]]$X20,gammasoutput[[112]]$X20,gammasoutput[[113]]$X20,gammasoutput[[114]]$X20,gammasoutput[[115]]$X20,gammasoutput[[116]]$X20,gammasoutput[[117]]$X20,gammasoutput[[118]]$X20,gammasoutput[[119]]$X20,gammasoutput[[120]]$X20,gammasoutput[[121]]$X20)
gammas$Y2 <- c(gammasoutput[[1]]$Y20,gammasoutput[[2]]$Y20,gammasoutput[[3]]$Y20,gammasoutput[[4]]$Y20,gammasoutput[[5]]$Y20,gammasoutput[[6]]$Y20,gammasoutput[[7]]$Y20,gammasoutput[[8]]$Y20,gammasoutput[[9]]$Y20,gammasoutput[[10]]$Y20,gammasoutput[[11]]$Y20,gammasoutput[[12]]$Y20,gammasoutput[[13]]$Y20,gammasoutput[[14]]$Y20,gammasoutput[[15]]$Y20,gammasoutput[[16]]$Y20,gammasoutput[[17]]$Y20,gammasoutput[[18]]$Y20,gammasoutput[[19]]$Y20,gammasoutput[[20]]$Y20,gammasoutput[[21]]$Y20,gammasoutput[[22]]$Y20,gammasoutput[[23]]$Y20,gammasoutput[[24]]$Y20,gammasoutput[[25]]$Y20,gammasoutput[[26]]$Y20,gammasoutput[[27]]$Y20,gammasoutput[[28]]$Y20,gammasoutput[[29]]$Y20,gammasoutput[[30]]$Y20,gammasoutput[[31]]$Y20,gammasoutput[[32]]$Y20,gammasoutput[[33]]$Y20,gammasoutput[[34]]$Y20,gammasoutput[[35]]$Y20,gammasoutput[[36]]$Y20,gammasoutput[[37]]$Y20,gammasoutput[[38]]$Y20,gammasoutput[[39]]$Y20,gammasoutput[[40]]$Y20,gammasoutput[[41]]$Y20,gammasoutput[[42]]$Y20,gammasoutput[[43]]$Y20,gammasoutput[[44]]$Y20,gammasoutput[[45]]$Y20,gammasoutput[[46]]$Y20,gammasoutput[[47]]$Y20,gammasoutput[[48]]$Y20,gammasoutput[[49]]$Y20,gammasoutput[[50]]$Y20,gammasoutput[[51]]$Y20,gammasoutput[[52]]$Y20,gammasoutput[[53]]$Y20,gammasoutput[[54]]$Y20,gammasoutput[[55]]$Y20,gammasoutput[[56]]$Y20,gammasoutput[[57]]$Y20,gammasoutput[[58]]$Y20,gammasoutput[[59]]$Y20,gammasoutput[[60]]$Y20,gammasoutput[[61]]$Y20,gammasoutput[[62]]$Y20,gammasoutput[[63]]$Y20,gammasoutput[[64]]$Y20,gammasoutput[[65]]$Y20,gammasoutput[[66]]$Y20,gammasoutput[[67]]$Y20,gammasoutput[[68]]$Y20,gammasoutput[[69]]$Y20,gammasoutput[[70]]$Y20,gammasoutput[[71]]$Y20,gammasoutput[[72]]$Y20,gammasoutput[[73]]$Y20,gammasoutput[[74]]$Y20,gammasoutput[[75]]$Y20,gammasoutput[[76]]$Y20,gammasoutput[[77]]$Y20,gammasoutput[[78]]$Y20,gammasoutput[[79]]$Y20,gammasoutput[[80]]$Y20,gammasoutput[[81]]$Y20,gammasoutput[[82]]$Y20,gammasoutput[[83]]$Y20,gammasoutput[[84]]$Y20,gammasoutput[[85]]$Y20,gammasoutput[[86]]$Y20,gammasoutput[[87]]$Y20,gammasoutput[[88]]$Y20,gammasoutput[[89]]$Y20,gammasoutput[[90]]$Y20,gammasoutput[[91]]$Y20,gammasoutput[[92]]$Y20,gammasoutput[[93]]$Y20,gammasoutput[[94]]$Y20,gammasoutput[[95]]$Y20,gammasoutput[[96]]$Y20,gammasoutput[[97]]$Y20,gammasoutput[[98]]$Y20,gammasoutput[[99]]$Y20,gammasoutput[[100]]$Y20,gammasoutput[[101]]$Y20,gammasoutput[[102]]$Y20,gammasoutput[[103]]$Y20,gammasoutput[[104]]$Y20,gammasoutput[[105]]$Y20,gammasoutput[[106]]$Y20,gammasoutput[[107]]$Y20,gammasoutput[[108]]$Y20,gammasoutput[[109]]$Y20,gammasoutput[[110]]$Y20,gammasoutput[[111]]$Y20,gammasoutput[[112]]$Y20,gammasoutput[[113]]$Y20,gammasoutput[[114]]$Y20,gammasoutput[[115]]$Y20,gammasoutput[[116]]$Y20,gammasoutput[[117]]$Y20,gammasoutput[[118]]$Y20,gammasoutput[[119]]$Y20,gammasoutput[[120]]$Y20,gammasoutput[[121]]$Y20)
colnames(gammas)[colnames(gammas)=="value"] <- "X1"
colnames(gammas)[colnames(gammas)=="variable"] <- "gammaA"
gammas$gammaA <- gammaAoutput
gammas$gammaB <-gammaBoutput
gammas

df_gammas <- as.data.frame(lapply(gammas, unlist))
df_gammas <- subset(df_gammas, select = -c(population))
df_gammas<- mutate(df_gammas, M = (X2+Y2)/(X1+Y1+X2+Y2),R = (X1+Y1)/(X1+Y1+X2+Y2)) # create new columns with Migrant and Resident ratios
df_gammas