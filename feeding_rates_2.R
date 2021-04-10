#############################################################################################
#### code to compute feeding rates with fixed number of predators and prey
#############################################################################################


#############################################################################################
#### set language seed and other settings
#### load packages
#############################################################################################

if ( !require("rlist") )   { install.packages("rlist"); library("rlist") }
if ( !require("deSolve") ) { install.packages("deSolve"); library("deSolve") }

rm(list=ls())
Sys.setenv(LANG = "en")
set.seed(98765678)

PLOT_IN_TERMINAL <- F
COMPUTE_AVERAGES <- F
CHECK_MEAN       <- F

#############################################################################################
#### paramteres of the simulations
#############################################################################################



#############################################################################################
#### set spatial lattice, cells and individuaks
#############################################################################################

#### size of the lattice
L <- 1
#### number of cells
M <- L*L
#### initial number of resource individuals
N_R  <- 1000
#### initial number of consumer individuals
N_A  <- 100
#### initial number of compounds
N_C  <- 100

#############################################################################################
#### model parameters: Rates of the process 
#############################################################################################

#### rate of movement of resources
mu_R     <- 0
#### rate of immigration of resources
lambda_R <- 0
#### death rate of resources
delta_R  <- 0.5
### birth rate of resources
beta_R   <- 1.5
#### rate of movement of consumer A
mu_A     <- 0
#### rate of immigration of consumer A
lambda_A <- 0
#### death rate of consumer A
delta_A  <- 1
#### encounter rate between consumer and resource
alpha    <- 2.5
#### handling time
h        <- 1
#### "decomposition rate" of consumer-resource pairs (compounds)
nu       <- 1/h
#### system's Volume (max number of resource individuals per cell)
V        <- 10000

#############################################################################################
#### stochastic simulation settings
#############################################################################################

#### number of events
NEV    <- 10000
### number of replicates
NREP   <- 2000

### transient
TRANS  <- ceiling(50*NEV/100)

#############################################################################################
#### start dynamics
#############################################################################################

#### definelist to put results of simulations
list_simu <- list()
####loop in replicates
kk <- 1
for(kk in 1:NREP){
  list_simu[[kk]] <- list()
  #### initial time
  tt <- 0
  #### define a matrix to describe resources in the lattice (with all zeros)
  spatial_grid_R <-  matrix(nrow = L, ncol = L, 0)
  #### define a matrix to describe consumers in the lattice (with all zeros)
  spatial_grid_A <-  matrix(nrow = L, ncol = L, 0)
  #### define a matrix to describe compounds in the lattice (with all zeros)
  spatial_grid_C <-  matrix(nrow = L, ncol = L, 0)
  #### set initial conditions
  #### put all individuals at the center of the lattice
  spatial_grid_R[ceiling(L/2),ceiling(L/2)] <- V - kk*5
  spatial_grid_A[ceiling(L/2),ceiling(L/2)] <- N_A
  spatial_grid_C[ceiling(L/2),ceiling(L/2)] <- N_C
  
  ##### transition rate at beginning
  # R_movement   <- (mu_R+delta_R)*sum(spatial_grid_R)+(mu_A+delta_A)*sum(spatial_grid_A)
  # R_immigration<- (lambda_R+lambda_A)*M
  # R_feeding    <- (alpha/V)*sum(spatial_grid_A*spatial_grid_R) + nu*sum(spatial_grid_C)
  # 
  # R <- R_movement+R_immigration+R_feeding
  
  ### set counters for feeding rates
  count_r_dd <- 0
  count_r_bb <- 0
  count_r_pp <- 0
  
  count_pp <- 0  
  count_dd <- 0
  count_hh <- 0
  #### loop in events
  jj <- 1
  for(jj in 1:NEV){
    
    #### calculate the total rate of events
    R_movement   <- (mu_R+delta_R)*sum(spatial_grid_R)+(mu_A+delta_A)*sum(spatial_grid_A)
    R_birth_R    <- beta_R*sum(spatial_grid_R*(V - spatial_grid_R)/V)      
    R_immigration<- lambda_R*(M*V - sum(spatial_grid_R)) + lambda_A*M
    R_feeding    <- alpha*sum(spatial_grid_A*spatial_grid_R/V) + nu*sum(spatial_grid_C)
    
    R <- R_movement + R_birth_R + R_immigration + R_feeding
    
    #### draw time step for next event
    pp      <- runif(1,0,1)
    delta_t <- -log(pp)/R
    tt      <- tt + delta_t
    
    #### calculate matrix of rates
    spatial_grid_rates <- (mu_R+delta_R)*spatial_grid_R+(mu_A+delta_A)*spatial_grid_A
    spatial_grid_rates <- spatial_grid_rates + lambda_R*(V-spatial_grid_R) + lambda_A
    spatial_grid_rates <- spatial_grid_rates + beta_R*spatial_grid_R*(V - spatial_grid_R)/V
    spatial_grid_rates <- spatial_grid_rates + alpha*spatial_grid_A*spatial_grid_R/V + nu*spatial_grid_C
    
    #### normalize to obtain the matrix of probabilities
    spatial_grid_pp    <- spatial_grid_rates/sum(spatial_grid_rates)
    
    #### pick a site where the event happens
    repeat{
      rand_col <- sample(1:ncol(spatial_grid_pp),1)
      rand_row <- sample(1:nrow(spatial_grid_pp),1)
      pp <- runif(1,0,1)
      if(spatial_grid_pp[rand_row,rand_col] > pp)break
    }
    #### check the same thing with sample (probabilities from matrix)
    #### site has been selected
    
    #### probabilities of the different events 
    #### movement
    p_m_R    <- mu_R*spatial_grid_R[rand_row,rand_col]/spatial_grid_rates[rand_row,rand_col]
    p_m_A    <- mu_A*spatial_grid_A[rand_row,rand_col]/spatial_grid_rates[rand_row,rand_col]
    #### immigration
    p_i_R    <- lambda_R*(V-spatial_grid_R[rand_row,rand_col])/spatial_grid_rates[rand_row,rand_col]
    p_i_A    <- lambda_A/spatial_grid_rates[rand_row,rand_col]
    #### deaths
    p_d_R    <- delta_R*spatial_grid_R[rand_row,rand_col]/spatial_grid_rates[rand_row,rand_col]
    p_d_A    <- delta_A*spatial_grid_A[rand_row,rand_col]/spatial_grid_rates[rand_row,rand_col]
    #### birth of resources
    p_b_R    <- (beta_R*spatial_grid_R[rand_row,rand_col]*(V - spatial_grid_R[rand_row,rand_col])/V)/spatial_grid_rates[rand_row,rand_col]
    #### feeding
    p_e      <- (alpha/V)*spatial_grid_R[rand_row,rand_col]*spatial_grid_A[rand_row,rand_col]/spatial_grid_rates[rand_row,rand_col]
    p_d      <- nu*spatial_grid_C[rand_row,rand_col]/spatial_grid_rates[rand_row,rand_col]
    #### check normalization to 1
    p_m_R+p_m_A+p_i_R+p_i_A+p_d_R+p_d_A+p_b_R+p_e+p_d
    
    #### now select the event that happens in the site
    pp <- runif(1,0,1)
    
    #### immigration of resources
    p1 <- p_i_R
    if(pp<p1){
      spatial_grid_R[rand_row,rand_col]<- spatial_grid_R[rand_row,rand_col]+1
    }
    #### immigration of consumers
    p2 <- p1 + p_i_A
    if((pp>=p1)&&(pp<p2)){
      spatial_grid_A[rand_row,rand_col]<- spatial_grid_A[rand_row,rand_col]+1
    }
    #### death of resources
    p1 <- p2
    p2 <- p1 + p_d_R
    if((pp>=p1)&&(pp<p2)){
      #spatial_grid_R[rand_row,rand_col]<- spatial_grid_R[rand_row,rand_col]-1
      count_r_dd <- count_r_dd + 1
    }
    #### death of consumers
    p1 <- p2
    p2 <- p1 + p_d_A
    if((pp>=p1)&&(pp<p2)){
      #spatial_grid_A[rand_row,rand_col]<- spatial_grid_A[rand_row,rand_col]-1
      count_dd <- count_dd + 1
    }
    #### birth of resources
    p1 <- p2
    p2 <- p1 + p_b_R
    if((pp>=p1)&&(pp<p2)){
      #spatial_grid_R[rand_row,rand_col]<- spatial_grid_R[rand_row,rand_col]+1
      count_r_bb <- count_r_bb + 1
    }
    #### formation of compounds through encounters
    p1 <- p2
    p2 <- p1 + p_e
    if((pp>=p1)&&(pp<p2)){
      #spatial_grid_R[rand_row,rand_col]<- spatial_grid_R[rand_row,rand_col]-1
      spatial_grid_A[rand_row,rand_col]<- spatial_grid_A[rand_row,rand_col]-1
      spatial_grid_C[rand_row,rand_col]<- spatial_grid_C[rand_row,rand_col]+1
      count_pp   <- count_pp + 1
      count_r_pp <- count_r_pp + 1
    }
    #### decomposition or handling
    p1 <- p2
    p2 <- p1 + p_d
    if((pp>=p1)&&(pp<p2)){
      spatial_grid_C[rand_row,rand_col]<- spatial_grid_C[rand_row,rand_col]-1
      spatial_grid_A[rand_row,rand_col]<- spatial_grid_A[rand_row,rand_col]+1
      count_hh <- count_hh + 1
    }
    #### movement of resources
    p1 <- p2
    p2 <- p1 + p_m_R
    if((pp>=p1)&&(pp<p2)){
      
      ##### remove one individual from the selected site
      spatial_grid_R[rand_row,rand_col] <- spatial_grid_R[rand_row,rand_col]-1
      
      #### select a new site where the individual moves
      pp <- runif(1,0,1)
      
      #### go down
      if(pp<0.25){
        next_row <- rand_row - 1
        next_col <- rand_col
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_R[next_row,next_col] <- spatial_grid_R[next_row,next_col]+1
      }
      #### go up
      if((pp>=0.25)&&(pp<0.5)){
        next_row <- rand_row + 1
        next_col <- rand_col
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_R[next_row,next_col] <- spatial_grid_R[next_row,next_col]+1
      }
      #### go left
      if((pp>=0.5)&&(pp<0.75)){
        next_row <- rand_row
        next_col <- rand_col - 1
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_R[next_row,next_col] <- spatial_grid_R[next_row,next_col]+1
      }
      #### go right
      if(pp>=0.75){
        next_row <- rand_row
        next_col <- rand_col+1
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_R[next_row,next_col] <- spatial_grid_R[next_row,next_col]+1
      }
    }
    #### movement of consumers
    p1 <- p2
    p2 <- p1 + p_m_A
    if((pp>=p1)&&(pp<p2)){
      
      ##### remove one individual from the selected site
      spatial_grid_A[rand_row,rand_col] <- spatial_grid_A[rand_row,rand_col]-1
      
      #### select a new site where the individual moves
      pp <- runif(1,0,1)
      
      #### go down
      if(pp<0.25){
        next_row <- rand_row - 1
        next_col <- rand_col
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_A[next_row,next_col] <- spatial_grid_A[next_row,next_col]+1
      }
      #### go up
      if((pp>=0.25)&&(pp<0.5)){
        next_row <- rand_row + 1
        next_col <- rand_col
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_A[next_row,next_col] <- spatial_grid_A[next_row,next_col]+1
      }
      #### go left
      if((pp>=0.5)&&(pp<0.75)){
        next_row <- rand_row
        next_col <- rand_col - 1
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_A[next_row,next_col] <- spatial_grid_A[next_row,next_col]+1
      }
      #### go right
      if(pp>=0.75){
        next_row <- rand_row
        next_col <- rand_col+1
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_A[next_row,next_col] <- spatial_grid_A[next_row,next_col]+1
      }
    }
    
    #############################################################################################
    #### save results in a list
    list_simu[[kk]][[jj]] <- list(tt,spatial_grid_R,spatial_grid_A,spatial_grid_C,
                                  c(count_r_dd,count_r_bb,count_r_pp,count_pp,count_dd,count_hh))
    
    #### plot results in terminal if needed
    if(PLOT_IN_TERMINAL){
      print(paste("replica ",kk,sep=""))
      print(paste("event "  ,jj,sep=""))
      print(paste("time "   ,tt,sep=""))
      print("Population Matrix Resources")
      print(spatial_grid_R)
      print(paste("Total number of resources ", sum(spatial_grid_R),sep=""))
      print("Population Matrix Consumers")
      print(spatial_grid_A)
      print(paste("Total number of consumers ", sum(spatial_grid_A),sep=""))
      print("Population Matrix Compounds")
      print(spatial_grid_C)
      print(paste("Total number of compounds ", sum(spatial_grid_C),sep=""))
    }
  }###   close loop in events
  print(kk)
}### close loop in replicates

#############################################################################################
### sampling of the results of stochastic simulation 
#############################################################################################

### time step
delta_t <- 0.1

### load list with simulation results
list_results <- list_simu

FF_R           <- data.frame(matrix(nrow = NREP,ncol = 2,0))
colnames(FF_R) <- c("Prey_Density","N_Feeding_Events")

kk <- 1
for(kk in 1:NREP){
  FF_R[kk,1]  <- as.numeric(list_results[[kk]][[1]][[2]])
  N_FEEDING_EVENTS <- as.numeric(list_results[[kk]][[NEV]][[5]][4])-as.numeric(list_results[[kk]][[NEV-TRANS]][[5]][4])
  DELTA_T_FEEDING  <- as.numeric(list_results[[kk]][[NEV]][[1]])-as.numeric(list_results[[kk]][[NEV-TRANS]][[1]])
  FF_R[kk,2]       <- N_FEEDING_EVENTS/(DELTA_T_FEEDING*(N_A+N_C))
}

#############################################################################
### plotting

PP  <- 12

b_b <- PP
l_l <- PP
t_t <- PP
r_r <- PP

X11(width=300,height=300)
par(oma=c(6,3,8,3))
par(oma=c(2,2,2,2))

layout.matrix <- matrix(seq(1,3,1),
                        nrow = 1, ncol = 1,
                        byrow = F)

layout(mat = layout.matrix,
       heights = rep(1,1),   # Heights of rows
       widths  = rep(1,1))   # Widths of columns

layout.show(1)

par(mar = c(b_b,l_l,t_t,r_r))
#############################################################################

plot(FF_R[,1],FF_R[,2],
     type='p',lwd=2,
     xlab="RESOURCES",
     ylab="PER CAPITA RATE",
     main="",
     cex.lab=3,
     cex.axis=2,
     ylim=c(0,max(FF_R[,2])),
     xlim=c(0,max(FF_R[,1])+1)
     #ylim=c(0,max(ff_R[,(NREP+2):(2*NREP+1)]))
)


lines(FF_R[,1],
      FF_R[,1]*(alpha/V)/(1+(alpha/V)*h*FF_R[,1]),
      type='l',col="red",lwd=3)
legend("bottomright",legend=c("STOCHASTIC DYNAMICS", "HOLLING TYPE II"),col=c("black","red"),lwd=3,cex=2)

mtext(paste("PARAMETERS\n",
            " Lambda_R = ",lambda_R," Delta_R = ",delta_R," Beta_R = ", beta_R,"\n",
            " Lambda_A = ",lambda_A," Delta_A = ",delta_A,"\n",
            " Alpha    = ",alpha, "   Nu      = ",signif(nu,digits=3),"\n","V = ",V,"\n",
            " NEVENTS  = ",NEV,"      NREP    = ",NREP,"\n",
            "TRANSIENT = ",TRANS,
            sep=""),
      side=3,outer=T,at=0.5,cex=1,line=0)


