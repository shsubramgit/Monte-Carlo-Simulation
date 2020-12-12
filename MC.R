#*****************************************************************************************
# Contents:  
#
#    1.  Monte-Carlo Simulation - Generate and plot a Random Walk path
#    2.  Time-Series Analysis:  Generate AR(p) process and plot their ACF and PACF
#*****************************************************************************************

#-----------------------------------------------------------------------------------------
#                                 1. Monte Carlo Simulation
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#                           Generate sample Brownian path ( Random Walk )
#-----------------------------------------------------------------------------------------

# set Brownian path time interval
t_initial <- 0.0  ;
t_final   <- 1.0  ;

nx        <- 1000 ;  # number of points in the time interval

# time step
bm_dt <- (t_final - t_initial)/(nx-1.0) ;
bm_t  <- seq(t_initial, t_final, by = bm_dt ) ;

# generate standard normal random numbers
bm_rnorm <- rnorm(nx, mean = 0.0, sd = 1.0) ;
hist(bm_rnorm, breaks = 100) ;

# set mean and volatility of Brownian path
bm_mean       <- 0.0 ;
bm_volatility <- 1.0*sqrt(bm_dt) ;

# Generate Brownian path
x    <- 0.0*c(1:nx)   ;  # Initialize
x[1] <- 0.0           ;

for (i in {1:(nx-1)} )
{
  x[i+1] <- x[i] + bm_mean*bm_dt + bm_volatility*bm_rnorm[i]    ; 
} ;

# plot the Brownian path
plot(
  bm_t,                                    # time
  x,                                       # Brownian path
  type = "l",                              # plot as line
  col  = "blue",                           # Blue color line
  xlab = "Time, t",                        # x-label
  ylab = "S(t)",                           # y-label
  main = "Brownian Sample Path, S(t)"      # Main title
);

#------------------------------------------------------------------------------------------------------
#                                2. Time Series Analysis  
#
#                  Auto-Correlation and Partial Auto-Correlaiion of AR(p) process
#------------------------------------------------------------------------------------------------------
# Reference:
#
#  Example 2:                                ARIMA
#                                 
#  Chapter 4:  Stochastic Models
#
#  Ref:
#  1. "Introductory time series with R",
#      Paul S.P. Cowpertwait & Andrew V. Metcalfe ; Springer, 2009
#
#  1. Play with AR coeff size and sign ; and observe ACF & PACF
#  2. For AR(P) to be stationary, ALL abs(root) of AR polynomial > 1.0
#     If any abs(root) < 1, the process explodes.
#     if one root = 1 , the process is non-stationary
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# 1. ACF of Brownian process ( Random Walk )
#------------------------------------------------------------------------------------------------------
nb <- 1000           ; # number of random numbers for Brownian process
w <- rnorm(nb)       ; # generate std gaussian random numbers

x <- array(0, dim = c(nb,1)) ;
x[1] <- w[1]         ; # initialize Brownian process
for (i in 2:nb) { x[i] <- x[i - 1] + w[i] } ; # generate Brownian process
plot(w, type = "l" ) ; # plot Gaussian random numbers
plot(x, type = "l" ) ; # plot Brownian path

acf(x)               ; # ACF of Brownian process ; Note the slow decay of acf coeff. == pg 72, eqn 4.10
pacf(x)              ; # PACF is a delta function === spike 1 at lag 0, and 0 elsewhere

acf(diff(x))         ; # ACF of 1-st diff. is a delta function === spike 1 at lag 0, and 0 elsewhere
pacf(diff(x))        ; # zero at all lags. why ?


#------------------------------------------------------------------------------------------------------
# 2. ACF of AR(1)
#------------------------------------------------------------------------------------------------------

nb   <- 1000     ; # number of random numbers for Brownian path
ar_1 <- 0.7      ; # coeff of AR(1) process

polyroot(c(1,-ar_1)) ;  # Find root of AR(1) process

w <- rnorm(nb) ; # generate std gaussian random numbers
x <- array(0, dim = c(nb,1)) ;
x[1] <- w[1]   ; # initialize Brownian path
for (i in 2:nb) { x[i] <- ar_1*x[i - 1] + w[i] } ; # generate Brownian path
plot(w, type = "l" ) ; # plot Gaussian random numbers
plot(x, type = "l" ) ; # plot AR(1) process
acf(x)               ; # ACF of AR(1); decays faster than brownian
pacf(x)              ; # Lag 1 spike is significant, others are not


#------------------------------------------------------------------------------------------------------
# 3. ACF of AR(2)
#------------------------------------------------------------------------------------------------------

nb   <- 1000     ; # number of random numbers for Brownian path
ar_1 <- 0.7      ; # coeff of AR(2) process
ar_2 <- -0.4      ; # coeff of AR(2) process

polyroot(c(1,-ar_1,-ar_2)) ;   # Find root of AR(2) process

w <- rnorm(nb) ; # generate std gaussian random numbers
x <- array(0, dim = c(nb,1)) ;
x[1] <- w[1]   ; # initialize Brownian path
x[2] <- w[2]   ; # initialize Brownian path
for (i in 3:nb) { x[i] <- ar_1*x[i - 1] + ar_2*x[i - 2] + w[i] } ; # generate Brownian path
plot(w, type = "l" ) ; # plot Gaussian random numbers
plot(x, type = "l" ) ; # plot AR(2) process
acf(x)               ; # ACF decays slow
pacf(x)              ; # Lag 1,2 spikes is significant, others are not


#------------------------------------------------------------------------------------------------------
# 4. ACF of AR(3)
#------------------------------------------------------------------------------------------------------

nb   <- 1000     ; # number of random numbers for Brownian path
ar_1 <- 0.7      ; # coeff of AR(3) process
ar_2 <- -0.4      ; # coeff of AR(3) process
ar_3 <- 0.2      ; # coeff of AR(3) process

polyroot(c(1, -ar_1, -ar_2, -ar_3)) ;   # Find root of AR(3) process

w <- rnorm(nb) ; # generate std gaussian random numbers
x <- array(0, dim = c(nb,1)) ;
x[1] <- w[1]   ; # initialize Brownian path
x[2] <- w[2]   ; # initialize Brownian path
x[3] <- w[3]   ; # initialize Brownian path
for (i in 4:nb) { x[i] <- ar_1*x[i - 1] + ar_2*x[i - 2] + ar_3*x[i - 3] + w[i] } ; # generate Brownian path
plot(w, type = "l" ) ; # plot Gaussian random numbers
plot(x, type = "l" ) ; # plot AR(3) process
acf(x)               ; # ACF decays slow
pacf(x)              ; # Lag 1,2,3 spikes is significant, others are not




