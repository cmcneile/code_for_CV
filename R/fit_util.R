#  source(file="fit_mvec.R")
#
# Try a linear fit to the mass of the vector meson
#
# http://wiki.r-project.org/rwiki/doku.php?id=rdoc:stats:optim
#
#

source("adalaide_util.R")

#
# linear fit model
#
linear_rho <- function(par,x)
{
  ss <- length(par)
  if ( ss != 2 )
  { stop("number of parameters out of range")   }
  return (par[1] + x * par[2] ) 
}




#
#  fit model from down under 
#
adalaide_rho <- function(par,x)
{
  ss <- length(par)
  if ( ss != 3 )
  { stop("number of parameters out of range")   }
 mPS.a <- sqrt( x ) 
 Lamb <- par[3]  


integrand_Piw <- function(k) 
{

# Lamb and mPS are global variables
aaa <- sigmaPiw(Lamb,k,mPS.a)

return (aaa)
}



integrand_PiPi <- function(k) 
{
##cat("DEBUG mPS= ",mPS.a, " Lamb=",Lamb, " k= ", k,  " \n" ) 

Lamb.pipi <- 2.0 * sqrt( Lamb^2 - 135^2 ) 

# Lamb and mPS are global variables
aaa <- sigmaPiPi(Lamb.pipi,k,mPS.a)
#cat("DEBUG mPS= ",mPS.a, " Lamb=",Lamb, " k= ", k,   " aaa=",aaa, " \n" ) 
return (aaa)
}


# pi -pi self energy
pole <- 770^2/4 - mPS.a^2
epsilon <- 4

if( pole < 0 ) 
{
 sig2 <- integrate(integrand_PiPi , lower = 0, upper = Inf)

 sig2.ans <- sig2$value 
}
else
{
 b.l <- sqrt(pole) - epsilon 
 b.u <- sqrt(pole) + epsilon 

 sig2.l <- integrate(integrand_PiPi , lower = 0, upper = b.l)
 sig2.u <- integrate(integrand_PiPi , lower = b.u, upper = Inf)

 sig2.ans <- sig2.l$value + sig2.u$value 

}

  sig.pipi <-  ( - 6.028^2  * sig2.ans ) / (6 * pi^2) 

# w -- pi integral
 sig1 <- integrate(integrand_Piw , lower = 0, upper = Inf)
 sig.piw  <- -1.0 * sig1$value * ( (16/1000)^2 * 770 ) / ( 12 * pi^2 )

  ss <- par[1] + x * par[2]
  ee <- (sig.pipi + sig.piw)/(2.0*( par[1] + x * par[2] ))

##  cat("DEBUG ss = " , ss , "  ee= ", ee , "\n")

 ans <- par[1] + x * par[2] + (sig.pipi + sig.piw)/(2.0*( par[1] + x * par[2] ))

 return ( ans ) 
}



#
# next order
#
linear_quad_rho <- function(par,x)
{
  ss <- length(par)
  if ( ss != 3 )
  { stop("number of parameters out of range")   }
  return (par[1] + x * par[2] + x**(3/2) * par[3]  ) 
}


#
# function to be minimized
#
ChiSqr.1mass <- function(par, x, y, err) {
  # index of mass
  # l <- length(par)/no.masses
  tr <- 4 ; 
##  ii <- c(1:tr)

  mass  <- par[1] 
  slope <- par[2]

##  Sumall = Sumall + sum(((y[ii] - mass + slope * x[ii] )/err[ii])^2)
##  Sumall = Sumall + sum(((y - mass + slope * x)/err)^2)
##  Sumall <- sum(((y - mass + slope * x)/err)^2)
##    Sumall <- sum(((y - linear_rho(par,x))/err)^2)


##    Sumall <- sum(((y[ii] - adalaide_rho(par,x[ii]))/err[ii])^2)

 Sumall <- 0.
  for(ii in 1:length(x)) {
    Sumall <- Sumall + ((y[ii] - adalaide_rho(par,x[ii]))/err[ii])^2
#     Sumall <- Sumall + ((y[ii] - linear_rho(par,x[ii]))/err[ii])^2
}

##   Sumall <- sum(((y - linear_quad_rho(par,x))/err)^2)
##  Sumall <- sum(((y - f(par,x))/err)^2)

## cat("DEBUG par = ", par[1] , " " , par[2] , " " , par[3] , "\n" ) 

  return(Sumall)
}

