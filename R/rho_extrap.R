#
#  Load the bootstrap data
#
#
#

source("rho_data.R") 


mps.sq <- mps^2


ainv.MeV <- 1000 * 0.1973 / a.fm
ainv.GeV <- 0.1973 / a.fm

cat("Data in lattice units\n") 
for (ii in seq(1 :length(rhomass.full) ))
cat(mu[ii], "  " , mps.sq[ii] , " " , rhomass.full[ii] , "  " , rhomasserr[ii] , "\n") 

cat ("a-1 MeV = ", ainv.MeV , "\n") 

##mps.sq.GeVsq <- ainv.GeV^2 *  mps.sq

mps.sq.MeVsq <- ainv.MeV^2 *  mps.sq
rhomass.boot.MeV <-   ainv.MeV * rhomass.boot
rhomasserr.MeV <-   ainv.MeV *rhomasserr
rhomass.full.MeV <-  ainv.MeV *rhomass.full

cat("Data in GeV^2 or MeV\n") 
for (ii in seq(1 :length(rhomass.full) ))
cat(mu[ii], "  " , mps.sq.MeVsq[ii]/1000^2 , " " , rhomass.full.MeV[ii] , "  " , rhomasserr.MeV[ii] , "\n") 

cat("Start to do some fits") 
source("fit_util.R") 



# starting guess in MeV
##par <- c(0.770,0.16, 0.8 )
par <- c(770,0.001 , 900 )
##par <- c(770,0.001  )


ans <- optim(par=par, fn=ChiSqr.1mass, method="Nelder-Mead",x=mps.sq.MeVsq ,y=rhomass.full.MeV,err=rhomasserr.MeV)

##ans <- optim(par=par, fn=ChiSqr.1mass, method="BFGS",x=mps.sq.MeVsq ,y=rhomass.full.MeV,err=rhomasserr.MeV)


print(ans$par[1])
print(ans$par[2])
print(ans$par[3])

full.par <- ans$par 

dof <- length(mps.sq.MeVsq) -2
cat("chi^2 / dof " , ans$value , " / "  , dof , "\n") 

cat("Length =  " , length(mps.sq.MeVsq) , "\n") 

mpssq.phys <- 135^2
m_rho_phys_full <-   adalaide_rho(ans$par,mpssq.phys)
cat("m_rho (physical) = " , m_rho_phys_full , "\n")


##stop("early end")

## plot the fit 

cat("Fit + Data in GeV^2 or MeV\n") 
for (ii in seq(1 :length(rhomass.boot) ))
#cat(mu[ii], "  " , mps.sq.MeVsq[ii]/1000^2 , " " , rhomass.full.MeV[ii] , "  " , rhomasserr.MeV[ii] , " from fit " , adalaide_rho(ans$par,mps.sq.MeVsq[ii])  ,  "\n") 

#cat(mu[ii], "  " , mps.sq.MeVsq[ii]/1000^2 , " " , rhomass.full.MeV[ii] , "  " , rhomasserr.MeV[ii] , " from fit " , linear_rho(ans$par,mps.sq.MeVsq[ii])  ,  "\n") 


xx_start <- 0.01
#xx_start <- 0.2

xx_end   <- 0.3
#tot <- 1
tot <- 500
delta <- ( xx_end - xx_start ) / tot

xx <- xx_start

for ( i in 1:tot ) {

mPS.sq <- xx * 1000^2
mrho <-   adalaide_rho(ans$par,mPS.sq)
##mrho <- linear_rho(ans$par,mPS.sq)

cat(xx, "  " , mrho , "\n"  )

xx <- xx + delta
}



#
# bootstrap analysis
#

cat("Starting bootstrap analysis \n")

m_rho_phys_boot <- c( 0 , 0 , 0 )
par1.boot <- c( 0 , 0 , 0 )
par2.boot <- c( 0 , 0 , 0 )
par3.boot <- c( 0 , 0 , 0 )

for (ii in seq(1 :length(myDataFrame4$mass) ))
{
 yy <- c( myDataFrame4$mass[ii] , myDataFrame64$mass[ii], myDataFrame8$mass[ii] , myDataFrame10$mass[ii] )
## cat( yy[1] , " " , yy[2] , " " , yy[3] , " " , yy[4] , "\n" )
 yy_MeV <- ainv.MeV * yy 

##ans <- optim(par, ChiSqr.1mass, method="BFGS",x=mps.sq.MeVsq ,y=yy_MeV,err=rhomasserr.MeV)
##ans <- optim(par, ChiSqr.1mass, method="Nelder-Mead",x=mps.sq.MeVsq ,y=yy_MeV,err=rhomasserr.MeV)
ans <- optim(par, ChiSqr.1mass, x=mps.sq.MeVsq ,y=yy_MeV,err=rhomasserr.MeV)

m_rho_phys_boot[ii] <- adalaide_rho(ans$par,mpssq.phys)
par1.boot[ii] <- ans$par[1] 
par2.boot[ii] <- ans$par[2] 
par3.boot[ii] <- ans$par[3] 

}

cat("chi^2 / dof " , ans$value , " / "  , dof , "\n") 
m_rho_phys_boot_err <- sd(m_rho_phys_boot ) 

part1.err <- sd(par1.boot)
part2.err <- sd(par2.boot)
part3.err <- sd(par3.boot)

cat(" m_rho = ", m_rho_phys_full, " +/- " , m_rho_phys_boot_err , "\n")

cat("m_rho^0 = " , full.par[1], " +/- " , part1.err )

cat("m_rho^0 + x m_ps. x = " , full.par[2], " +/- " , part2.err )

cat("lambda_piw =  = " , full.par[3], " +/- " , part3.err )





