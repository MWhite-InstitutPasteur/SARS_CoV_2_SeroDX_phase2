
ant_names <- "NP_NA_IgG"

##########################################################################
##########################################################################
##                                                                      ## 
##  #####   ####  #####  ##  ## ##     ####  ###### ####  ####  #   ##  ##
##  ##  ## ##  ## ##  ## ##  ## ##    ##  ##   ##    ##  ##  ## ##  ##  ##
##  #####  ##  ## #####  ##  ## ##    ######   ##    ##  ##  ## ### ##  ##
##  ##     ##  ## ##     ##  ## ##    ##  ##   ##    ##  ##  ## ## ###  ##
##  ##      ####  ##      ####  ##### ##  ##   ##   ####  ####  ##  ##  ##
##                                                                      ## 
##########################################################################
##########################################################################

pop_file <- paste( "C:/U/CoronaVirus/SeroSig_Phase2/CoV_kin/IgG_mod_1/OUTPUT/", ant_names, "_dil_global_med.txt", sep="" )

MCMC_pop <- read.table( pop_file )
colnames(MCMC_pop) <- c( "AB_bg", "beta", "tau", "t_delta",  "t_short",     "t_long",     "t_IgG",     "rho",   
                         "sig_Ab_bg", "sig_beta", "sig_tau", "sig_t_delta", "sig_t_short", "sig_t_long", "sig_t_IgG", "sig_rho",  
                         "sig_obs", "likelihood", "prior")



N_par = ncol(MCMC_pop) - 2


for(k in 9:15)
{
	MCMC_pop[,k] = 1/sqrt(MCMC_pop[,k])
}

for(k in 1:7)
{
	MCMC_pop[,k]   = exp( MCMC_pop[,k] + 0.5*MCMC_pop[,8+k]^2 )
	MCMC_pop[,8+k] = sqrt( exp(MCMC_pop[,8+k]^2) -1 )*MCMC_pop[,k]
}


for(k in 1:nrow(MCMC_pop))
{
	f_m <- function(x)
	{ 
		0.3989423*exp( -0.5*( ( log(x/(1-x))-MCMC_pop[k,8] )/MCMC_pop[k,16] )^2 )/( MCMC_pop[k,16]*(1-x) )
	}

	f_m2 <- function(x)
	{ 
		x*0.3989423*exp( -0.5*( (log(x/(1-x))-MCMC_pop[k,8])/MCMC_pop[k,16] )^2 )/( MCMC_pop[k,16]*(1-x) )
	}
	
	tryCatch(
	{
		moment_1 <- integrate(f_m, lower=0, upper=1)$value
		moment_2 <- integrate(f_m2, lower=0, upper=1)$value
	}, error=function(e){ NULL }
	)

	MCMC_pop[k,8] <- moment_1
	MCMC_pop[k,16] <- sqrt( moment_2 - moment_1^2 )
}


MCMC_pop_burn <- MCMC_pop[floor(0.1*nrow(MCMC_pop)):(nrow(MCMC_pop)-1),]



#############
## Frame 1 ##
#############

##############################################
## First make some trace plots

par(ask=TRUE)

par(mfrow=c(4,5))

for(i in 1:N_par)
{
	plot(x=1:nrow(MCMC_pop), y=MCMC_pop[,i], 
           pch=19, cex=0.01, col="grey", 
	     ##  ylim=c( min(MCMC_pop[,i], na.rm=TRUE), max(MCMC_pop[,i], na.rm=TRUE) ),
	     ylim=quantile(MCMC_pop[,i], prob=c(0.01,0.99), na.rm=TRUE),
	     main=paste( colnames(MCMC_pop)[i] ), ylab="", xlab="")
}


plot(x=1:nrow(MCMC_pop), y=MCMC_pop[,N_par+1], 
pch=19, cex=0.01, col="grey", 
ylim=quantile(MCMC_pop[,N_par+1], prob=c(0.0001,1), na.rm=TRUE),
main=paste( colnames(MCMC_pop)[N_par+1] ), ylab="", xlab="")

##points(x=1:nrow(MCMC_pop_2), y=MCMC_pop_2[,N_par+1], 
##pch=19, cex=0.01, col="yellowgreen") 


#############
## Frame 2 ##
#############

par(mfrow=c(4,5))


for(i in 1:N_par)
{
	DEN_1 = density( log(MCMC_pop_burn[,i]), na.rm=TRUE )

	DEN_1$x <- exp(DEN_1$x)

	QUANT_1 = quantile( MCMC_pop_burn[,i], prob=c(0.025, 0.5, 0.975), na.rm=TRUE )



	###############################
	## MCMC chain 1

	plot(x=DEN_1$x, y=DEN_1$y, type='l', log="x",
	xlim=quantile(MCMC_pop_burn[,i], prob=c(0.0,0.995), na.rm=TRUE),
	main=paste( colnames(MCMC_pop_burn)[i] ), ylab=NULL, xlab=NULL)

	low_index  = which(DEN_1$x<QUANT_1[1])
	mid_index  = intersect( which(DEN_1$x>=QUANT_1[1]), which(DEN_1$x<=QUANT_1[3]) )
	high_index = which(DEN_1$x>QUANT_1[3])

	polygon( x=c( DEN_1$x[low_index], rev(DEN_1$x[low_index]) ),
		   y=c( rep(0,length(low_index)), rev(DEN_1$y[low_index]) ), 
               col="pink")

	polygon( x=c( DEN_1$x[mid_index], rev(DEN_1$x[mid_index]) ),
		   y=c( rep(0,length(mid_index)), rev(DEN_1$y[mid_index]) ), 
               col="grey")

	polygon( x=c( DEN_1$x[high_index], rev(DEN_1$x[high_index]) ),
		   y=c( rep(0,length(high_index)), rev(DEN_1$y[high_index]) ), 
               col="pink")

	points(x=rep(QUANT_1[2],2), y=c(0,max(DEN_1$y)), type='l', lty="dashed", lwd=2)
}




for(i in 1:(N_par/2))
{
	print( colnames(MCMC_pop_burn)[i] )

	print(paste( quantile( MCMC_pop_burn[,i], prob=c(0.5, 0.025, 0.975), na.rm=TRUE ) ))
}



for(i in (1 + N_par/2):(N_par))
{
	print( colnames(MCMC_pop_burn)[i] )

	print(paste( quantile( MCMC_pop_burn[,i], prob=c(0.5, 0.025, 0.975), na.rm=TRUE ) ))
}

quantile(MCMC_pop_burn$sig_obs, prob=c(0.5, 0.025, 0.975) )





#######################################################################
#######################################################################
##                                                                   ##
##  #### #   ## ####   #### ##   ## #### ####   ##  ##  ####  ##     ##
##   ##  ##  ## ## ##   ##  ##   ##  ##  ## ##  ##  ## ##  ## ##     ##
##   ##  ### ## ##  ##  ##   ## ##   ##  ##  ## ##  ## ###### ##     ##
##   ##  ## ### ## ##   ##    ###    ##  ## ##  ##  ## ##  ## ##     ## 
##  #### ##  ## ####   ####    #    #### ####    ####  ##  ## #####  ##
##                                                                   ## 
#######################################################################
#######################################################################


AB_data = read.table( paste("C:/U/CoronaVirus/SeroSig_Phase2/DATA/ABkin_input/", ant_names, "_dil.txt", sep="") )

N_tt <- 15

N_part <- nrow(AB_data)


AB_max = 0.2
tt_max = max( AB_data[,4:(3+N_tt)] )


AB_min = 1.95e-5


tt_plot <- seq(from=-10, to=50, by=1)
log2 <- log(2)


ind_file <-  paste( "C:/U/CoronaVirus/SeroSig_Phase2/CoV_kin/IgG_mod_1/OUTPUT/", ant_names, "_dil_local_med.txt", sep="" )




MCMC_ind <- read.table( ind_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]




par(ask=TRUE)

lay.mat <- rbind( c( 1, 1, 2, 2, 3, 4, 5, 6, 7),
                  c( 1, 1, 2, 2, 8, 9,10,11,12),
	            c(13,13,14,14,15,16,17,18,19),
                  c(13,13,14,14,20,21,22,23,24) )
layout(lay.mat)
layout.show(24)


for(n in 1:194 )
{
	###################################
	## Data for participant n

	N_sam = length(which(AB_data[n,4:(3+N_tt)] > -0.5))

	tt = as.numeric(AB_data[n,4:(3+N_sam)])
	AB = as.numeric(AB_data[n,(4+N_tt):(3+N_tt+N_sam)])

	

	###################################
	## Model prediction for participant i
	## Posterior projections

	AB_mod <- matrix(NA, nrow=nrow(MCMC_ind), ncol=length(tt_plot))
	
	for(i in 1:nrow(AB_mod))
	{
		Ab_0    = MCMC_ind[i,(n-1)*9+1]
		beta    = MCMC_ind[i,(n-1)*9+2]
		tau     = MCMC_ind[i,(n-1)*9+3]		
		t_delta = MCMC_ind[i,(n-1)*9+4]
		t_short = MCMC_ind[i,(n-1)*9+5]
		t_long  = MCMC_ind[i,(n-1)*9+6]
		t_IgG   = MCMC_ind[i,(n-1)*9+7]
		rho     = MCMC_ind[i,(n-1)*9+8]

		r_delta = log(2)/t_delta
		r_cs = log(2)/t_short
		r_cl = log(2)/t_long
		r_a  = log(2)/t_IgG

		part_1 = rho / ((r_cs - r_a)*(r_cs - r_delta))
		part_2 = (1 - rho) / ((r_cl - r_a)*(r_cl - r_delta))
		part_3 = (r_a + r_cs*(rho-1) - r_cl*rho) / ((r_cs - r_a)*(r_cl - r_a)*(r_a - r_delta))
		part_4 = - (r_delta + r_cs*(rho - 1) - r_cl*rho) / ((r_cs - r_delta)*(r_cl - r_delta)*(r_a - r_delta))

		##AB_tt = Ab_0*exp(-r_cl*tt_plot)
		AB_tt =  rep( Ab_0, length(tt_plot) )

		tt_temp = (tt_plot - tau)[which(tt_plot >= tau )]
		AB_tt[which(tt_plot >= tau)] = AB_tt[which(tt_plot >= tau)] +
                                                  beta*( part_1*exp(-r_cs*tt_temp) + part_2*exp(-r_cl*tt_temp) +  
                                                         part_3*exp(-r_a*tt_temp) + part_4*exp(-r_delta*tt_temp) ) 

		AB_mod[i,] = AB_tt
	}

	AB_quant <- matrix(NA, nrow=3, ncol=length(tt_plot))
	for(j in 1:length(tt_plot))
	{
		AB_quant[,j] <- quantile( AB_mod[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}


	########################
	## PANEL 1

	if( AB_data[n,1] == 1 )
	{
		plot_title <- paste("Bichat case #", n, sep="")
	}
	if( AB_data[n,1] == 2 )
	{
		plot_title <- paste("Strasbourg case #", n - 4, sep="")
	}
	if( AB_data[n,1] == 3 )
	{
		plot_title <- paste("Cochin case #", n - 4 - 157, sep="")
	}

	plot( x=tt, y=AB, pch=19, cex=1, col="red", 
            xlim=c(-0,50), ylim=c(0,0.025),
		xlab="time (days)", ylab="antibody titre",
            main=plot_title )

	points(x=tt_plot, y=AB_quant[2,], type='l')

	polygon(x=c(tt_plot, rev(tt_plot)), 
		  y=c( AB_quant[1,], rev(AB_quant[3,]) ),
	        col=rgb(190/256,190/256,190/256,0.4), border=NA)

	points(x=c(-100,1000), y=c(0.02,0.02), type='l', lty="dashed")


	########################
	## PANEL 2

	plot( x=tt, y=AB, pch=19, cex=1, col="red", 
            xlim=c(-0,50), ylim=c(1e-5,0.025), log="y",
		xlab="time (days)", ylab="antibody titre",
            main=plot_title )


	points(x=tt_plot, y=AB_quant[2,], type='l')


	polygon(x=c(tt_plot, rev(tt_plot)), 
		y=c( AB_quant[1,], rev(AB_quant[3,]) ),
		col=rgb(190/256,190/256,190/256,0.4), border=NA)

	points(x=c(-100,1000), y=c(1.95e-5,1.95e-5), type='l', lty="dashed")
	points(x=c(-100,1000), y=c(0.02,0.02), type='l', lty="dashed")



	#################################
	## PANEL 3: Ab_0 chain

	plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*9+1], 
           pch=19, cex=0.01, col="grey", 
           xlab="MCMC iteration", ylab="", main="Ab_0" )


	#################################
	## PANEL 4: beta chain

	plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*9+2], 
           pch=19, cex=0.01, col="grey", 
           xlab="MCMC iteration", ylab="", main="beta" )


	#################################
	## PANEL 5: tau chain

	plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*9+3], 
           pch=19, cex=0.01, col="grey", 
           xlab="MCMC iteration", ylab="", main="tau" )


	#################################
	## PANEL 6: t_delta chain

	plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*9+4], 
           pch=19, cex=0.01, col="grey", 
           xlab="MCMC iteration", ylab="", main="t_delta" )


	#################################
	## PANEL 7: t_short chain

	plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*9+5], 
           pch=19, cex=0.01, col="grey", 
           xlab="MCMC iteration", ylab="", main="t_short" )


	#################################
	## PANEL 8: t_long chain

	plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*9+6], 
           pch=19, cex=0.01, col="grey", log="y",
           xlab="MCMC iteration", ylab="", main="t_long" )


	#################################
	## PANEL 9: t_IgG chain

	plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*9+7], 
           pch=19, cex=0.01, col="grey", 
           xlab="MCMC iteration", ylab="", main="t_IgG" )


	#################################
	## PANEL 10: rho chain

	plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*9+8], 
           pch=19, cex=0.01, col="grey", ylim=c(0,1),
           xlab="MCMC iteration", ylab="", main="rho" )


	#################################
	## PANEL 11: lhood chain

	plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*9+9], 
           pch=19, cex=0.01, col="grey", 
           xlab="MCMC iteration", ylab="", main="likelihood" )

	#################################
	## PANEL 12: empty

	plot.new()
}


##################################################
##################################################
##                                              ##
##  #   ## #####  ####     ####   ####  #   ##  ##
##  ##  ## ##    ##       ##  ## ##  ## ##  ##  ##
##  ### ## ####  ## ###   ##     ##  ## ### ##  ##
##  ## ### ##    ##  ##   ##  ## ##  ## ## ###  ##  
##  ##  ## #####  ####     ####   ####  ##  ##  ##
##                                              ##
################################################## 
##################################################

N_negcon <- 270



par(mfrow=c(3,6))


for(n in 195:464)
{
	###################################
	## Data for participant n

	Ab_bg   = MCMC_ind[,n*9+1]


	########################
	## PANEL 1

	plot( x=0, y=AB_data[n,19], pch=19, cex=1, col="red", 
            xlim=c(-50,50), ylim=c(1e-6,0.3), log="y",
		xlab="time (days)", ylab="antibody titre",
            main=paste("negative control #", n , sep="") )


	points(x=c(-50,450), y=rep(median(Ab_bg),2), type='l')

	points(x=c(-50,450), y=rep(1.95e-5,2), type='l', lty="dashed")


	polygon(x=c(-50, 450, 450, -50), 
		  y=quantile( Ab_bg, prob=c(0.025, 0.025, 0.975, 0.975)),
	        col=rgb(190/256,190/256,190/256,0.4), border=NA)


	#################################
	## PANEL 2: Ab_bg chain

	plot(x=1:nrow(MCMC_ind), y=MCMC_ind[,(n-1)*9+1], 
           pch=19, cex=0.01, col="grey", 
           xlab="MCMC iteration", ylab="", main="Ab_bg" )

}





