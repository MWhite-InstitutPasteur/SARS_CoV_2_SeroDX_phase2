library(binom)
library(MASS)
library(randomForest)
library(pROC)


set.seed(1234)



###################################
###################################
##                               ##
##  ####    ####  ######  ####   ##
##  ## ##  ##  ##   ##   ##  ##  ##
##  ##  ## ######   ##   ######  ##
##  ## ##  ##  ##   ##   ##  ##  ##
##  ####   ##  ##   ##   ##  ##  ##
##                               ##
###################################
###################################


AB_data <- read.csv("C:\\U\\CoronaVirus\\SeroSig_Phase2\\DATA\\SARS-CoV-2_serology.csv")

AB_data <- AB_data[which(is.na(AB_data$S1RBD_NA_IgG_dil)==FALSE),]

N_positive <- length(which(AB_data$status == "positive"))
N_negative <- length(which(AB_data$status == "negative"))


N_pos_kin <- 194

SP_cut <- rep(NA, 7)
names(SP_cut) <- colnames(AB_data)[10:16]

for(k in 1:7)
{
	SP_cut[k] <- quantile( AB_data[which(AB_data$status=="negative"),9+k], prob=0.99 ) 
}


t_bin_widths <- c(0, 10, 15, 20, 25, 30, 40)
t_bin_mids <- 0.5*( t_bin_widths[-1] + t_bin_widths[-length(t_bin_widths)] )
N_bins <- length(t_bin_mids)

IgG_bins <- array(NA, dim=c(7, N_bins, 3) )

for( k in 1:7 )
{
	ab_temp <- AB_data[which(AB_data$status == "positive"),9+k]
	tt_temp <- AB_data$days_post[which(AB_data$status == "positive")]

	for(g in 1:N_bins)
	{
		index <- which( tt_temp > t_bin_widths[g] & tt_temp <= t_bin_widths[g+1] )
		
		IgG_bins[k,g,2] <- exp(mean(log(ab_temp[index])))
		IgG_bins[k,g,c(1,3)] <- quantile( ab_temp[index], prob=c(0.025,0.975) )
	}
}




IgG_NegCon_bins <- matrix(NA, nrow=7, ncol=3 )

for( k in 1:7 )
{
	ab_temp <- AB_data[which(AB_data$status == "negative"),9+k]

	IgG_NegCon_bins[k,2]      <- exp(mean(log(ab_temp)))
	IgG_NegCon_bins[k,c(1,3)] <- quantile( ab_temp, prob=c(0.025,0.975) )
}




IgM_bins <- array(NA, dim=c(7, N_bins, 3) )

for( k in 1:7 )
{
	ab_temp <- AB_data[which(AB_data$status == "positive"),22+k]
	tt_temp <- AB_data$days_post[which(AB_data$status == "positive")]

	for(g in 1:N_bins)
	{
		index <- which( tt_temp > t_bin_widths[g] & tt_temp <= t_bin_widths[g+1] )
		
		IgM_bins[k,g,2] <- exp(mean(log(ab_temp[index])))
		IgM_bins[k,g,c(1,3)] <- quantile( ab_temp[index], prob=c(0.025,0.975) )
	}
}



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

ant_names <- c("Stri_IPP_IgG_dil", "S1RBD_NA_IgG_dil", "RBD_IPP_IgG_dil", "S1_NA_IgG_dil",
              "S2_NA_IgG_dil", "NP_IPP_IgG_dil", "NP_NA_IgG_dil")

N_tt_plot <- 100

tt_plot <- exp(seq(from=log(1), to=log(500), length=N_tt_plot))
log2 <- log(2)



N_posterior <- 100

AB_mod <- array(NA, dim=c(7, N_pos_kin, N_posterior, N_tt_plot ))


for(k in 1:7)
{
	Ant_file <- paste( "C:/U/CoronaVirus/SeroSig_Phase2/CoV_kin/IgG_mod_1/OUTPUT/", ant_names[k], "_local_med.txt", sep="" )

	MCMC_ind <- read.table( Ant_file )

	MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]
	MCMC_ind <- MCMC_ind[seq(from=1, to=nrow(MCMC_ind), length=N_posterior),]


	for(n in 1:N_pos_kin)
	{
		###################################
		## Model prediction for participant n
		## Posterior projections
	
		for(i in 1:N_posterior)
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

			AB_tt =  rep( Ab_0, N_tt_plot )

			tt_temp = (tt_plot - tau)[which(tt_plot >= tau )]
			AB_tt[which(tt_plot >= tau)] = AB_tt[which(tt_plot >= tau)] +
            	                                      beta*( part_1*exp(-r_cs*tt_temp) + part_2*exp(-r_cl*tt_temp) +  
            	                                             part_3*exp(-r_a*tt_temp) + part_4*exp(-r_delta*tt_temp) ) 

			AB_mod[k,n,i,] = AB_tt
		}
	}
}


AB_ind_mod_quant <- array(NA, dim=c(7, N_pos_kin, 3, N_tt_plot ))


for(k in 1:7)
{
	for(n in 1:N_pos_kin)
	{
		for(j in 1:N_tt_plot)
		{
			AB_ind_mod_quant[k,n,,j] <- quantile( AB_mod[k,n,,j], prob=c(0.025, 0.5, 0.975) )
		}
	}
}

AB_pop_mod_quant <- array(NA, dim=c(7, 5, N_tt_plot ))


for(k in 1:7)
{
	for(j in 1:N_tt_plot)
	{
		AB_pop_mod_quant[k,,j] <- quantile( AB_mod[k,,,j], prob=c(0.025, 0.25, 0.5, 0.75, 0.975) )
	}
}




###########################################################
###########################################################
##                                                       ##  
##   ####  ##     ####   ####   ####  #### ##### ##  ##  ## 
##  ##  ## ##    ##  ## ##     ##      ##  ##    ##  ##  ##
##  ##     ##    ######  ####   ####   ##  ####   ####   ## 
##  ##  ## ##    ##  ##     ##     ##  ##  ##      ##    ##
##   ####  ##### ##  ##  ####   ####  #### ##      ##    ##
##                                                       ##
###########################################################
###########################################################


SP_mod <- array(NA, dim=c(7, N_posterior, N_tt_plot ))

for(k in 1:7)
{
	for(i in 1:N_posterior)
	{
		for(j in 1:N_tt_plot)
		{
			SP_mod[k,i,j] <- length(which(AB_mod[k,,i,j] > SP_cut[k]))/N_pos_kin
		}
	}
}			


SP_mod_quant <- array(NA, dim=c(7, 3, N_tt_plot ))

for(k in 1:7)
{
	for(j in 1:N_tt_plot)
	{
		SP_mod_quant[k,,j] <- quantile( SP_mod[k,,j], prob=c(0.025, 0.5, 0.975) )
	}
}
			






##################################
##################################
##                              ##  
##  #####  ##     ####  ######  ##
##  ##  ## ##    ##  ##   ##    ##
##  #####  ##    ##  ##   ##    ##
##  ##     ##    ##  ##   ##    ## 
##  ##     #####  ####    ##    ##
##                              ##
##################################
##################################


ind_plot <- 4
ind_data_id <- "B01-S002"

ant_cols <- c("red", "brown", "orangered", "gold", "yellow",
              "forestgreen", "limegreen", "mediumblue", "royalblue", 
              "burlywood", "darkkhaki", "cornsilk4",
              "red", "brown", "orangered", "gold", "yellow",
              "forestgreen", "limegreen", "mediumblue", "royalblue", 
              "burlywood", "darkkhaki", "cornsilk4")	

plot_names <- c("anti-Stri IgG", "anti-RBDv1 IgG", "anti-RBDv2 IgG",  "anti-S1 IgG", "anti-S2 IgG", 
                "anti-NPv1 IgG", "anti-NPv2 IgG")


tiff(file="Fig3_AB_kinetics.tif", width=40, height=25, units="cm", res=500)

lay.mat <- rbind( c( 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2 ), 
                  c( 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9 ),
                  c(10,10,10,10,10,10,10,11,11,11,11,11,11,11 ), 
                  c(12,12,13,13,14,14,15,15,16,16,17,17,18,18 ),
                  c(19,19,19,19,19,19,19,20,20,20,20,20,20,20 ),
                  c(21,21,22,22,23,23,24,24,25,25,26,26,27,27 ) )
layout(lay.mat, heights=c(1,8,1,8,1,8), widths=c(1.1,1.1,1,1,1,1,1,1,1,1,1,1,1,1))
layout.show(27)


############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))


plot.new()
title( "(A) Individual-level IgG kinetics (Bichat patient)", 
        cex.main=3.0, line=-2)


###############	
##           ##
##  LEGEND   ##
##           ##
###############

plot.new()

legend(x='center', 
       legend = c("IgG antibody dilution  ", "IgM antibody dilution  ", "IgG model prediction"), 
       col = c("black", "black", "black"), 
       pch=c(19,8,NA), lty=c(0,0,1), lwd=2,
       ncol=3, cex=1.5, bty="n" )




par(mar=c(4,6.5,3,1))
par(mgp=c(2.5, 1, 0))

point.size = 1.5
lab.size   = 1.5
axis.size  = 1.0
main.size  = 1.5
line.size  = 2

############################
## PANELS 1:7 Individual kinetics

line_seq_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3)
line_seq_y <- c(1, 3, 5, 10, 30, 50, 100, 300, 500)

for( k in 1:7 )
{
	if( k == 1 )
	{
		par(mar=c(4,5.5,3,1))
	}
	if( k %in% 2:7 )
	{
		par(mar=c(4,3,3,1))
	}

	plot(x = AB_data$days_post[which(AB_data$part_id == ind_data_id)], 
           y = AB_data[which(AB_data$part_id == ind_data_id),9+k],
	     col=ant_cols[k], pch=19, cex=point.size,
	     xlim=c(1,500), ylim=c(1e-5,0.02), log="xy",
	     yaxt='n', xaxt='n', bty='n',
	     ylab="", xlab="days after symptoms",
	     main=plot_names[k],
	     cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

	for(i in 1:length(line_seq_x))
	{
		points(x=c(0.1,1000), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
	}
	for(i in 1:length(line_seq_y))
	{
		points(x=rep(line_seq_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
	}

	points(x=c(0.1,1000), y=rep(1.95e-5,2), type='l', col="black", lty="dashed")
	points(x=c(0.1,1000), y=rep(0.02,2), type='l', col="black", lty="dashed")


	points(x=tt_plot, y=AB_ind_mod_quant[k,ind_plot,2,], 
             type='l', col="black", lwd=line.size)

	polygon(x=c(tt_plot, rev(tt_plot)), 
		  y=c( AB_ind_mod_quant[k,ind_plot,1,], rev(AB_ind_mod_quant[k,ind_plot,3,]) ),
	        col=rgb(190/256,190/256,190/256,0.4), border=NA)


	points( x = AB_data$days_post[which(AB_data$part_id == ind_data_id)], 
              y = AB_data[which(AB_data$part_id == ind_data_id),9+k],
	        col=ant_cols[k], pch=19, cex=point.size )

	points( x = AB_data$days_post[which(AB_data$part_id == ind_data_id)], 
              y = AB_data[which(AB_data$part_id == ind_data_id),22+k],
	        col=ant_cols[k], pch=8, cex=point.size )

	points( x=tt_plot, y=rep(SP_cut[k],length(tt_plot)), 
	        type='l', lty="dashed", col=ant_cols[k], lwd=line.size)

	if( k == 1 )
	{
		mtext(side = 2, line = 4, 
		cex=axis.size, 
		text="antibody dilution")
	}

	axis(1, at=c(1, 3, 5, 10, 30, 50, 100, 300, 500), 
	        label=c("1", "3", "5", "10", "30", "50", "100", "300", "500"), 
	        cex.axis=0.8*axis.size )

	axis(2, at = c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03), 
	        label = c("0.00001", "", "0.0001", "", "0.001", "", "0.01", ""),
	        las=2, cex.axis=axis.size )
}




############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))


plot.new()
title( "(B) Population-level IgG kinetics (all participants)", 
        cex.main=3.0, line=-2)


###############	
##           ##
##  LEGEND   ##
##           ##
###############

plot.new()

legend(x='center', 
       legend = c("IgG GMT (positive)", "IgG GMT (negative controls)", "IgG model prediction"), 
       col = c("black", "black", "black"), 
       pch=c(19,2,NA), lty=c(0,0,1), lwd=2,
       ncol=3, cex=1.5, bty="n" )




############################
## PANELS 8:14 Population kinetics

par(mar=c(4,6.5,3,1))
par(mgp=c(2.5, 1, 0))



line_seq_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3)
line_seq_y <- c(1, 3, 5, 10, 30, 50, 100, 300, 500)

for( k in 1:7 )
{
	if( k == 1 )
	{
		par(mar=c(4,5.5,3,1))
	}
	if( k %in% 2:7 )
	{
		par(mar=c(4,3,3,1))
	}

	plot(x = 100, 
           y =100,
	     col=ant_cols[k], pch=19, cex=point.size,
	     xlim=c(1,500), ylim=c(1e-5,0.02), log="xy",
	     yaxt='n', xaxt='n', bty='n',
	     ylab="", xlab="days after symptoms",
	     main=plot_names[k],
	     cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

	for(i in 1:length(line_seq_x))
	{
		points(x=c(0.1,1000), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
	}
	for(i in 1:length(line_seq_y))
	{
		points(x=rep(line_seq_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
	}

	points(x=c(0.1,1000), y=rep(1.95e-5,2), type='l', col="black", lty="dashed")
	points(x=c(0.1,1000), y=rep(0.02,2), type='l', col="black", lty="dashed")





	polygon(x=c(tt_plot, rev(tt_plot)), 
		  y=c( AB_pop_mod_quant[k,1,], rev(AB_pop_mod_quant[k,5,]) ),
	        col=rgb(190/256,190/256,190/256,0.4), border=NA)


	#points( x = AB_data$days_post[which(AB_data$status == "positive")], 
      #     y = AB_data[which(AB_data$status == "positive"),9+k],
	#        col=ant_cols[k], pch=19, cex=point.size )

	for(g in 1:N_bins)
	{
		points( x=t_bin_mids[g], y=IgG_bins[k,g,2],
                    col=ant_cols[k], pch=19, cex=point.size )

		arrows( x0=t_bin_mids[g], x1=t_bin_mids[g], y0=IgG_bins[k,g,1], y1=IgG_bins[k,g,3], 
              code=3, col=ant_cols[k], angle=90, lwd=1.5, length=0.1)  

		#points( x=t_bin_mids[g], y=IgM_bins[k,g,2],
            #        col=ant_cols[k], pch=8, cex=point.size )
	}


	points( x=tt_plot, y=rep(SP_cut[k],length(tt_plot)), 
	        type='l', lty="dashed", col=ant_cols[k], lwd=line.size)

	points(x=tt_plot, y=AB_pop_mod_quant[k,3,], 
             type='l', col="black", lwd=line.size)


	#points( x = exp(runif(n=N_negative, log(1), log(3))),  
      #        y = AB_data[which(AB_data$status == "negative"),9+k],
	#        col=ant_cols[k], pch=2, cex=0.5*point.size )

	points( x=2, y=IgG_NegCon_bins[k,2],
                col=ant_cols[k], pch=2, cex=point.size )

	arrows( x0=2, x1=2, y0=IgG_NegCon_bins[k,1], y1=IgG_NegCon_bins[k,3], 
              code=3, col=ant_cols[k], angle=90, lwd=1.5, length=0.1) 


	if( k == 1 )
	{
		mtext(side = 2, line = 4, 
		cex=axis.size, 
		text="antibody dilution")
	}


	axis(1, at=c(1, 3, 5, 10, 30, 50, 100, 300, 500), 
	        label=c("1", "3", "5", "10", "30", "50", "100", "300", "500"), 
	        cex.axis=0.8*axis.size )

	axis(2, at = c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03), 
	        label = c("0.00001", "", "0.0001", "", "0.001", "", "0.01", ""),
	        las=2, cex.axis=axis.size )
}



############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))

plot.new()
title( "(C) Model-predicted sensitivity", 
        cex.main=3.0, line=-2)

###############	
##           ##
##  LEGEND   ##
##           ##
###############

plot.new()



############################
## PANELS 15:21 Seropositivity

par(mar=c(4,6.5,3,1))
par(mgp=c(2.5, 1, 0))


line_seq_x <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1)
line_seq_y <- c(1, 3, 5, 10, 30, 50, 100, 300, 500)

for( k in 1:7 )
{
	if( k == 1 )
	{
		par(mar=c(4,5.5,3,1))
	}
	if( k %in% 2:7 )
	{
		par(mar=c(4,3,3,1))
	}

	plot(x = 100, 
           y = 100,
	     col=ant_cols[k], pch=19, cex=point.size,
	     xlim=c(1,500), ylim=c(0,1), log="x",
	     yaxt='n', xaxt='n', bty='n',
	     ylab="", xlab="days after symptoms",
	     main=plot_names[k],
	     cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

	for(i in 1:length(line_seq_x))
	{
		points(x=c(0.1,1000), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
	}
	for(i in 1:length(line_seq_y))
	{
		points(x=rep(line_seq_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
	}


	polygon(x=c(tt_plot, rev(tt_plot)), 
		  y=c( SP_mod_quant[k,1,], rev(SP_mod_quant[k,3,]) ),
	        col=rgb(190/256,190/256,190/256,0.4), border=NA)

	points(x=tt_plot, y=SP_mod_quant[k,2,], 
             type='l', col="black", lwd=line.size)

	if( k == 1 )
	{
		mtext(side = 2, line = 4, 
		cex=axis.size, 
		text="sensitivity")
	}


	axis(1, at=c(1, 3, 5, 10, 30, 50, 100, 300, 500), 
	        label=c("1", "3", "5", "10", "30", "50", "100", "300", "500"), 
	        cex.axis=0.8*axis.size )

	axis(2, at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
	        label = c("0%", "20%", "40%", "60%", "80%", "100%"),
	        las=2, cex.axis=axis.size )
}

dev.off()




#########################################
#########################################
##                                     ##
##  ####   #####  ####   ####  ##  ##  ##
##  ## ##  ##    ##  ## ##  ## ##  ##  ##
##  ##  ## ####  ##     ######  ####   ## 
##  ## ##  ##    ##  ## ##  ##   ##    ##
##  ####   #####  ####  ##  ##   ##    ##
##                                     ##  
#########################################
#########################################

index_6m <- which.min(abs(tt_plot - 180))
index_12m <- which.min(abs(tt_plot - 360))



AB_reduc_6m <- array(NA, dim=c(7,N_pos_kin,N_posterior))

for(k in 1:7)
{
	for(n in 1:N_pos_kin)
	{
		for(i in 1:N_posterior)
		{
			AB_reduc_6m[k,n,i] <- 1 - AB_mod[k,n,i,index_6m]/max(AB_mod[k,n,i,])
		}
	}
}

AB_reduc_6m_quant <- matrix(NA, nrow=7, ncol=3)

for(k in 1:7)
{
	AB_reduc_6m_quant[k,] <- quantile( AB_reduc_6m[k,,], prob=c(0.5, 0.025, 0.975) )
}		





AB_reduc_12m <- array(NA, dim=c(7,N_pos_kin,N_posterior))

for(k in 1:7)
{
	for(n in 1:N_pos_kin)
	{
		for(i in 1:N_posterior)
		{
			AB_reduc_12m[k,n,i] <- 1 - AB_mod[k,n,i,index_12m]/max(AB_mod[k,n,i,])
		}
	}
}

AB_reduc_12m_quant <- matrix(NA, nrow=7, ncol=3)

for(k in 1:7)
{
	AB_reduc_12m_quant[k,] <- quantile( AB_reduc_12m[k,,], prob=c(0.5, 0.025, 0.975) )
}		





SP_mod_quant[,,index_6m]

