library(binom)
library(MASS)
library(randomForest)
library(pROC)


N_tree <- 5000

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

AB_data <- AB_data[-which(AB_data$days_post <= 10),]

N_positive <- length(which(AB_data$status == "positive"))
N_negative <- length(which(AB_data$status == "negative"))


N_pos_kin <- 194


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

ant_names <- c("Stri_IPP_IgG_dil", "RBD_IPP_IgG_dil", "NP_IPP_IgG_dil", "S2_NA_IgG_dil")

N_tt_plot <- 200

tt_plot <- exp(seq(from=log(1), to=log(500), length=N_tt_plot))
log2 <- log(2)



N_posterior <- 150

AB_mod <- array(NA, dim=c(4, N_pos_kin, N_posterior, N_tt_plot ))


for(k in 1:4)
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


######################################################################
######################################################################
##                                                                  ##  
##   ####  ##     ####   ####   ####  #### ##### #### ##### #####   ## 
##  ##  ## ##    ##  ## ##     ##      ##  ##     ##  ##    ##  ##  ##
##  ##     ##    ######  ####   ####   ##  ####   ##  ####  #####   ## 
##  ##  ## ##    ##  ##     ##     ##  ##  ##     ##  ##    ## ##   ## 
##   ####  ##### ##  ##  ####   ####  #### ##    #### ##### ##  ##  ##
##                                                                  ##
######################################################################
######################################################################

######################################################################
##                                                                  ##
##  1 antigen                                                       ##
##                                                                  ##
######################################################################

Stri_cut <- quantile( AB_data[which(AB_data$status=="negative"),10], prob=0.99 ) 


Stri_class <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot )


for(i in 1:N_posterior)
{
	for(j in 1:N_tt_plot)
	{
		Stri_class[i,j] <- length(which(AB_mod[1,,i,j] > Stri_cut))/N_pos_kin
	}
}
			


Stri_class_quant <- matrix(NA, nrow=3, ncol=N_tt_plot )

for(j in 1:N_tt_plot)
{
	Stri_class_quant[,j] <- quantile( Stri_class[,j], prob=c(0.025, 0.5, 0.975) )
}

			



######################################################################
##                                                                  ##
##  2 antigens                                                      ##
##                                                                  ##
######################################################################


ant_seq2 <- c(10,12)

status_train <- AB_data$status
AB_train <- log(AB_data[,ant_seq2])

RF_2ant = randomForest( status_train ~ ., data=AB_train, 
                                      importance=TRUE, ntree=N_tree)

RF_2ant_roc <- roc(status_train, RF_2ant$votes[,2])

alpha2 <- RF_2ant_roc$thresholds[min(which(RF_2ant_roc$specificities > 0.99))]


RF_2score <- array(NA, dim=c(N_pos_kin, N_posterior, N_tt_plot ))

for(n in 1:N_pos_kin)
{
	for(i in 1:N_posterior)
	{
		AB_test <- log(AB_mod[1:2,n,i,])
		AB_test <- as.data.frame( t(AB_test) )
		colnames(AB_test) <- colnames(AB_train)

		RF_2ant_pred_obj = predict( RF_2ant, newdata=AB_test, predict.all=TRUE)

		RF_2ant_votes = rowSums(RF_2ant_pred_obj$individual=="negative")/N_tree

		RF_2score[n,i,] <- RF_2ant_votes
	}	
}

RF_2class <- array(0, dim=c(N_pos_kin, N_posterior, N_tt_plot ))

for(n in 1:N_pos_kin)
{
	for(i in 1:N_posterior)
	{
		index <- which(RF_2score[n,i,] < alpha2)
	
		if( length(index) > 0 )
		{
			RF_2class[n,i,index] <- 1
		}
	}
}

RF_2class_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	for(j in 1:N_tt_plot)
	{
		RF_2class_sens[i,j] <- mean(RF_2class[,i,j])
	}
}

RF_2class_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)

for(j in 1:N_tt_plot)
{
	RF_2class_sens_quant[,j] <- quantile( RF_2class_sens[,j], prob=c(0.025, 0.5, 0.975) )
}



######################################################################
##                                                                  ##
##  3 antigens                                                      ##
##                                                                  ##
######################################################################


ant_seq3 <- c(10,12,15)

status_train <- AB_data$status
AB_train <- log(AB_data[,ant_seq3])

index_rm <- unique(which(is.na(AB_train)==TRUE, arr.ind=TRUE)[,1])

status_train <- status_train[-index_rm]
AB_train <- AB_train[-index_rm,]




RF_3ant = randomForest( status_train ~ ., data=AB_train, 
                                      importance=TRUE, ntree=N_tree)

RF_3ant_roc <- roc(status_train, RF_3ant$votes[,2])

alpha3 <- RF_3ant_roc$thresholds[min(which(RF_3ant_roc$specificities > 0.99))]


RF_3score <- array(NA, dim=c(N_pos_kin, N_posterior, N_tt_plot ))

for(n in 1:N_pos_kin)
{
	for(i in 1:N_posterior)
	{
		AB_test <- log(AB_mod[1:3,n,i,])
		AB_test <- as.data.frame( t(AB_test) )
		colnames(AB_test) <- colnames(AB_train)

		RF_3ant_pred_obj = predict( RF_3ant, newdata=AB_test, predict.all=TRUE)

		RF_3ant_votes = rowSums(RF_3ant_pred_obj$individual=="negative")/N_tree

		RF_3score[n,i,] <- RF_3ant_votes
	}	
}

RF_3class <- array(0, dim=c(N_pos_kin, N_posterior, N_tt_plot ))

for(n in 1:N_pos_kin)
{
	for(i in 1:N_posterior)
	{
		index <- which(RF_3score[n,i,] < alpha3)
	
		if( length(index) > 0 )
		{
			RF_3class[n,i,index] <- 1
		}
	}
}

RF_3class_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	for(j in 1:N_tt_plot)
	{
		RF_3class_sens[i,j] <- mean(RF_3class[,i,j])
	}
}

RF_3class_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)

for(j in 1:N_tt_plot)
{
	RF_3class_sens_quant[,j] <- quantile( RF_3class_sens[,j], prob=c(0.025, 0.5, 0.975) )
}







######################################################################
##                                                                  ##
##  4 antigens                                                      ##
##                                                                  ##
######################################################################


ant_seq4 <- c(10,12,15,14)

status_train <- AB_data$status
AB_train <- log(AB_data[,ant_seq4])

index_rm <- unique(which(is.na(AB_train)==TRUE, arr.ind=TRUE)[,1])

status_train <- status_train[-index_rm]
AB_train <- AB_train[-index_rm,]


RF_4ant = randomForest( status_train ~ ., data=AB_train, 
                                      importance=TRUE, ntree=N_tree)

RF_4ant_roc <- roc(status_train, RF_4ant$votes[,2])

alpha4 <- RF_4ant_roc$thresholds[min(which(RF_4ant_roc$specificities > 0.99))]


RF_4score <- array(NA, dim=c(N_pos_kin, N_posterior, N_tt_plot ))

for(n in 1:N_pos_kin)
{
	for(i in 1:N_posterior)
	{
		AB_test <- log(AB_mod[1:4,n,i,])
		AB_test <- as.data.frame( t(AB_test) )
		colnames(AB_test) <- colnames(AB_train)

		RF_4ant_pred_obj = predict( RF_4ant, newdata=AB_test, predict.all=TRUE)

		RF_4ant_votes = rowSums(RF_4ant_pred_obj$individual=="negative")/N_tree

		RF_4score[n,i,] <- RF_4ant_votes
	}	
}

RF_4class <- array(0, dim=c(N_pos_kin, N_posterior, N_tt_plot ))

for(n in 1:N_pos_kin)
{
	for(i in 1:N_posterior)
	{
		index <- which(RF_4score[n,i,] < alpha4)
	
		if( length(index) > 0 )
		{
			RF_4class[n,i,index] <- 1
		}
	}
}

RF_4class_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	for(j in 1:N_tt_plot)
	{
		RF_4class_sens[i,j] <- mean(RF_4class[,i,j])
	}
}

RF_4class_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)

for(j in 1:N_tt_plot)
{
	RF_4class_sens_quant[,j] <- quantile( RF_4class_sens[,j], prob=c(0.025, 0.5, 0.975) )
}


save.image("Fig4_kinetic_classifier.RData")

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


tiff(file="Fig4_kinetic_classifier.tif", width=12, height=10, units="cm", res=300)

PP <- 10


point.size = 1
lab.size   = 1
axis.size  = 1
main.size  = 1


par(mar=c(4, 4, 1, 1))
par(mgp=c(1.5, 0.5, 0))



line_seq_x <- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)
line_seq_y <- c(1, 3, 5, 10, 30, 50, 100, 300, 500)


plot(x = 100, 
     y = 100,
     col=ant_cols[k], pch=19, cex=point.size,
     xlim=c(1,500), ylim=c(0,1), log="x",
     yaxt='n', xaxt='n', bty='n',
     ylab="", xlab="days after symptoms",
     main="",
     cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_x))
{
	points(x=c(0.1,500), y=rep(line_seq_x[i],2)^PP, type='l', col="grey", lty="dashed")
}
for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(-1,1), type='l', col="grey", lty="dashed")
}


polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( RF_4class_sens_quant[1,]^PP, rev(RF_4class_sens_quant[3,]^PP) ),
        col=rgb(190/256,190/256,190/256,0.4), border=NA)


points(x=tt_plot, y=RF_4class_sens_quant[2,]^PP, type='l', col="black")

points(x=tt_plot, y=RF_3class_sens_quant[2,]^PP, type='l', col="forestgreen")

points(x=tt_plot, y=RF_2class_sens_quant[2,]^PP, type='l', col="orangered")

points(x=tt_plot, y=Stri_class_quant[2,]^PP, type='l', col="red")

mtext(side = 2, line = 2.5, 
	cex=axis.size, 
	text="sensitivity")

mtext(side = 1, line = 1.5, 
	cex=axis.size, 
	text="days after symptom onset")






legend(x="topleft", 
cex=0.9, 
bg="white", box.col="white",
fill = c("red", "orangered", "forestgreen", "black"),
border = c("red", "orangered", "forestgreen", "black"), 
legend = c("Stri IgG", "+ RBDv2 IgG", "+ NPv1 IgG", "+ S2 IgG") )  

axis(1, at=c(1, 3, 5, 10, 30, 50, 100, 300, 500), 
        label=c("1", "3", "5", "10", "30", "50", "100", "300", "500"), 
        cex.axis=1 )

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=1 ) 


dev.off()




Stri_class_quant[,which.min(abs(tt_plot - 365))]



RF_2class_sens_quant[,which.min(abs(tt_plot - 365))]

RF_3class_sens_quant[,which.min(abs(tt_plot - 365))]

RF_4class_sens_quant[,which.min(abs(tt_plot - 365))]






RF_2class_sens_quant[,which.min(abs(tt_plot - 180))]

RF_3class_sens_quant[,which.min(abs(tt_plot - 180))]

RF_4class_sens_quant[,which.min(abs(tt_plot - 180))]





