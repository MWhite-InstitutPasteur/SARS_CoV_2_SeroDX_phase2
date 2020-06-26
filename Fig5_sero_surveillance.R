library(binom)
library(MASS)
library(randomForest)
library(pROC)


N_tree <- 10000 # 5000

N_rep <- 1000 # 1000

N_SS_cut <- 5000  #  5000

SS_seq <- seq(from=0, to=1, length=N_SS_cut)

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


site_update <- as.vector(AB_data$site)

site_update[which( AB_data$site == "Bichat" & AB_data$days_post <= 10 )] <- "Bichat_early"
site_update[which( AB_data$site == "Bichat" & AB_data$days_post > 10 )] <- "Bichat_late"

AB_data$site <- site_update

###################################
###################################
##                               ##
##   ##     ####  #   ## ######  ##
##  ###    ##  ## ##  ##   ##    ##
##   ##    ###### ######   ##    ##
##   ##    ##  ## ## ###   ##    ##
##  ####   ##  ## ##  ##   ##    ##
##                               ##
###################################
###################################

##N_dil_cut <- 10000 
##
#3dil_cut <- exp(seq(from=log(1e-5), to=log(0.3), length=N_dil_cut))



dil_cut <- sort(unique(AB_data$Stri_IPP_IgG_dil))

dil_cut <- c( 1e-5, dil_cut, 0.5*(dil_cut[-1] + dil_cut[-length(dil_cut)]), 0.1 )

dil_cut <- sort(dil_cut)

N_dil_cut <- length(dil_cut)


#########################################################
## Calculation of 1 antigen sensitivity and specificity
## with uncertainty calculated with binomial distribution

sens_1ant_binom <- matrix(NA, ncol=3, nrow=N_dil_cut)
spec_1ant_binom <- matrix(NA, ncol=3, nrow=N_dil_cut)

colnames(sens_1ant_binom) <- c( "med", "lwr", "upr" )
colnames(spec_1ant_binom) <- c( "med", "lwr", "upr" )

for(i in 1:N_dil_cut)
{
	N_pos <- length(which(is.na(AB_data$Stri_IPP_IgG_dil[which(AB_data$status == "positive")])==FALSE))
	N_neg <- length(which(is.na(AB_data$Stri_IPP_IgG_dil[which(AB_data$status == "negative")])==FALSE))

	sens_1ant_binom[i,1] <- length(which( AB_data$Stri_IPP_IgG_dil[which(AB_data$status == "positive")] > dil_cut[i] ))/N_pos
	spec_1ant_binom[i,1] <- length(which( AB_data$Stri_IPP_IgG_dil[which(AB_data$status == "negative")] < dil_cut[i] ))/N_neg

	sens_1ant_binom[i,2:3] <- as.numeric(as.vector(
                                  binom.confint( sens_1ant_binom[i,1]*N_positive, N_positive, method="wilson")[1,5:6] 
                                                   ))

	spec_1ant_binom[i,2:3] <- as.numeric(as.vector(
                                  binom.confint( spec_1ant_binom[i,1]*N_negative, N_negative, method="wilson")[1,5:6] 
                                                   ))
}


#########################################################
## Calculation of 1 antigen sensitivity and specificity
## with cross-validated uncertainty
## (essentially bootstrapping)



Stri_IgG_pos <- AB_data$Stri_IPP_IgG_dil[which(AB_data$status == "positive")]
Stri_IgG_neg <- AB_data$Stri_IPP_IgG_dil[which(AB_data$status == "negative")]


sens_1ant_xval_mat <- matrix(NA, nrow=N_dil_cut, ncol=N_rep)
spec_1ant_xval_mat <- matrix(NA, nrow=N_dil_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	Stri_IgG_pos_sub <- Stri_IgG_pos[sample( length(Stri_IgG_pos), (2/3)*length(Stri_IgG_pos) )]		
	Stri_IgG_neg_sub <- Stri_IgG_neg[sample( length(Stri_IgG_neg), (2/3)*length(Stri_IgG_neg) )]		

	for(i in 1:N_dil_cut)
	{
		sens_1ant_xval_mat[i,n] <- length(which( Stri_IgG_pos_sub > dil_cut[i] ))/length(Stri_IgG_pos_sub)
 		spec_1ant_xval_mat[i,n] <- length(which( Stri_IgG_neg_sub < dil_cut[i] ))/length(Stri_IgG_neg_sub)
	}
}


sens_1ant_xval_quant <- matrix(NA, nrow=N_dil_cut, ncol=3)
spec_1ant_xval_quant <- matrix(NA, nrow=N_dil_cut, ncol=3)

colnames(sens_1ant_xval_quant) <- c( "med", "lwr", "upr" )
colnames(spec_1ant_xval_quant) <- c( "med", "lwr", "upr" )

for(i in 1:N_dil_cut)
{
	sens_1ant_xval_quant[i,] <- quantile( sens_1ant_xval_mat[i,], prob=c(0.5, 0.025, 0.975) )
	spec_1ant_xval_quant[i,] <- quantile( spec_1ant_xval_mat[i,], prob=c(0.5, 0.025, 0.975) )
}




##############################################################
##############################################################
##                                                          ##
##  ##    # ##  ## ##   ###### ####    ####  #   ## ######  ##
##  ##   ## ##  ## ##     ##    ##    ##  ## ##  ##   ##    ##
##  ####### ##  ## ##     ##    ##    ###### ### ##   ##    ##
##  ## # ## ##  ## ##     ##    ##    ##  ## ## ###   ##    ## 
##  ##   ##  ##### #####  ##   ####   ##  ## ##  ##   ##    ##
##                                                          ##
##############################################################
##############################################################

status <- AB_data$status

AB <- log(AB_data[,c(10:16,43:47,23:29,56:60)])

AB <- AB[,c(1,3,6,5,14,18)]

index_trim <- unique(which(is.na(AB)==TRUE, arr.ind=TRUE)[,1])

AB_trim <- AB[-index_trim,]
status_trim <- status[-index_trim]



##############################################
## Uncertainty in sensitivity and specificity 
## for algorithm fitted to all the data

RF_multi_ant = randomForest( status_trim ~ ., data=AB_trim, 
                                        importance=TRUE, ntree=N_tree)

RF_multi_ant_ROC <- roc(status_trim, RF_multi_ant$votes[,2])

varImpPlot( RF_multi_ant )

N_SS_multi_ant <- length(RF_multi_ant_ROC$sensitivities)



##############################################
## Calculate uncertainty in sensitivity 

sens_multi_ant_binom <- matrix(NA, nrow=length(RF_multi_ant_ROC$sensitivities), ncol=3)
colnames(sens_multi_ant_binom) <- c("med", "lwr", "upr")

for(i in 1:nrow(sens_multi_ant_binom))
{
	N_pos <- length(which(status_trim=="positive"))

	sens_multi_ant_binom[i,] <- as.numeric(as.vector(
                                          binom.confint( RF_multi_ant_ROC$sensitivities[i]*N_pos, 
                                                         N_pos, method="wilson")[1,4:6] ))
}

##############################################
## Calculate uncertainty in sensitivity 

spec_multi_ant_binom <- matrix(NA, nrow=length(RF_multi_ant_ROC$specificities), ncol=3)
colnames(spec_multi_ant_binom) <- c("med", "lwr", "upr")

for(i in 1:nrow(spec_multi_ant_binom))
{
	N_neg <- length(which(status_trim=="negative"))

	spec_multi_ant_binom[i,] <- as.numeric(as.vector(
                                       binom.confint( RF_multi_ant_ROC$specificities[i]*N_neg, 
                                                        N_neg, method="wilson")[1,4:6] ))
}



####################################
####################################
##                                ##
##  ##   ## ##   ##  ####  ##     ##
##   ## ##  ##   ## ##  ## ##     ##
##    ###    ## ##  ###### ##     ##
##   ## ##    ###   ##  ## ##     ##
##  ##   ##    #    ##  ## #####  ##
##                                ##
####################################
####################################

####################################
## Uncertainty in sensitivity and specificity 
## for algorithmvia cross-validation - randomly
## leaving 1/3 out

#############################################################################
## Random forests classifier for disjoint training and testing data

RF_multi_ant_SS = function( AB_train, status_train, AB_test, status_test, antigen_list )
{
	AB_train = as.matrix(AB_train[,antigen_list])

	if( length(which(is.na(AB_train))) > 0 )
	{
		status_train = status_train[-which(is.na(AB_train), arr.ind=TRUE)[,1]]

		AB_train = AB_train[-which(is.na(AB_train), arr.ind=TRUE)[,1],]
	}


	AB_test = as.matrix(AB_test[,antigen_list])

	if( length(which(is.na(AB_test))) > 0 )
	{
		status_test = status_test[-which(is.na(AB_test), arr.ind=TRUE)[,1]]

		AB_test = AB_test[-which(is.na(AB_test), arr.ind=TRUE)[,1],]
	}

	############################################
	## Fit Rforest and create prediction object

	tryCatch(
	{
		RF_multi_ant = randomForest( as.factor(status_train) ~ ., data=as.data.frame(AB_train), 
                                      importance=TRUE, ntree=N_tree )

		RF_multi_ant_pred_obj = predict( RF_multi_ant, newdata=AB_test, predict.all=TRUE)

		RF_multi_ant_votes = rowSums(RF_multi_ant_pred_obj$individual=="negative")/N_tree

		RF_multi_ant_roc <- roc(status_test, RF_multi_ant_votes)

	}, error=function(e){ NULL }
	)

	SS_mat = cbind( RF_multi_ant_roc$sensitivities, RF_multi_ant_roc$specificities  )
	colnames(SS_mat) = c("sens", "spec" )

	SS_mat
}


#############################################################################
## Cross-validated random forests classification algorithm

RF_multi_ant_SS_xval = function( AB_dd, status_data, antigen_list, train_prop, N_rep )
{
	AB_class <- list()	

	for(n in 1:N_rep)
	{
		############################################
		## Prepare testing and training data

		index_train = sample( nrow(AB_dd), train_prop*nrow(AB_dd) )
		index_test  = setdiff( 1:nrow(AB_dd), index_train )

		status_train = status_data[index_train]
		status_test  = status_data[index_test]

		AB_train = AB_dd[index_train,]
		AB_test  = AB_dd[index_test,]

		
		AB_class[[n]] = RF_multi_ant_SS( AB_train, status_train, AB_test, status_test, antigen_list )
	}

	AB_class
}

#########################################
## Cross-validated algorithm (mixed)

RF_multi_ant_xval_rep <- RF_multi_ant_SS_xval( AB, status, antigen_list=c(1,2,3,4,5,6), train_prop=2/3, N_rep=N_rep )



#########################################
## Cross-validated variation in sensitivity

RF_multi_ant_xval_se_vary_sp_fix <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:N_SS_cut)
	{
		RF_multi_ant_xval_se_vary_sp_fix[i,n] <- RF_multi_ant_xval_rep[[n]][min(which(RF_multi_ant_xval_rep[[n]][,2] > SS_seq[i])),1]
	}
}

RF_multi_ant_xval_se_vary_sp_fix[N_SS_cut,] <- 0


RF_multi_ant_xval_se_vary_sp_fix_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
colnames(RF_multi_ant_xval_se_vary_sp_fix_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_cut)
{
	RF_multi_ant_xval_se_vary_sp_fix_quant[i,] <- quantile( RF_multi_ant_xval_se_vary_sp_fix[i,], prob=c(0.5, 0.025, 0.975) )
}


#########################################
## Cross-validated variation in specificity

RF_multi_ant_xval_se_fix_sp_vary <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:N_SS_cut)
	{
		RF_multi_ant_xval_se_fix_sp_vary[i,n] <- RF_multi_ant_xval_rep[[n]][min(which(RF_multi_ant_xval_rep[[n]][,1] < SS_seq[i])),2]
	}
}

RF_multi_ant_xval_se_fix_sp_vary[1,] <- 1

RF_multi_ant_xval_se_fix_sp_vary_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
colnames(RF_multi_ant_xval_se_fix_sp_vary_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_cut)
{
	RF_multi_ant_xval_se_fix_sp_vary_quant[i,] <- quantile( RF_multi_ant_xval_se_fix_sp_vary[i,], prob=c(0.5, 0.025, 0.975) )
}


####################################
####################################
##                                ## 
##   ####  ##  ## #####  ##   ##  ##
##  ##     ##  ## ##  ## ##   ##  ## 
##   ####  ##  ## #####   ## ##   ## 
##      ## ##  ## ## ##    ###    ##
##   ####   ####  ##  ##    #     ##
##                                ##
####################################
####################################

####################################
##                                ##
##   Multiple antigens            ## 
##                                ##
####################################

SS_multi_ant <- matrix(NA, nrow=N_SS_multi_ant, ncol=2+2*N_rep)
colnames(SS_multi_ant) <- c("sens_hat", "spec_hat", 1:(2*N_rep))

for(j in 1:N_rep)
{
	colnames(SS_multi_ant)[2+j] <- paste("sens_rep_", j, sep="")
	colnames(SS_multi_ant)[N_rep+2+j] <- paste("spec_rep_", j, sep="")
}

SS_multi_ant[,1] <- RF_multi_ant_ROC$sensitivities
SS_multi_ant[,2] <- RF_multi_ant_ROC$specificities

for(i in 1:N_SS_multi_ant)
{
	for(j in 1:N_rep)
	{
		sens_index <- max(which(RF_multi_ant_ROC$sensitivities[i] < RF_multi_ant_xval_rep[[j]][,1]))

		SS_multi_ant[i,2+N_rep+j] <- RF_multi_ant_xval_rep[[j]][sens_index,2]

		spec_index <- max(which(RF_multi_ant_ROC$specificities[i] > RF_multi_ant_xval_rep[[j]][,2]))+1
		
		SS_multi_ant[i,2+j] <- RF_multi_ant_xval_rep[[j]][spec_index,1]
	}
}
	
SS_multi_ant[which(SS_multi_ant[,1]==1),3:(2+N_rep)] <- 1
SS_multi_ant[which(SS_multi_ant[,1]==1),(3+N_rep):(2+2*N_rep)] <- 0


####################################
##                                ##
##   Single antigen               ## 
##                                ##
####################################

SS_one_ant <- matrix(NA, nrow=N_dil_cut, ncol=2+2*N_rep)
colnames(SS_one_ant) <- c("sens_hat", "spec_hat", 1:(2*N_rep))

for(j in 1:N_rep)
{
	colnames(SS_one_ant)[2+j] <- paste("sens_rep_", j, sep="")
	colnames(SS_one_ant)[N_rep+2+j] <- paste("spec_rep_", j, sep="")
}

SS_one_ant[,1] <- sens_1ant_binom[,1]
SS_one_ant[,2] <- spec_1ant_binom[,1]

for(i in 1:N_dil_cut)
{
	for(j in 1:N_rep)
	{
		sens_index <- max(which(sens_1ant_binom[i,1] < sens_1ant_xval_mat[,j]))

		SS_one_ant[i,2+N_rep+j] <- spec_1ant_xval_mat[sens_index,j]

		spec_index <- max(which(spec_1ant_binom[i,1] > spec_1ant_xval_mat[,j]))+1
		
		SS_one_ant[i,2+j] <- sens_1ant_xval_mat[spec_index,j]
	}
}
	
SS_one_ant[which(SS_one_ant[,1]==1),3:(2+N_rep)] <- 1
SS_one_ant[which(SS_one_ant[,2]==0),3:(2+N_rep)] <- 1
SS_one_ant[which(SS_one_ant[,1]==1),(3+N_rep):(2+2*N_rep)] <- 0





#####################
#####################
##                 ##
##  #####  ##  ##  ##
##  ##        ##   ## 
##  ####     ##    ##
##     ##   ##     ##
##  ####   ##  ##  ##
##                 ##
#####################
#####################


T_prev <- 0.05


########################################################
##                                                    ## 
##  Multi-antigen: 5% true prevalence                 ## 
##                                                    ##
########################################################

M_multi_ant_5prev <- T_prev*SS_multi_ant[,(2+1):(2+N_rep)] + (1 - T_prev)*(1 - SS_multi_ant[,(2+N_rep+1):(2+2*N_rep)])

M_multi_ant_5prev_quant <- matrix(NA, nrow=N_SS_multi_ant, ncol=3)
colnames(M_multi_ant_5prev_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_multi_ant)
{
	M_multi_ant_5prev_quant[i,] <- quantile( M_multi_ant_5prev[i,], prob=c(0.5, 0.025, 0.975) )
}


##################################################
## Variation in sensitivity and specificity

spec_mat <- matrix( rep(RF_multi_ant_ROC$specificities, N_rep), ncol=N_rep)
sens_mat <- matrix( rep(RF_multi_ant_ROC$sensitivities, N_rep), ncol=N_rep)


##################################################
## Adjusted prevalence

A_multi_ant_5prev <- (M_multi_ant_5prev + spec_mat - 1)/( sens_mat + spec_mat - 1)

A_multi_ant_5prev[which(A_multi_ant_5prev < 0, arr.ind=TRUE)] <- 0
A_multi_ant_5prev[which(is.na(A_multi_ant_5prev) == TRUE, arr.ind=TRUE)] <- 0

A_multi_ant_5prev_quant <- matrix(NA, nrow=N_SS_multi_ant, ncol=3)
colnames(A_multi_ant_5prev_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_multi_ant)
{
	A_multi_ant_5prev_quant[i,] <- quantile( A_multi_ant_5prev[i,], prob=c(0.5, 0.025, 0.975) )
}


##################################################
## Relative error in adjusted prevalence

##A_multi_ant_5prev_rel_error <- abs(A_multi_ant_5prev/T_prev)
	
A_multi_ant_5prev_rel_error <- abs(A_multi_ant_5prev - T_prev)/T_prev

A_multi_ant_5prev_rel_error <- rowSums(A_multi_ant_5prev_rel_error)/N_rep


##################################################
## Absolute error in adjusted prevalence

A_multi_ant_5prev_abs_error <- abs(A_multi_ant_5prev - T_prev)

A_multi_ant_5prev_abs_error <- rowSums(A_multi_ant_5prev_abs_error)/N_rep





########################################################
##                                                    ## 
##  Single antigen: 5% true prevalence                 ## 
##                                                    ##
########################################################

M_one_ant_5prev <- T_prev*SS_one_ant[,(2+1):(2+N_rep)] + (1 - T_prev)*(1 - SS_one_ant[,(2+N_rep+1):(2+2*N_rep)])

M_one_ant_5prev_quant <- matrix(NA, nrow=N_dil_cut, ncol=3)
colnames(M_one_ant_5prev_quant) <- c("med", "lwr", "upr")

for(i in 1:N_dil_cut)
{
	M_one_ant_5prev_quant[i,] <- quantile( M_one_ant_5prev[i,], prob=c(0.5, 0.025, 0.975) )
}


##################################################
## Variation in sensitivity and specificity

spec_mat <- matrix( rep(spec_1ant_binom[,1], N_rep), ncol=N_rep)
sens_mat <- matrix( rep(sens_1ant_binom[,1], N_rep), ncol=N_rep)


##################################################
## Adjusted prevalence

A_one_ant_5prev <- (M_one_ant_5prev + spec_mat - 1)/( sens_mat + spec_mat - 1)

A_one_ant_5prev[which(A_one_ant_5prev < 0, arr.ind=TRUE)] <- 0
A_one_ant_5prev[which(is.na(A_one_ant_5prev) == TRUE, arr.ind=TRUE)] <- 0

A_one_ant_5prev_quant <- matrix(NA, nrow=N_dil_cut, ncol=3)
colnames(A_one_ant_5prev_quant) <- c("med", "lwr", "upr")

for(i in 1:N_dil_cut)
{
	A_one_ant_5prev_quant[i,] <- quantile( A_one_ant_5prev[i,], prob=c(0.5, 0.025, 0.975) )
}


##################################################
## Relative error in adjusted prevalence

##A_one_ant_5prev_rel_error <- abs(A_one_ant_5prev/T_prev)
	
A_one_ant_5prev_rel_error <- abs(A_one_ant_5prev - T_prev)/T_prev

A_one_ant_5prev_rel_error <- rowSums(A_one_ant_5prev_rel_error)/N_rep


##################################################
## Absolute error in adjusted prevalence

A_one_ant_5prev_abs_error <- abs(A_one_ant_5prev - T_prev)

A_one_ant_5prev_abs_error <- rowSums(A_one_ant_5prev_abs_error)/N_rep



########################################################
########################################################
##                                                    ##
##   ####  #####  ######   ####  #####  #####  ####   ## 
##  ##  ## ##  ##   ##    ##     ##  ## ##    ##  ##  ##
##  ##  ## #####    ##     ####  #####  ####  ##      ##
##  ##  ## ##       ##        ## ##     ##    ##  ##  ##
##   ####  ##       ##     ####  ##     #####  ####   ## 
##                                                    ##
########################################################
########################################################

drop_sens <- which(RF_multi_ant_ROC$sensitivities == 1)
drop_sens <- drop_sens[-length(drop_sens)]

drop_spec <- which(RF_multi_ant_ROC$specificities==1)
drop_spec <- drop_spec[-1]

drop_multi <- c(drop_sens, drop_spec)



########################################################
##                                                    ## 
##  Multi-antigen prevalence estimating function      ## 
##                                                    ##
########################################################

T_prev_multi_ant_func <- function( T_prevv, plot_out=TRUE )
{
	M_prev <- T_prevv*SS_multi_ant[,(2+1):(2+N_rep)] + (1 - T_prevv)*(1 - SS_multi_ant[,(2+N_rep+1):(2+2*N_rep)])

	M_prev_quant <- matrix(NA, nrow=N_SS_multi_ant, ncol=3)
	colnames(M_prev_quant) <- c("med", "lwr", "upr")

	for(i in 1:N_SS_multi_ant)
	{
		M_prev_quant[i,] <- quantile( M_prev[i,], prob=c(0.5, 0.025, 0.975) )
	}


	##################################################
	## Variation in sensitivity and specificity
	
	spec_mat <- matrix( rep(RF_multi_ant_ROC$specificities, N_rep), ncol=N_rep)
	sens_mat <- matrix( rep(RF_multi_ant_ROC$sensitivities, N_rep), ncol=N_rep)


	##################################################
	## Adjusted prevalence

	A_prev <- (M_prev + spec_mat - 1)/( sens_mat + spec_mat - 1)

	A_prev[which(A_prev < 0, arr.ind=TRUE)] <- 0
	A_prev[which(is.na(A_prev) == TRUE, arr.ind=TRUE)] <- 0

	A_prev_quant <- matrix(NA, nrow=N_SS_multi_ant, ncol=3)
	colnames(A_prev_quant) <- c("med", "lwr", "upr")

	for(i in 1:N_SS_multi_ant)
	{
		A_prev_quant[i,] <- quantile( A_prev[i,], prob=c(0.5, 0.025, 0.975) )
	}


	##################################################
	## Relative error in adjusted prevalence

	##A_prev_rel_error <- abs(A_prev - T_prevv)/T_prevv

	##A_prev_rel_error <- rowSums(A_prev_rel_error)/N_rep

	A_prev_rel_error <- abs(A_prev - T_prevv)/T_prevv

	A_prev_rel_error <- rowSums(A_prev_rel_error)/N_rep


	##################################################
	## Absolute error in adjusted prevalence

	##A_prev_abs_error <- abs(A_prev - T_prevv)

	##A_prev_abs_error <- rowSums(A_prev_abs_error)/N_rep

	A_prev_abs_error <- abs(A_prev - T_prevv)

	A_prev_abs_error <- rowSums(A_prev_abs_error)/N_rep



	##################################################
	##################################################
	## Summarise findings in a plot

	if( plot_out==TRUE )
	{
		par(mfrow=c(2,2))
		par(mar=c(3,3,2,1))

		##################################################
		## Panel A: measured prevalence

		plot(x=100, y=100, type='l', lty="dashed",
		xlim=c(0,0.1), ylim=c(0,0.5),
		#xaxs='i', yaxs='i', 
		xlab="1 - specificity", ylab="prevalence",
		main="(A) Measured prevalence")

		points( x=1 - RF_multi_ant_ROC$specificities, y=M_prev_quant[,1], 
    			  type='l', lwd=2, col="black" )

		for(i in 1:N_SS_multi_ant)
		{
			polygon(x=1 - RF_multi_ant_ROC$specificities[c(i-1,i,i,i-1)], 
			  y=c( M_prev_quant[c(i,i),2], M_prev_quant[c(i-1,i-1),3] ), 
      	  	col=rgb(190/256,190/256,190/256,0.4), border=NA)
		}

		points(x=c(0,1), y=rep(T_prevv,2), type='l', col="red")


		##################################################
		## Panel B: adjusted prevalence

		plot(x=100, y=100, type='l', lty="dashed",
		xlim=c(0,0.1), ylim=c(0,0.5), 
		#xaxs='i', yaxs='i', 
		xlab="1 - specificity", ylab="prevalence",
		main="(B) Adjusted prevalence")

		points( x=1 - RF_multi_ant_ROC$specificities, y=A_prev_quant[,1], 
    		  	  type='l', lwd=2, col="black" )
	
		for(i in 1:N_SS_multi_ant)
		{
			polygon(x=1 - RF_multi_ant_ROC$specificities[c(i-1,i,i,i-1)], 
			  y=c( A_prev_quant[c(i,i),2], A_prev_quant[c(i-1,i-1),3] ), 
      	  	col=rgb(190/256,190/256,190/256,0.4), border=NA)
		}

		points(x=c(0,1), y=rep(T_prevv,2), type='l', col="red")


		##################################################
		## Panel C: Expected relative error

		plot(x=1 - RF_multi_ant_ROC$specificities, y=A_prev_rel_error, 
		    	  type='l', lwd=2, col="black",
		xlim=c(0,0.1), log="y",  
		#xaxs='i', yaxs='i', 
		xlab="1 - specificity",
		main="(C) Expected relative error")

	

		##################################################
		## Panel D: Expected relative error

		plot(x=1 - RF_multi_ant_ROC$specificities, y=A_prev_abs_error, 
	    		  type='l', lwd=2, col="black",
		xlim=c(0,0.1), 
		#xaxs='i', yaxs='i', 
		xlab="1 - specificity",
		main="(D) Expected absolute error")

	}	

	##################################################
	##################################################
	## Summarise findings in a plot

	sens_out <- RF_multi_ant_ROC$sensitivities#[which( RF_multi_ant_ROC$specificities > 0 & RF_multi_ant_ROC$specificities < 1 )]
	spec_out <- RF_multi_ant_ROC$specificities#[which( RF_multi_ant_ROC$specificities > 0 & RF_multi_ant_ROC$specificities < 1 )]
	A_prev_rel_error_out <- A_prev_rel_error#[which( RF_multi_ant_ROC$specificities > 0 & RF_multi_ant_ROC$specificities < 1 )]
	A_prev_abs_error_out <- A_prev_abs_error#[which( RF_multi_ant_ROC$specificities > 0 & RF_multi_ant_ROC$specificities < 1 )]


	OUTP <- c( sens_out[which.min(A_prev_rel_error_out)],
		     spec_out[which.min(A_prev_rel_error_out)],
                 min(A_prev_rel_error_out),
                 min(A_prev_abs_error_out) )

	OUTP	

}


########################################################
##                                                    ## 
##  Single antigen prevalence estimating function     ## 
##                                                    ##
########################################################

T_prev_1ant_func <- function( T_prevv, plot_out=TRUE )
{
	##################################################
	## Measured prevalence  

	M_prev <-  T_prevv*SS_one_ant[,(2+1):(2+N_rep)] + (1 - T_prevv)*(1 - SS_one_ant[,(2+N_rep+1):(2+2*N_rep)])

	M_prev_quant <- matrix(NA, nrow=N_dil_cut, ncol=3)
	colnames(M_prev_quant) <- c("med", "lwr", "upr")

	for(i in 1:N_dil_cut)
	{
		M_prev_quant[i,] <- quantile( M_prev[i,], prob=c(0.5, 0.025, 0.975) )
	}

	##################################################
	## Variation in sensitivity and specificity

	#spec_mat <- matrix( rep(spec_1ant_xval_quant[,1], N_rep), ncol=N_rep)
	#sens_mat <- matrix( rep(sens_1ant_xval_quant[,1], N_rep), ncol=N_rep)

	spec_mat <- matrix( rep(spec_1ant_binom[,1], N_rep), ncol=N_rep)
	sens_mat <- matrix( rep(sens_1ant_binom[,1], N_rep), ncol=N_rep)


	##################################################
	## Adjusted prevalence

	A_prev <- (M_prev + spec_mat - 1)/( sens_mat + spec_mat - 1)

	A_prev[which(A_prev < 0, arr.ind=TRUE)] <- 0
	A_prev[which(is.na(A_prev) == TRUE, arr.ind=TRUE)] <- 0

	A_prev_quant <- matrix(NA, nrow=N_dil_cut, ncol=3)
	colnames(A_prev_quant) <- c("med", "lwr", "upr")

	for(i in 1:N_dil_cut)
	{
		A_prev_quant[i,] <- quantile( A_prev[i,], prob=c(0.5, 0.025, 0.975) )
	}


	##################################################
	## Relative error in adjusted prevalence

	##A_prev_rel_error <- abs(A_prev/T_prevv)
		
	A_prev_rel_error <- abs(A_prev - T_prevv)/T_prevv

	A_prev_rel_error <- rowSums(A_prev_rel_error)/N_rep


	##################################################
	## Absolute error in adjusted prevalence

	A_prev_abs_error <- abs(A_prev - T_prevv)

	A_prev_abs_error <- rowSums(A_prev_abs_error)/N_rep


	##################################################
	##################################################
	## Summarise findings in a plot

	if( plot_out==TRUE )
	{
		par(mfrow=c(2,2))

		##################################################
		## Panel A: measured prevalence

		plot(x=100, y=100, type='l', lty="dashed",
		xlim=c(0,0.1), ylim=c(0,0.5),
		#xaxs='i', yaxs='i', 
		xlab="1 - specificity", ylab="prevalence",
		main="(A) Measured prevalence")

		points( x=1 - spec_1ant_binom[,1], y=M_prev_quant[,1], 
    			  type='l', lwd=2, col="black" )

		for(i in 1:N_dil_cut)
		{
			polygon(x=1 - spec_1ant_binom[c(i-1,i,i,i-1),1], 
			  y=c( M_prev_quant[c(i,i),2], M_prev_quant[c(i-1,i-1),3] ), 
      	  	col=rgb(190/256,190/256,190/256,0.4), border=NA)
		}

		points(x=c(0,1), y=rep(T_prevv,2), type='l', col="red")


		##################################################
		## Panel B: adjusted prevalence

		plot(x=100, y=100, type='l', lty="dashed",
		xlim=c(0,0.1), ylim=c(0,0.5), 
		#xaxs='i', yaxs='i', 
		xlab="1 - specificity", ylab="prevalence",
		main="(B) Adjusted prevalence")

		points( x=1 - spec_1ant_binom[,1], y=A_prev_quant[,1], 
    		  	  type='l', lwd=2, col="black" )
	
		for(i in 1:N_dil_cut)
		{
			polygon(x=1 - spec_1ant_binom[c(i-1,i,i,i-1),1], 
			  y=c( A_prev_quant[c(i,i),2], A_prev_quant[c(i-1,i-1),3] ), 
      	  	col=rgb(190/256,190/256,190/256,0.4), border=NA)
		}

		points(x=c(0,1), y=rep(T_prevv,2), type='l', col="red")


		##################################################
		## Panel C: Expected relative error

		plot(x=1 - spec_1ant_binom[,1], y=A_prev_rel_error, 
		    	  type='l', lwd=2, col="black",
		xlim=c(0,0.1),  log="y",  
		#xaxs='i', yaxs='i', 
		xlab="1 - specificity",
		main="(C) Expected relative error")



		##################################################
		## Panel D: Expected relative error

		plot(x=1 - spec_1ant_binom[,1], y=A_prev_abs_error, 
	    		  type='l', lwd=2, col="black",
		xlim=c(0,0.1),
		#xaxs='i', yaxs='i', 
		xlab="1 - specificity",
		main="(D) Expected absolute error")


	}	

	##################################################
	##################################################
	## Summarise findings in a plot

	sens_out <- sens_1ant_binom#[which( spec_1ant_binom[,1] > 0 & spec_1ant_binom[,1] < 1 ),1]
	spec_out <- spec_1ant_binom#[which( spec_1ant_binom[,1] > 0 & spec_1ant_binom[,1] < 1 ),1]
	A_prev_rel_error_out <- A_prev_rel_error#[which( spec_1ant_binom[,1] > 0 & spec_1ant_binom[,1] < 1 )]
	A_prev_abs_error_out <- A_prev_abs_error#[which( spec_1ant_binom[,1] > 0 & spec_1ant_binom[,1] < 1 )]


	OUTP <- c( sens_out[which.min(A_prev_rel_error_out)],
		     spec_out[which.min(A_prev_rel_error_out)],
                 min(A_prev_rel_error_out),
                 min(A_prev_abs_error_out) )

	OUTP	

}


########################################################
##                                                    ## 
##  Calculate optimum sensitivity and specificity     ##
##  as a function of prevalence                       ## 
##                                                    ##
########################################################



N_T_prev <- 3000

T_prev_seq <- exp(seq(from=log(0.0001), to=log(1), length=N_T_prev))

T_prev_multi_ant <- matrix(NA, nrow=N_T_prev, ncol=4)
colnames(T_prev_multi_ant) <- c("sensitivity", "specificity", "rel_error", "abs_error")

for(k in 1:N_T_prev)
{
	T_prev_multi_ant[k,] <- T_prev_multi_ant_func( T_prev_seq[k], plot_out=FALSE )
}


T_prev_1ant <- matrix(NA, nrow=N_T_prev, ncol=4)
colnames(T_prev_1ant) <- c("sensitivity", "specificity", "rel_error", "abs_error")

for(k in 1:N_T_prev)
{
	T_prev_1ant[k,] <- T_prev_1ant_func( T_prev_seq[k], plot_out=FALSE )
}


save.image("Fig5_sero_surveillance.RData")

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



PP <- 10

tiff(file="Fig5_sero_surveillance.tif", width=25, height=40, units="cm", res=500)

lay.mat <- rbind( c( 1, 1 ),
                  c( 2, 2 ),
                  c( 3, 4 ),
                  c( 5, 5 ),
                  c( 6, 7 ),
                  c( 8, 8 ),
                  c( 9,10 ),
                  c(11,11 ) )
layout(lay.mat, heights=c(20,1.5,10,1.5,10,1.5,10,11), widths=c(1,1))
layout.show(11)







#####################################
#####################################
##                                 ## 
##  PANEL 1                        ##
##  Approach 2: cross-validated    ##
##  Variation in specificity       ##
##                                 ##
#####################################
##################################### 


lab.size   = 2
axis.size  = 1.75
main.size  = 3
line.size  = 3

par(mar=c(5,7,2,1.5))
par(mgp=c(2.5, 1.25,0))


line_seq_x <- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)
line_seq_y <- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="(A) ROC analysis with cross-validation",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3.5, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 4.5, 
cex=lab.size, 
text="sensitivity")

for(i in 1:length(line_seq_y))
{
	points(x=c(0,1), y=rep(line_seq_y[i]^PP,2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_x))
{
	points(x=rep(1 - line_seq_x[i]^PP,2), y=c(0,1), type='l', col="grey", lty="dashed")
}


points( x=1 - spec_1ant_binom[,1]^PP, y=sens_1ant_binom[,1]^PP, 
    	  type='S', lwd=line.size, col="orangered" )

for(i in 1:(N_dil_cut-1))
{
	polygon(x=c(1 - spec_1ant_xval_quant[c(i,i+1),2]^PP, rev(1 - spec_1ant_xval_quant[c(i,i+1),3]^PP)), 
		  y=c( sens_1ant_binom[c(i,i+1),1]^PP, rev(sens_1ant_binom[c(i,i+1),1]^PP) ),
	        col=rgb(190/256,90/256,90/256,0.2), border=NA)
}


points( x=1 - RF_multi_ant_xval_se_fix_sp_vary_quant[,1]^PP, y=SS_seq^PP, 
    	  type='S', lwd=line.size, col="aquamarine4" )


for(i in 1:(N_SS_cut-1))
{
	polygon(x=c(1 - RF_multi_ant_xval_se_fix_sp_vary_quant[c(i,i+1),2]^PP, rev(1 - RF_multi_ant_xval_se_fix_sp_vary_quant[c(i,i+1),3]^PP)), 
		  y=c( SS_seq[c(i,i+1)]^PP, rev(SS_seq[c(i,i+1)]^PP) ),
      	  col=rgb(69/256,139/256,116/256,0.2), border=NA)
}



legend(x="bottomright", 
cex=3, 
bg="white", box.col="white",
fill = c("orangered", "aquamarine4"),
border = c("orangered", "aquamarine4"), 
legend = c("Stri IgG", "multiplex") )  

axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 


############################
## Label (B)              ## 
############################

par(mar = c(0,0,0,0))


plot.new()
title( "(B) Measured seroprevalence", 
        cex.main=3.0, line=-2)



#####################################
#####################################
##                                 ## 
##  PANEL 2                        ##
##  Measured prevalence            ##
##  single antigen                 ##
##                                 ##
#####################################
##################################### 

lab.size   = 1.5
axis.size  = 1.75
main.size  = 3
line.size  = 3

par(mar=c(5,8,1.5,1.5))
par(mgp=c(2.5, 1.25,0))


line_seq_x <- c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1)
line_seq_y <- c(0.0, 0.05, 0.1, 0.15, 0.2)

plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0.0,0.1), ylim=c(0,0.2), 
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3.4, 
cex=lab.size, 
text="1 - specificity (= false positive rate)")

mtext(side = 2, line = 4.5, 
cex=lab.size, 
text="seroprevalence")

for(i in 1:length(line_seq_y))
{
	points(x=c(1e-5,1), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}

pos_index <- which( 1 - spec_1ant_binom[,1] > -1 )

points( x=(1 - spec_1ant_binom[pos_index,1]), y=M_one_ant_5prev_quant[pos_index,1], 
	  type='l', lwd=line.size, col="orangered" )

polygon(x=c( (1 - spec_1ant_binom[pos_index,1]), rev( (1 - spec_1ant_binom[pos_index,1]) ) ),
	  y=c( M_one_ant_5prev_quant[pos_index,2], rev(M_one_ant_5prev_quant[pos_index,3]) ), 
     	  col=rgb(190/256,90/256,90/256,0.2), border=NA)


points(x=c(1e-5,1), y=rep(0.05,2), 
type='l', col="black", lty="dashed", lwd=line.size)

text(x=0.07, y=0.035,
labels="true seroprevalence = 5%", cex=1.75)

axis(1, at=c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1), 
        labels=c("0%", "2%", "4%", "6%", "8%", "10%"), 
        cex.axis=axis.size) 

axis(2, at=c(0, 0.05, 0.10, 0.15, 0.20), 
        labels=c("0%", "5%", "10%", "15%", "20%"), 
        las=2, cex.axis=axis.size ) 



#####################################
#####################################
##                                 ## 
##  PANEL 3                        ##
##  Measured prevalence            ##
##  multipe antigens               ##
##                                 ##
#####################################
##################################### 




line_seq_x <- c(0.0, 0.01, 0.05, 0.1)
line_seq_y <- c(0.0, 0.05, 0.1, 0.15, 0.2)



plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0.00,0.1), ylim=c(0,0.2), 
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3.4, 
cex=lab.size, 
text="1 - specificity (= false positive rate)")

mtext(side = 2, line = 4.5, 
cex=lab.size, 
text="seroprevalence")

for(i in 1:length(line_seq_y))
{
	points(x=c(1e-5,1), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}

pos_index <- which( 1 - RF_multi_ant_ROC$specificities > -1 )


points( x=(1 - RF_multi_ant_ROC$specificities[pos_index]), y=M_multi_ant_5prev_quant[pos_index,1], 
	  type='l', lwd=line.size, col="aquamarine4" )


polygon(x=c( (1 - RF_multi_ant_ROC$specificities[pos_index]), rev( (1 - RF_multi_ant_ROC$specificities[pos_index])) ),
	  y=c( M_multi_ant_5prev_quant[pos_index,2], rev(M_multi_ant_5prev_quant[pos_index,3]) ), 
     	  col=rgb(69/256,139/256,116/256,0.2), border=NA)



points(x=c(1e-5,1), y=rep(0.05,2), 
type='l', col="black", lty="dashed", lwd=line.size)

text(x=0.07, y=0.035,
labels="true seroprevalence = 5%", cex=1.75)

axis(1, at=c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1), 
        labels=c("0%", "2%", "4%", "6%", "8%", "10%"), 
        cex.axis=axis.size) 

axis(2, at=c(0, 0.05, 0.10, 0.15, 0.20), 
        labels=c("0%", "5%", "10%", "15%", "20%"), 
        las=2, cex.axis=axis.size ) 




############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))


plot.new()
title( "(C) Adjusted seroprevalence", 
        cex.main=3.0, line=-2)



#####################################
#####################################
##                                 ## 
##  PANEL 4                        ##
##  Adjusted prevalence            ##
##  single antigen                 ##
##                                 ##
#####################################
##################################### 


par(mar=c(5,8,1.5,1.5))
par(mgp=c(2.5, 1.25,0))


line_seq_x <- c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1)
line_seq_y <- c(0.0, 0.05, 0.1, 0.15, 0.2)


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0.00,0.1), ylim=c(0,0.2), 
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3.4, 
cex=lab.size, 
text="1 - specificity (= false positive rate)")

mtext(side = 2, line = 4.5, 
cex=lab.size, 
text="seroprevalence")

for(i in 1:length(line_seq_y))
{
	points(x=c(1e-5,1), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}

pos_index <- which( 1 - spec_1ant_binom[,1] > 0 )

pos_index <- c(pos_index, max(pos_index)+1)

points( x=1 - spec_1ant_binom[pos_index,1], y=A_one_ant_5prev_quant[pos_index,1], 
	  type='l', lwd=line.size, col="orangered" )

polygon(x=c( 1 - spec_1ant_binom[pos_index,1], rev( 1 - spec_1ant_binom[pos_index,1]) ),
	  y=c( A_one_ant_5prev_quant[pos_index,2], rev(A_one_ant_5prev_quant[pos_index,3]) ), 
     	  col=rgb(190/256,90/256,90/256,0.2), border=NA)


points(x=c(1e-5,1), y=rep(0.05,2), 
type='l', col="black", lty="dashed", lwd=line.size)

text(x=0.07, y=0.03,
labels="true seroprevalence = 5%", cex=1.75)

axis(1, at=c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1), 
        labels=c("0%", "2%", "4%", "6%", "8%", "10%"), 
        cex.axis=axis.size) 

axis(2, at=c(0, 0.05, 0.10, 0.15, 0.20), 
        labels=c("0%", "5%", "10%", "15%", "20%"), 
        las=2, cex.axis=axis.size ) 
 



#####################################
#####################################
##                                 ## 
##  PANEL 5                        ##
##  Measured prevalence            ##
##  multipe antigens               ##
##                                 ##
#####################################
##################################### 

line_seq_x <- c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1)
line_seq_y <- c(0.0, 0.05, 0.1, 0.15, 0.2)

plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0.00,0.1), ylim=c(0,0.2), #log="x",
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3.4, 
cex=lab.size, 
text="1 - specificity (= false positive rate)")

mtext(side = 2, line = 4.5, 
cex=lab.size, 
text="seroprevalence")

for(i in 1:length(line_seq_y))
{
	points(x=c(1e-5,1), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}

pos_index <- which( 1 - RF_multi_ant_ROC$specificities[] > 0 )
pos_index <- c(pos_index, max(pos_index)+1)

points( x=1 - RF_multi_ant_ROC$specificities[pos_index], y=A_multi_ant_5prev_quant[pos_index,1], 
	  type='l', lwd=line.size, col="aquamarine4" )


polygon(x=c( 1 - RF_multi_ant_ROC$specificities[pos_index], rev( 1 - RF_multi_ant_ROC$specificities[pos_index]) ),
	  y=c( A_multi_ant_5prev_quant[pos_index,2], rev(A_multi_ant_5prev_quant[pos_index,3]) ), 
     	  col=rgb(69/256,139/256,116/256,0.2), border=NA)



points(x=c(1e-5,1), y=rep(0.05,2), 
type='l', col="black", lty="dashed", lwd=line.size)

text(x=0.07, y=0.075,
labels="true seroprevalence = 5%", cex=1.75)

axis(1, at=c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1), 
        labels=c("0%", "2%", "4%", "6%", "8%", "10%"), 
        cex.axis=axis.size) 

axis(2, at=c(0, 0.05, 0.10, 0.15, 0.20), 
        labels=c("0%", "5%", "10%", "15%", "20%"), 
        las=2, cex.axis=axis.size ) 







############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))


plot.new()
title( "(D) Optimal sensitivity and specificity", 
        cex.main=3.0, line=-2)



###########################################
###########################################
##                                       ## 
##  PANEL 6                              ##
##  Optimal sensitivity and specificity  ##
##  single antigen                       ##
##                                       ##
###########################################
########################################### 



par(mar=c(5,8,1.5,1.5))
par(mgp=c(2.5, 1.25,0))



line_seq_x <- c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1)
line_seq_y <- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)


plot(x=100, y=100, 
xlim=c(0.001,1), ylim=c(0,1), log="x",
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)



for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=c(1e-5,1), y=rep(line_seq_y[i]^PP,2), type='l', col="grey", lty="dashed")
}


mtext(side = 1, line = 3.4, 
cex=lab.size, 
text="true seroprevalence")

mtext(side = 2, line = 4.5, 
cex=lab.size, 
text="sensitivity \n & specificity")

points( x=T_prev_seq, y=T_prev_1ant[,1]^PP, 
    		  type='l', lwd=line.size, col="dodgerblue" )

points( x=T_prev_seq, y=T_prev_1ant[,2]^PP, 
    		  type='l', lwd=line.size, col="black" )


text(x=0.005, y=0.99^PP,
labels="specificity", 
col="black", cex=2.5)

text(x=0.005, y=0.85^PP,
labels="sensitivity", 
col="dodgerblue", cex=2.5)

axis(1, at=c(0.001, 0.01, 0.1,  1), 
        labels=c("0.1%", "1%", "10%", "100%"), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 


###########################################
###########################################
##                                       ## 
##  PANEL 7                              ##
##  Optimal sensitivity and specificity  ##
##  multiple antigens                    ##
##                                       ##
###########################################
########################################### 

line_seq_x <- c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1)
line_seq_y <- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)


plot(x=100, y=100, 
xlim=c(0.001,1), ylim=c(0,1), log="x",
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=c(1e-5,1), y=rep(line_seq_y[i]^PP,2), type='l', col="grey", lty="dashed")
}


mtext(side = 1, line = 3.4, 
cex=lab.size, 
text="true seroprevalence")

mtext(side = 2, line = 4.5, 
cex=lab.size, 
text="sensitivity \n & specificity")



points( x=T_prev_seq, y=T_prev_multi_ant[,1]^PP, 
    		  type='l', lwd=2, col="dodgerblue" )

points( x=T_prev_seq, y=T_prev_multi_ant[,2]^PP, 
    		  type='l', lwd=2, col="black" )

text(x=0.005, y=0.99^PP,
labels="specificity", 
col="black", cex=2.5)

text(x=0.005, y=0.88^PP,
labels="sensitivity", 
col="dodgerblue", cex=2.5)

axis(1, at=c(0.0001, 0.001, 0.01, 0.1,  1), 
        labels=c("0.01%", "0.1%", "1%", "10%", "100%"), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 



###########################################
###########################################
##                                       ## 
##  PANEL 8                              ##
##  Expected relative error              ##
##  single antigen                       ##
##                                       ##
###########################################
########################################### 

par(mar=c(5,8,2.5,1.5))
par(mgp=c(2.5, 1.25,0))

line_seq_x <- c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1)
line_seq_y <- c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12)


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0.001,1), ylim=c(0.0,0.12), log="x",
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="(E) Expected relative error",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_x))
{
	points(x=rep(line_seq_x[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=c(1e-5,1), y=rep(line_seq_y[i],2), type='l', col="grey", lty="dashed")
}


mtext(side = 1, line = 3.4, 
cex=lab.size, 
text="true seroprevalence")

mtext(side = 2, line = 4.5, 
cex=lab.size, 
text="expected error")


points( x=T_prev_seq, y=T_prev_1ant[,3], 
    		  type='l', lwd=line.size, col="orangered" )

points( x=T_prev_seq, y=T_prev_multi_ant[,3], 
    		  type='l', lwd=line.size, col="aquamarine4" )

axis(1, at=c(0.001, 0.01, 0.1,  1), 
        labels=c("0.1%", "1%", "10%", "100%"), 
        cex.axis=axis.size)

axis(2, at=c(0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12), 
        labels=c("0%", "2%", "4%", "6%", "8%", "10%", "12%"), 
        las=2, cex.axis=axis.size) 







dev.off()




