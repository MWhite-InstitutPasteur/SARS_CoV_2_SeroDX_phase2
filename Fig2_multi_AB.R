library(binom)

AB_data <- read.csv("C:\\U\\CoronaVirus\\SeroSig_Phase2\\DATA\\SARS-CoV-2_serology.csv")

AB_data <- AB_data[-which(AB_data$days_post <= 10),]

N_positive <- length(which(AB_data$status == "positive"))
N_negative <- length(which(AB_data$status == "negative"))



###################################################
###################################################
##                                               ##
##   ####  #   ## #####     ####  #   ## ######  ##
##  ##  ## ##  ## ##       ##  ## ##  ##   ##    ##
##  ##  ## ### ## ####     ###### ### ##   ##    ##
##  ##  ## ## ### ##       ##  ## ## ###   ##    ## 
##   ####  ##  ## #####    ##  ## ##  ##   ##    ##
##                                               ##
###################################################
###################################################


N_ant <- 13
N_SS <- 50000

dil_cut <- exp(seq(from=log(1e-5), to=log(0.3), length=N_SS))
MFI_cut <- exp(seq(from=log(1), to=log(100000), length=N_SS))


sens_mat <- matrix(NA, ncol=2*N_ant, nrow=N_SS)
spec_mat <- matrix(NA, ncol=2*N_ant, nrow=N_SS)

colnames(sens_mat) <- c( "Stri_IPP_IgG_dil","S1RBD_NA_IgG_dil", "RBD_IPP_IgG_dil",  "S1_NA_IgG_dil", "S2_NA_IgG_dil",
                         "NP_IPP_IgG_dil", "NP_NA_IgG_dil", "X229E_NP_NA_IgG_MFI", "NL63_NP_NA_IgG_MFI", "fluA_NA_IgG_MFI",
                         "ade40_NA_IgG_MFI", "rub_NA_IgG_MFI", "mumps_NA_IgG_MFI",
                         "Stri_IPP_IgM_dil", "S1RBD_NA_IgM_dil", "RBD_IPP_IgM_dil",  "S1_NA_IgM_dil", "S2_NA_IgM_dil",
                         "NP_IPP_IgM_dil", "NP_NA_IgM_dil", "X229E_NP_NA_IgM_MFI", "NL63_NP_NA_IgM_MFI", "fluA_NA_IgM_MFI",    
                         "ade40_NA_IgM_MFI", "rub_NA_IgM_MFI", "mumps_NA_IgM_MFI" ) 

colnames(spec_mat) <- c( "Stri_IPP_IgG_dil", "S1RBD_NA_IgG_dil", "RBD_IPP_IgG_dil", "S1_NA_IgG_dil", "S2_NA_IgG_dil",
                         "NP_IPP_IgG_dil", "NP_NA_IgG_dil", "X229E_NP_NA_IgG_MFI", "NL63_NP_NA_IgG_MFI", "fluA_NA_IgG_MFI",
                         "ade40_NA_IgG_MFI", "rub_NA_IgG_MFI", "mumps_NA_IgG_MFI",
                         "Stri_IPP_IgM_dil", "S1RBD_NA_IgM_dil", "RBD_IPP_IgM_dil", "S1_NA_IgM_dil", "S2_NA_IgM_dil",
                         "NP_IPP_IgM_dil", "NP_NA_IgM_dil", "X229E_NP_NA_IgM_MFI", "NL63_NP_NA_IgM_MFI", "fluA_NA_IgM_MFI",    
                         "ade40_NA_IgM_MFI", "rub_NA_IgM_MFI", "mumps_NA_IgM_MFI" )

for(i in 1:N_SS)
{
	sens_mat[i,1] <- length(which( AB_data$Stri_IPP_IgG_dil[which(AB_data$status == "positive")] > dil_cut[i] ))/length(which(is.na(AB_data$Stri_IPP_IgG_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,1] <- length(which( AB_data$Stri_IPP_IgG_dil[which(AB_data$status == "negative")] < dil_cut[i] ))/length(which(is.na(AB_data$Stri_IPP_IgG_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,2] <- length(which( AB_data$S1RBD_NA_IgG_dil[which(AB_data$status == "positive")] > dil_cut[i] ))/length(which(is.na(AB_data$S1RBD_NA_IgG_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,2] <- length(which( AB_data$S1RBD_NA_IgG_dil[which(AB_data$status == "negative")] < dil_cut[i] ))/length(which(is.na(AB_data$S1RBD_NA_IgG_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,3] <- length(which( AB_data$RBD_IPP_IgG_dil[which(AB_data$status == "positive")] > dil_cut[i] ))/length(which(is.na(AB_data$RBD_IPP_IgG_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,3] <- length(which( AB_data$RBD_IPP_IgG_dil[which(AB_data$status == "negative")] < dil_cut[i] ))/length(which(is.na(AB_data$RBD_IPP_IgG_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,4] <- length(which( AB_data$S1_NA_IgG_dil[which(AB_data$status == "positive")] > dil_cut[i] ))/length(which(is.na(AB_data$S1_NA_IgG_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,4] <- length(which( AB_data$S1_NA_IgG_dil[which(AB_data$status == "negative")] < dil_cut[i] ))/length(which(is.na(AB_data$S1_NA_IgG_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,5] <- length(which( AB_data$S2_NA_IgG_dil[which(AB_data$status == "positive")] > dil_cut[i] ))/length(which(is.na(AB_data$S2_NA_IgG_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,5] <- length(which( AB_data$S2_NA_IgG_dil[which(AB_data$status == "negative")] < dil_cut[i] ))/length(which(is.na(AB_data$S2_NA_IgG_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,6] <- length(which( AB_data$NP_IPP_IgG_dil[which(AB_data$status == "positive")] > dil_cut[i] ))/length(which(is.na(AB_data$NP_IPP_IgG_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,6] <- length(which( AB_data$NP_IPP_IgG_dil[which(AB_data$status == "negative")] < dil_cut[i] ))/length(which(is.na(AB_data$NP_IPP_IgG_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,7] <- length(which( AB_data$NP_NA_IgG_dil[which(AB_data$status == "positive")] > dil_cut[i] ))/length(which(is.na(AB_data$NP_NA_IgG_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,7] <- length(which( AB_data$NP_NA_IgG_dil[which(AB_data$status == "negative")] < dil_cut[i] ))/length(which(is.na(AB_data$NP_NA_IgG_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,8] <- length(which( AB_data$X229E_NP_NA_IgG_MFI[which(AB_data$status == "positive")] > MFI_cut[i] ))/length(which(is.na(AB_data$X229E_NP_NA_IgG_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,8] <- length(which( AB_data$X229E_NP_NA_IgG_MFI[which(AB_data$status == "negative")] < MFI_cut[i] ))/length(which(is.na(AB_data$X229E_NP_NA_IgG_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,9] <- length(which( AB_data$NL63_NP_NA_IgG_MFI[which(AB_data$status == "positive")] > MFI_cut[i] ))/length(which(is.na(AB_data$NL63_NP_NA_IgG_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,9] <- length(which( AB_data$NL63_NP_NA_IgG_MFI[which(AB_data$status == "negative")] < MFI_cut[i] ))/length(which(is.na(AB_data$NL63_NP_NA_IgG_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,10] <- length(which( AB_data$fluA_NA_IgG_MFI[which(AB_data$status == "positive")] > MFI_cut[i] ))/length(which(is.na(AB_data$fluA_NA_IgG_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,10] <- length(which( AB_data$fluA_NA_IgG_MFI[which(AB_data$status == "negative")] < MFI_cut[i] ))/length(which(is.na(AB_data$fluA_NA_IgG_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,11] <- length(which( AB_data$ade40_NA_IgG_MFI[which(AB_data$status == "positive")] > MFI_cut[i] ))/length(which(is.na(AB_data$ade40_NA_IgG_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,11] <- length(which( AB_data$ade40_NA_IgG_MFI[which(AB_data$status == "negative")] < MFI_cut[i] ))/length(which(is.na(AB_data$ade40_NA_IgG_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,12] <- length(which( AB_data$rub_NA_IgG_MFI[which(AB_data$status == "positive")] > MFI_cut[i] ))/length(which(is.na(AB_data$rub_NA_IgG_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,12] <- length(which( AB_data$rub_NA_IgG_MFI[which(AB_data$status == "negative")] < MFI_cut[i] ))/length(which(is.na(AB_data$rub_NA_IgG_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,13] <- length(which( AB_data$mumps_NA_IgG_MFI[which(AB_data$status == "positive")] > MFI_cut[i] ))/length(which(is.na(AB_data$mumps_NA_IgG_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,13] <- length(which( AB_data$mumps_NA_IgG_MFI[which(AB_data$status == "negative")] < MFI_cut[i] ))/length(which(is.na(AB_data$mumps_NA_IgG_dil[which(AB_data$status == "negative")])==FALSE))
}

for(i in 1:N_SS)
{
	sens_mat[i,14] <- length(which( AB_data$Stri_IPP_IgM_dil[which(AB_data$status == "positive")] > dil_cut[i] ))/length(which(is.na(AB_data$Stri_IPP_IgM_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,14] <- length(which( AB_data$Stri_IPP_IgM_dil[which(AB_data$status == "negative")] < dil_cut[i] ))/length(which(is.na(AB_data$Stri_IPP_IgM_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,15] <- length(which( AB_data$S1RBD_NA_IgM_dil[which(AB_data$status == "positive")] > dil_cut[i] ))/length(which(is.na(AB_data$S1RBD_NA_IgM_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,15] <- length(which( AB_data$S1RBD_NA_IgM_dil[which(AB_data$status == "negative")] < dil_cut[i] ))/length(which(is.na(AB_data$S1RBD_NA_IgM_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,16] <- length(which( AB_data$RBD_IPP_IgM_dil[which(AB_data$status == "positive")] > dil_cut[i] ))/length(which(is.na(AB_data$RBD_IPP_IgM_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,16] <- length(which( AB_data$RBD_IPP_IgM_dil[which(AB_data$status == "negative")] < dil_cut[i] ))/length(which(is.na(AB_data$RBD_IPP_IgM_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,17] <- length(which( AB_data$S1_NA_IgM_dil[which(AB_data$status == "positive")] > dil_cut[i] ))/length(which(is.na(AB_data$S1_NA_IgM_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,17] <- length(which( AB_data$S1_NA_IgM_dil[which(AB_data$status == "negative")] < dil_cut[i] ))/length(which(is.na(AB_data$S1_NA_IgM_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,18] <- length(which( AB_data$S2_NA_IgM_dil[which(AB_data$status == "positive")] > dil_cut[i] ))/length(which(is.na(AB_data$S2_NA_IgM_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,18] <- length(which( AB_data$S2_NA_IgM_dil[which(AB_data$status == "negative")] < dil_cut[i] ))/length(which(is.na(AB_data$S2_NA_IgM_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,19] <- length(which( AB_data$NP_IPP_IgM_dil[which(AB_data$status == "positive")] > dil_cut[i] ))/length(which(is.na(AB_data$NP_IPP_IgM_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,19] <- length(which( AB_data$NP_IPP_IgM_dil[which(AB_data$status == "negative")] < dil_cut[i] ))/length(which(is.na(AB_data$NP_IPP_IgM_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,20] <- length(which( AB_data$NP_NA_IgM_dil[which(AB_data$status == "positive")] > dil_cut[i] ))/length(which(is.na(AB_data$NP_NA_IgM_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,20] <- length(which( AB_data$NP_NA_IgM_dil[which(AB_data$status == "negative")] < dil_cut[i] ))/length(which(is.na(AB_data$NP_NA_IgM_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,21] <- length(which( AB_data$X229E_NP_NA_IgM_MFI[which(AB_data$status == "positive")] > MFI_cut[i] ))/length(which(is.na(AB_data$X229E_NP_NA_IgM_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,21] <- length(which( AB_data$X229E_NP_NA_IgM_MFI[which(AB_data$status == "negative")] < MFI_cut[i] ))/length(which(is.na(AB_data$X229E_NP_NA_IgM_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,22] <- length(which( AB_data$NL63_NP_NA_IgM_MFI[which(AB_data$status == "positive")] > MFI_cut[i] ))/length(which(is.na(AB_data$NL63_NP_NA_IgM_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,22] <- length(which( AB_data$NL63_NP_NA_IgM_MFI[which(AB_data$status == "negative")] < MFI_cut[i] ))/length(which(is.na(AB_data$NL63_NP_NA_IgM_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,23] <- length(which( AB_data$fluA_NA_IgM_MFI[which(AB_data$status == "positive")] > MFI_cut[i] ))/length(which(is.na(AB_data$fluA_NA_IgM_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,23] <- length(which( AB_data$fluA_NA_IgM_MFI[which(AB_data$status == "negative")] < MFI_cut[i] ))/length(which(is.na(AB_data$fluA_NA_IgM_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,24] <- length(which( AB_data$ade40_NA_IgM_MFI[which(AB_data$status == "positive")] > MFI_cut[i] ))/length(which(is.na(AB_data$ade40_NA_IgM_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,24] <- length(which( AB_data$ade40_NA_IgM_MFI[which(AB_data$status == "negative")] < MFI_cut[i] ))/length(which(is.na(AB_data$ade40_NA_IgM_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,25] <- length(which( AB_data$rub_NA_IgM_MFI[which(AB_data$status == "positive")] > MFI_cut[i] ))/length(which(is.na(AB_data$rub_NA_IgM_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,25] <- length(which( AB_data$rub_NA_IgM_MFI[which(AB_data$status == "negative")] < MFI_cut[i] ))/length(which(is.na(AB_data$rub_NA_IgM_dil[which(AB_data$status == "negative")])==FALSE))

	sens_mat[i,26] <- length(which( AB_data$mumps_NA_IgM_MFI[which(AB_data$status == "positive")] > MFI_cut[i] ))/length(which(is.na(AB_data$mumps_NA_IgM_dil[which(AB_data$status == "positive")])==FALSE))
	spec_mat[i,26] <- length(which( AB_data$mumps_NA_IgM_MFI[which(AB_data$status == "negative")] < MFI_cut[i] ))/length(which(is.na(AB_data$mumps_NA_IgM_dil[which(AB_data$status == "negative")])==FALSE))
}





######################################
## Calculate Area Under Curve (AUC)


AUC_one_ant <- rep(NA, 2*N_ant)

for(j in 1:(2*N_ant))
{
	AUC_one_ant[j] <- sum( (sens_mat[1:(N_SS-1),j] - sens_mat[2:N_SS,j])*
                             0.5*(spec_mat[1:(N_SS-1),j] + spec_mat[2:N_SS,j]) )
}





sens_target <- matrix(NA, nrow=3, ncol=2*N_ant)
colnames(sens_target) <- colnames(sens_mat)
rownames(sens_target) <- c("sens_99", "sens_spec", "spec_99")

spec_target <- matrix(NA, nrow=3, ncol=2*N_ant)
colnames(spec_target) <- colnames(spec_mat)
rownames(spec_target) <- c("sens_99", "sens_spec", "spec_99")

for(j in 1:(2*N_ant))
{
	sens_target[1,j] <- sens_mat[max(which(sens_mat[,j] > 0.99)),j]
	spec_target[1,j] <- spec_mat[max(which(sens_mat[,j] > 0.99)),j]

	sens_target[2,j] <- sens_mat[which.min(abs(sens_mat[,j] - spec_mat[,j])),j]
	spec_target[2,j] <- spec_mat[which.min(abs(sens_mat[,j] - spec_mat[,j])),j]

	sens_target[3,j] <- sens_mat[min(which(spec_mat[,j] > 0.99)),j]
	spec_target[3,j] <- spec_mat[min(which(spec_mat[,j] > 0.99)),j]
}




sens_target_lwr <- sens_target

for(i in 1:nrow(sens_target_lwr))
{
	for(j in 1:ncol(sens_target_lwr))
	{
		sens_target_lwr[i,j] <- binom.confint( sens_target_lwr[i,j]*N_positive, N_positive, method="wilson")[1,5]
		sens_target_lwr[i,j] <- round( 100*sens_target_lwr[i,j], 1)  
	}

}

sens_target_upr <- sens_target

for(i in 1:nrow(sens_target_upr))
{
	for(j in 1:ncol(sens_target_upr))
	{
		sens_target_upr[i,j] <- binom.confint( sens_target_upr[i,j]*N_positive, N_positive, method="wilson")[1,6]
		sens_target_upr[i,j] <- round( 100*sens_target_upr[i,j], 1)  
	}
}


spec_target_lwr <- spec_target

for(i in 1:nrow(spec_target_lwr))
{
	for(j in 1:ncol(sens_target_upr))
	{
		spec_target_lwr[i,j] <- binom.confint( spec_target_lwr[i,j]*N_negative, N_negative, method="wilson")[1,5]
		spec_target_lwr[i,j] <- round( 100*spec_target_lwr[i,j], 1)  
	}
}

spec_target_upr <- spec_target

for(i in 1:nrow(spec_target_upr))
{
	for(j in 1:ncol(sens_target_upr))
	{
		spec_target_upr[i,j] <- binom.confint( spec_target_upr[i,j]*N_negative, N_negative, method="wilson")[1,6]
		spec_target_upr[i,j] <- round( 100*spec_target_upr[i,j], 1)  
	}
}


sens_target <- round( 100*sens_target, 1)



spec_target <- round( 100*spec_target, 1)





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

set.seed(1234)

N_tree <- 50000


library(MASS)
library(ROCR)
library(randomForest)
library(pROC)

status <- AB_data$status

AB <- log(AB_data[,c(10:16,43:47,23:29,56:60)])


####################################
## Variable Importance
## (all antigens)

index_trim <- unique(which(is.na(AB)==TRUE, arr.ind=TRUE)[,1])
AB_trim <- AB[-index_trim,]
status_trim <- status[-index_trim]

RF_all = randomForest( status_trim ~ ., data=AB_trim, 
                                      importance=TRUE, ntree=N_tree)

rf.roc_all <- roc(status, RF_all$votes[,2])

varImpPlot( RF_all )





RF_summary <- function( antigen_list )
{
	AB_trim <- AB[,antigen_list]
	status_trim <- status

	index_trim <- unique(which(is.na(AB_trim)==TRUE, arr.ind=TRUE)[,1])
	if( length(index_trim) > 0 )
	{
		AB_trim <- AB_trim[-index_trim,]
		status_trim <- status[-index_trim]
	}


	RF_ant = randomForest( status_trim ~ ., data=AB_trim, 
                                      importance=TRUE, ntree=N_tree)

	RF_ant_roc <- roc(status_trim, RF_ant$votes[,2])

	c( RF_ant_roc$sensitivities[min(which(RF_ant_roc$specificities > 0.99))],
         RF_ant_roc$auc,
         RF_ant_roc$specificities[max(which(RF_ant_roc$sensitivities > 0.99))] )
}


####################################
## 2 antigens: Spike IgG & other antigens

RF_2ant_summary <- matrix(NA, nrow=23, ncol=3)
rownames(RF_2ant_summary) <- colnames(AB)[2:24]
colnames(RF_2ant_summary) <- c("high_spec", "auc", "high_sens")

for(i in 1:23)
{
	RF_2ant_summary[i,] <- RF_summary( c(1,1+i) )
}


## Choose RBD IgG (3)


RF_2ant = randomForest( status ~ ., data=AB[,c(1,3)], 
                                      importance=TRUE, ntree=N_tree)

RF_2ant_roc <- roc(status, RF_2ant$votes[,2])


####################################
## 3 antigens: Spike IgG, RBD IgG, & other antigens

RF_3ant_summary <- matrix(NA, nrow=22, ncol=3)
rownames(RF_3ant_summary) <- colnames( AB)[c(2,4:24)]
colnames(RF_3ant_summary) <- c("high_spec", "auc", "high_sens")

test_seq <- c(2,4:24)

for(i in 1:22)
{
	RF_3ant_summary[i,] <- RF_summary( c(1,3,test_seq[i]) )
}

## Choose NP IgG (6)


AB_trim <- AB[,c(1,3,6)]
index_trim <- unique(which(is.na(AB_trim)==TRUE, arr.ind=TRUE)[,1])
AB_trim <- AB_trim[-index_trim,]
status_trim <- status[-index_trim]

RF_3ant = randomForest( status_trim ~ ., data=AB_trim, 
                                      importance=TRUE, ntree=N_tree)


RF_3ant_roc <- roc(status_trim, RF_3ant$votes[,2])


####################################
## 4 antigens: Spike IgG, RBD IgG, NP IgG & other antigens

RF_4ant_summary <- matrix(NA, nrow=21, ncol=3)
rownames(RF_4ant_summary) <- colnames(AB)[c(2,4,5,7:24)]
colnames(RF_4ant_summary) <- c("high_spec", "auc", "high_sens")

test_seq <- c(2,4,5,7:24)

for(i in 1:21)
{
	RF_4ant_summary[i,] <- RF_summary( c(1,3,6,test_seq[i]) )
}

## Choose S2 IgG (5)


AB_trim <- AB[,c(1,3,6,5)]
index_trim <- unique(which(is.na(AB_trim)==TRUE, arr.ind=TRUE)[,1])
AB_trim <- AB_trim[-index_trim,]
status_trim <- status[-index_trim]

RF_4ant = randomForest( status_trim ~ ., data=AB_trim, 
                                      importance=TRUE, ntree=N_tree)

RF_4ant_roc <- roc(status_trim, RF_4ant$votes[,2])


####################################
## 5 antigens: Spike IgG, RBD IgG, NP IgG, S2 IgG & other antigens

RF_5ant_summary <- matrix(NA, nrow=20, ncol=3)
rownames(RF_5ant_summary) <- colnames( AB)[c(2,4,7:24)]
colnames(RF_5ant_summary) <- c("high_spec", "auc", "high_sens")

test_seq <- c(2,4,7:24)

for(i in 1:20)
{
	RF_5ant_summary[i,] <- RF_summary( c(1,3,6,5,test_seq[i]) )
}

## Choose S1RBD IgM (14)



AB_trim <- AB[,c(1,3,6,5,14)]
index_trim <- unique(which(is.na(AB_trim)==TRUE, arr.ind=TRUE)[,1])
AB_trim <- AB_trim[-index_trim,]
status_trim <- status[-index_trim]

RF_5ant = randomForest( status_trim ~ ., data=AB_trim, 
                                      importance=TRUE, ntree=N_tree)


RF_5ant_roc <- roc(status_trim, RF_5ant$votes[,2])


####################################
## 6 antigens: Spike IgG, RBD IgG, NP IgG, S1 IgG, S1RBD IgM & other antigens

RF_6ant_summary <- matrix(NA, nrow=19, ncol=3)
rownames(RF_6ant_summary) <- colnames( AB)[c(2,4,7:13,15:24)]
colnames(RF_6ant_summary) <- c("high_spec", "auc", "high_sens")

test_seq <- c(2,4,7:13,15:24)

for(i in 1:19)
{
	RF_6ant_summary[i,] <- RF_summary( c(1,3,6,5,14,test_seq[i]) )
}

## Choose NP IgM (18)


AB_trim <- AB[,c(1,3,6,5,14,18)]
index_trim <- unique(which(is.na(AB_trim)==TRUE, arr.ind=TRUE)[,1])
AB_trim <- AB_trim[-index_trim,]
status_trim <- status[-index_trim]

RF_6ant = randomForest( status_trim ~ ., data=AB_trim, 
                                      importance=TRUE, ntree=N_tree)



RF_6ant_roc <- roc(status_trim, RF_6ant$votes[,2])





####################################
##

Sens_Combo <- matrix(NA, nrow=6, ncol=3)

Sens_Combo[1,] <- c( spec_target[1,1], spec_target_lwr[1,1], spec_target_upr[1,1] )/100

Sens_Combo[2,1] <- RF_2ant_summary[2,1]

Sens_Combo[3,1] <- RF_3ant_summary[4,1]

Sens_Combo[4,1] <- RF_4ant_summary[3,1]

Sens_Combo[5,1] <- RF_5ant_summary[10,1]

Sens_Combo[6,1] <- RF_6ant_summary[13,1]

for(k in 2:6)
{
	Sens_Combo[k,2] <- binom.confint( Sens_Combo[k,1]*N_positive, N_positive, method="wilson")[1,5]
	Sens_Combo[k,3] <- binom.confint( Sens_Combo[k,1]*N_positive, N_positive, method="wilson")[1,6]
}

rownames(Sens_Combo) <- c("Stri IgG", "RBD IgG", "NP IgG", "S1 IgG",  "S1RBD IgM", "NP IgM")


Sens_Combo_col <- c("red", "orangered", "forestgreen", "yellow", "brown", "aquamarine4")


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


N_tree <- 10000

N_rep <- 1000

N_SS_cut <- 5000

SS_seq <- seq(from=0, to=1, length=N_SS_cut)




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
##                                     ##
##   1 antigens                        ##
##                                     ##
#########################################


Stri_IgG_pos <- AB_data$Stri_IPP_IgG_dil[which(AB_data$status == "positive")]
Stri_IgG_neg <- AB_data$Stri_IPP_IgG_dil[which(AB_data$status == "negative")]


sens_1ant_xval_mat <- matrix(NA, nrow=N_SS, ncol=N_rep)
spec_1ant_xval_mat <- matrix(NA, nrow=N_SS, ncol=N_rep)

for(n in 1:N_rep)
{
	Stri_IgG_pos_sub <- Stri_IgG_pos[sample( length(Stri_IgG_pos), (2/3)*length(Stri_IgG_pos) )]		
	Stri_IgG_neg_sub <- Stri_IgG_neg[sample( length(Stri_IgG_neg), (2/3)*length(Stri_IgG_neg) )]		

	for(i in 1:N_SS)
	{
		sens_1ant_xval_mat[i,n] <- length(which( Stri_IgG_pos_sub > dil_cut[i] ))/length(Stri_IgG_pos_sub)
 		spec_1ant_xval_mat[i,n] <- length(which( Stri_IgG_neg_sub < dil_cut[i] ))/length(Stri_IgG_neg_sub)
	}
}


sens_1ant_xval_quant <- matrix(NA, nrow=N_SS, ncol=3)
spec_1ant_xval_quant <- matrix(NA, nrow=N_SS, ncol=3)

colnames(sens_1ant_xval_quant) <- c( "med", "lwr", "upr" )
colnames(spec_1ant_xval_quant) <- c( "med", "lwr", "upr" )

for(i in 1:N_SS)
{
	sens_1ant_xval_quant[i,] <- quantile( sens_1ant_xval_mat[i,], prob=c(0.5, 0.025, 0.975) )
	spec_1ant_xval_quant[i,] <- quantile( spec_1ant_xval_mat[i,], prob=c(0.5, 0.025, 0.975) )
}

index_1 <- which.min(abs(dil_cut - quantile( Stri_IgG_neg, prob=0.99)))

sens_xval_1ant <- sens_1ant_xval_quant[index_1,]


#########################################
##                                     ##
##   2 antigens                        ##
##                                     ##
#########################################


#########################################
## Cross-validated algorithm (mixed)

RF_2_ant_xval_rep <- RF_multi_ant_SS_xval( AB, status, antigen_list=c(1,3), train_prop=2/3, N_rep=N_rep )



#########################################
## Cross-validated variation in sensitivity

RF_2_ant_xval_se_vary_sp_fix <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:N_SS_cut)
	{
		RF_2_ant_xval_se_vary_sp_fix[i,n] <- RF_2_ant_xval_rep[[n]][min(which(RF_2_ant_xval_rep[[n]][,2] > SS_seq[i])),1]
	}
}

RF_2_ant_xval_se_vary_sp_fix[N_SS_cut,] <- 0


RF_2_ant_xval_se_vary_sp_fix_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
colnames(RF_2_ant_xval_se_vary_sp_fix_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_cut)
{
	RF_2_ant_xval_se_vary_sp_fix_quant[i,] <- quantile( RF_2_ant_xval_se_vary_sp_fix[i,], prob=c(0.5, 0.025, 0.975) )
}


index_2 <- which.min(abs(SS_seq - 0.99))


sens_xval_2ant <- RF_2_ant_xval_se_vary_sp_fix_quant[index_2,]



#########################################
##                                     ##
##   3 antigens                        ##
##                                     ##
#########################################


#########################################
## Cross-validated algorithm (mixed)

RF_3_ant_xval_rep <- RF_multi_ant_SS_xval( AB, status, antigen_list=c(1,3,6), train_prop=2/3, N_rep=N_rep )



#########################################
## Cross-validated variation in sensitivity

RF_3_ant_xval_se_vary_sp_fix <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:N_SS_cut)
	{
		RF_3_ant_xval_se_vary_sp_fix[i,n] <- RF_3_ant_xval_rep[[n]][min(which(RF_3_ant_xval_rep[[n]][,2] > SS_seq[i])),1]
	}
}

RF_3_ant_xval_se_vary_sp_fix[N_SS_cut,] <- 0


RF_3_ant_xval_se_vary_sp_fix_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
colnames(RF_3_ant_xval_se_vary_sp_fix_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_cut)
{
	RF_3_ant_xval_se_vary_sp_fix_quant[i,] <- quantile( RF_3_ant_xval_se_vary_sp_fix[i,], prob=c(0.5, 0.025, 0.975) )
}

index_3 <- which.min(abs(SS_seq - 0.99))


sens_xval_3ant <- RF_3_ant_xval_se_vary_sp_fix_quant[index_3,]



#########################################
##                                     ##
##   4 antigens                        ##
##                                     ##
#########################################


#########################################
## Cross-validated algorithm (mixed)

RF_4_ant_xval_rep <- RF_multi_ant_SS_xval( AB, status, antigen_list=c(1,3,6,5), train_prop=2/3, N_rep=N_rep )



#########################################
## Cross-validated variation in sensitivity

RF_4_ant_xval_se_vary_sp_fix <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:N_SS_cut)
	{
		RF_4_ant_xval_se_vary_sp_fix[i,n] <- RF_4_ant_xval_rep[[n]][min(which(RF_4_ant_xval_rep[[n]][,2] > SS_seq[i])),1]
	}
}

RF_4_ant_xval_se_vary_sp_fix[N_SS_cut,] <- 0


RF_4_ant_xval_se_vary_sp_fix_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
colnames(RF_4_ant_xval_se_vary_sp_fix_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_cut)
{
	RF_4_ant_xval_se_vary_sp_fix_quant[i,] <- quantile( RF_4_ant_xval_se_vary_sp_fix[i,], prob=c(0.5, 0.025, 0.975) )
}


index_4 <- which.min(abs(SS_seq - 0.99))


sens_xval_4ant <- RF_4_ant_xval_se_vary_sp_fix_quant[index_4,]




#########################################
##                                     ##
##   5 antigens                        ##
##                                     ##
#########################################


#########################################
## Cross-validated algorithm (mixed)

RF_5_ant_xval_rep <- RF_multi_ant_SS_xval( AB, status, antigen_list=c(1,3,6,5,14), train_prop=2/3, N_rep=N_rep )



#########################################
## Cross-validated variation in sensitivity

RF_5_ant_xval_se_vary_sp_fix <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:N_SS_cut)
	{
		RF_5_ant_xval_se_vary_sp_fix[i,n] <- RF_5_ant_xval_rep[[n]][min(which(RF_5_ant_xval_rep[[n]][,2] > SS_seq[i])),1]
	}
}

RF_5_ant_xval_se_vary_sp_fix[N_SS_cut,] <- 0


RF_5_ant_xval_se_vary_sp_fix_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
colnames(RF_5_ant_xval_se_vary_sp_fix_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_cut)
{
	RF_5_ant_xval_se_vary_sp_fix_quant[i,] <- quantile( RF_5_ant_xval_se_vary_sp_fix[i,], prob=c(0.5, 0.025, 0.975) )
}


index_5 <- which.min(abs(SS_seq - 0.99))


sens_xval_5ant <- RF_5_ant_xval_se_vary_sp_fix_quant[index_5,]




#########################################
##                                     ##
##   6 antigens                        ##
##                                     ##
#########################################


#########################################
## Cross-validated algorithm (mixed)

RF_6_ant_xval_rep <- RF_multi_ant_SS_xval( AB, status, antigen_list=c(1,3,6,5,14,18), train_prop=2/3, N_rep=N_rep )



#########################################
## Cross-validated variation in sensitivity

RF_6_ant_xval_se_vary_sp_fix <- matrix(NA, nrow=N_SS_cut, ncol=N_rep)

for(n in 1:N_rep)
{
	for(i in 1:N_SS_cut)
	{
		RF_6_ant_xval_se_vary_sp_fix[i,n] <- RF_6_ant_xval_rep[[n]][min(which(RF_6_ant_xval_rep[[n]][,2] > SS_seq[i])),1]
	}
}

RF_6_ant_xval_se_vary_sp_fix[N_SS_cut,] <- 0


RF_6_ant_xval_se_vary_sp_fix_quant <- matrix(NA, nrow=N_SS_cut, ncol=3)
colnames(RF_6_ant_xval_se_vary_sp_fix_quant) <- c("med", "lwr", "upr")

for(i in 1:N_SS_cut)
{
	RF_6_ant_xval_se_vary_sp_fix_quant[i,] <- quantile( RF_6_ant_xval_se_vary_sp_fix[i,], prob=c(0.5, 0.025, 0.975) )
}



index_6 <- which.min(abs(SS_seq - 0.99))


sens_xval_6ant <- RF_6_ant_xval_se_vary_sp_fix_quant[index_6,]















sens_target_combi <- matrix(NA, nrow=5, ncol=3)
rownames(sens_target_combi) <- c("RBD", "NP", "S2",
                           "S1RBD_IgM", "NP_IgM")
colnames(sens_target_combi) <- c("spec_99", "sens_spec", "sens_99")

sens_target_combi_lwr <- sens_target_combi
sens_target_combi_upr <- sens_target_combi

spec_target_combi <- sens_target_combi
spec_target_combi_lwr <- sens_target_combi
spec_target_combi_upr <- sens_target_combi



for(i in 1:5)
{
	if( i == 1 )
	{
		RFX <- RF_2ant_roc
	}

	if( i == 2 )
	{
		RFX <- RF_3ant_roc
	}		

	if( i == 3 )
	{
		RFX <- RF_4ant_roc
	}		

	if( i == 4 )
	{
		RFX <- RF_5ant_roc
	}		

	if( i == 5 )
	{
		RFX <- RF_6ant_roc
	}		


	################################
	## High specificity target

	sens_target_combi[i,1]     <- binom.confint( RFX$sensitivities[min(which(RFX$specificities > 0.99))]*N_positive, N_positive, method="wilson")[1,4]
	sens_target_combi_lwr[i,1] <- binom.confint( RFX$sensitivities[min(which(RFX$specificities > 0.99))]*N_positive, N_positive, method="wilson")[1,5]
	sens_target_combi_upr[i,1] <- binom.confint( RFX$sensitivities[min(which(RFX$specificities > 0.99))]*N_positive, N_positive, method="wilson")[1,6]

	spec_target_combi[i,1]     <- binom.confint( RFX$specificities[min(which(RFX$specificities > 0.99))]*N_negative, N_negative, method="wilson")[1,4]
	spec_target_combi_lwr[i,1] <- binom.confint( RFX$specificities[min(which(RFX$specificities > 0.99))]*N_negative, N_negative, method="wilson")[1,5]
	spec_target_combi_upr[i,1] <- binom.confint( RFX$specificities[min(which(RFX$specificities > 0.99))]*N_negative, N_negative, method="wilson")[1,6]


	################################
	## Balanced sensitivty and specificity target

	sens_target_combi[i,2]     <- binom.confint( RFX$sensitivities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_positive, N_positive, method="wilson")[1,4]
	sens_target_combi_lwr[i,2] <- binom.confint( RFX$sensitivities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_positive, N_positive, method="wilson")[1,5]
	sens_target_combi_upr[i,2] <- binom.confint( RFX$sensitivities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_positive, N_positive, method="wilson")[1,6]

	spec_target_combi[i,2]     <- binom.confint( RFX$specificities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_negative, N_negative, method="wilson")[1,4]
	spec_target_combi_lwr[i,2] <- binom.confint( RFX$specificities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_negative, N_negative, method="wilson")[1,5]
	spec_target_combi_upr[i,2] <- binom.confint( RFX$specificities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_negative, N_negative, method="wilson")[1,6]


	################################
	## High sensitivity target

	sens_target_combi[i,3]     <- binom.confint( RFX$sensitivities[max(which(RFX$sensitivities > 0.99))]*N_positive, N_positive, method="wilson")[1,4]
	sens_target_combi_lwr[i,3] <- binom.confint( RFX$sensitivities[max(which(RFX$sensitivities > 0.99))]*N_positive, N_positive, method="wilson")[1,5]
	sens_target_combi_upr[i,3] <- binom.confint( RFX$sensitivities[max(which(RFX$sensitivities > 0.99))]*N_positive, N_positive, method="wilson")[1,6]

	spec_target_combi[i,3]     <- binom.confint( RFX$specificities[max(which(RFX$sensitivities > 0.99))]*N_negative, N_negative, method="wilson")[1,4]
	spec_target_combi_lwr[i,3] <- binom.confint( RFX$specificities[max(which(RFX$sensitivities > 0.99))]*N_negative, N_negative, method="wilson")[1,5]
	spec_target_combi_upr[i,3] <- binom.confint( RFX$specificities[max(which(RFX$sensitivities > 0.99))]*N_negative, N_negative, method="wilson")[1,6]

}



spec_target_combi     <- round(100*spec_target_combi,1)
spec_target_combi_lwr <- round(100*spec_target_combi_lwr,1)
spec_target_combi_upr <- round(100*spec_target_combi_upr,1)

sens_target_combi     <- round(100*sens_target_combi,1)
sens_target_combi_lwr <- round(100*sens_target_combi_lwr,1)
sens_target_combi_upr <- round(100*sens_target_combi_upr,1)



save.image("Fig2_multi_AB.RData")



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

line_seq_dil_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01)
line_seq_dil_y <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01)

line_seq_MFI_x <- c(10, 30, 100, 300, 1000, 3000, 10000, 30000)
line_seq_MFI_y <- c(10, 30, 100, 300, 1000, 3000, 10000, 30000)



tiff(file="Fig2_multi_AB.tif", width=40, height=30, units="cm", res=500)

lay.mat <- rbind( c( 1, 1, 1, 1, 1, 1 ),
                  c( 2, 2, 3, 3, 4, 4 ),
                  c( 5, 5, 6, 6, 7, 7 ),
                  c( 8, 8, 8, 8, 8, 8 ),
                  c( 9, 9, 9,10,10,10 ) )
layout(lay.mat, heights=c(1.5,10,10,1.5,10), widths=c(1,1,1,1,1,1))
layout.show(10)



############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))


plot.new()
title( "(A) Pairwise serological signatures", 
        cex.main=3.0, line=-2)





par(mar=c(5,7.5,3,1))
par(mgp=c(3, 1, 0))

point.size = 0.75
lab.size   = 2
axis.size  = 1.5
main.size  = 2.5


###################################
## Panel 1: Spike_IPP IgG vs S1RBD_NA IgG

plot( x = AB_data$RBD_IPP_IgG_dil[which(AB_data$site=="Bichat")], 
      y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="Bichat")],
      xlim=c(1e-5, 0.03), ylim=c(1e-5, 0.03), log="xy",
      pch=19, col="yellowgreen",
	yaxt='n', xaxt='n', bty='n', 
	xlab=expression(paste( "anti-RBDv2 IgG dilution", sep="" )), 
	ylab="", 
	main=expression(paste( "anti-Stri", " IgG vs anti-RBDv2 IgG", sep="" )),
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_dil_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_dil_x[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_dil_y))
{
	points(x=rep(line_seq_dil_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points(x=rep(1.95e-5,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")
points(x=rep(0.02,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")

points(x=c(1e-10,1e10), y=rep(1.95e-5,2),  type='l', col="black", lty="dashed")
points(x=c(1e-10,1e10), y=rep(0.02,2), type='l', col="black", lty="dashed")

mtext(side = 2, line = 4.5, 
cex=axis.size, 
text=expression(paste( "anti-Stri", " IgG dilution", sep="" )))

points( x = AB_data$RBD_IPP_IgG_dil[which(AB_data$site=="Strasbourg")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="Strasbourg")],
        pch=19, col="sienna1" )

points( x = AB_data$RBD_IPP_IgG_dil[which(AB_data$site=="Cochin")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="Cochin")],
        pch=19, col="darkmagenta" )

points( x = AB_data$RBD_IPP_IgG_dil[which(AB_data$site=="EFS")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="EFS")],
        pch=19, col="royalblue" )

points( x = AB_data$RBD_IPP_IgG_dil[which(AB_data$site=="TRC")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="TRC")],
        pch=19, col="cornflowerblue" )

points( x = AB_data$RBD_IPP_IgG_dil[which(AB_data$site=="PNC")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="PNC")],
        pch=19, col="dodgerblue" )

axis(1, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        cex.axis=axis.size )

axis(2, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        las=2, cex.axis=axis.size )


###################################
## Panel 2: Spike_IPP IgG vs NP_IPP IgG

plot( x = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="Bichat")], 
      y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="Bichat")],
      xlim=c(1e-5, 0.03), ylim=c(1e-5, 0.03), log="xy",
      pch=19, col="yellowgreen",
	yaxt='n', xaxt='n', bty='n', 
	xlab=expression(paste( "anti-NPv1 IgG dilution", sep="" )), 
	ylab="", 
	main=expression(paste( "anti-Stri", " IgG vs anti-NPv1 IgG", sep="" )),
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_dil_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_dil_x[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_dil_y))
{
	points(x=rep(line_seq_dil_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points(x=rep(1.95e-5,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")
points(x=rep(0.02,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")

points(x=c(1e-10,1e10), y=rep(1.95e-5,2),  type='l', col="black", lty="dashed")
points(x=c(1e-10,1e10), y=rep(0.02,2), type='l', col="black", lty="dashed")

mtext(side = 2, line = 4.5, 
cex=axis.size, 
text=expression(paste( "anti-Stri", " IgG dilution", sep="" )))

points( x = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="Strasbourg")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="Strasbourg")],
        pch=19, col="sienna1" )

points( x = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="Cochin")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="Cochin")],
        pch=19, col="darkmagenta" )

points( x = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="EFS")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="EFS")],
        pch=19, col="royalblue" )

points( x = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="TRC")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="TRC")],
        pch=19, col="cornflowerblue" )

points( x = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="PNC")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="PNC")],
        pch=19, col="dodgerblue" )

axis(1, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        cex.axis=axis.size )

axis(2, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        las=2, cex.axis=axis.size )



###################################
## Panel 3: Spike_IPP IgG vs S1RBD_NA IgG

plot( x = AB_data$S1RBD_NA_IgG_dil[which(AB_data$site=="Bichat")], 
      y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="Bichat")],
      xlim=c(1e-5, 0.03), ylim=c(1e-5, 0.03), log="xy",
      pch=19, col="yellowgreen",
	yaxt='n', xaxt='n', bty='n', 
	xlab=expression(paste( "anti-RBDv1 IgG dilution", sep="" )), 
	ylab="", 
	main=expression(paste( "anti-Stri", " IgG vs anti-RBDv1 IgG", sep="" )),
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_dil_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_dil_x[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_dil_y))
{
	points(x=rep(line_seq_dil_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points(x=rep(1.95e-5,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")
points(x=rep(0.02,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")

points(x=c(1e-10,1e10), y=rep(1.95e-5,2),  type='l', col="black", lty="dashed")
points(x=c(1e-10,1e10), y=rep(0.02,2), type='l', col="black", lty="dashed")

mtext(side = 2, line = 4.5, 
cex=axis.size, 
text=expression(paste( "anti-Stri", " IgG dilution", sep="" )))

points( x = AB_data$S1RBD_NA_IgG_dil[which(AB_data$site=="Strasbourg")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="Strasbourg")],
        pch=19, col="sienna1" )

points( x = AB_data$S1RBD_NA_IgG_dil[which(AB_data$site=="Cochin")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="Cochin")],
        pch=19, col="darkmagenta" )

points( x = AB_data$S1RBD_NA_IgG_dil[which(AB_data$site=="EFS")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="EFS")],
        pch=19, col="royalblue" )

points( x = AB_data$S1RBD_NA_IgG_dil[which(AB_data$site=="TRC")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="TRC")],
        pch=19, col="cornflowerblue" )

points( x = AB_data$S1RBD_NA_IgG_dil[which(AB_data$site=="PNC")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="PNC")],
        pch=19, col="dodgerblue" )

axis(1, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        cex.axis=axis.size )

axis(2, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        las=2, cex.axis=axis.size )



###################################
## Panel 4: Spike_IPP IgG vs Spike_IPP IgM

plot( x = AB_data$Stri_IPP_IgM_dil[which(AB_data$site=="Bichat")], 
      y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="Bichat")],
      xlim=c(1e-5, 0.03), ylim=c(1e-5, 0.03), log="xy",
      pch=19, col="yellowgreen",
	yaxt='n', xaxt='n', bty='n', 
	xlab=expression(paste( "anti-Stri", " IgM dilution", sep="" )), 
	ylab="", 
	main=expression(paste( "anti-Stri", " IgG vs anti-Stri", " IgM", sep="" )),
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_dil_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_dil_x[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_dil_y))
{
	points(x=rep(line_seq_dil_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points(x=rep(1.95e-5,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")
points(x=rep(0.02,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")

points(x=c(1e-10,1e10), y=rep(1.95e-5,2),  type='l', col="black", lty="dashed")
points(x=c(1e-10,1e10), y=rep(0.02,2), type='l', col="black", lty="dashed")

mtext(side = 2, line = 4.5, 
cex=axis.size, 
text=expression(paste( "anti-Stri", " IgG dilution", sep="" )))

points( x = AB_data$Stri_IPP_IgM_dil[which(AB_data$site=="Strasbourg")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="Strasbourg")],
        pch=19, col="sienna1" )

points( x = AB_data$Stri_IPP_IgM_dil[which(AB_data$site=="Cochin")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="Cochin")],
        pch=19, col="darkmagenta" )

points( x = AB_data$Stri_IPP_IgM_dil[which(AB_data$site=="EFS")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="EFS")],
        pch=19, col="royalblue" )

points( x = AB_data$Stri_IPP_IgM_dil[which(AB_data$site=="TRC")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="TRC")],
        pch=19, col="cornflowerblue" )

points( x = AB_data$Stri_IPP_IgM_dil[which(AB_data$site=="PNC")], 
        y = AB_data$Stri_IPP_IgG_dil[which(AB_data$site=="PNC")],
        pch=19, col="dodgerblue" )


axis(1, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
      	  cex.axis=axis.size )

axis(2, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        las=2, cex.axis=axis.size )



###################################
## Panel 5: NP_IPP IgG vs NP_IPP IgM

plot( x = AB_data$NP_IPP_IgM_dil[which(AB_data$site=="Bichat")], 
      y = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="Bichat")],
      xlim=c(1e-5, 0.03), ylim=c(1e-5, 0.03), log="xy",
      pch=19, col="yellowgreen",
	yaxt='n', xaxt='n', bty='n', 
	xlab=expression(paste( "anti-NPv1 IgM MFI", sep="" )), 
	ylab="", 
	main=expression(paste( "anti-NPv1 IgG vs anti-NPv1 IgG", sep="" )),
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_dil_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_dil_x[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_dil_y))
{
	points(x=rep(line_seq_dil_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points(x=rep(1.95e-5,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")
points(x=rep(0.02,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")

points(x=c(1e-10,1e10), y=rep(1.95e-5,2),  type='l', col="black", lty="dashed")
points(x=c(1e-10,1e10), y=rep(0.02,2), type='l', col="black", lty="dashed")

mtext(side = 2, line = 4.5, 
cex=axis.size, 
text=expression(paste( "anti-NPv1 IgG dilution", sep="" )))

points( x = AB_data$NP_IPP_IgM_dil[which(AB_data$site=="Strasbourg")], 
        y = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="Strasbourg")],
        pch=19, col="sienna1" )

points( x = AB_data$NP_IPP_IgM_dil[which(AB_data$site=="Cochin")], 
        y = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="Cochin")],
        pch=19, col="darkmagenta" )

points( x = AB_data$NP_IPP_IgM_dil[which(AB_data$site=="EFS")], 
        y = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="EFS")],
        pch=19, col="royalblue" )

points( x = AB_data$NP_IPP_IgM_dil[which(AB_data$site=="TRC")], 
        y = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="TRC")],
        pch=19, col="cornflowerblue" )

points( x = AB_data$NP_IPP_IgM_dil[which(AB_data$site=="PNC")], 
        y = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="PNC")],
        pch=19, col="dodgerblue" )


axis(1, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
      	  cex.axis=axis.size )

axis(2, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        las=2, cex.axis=axis.size )


###################################
## Panel 6: NP_IPP IgG vs 229E_NP IgM

plot( x = AB_data$X229E_NP_NA_IgG_MFI[which(AB_data$site=="Bichat")], 
      y = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="Bichat")],
      xlim=c(10, 30000), ylim=c(1e-5, 0.03), log="xy",
      pch=19, col="yellowgreen",
	yaxt='n', xaxt='n', bty='n', 
	xlab=expression(paste( "anti-229E NP IgG MFI", sep="" )), 
	ylab="", 
	main=expression(paste( "anti-NPv1 IgG vs anti-229E NP IgG", sep="" )),
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_dil_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_dil_x[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_MFI_y))
{
	points(x=rep(line_seq_MFI_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points(x=rep(1.95e-5,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")
points(x=rep(0.02,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")

points(x=c(1e-10,1e10), y=rep(1.95e-5,2),  type='l', col="black", lty="dashed")
points(x=c(1e-10,1e10), y=rep(0.02,2), type='l', col="black", lty="dashed")

mtext(side = 2, line = 4.5, 
cex=axis.size, 
text=expression(paste( "anti-NPv1 IgG dilution", sep="" )))

points( x = AB_data$X229E_NP_NA_IgG_MFI[which(AB_data$site=="Strasbourg")], 
        y = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="Strasbourg")],
        pch=19, col="sienna1" )

points( x = AB_data$X229E_NP_NA_IgG_MFI[which(AB_data$site=="Cochin")], 
        y = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="Cochin")],
        pch=19, col="darkmagenta" )

points( x = AB_data$X229E_NP_NA_IgG_MFI[which(AB_data$site=="EFS")], 
        y = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="EFS")],
        pch=19, col="royalblue" )

points( x = AB_data$X229E_NP_NA_IgG_MFI[which(AB_data$site=="TRC")], 
        y = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="TRC")],
        pch=19, col="cornflowerblue" )

points( x = AB_data$X229E_NP_NA_IgG_MFI[which(AB_data$site=="PNC")], 
        y = AB_data$NP_IPP_IgG_dil[which(AB_data$site=="PNC")],
        pch=19, col="dodgerblue" )

axis(1, at=c(10, 30, 100, 300, 1000, 3000, 10000, 30000), 
      	  label=c("10", "", "100", "", "1000", "", "10000", ""),
      	  cex.axis=axis.size )

axis(2, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        las=2, cex.axis=axis.size )



###############	
##           ##
##  LEGEND   ##
##           ##
###############

par(mar = c(0,0,0,0))
plot.new()

legend(x='center', 
       legend = c("Bichat", "Strasbourg", "Cochin",  
                  "France negative", "Thai negative", "Peru negative"), 
       fill = c("yellowgreen", "sienna1", "darkmagenta", "dodgerblue", "royalblue", "cornflowerblue"), 
       border = c("yellowgreen", "sienna1", "darkmagenta", "dodgerblue", "royalblue", "cornflowerblue"), 
       ncol=6, cex=2, bty="n" )



#####################################
#####################################
##                                 ## 
##  PANELS 7                       ##
##  ROC analysis


PP <- 10


line_seq_x <- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)
line_seq_y <- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)



par(mar=c(5,7,3,1.5))
par(mgp=c(2.5, 1.25,0))


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="(B) Multi-antigen ROC analysis",
cex.lab=3, cex.axis=2, cex.main=2.5)

mtext(side = 1, line = 3.4, 
cex=2, 
text="1 - specificity")

mtext(side = 2, line = 4.5, 
cex=2, 
text="sensitivity")

for(i in 1:length(line_seq_y))
{
	points(x=c(0,1), y=rep(line_seq_y[i]^PP,2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_x))
{
	points(x=rep(1 - line_seq_x[i]^PP,2), y=c(0,1), type='l', col="grey", lty="dashed")
}


points( x=1-spec_mat[,1]^PP, y=sens_mat[,1]^PP, 
    	  type='S', lwd=2, col=Sens_Combo_col[1] )


points( x=1-RF_2ant_roc$specificities^PP, y=RF_2ant_roc$sensitivities^PP, 
    	  type='S', lwd=2, col=Sens_Combo_col[2] )

points( x=1-RF_3ant_roc$specificities^PP, y=RF_3ant_roc$sensitivities^PP, 
    	  type='S', lwd=2, col=Sens_Combo_col[3] )

points( x=1-RF_4ant_roc$specificities^PP, y=RF_4ant_roc$sensitivities^PP, 
    	  type='S', lwd=2, col=Sens_Combo_col[4] )

points( x=1-RF_5ant_roc$specificities^PP, y=RF_5ant_roc$sensitivities^PP, 
    	  type='S', lwd=2, col=Sens_Combo_col[5] )

points( x=1-RF_6ant_roc$specificities^PP, y=RF_6ant_roc$sensitivities^PP, 
    	  type='S', lwd=2, col=Sens_Combo_col[6] )


legend(x="bottomright", 
cex=2, 
bg="white", box.col="white",
fill = Sens_Combo_col,
border = Sens_Combo_col, 
legend = c("Stri IgG", "+ RBDv2 IgG", "+ NPv1 IgG", "+ S2 IgG", "+ RBDv1 IgM",  "+ NPv1 IgM") )  

axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=1.5) 

axis(1, at=1- c(0.0, 0.99^PP), 
        labels=rev(c( "1%",  "100%")), 
        cex.axis=1.5) 


axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=1.5 ) 




###################################
## Panel 8: High Spec Target


ant_cols <- c("red", "brown", "orangered", "gold", "yellow",
              "forestgreen", "limegreen", "mediumblue", "royalblue", 
              "burlywood", "darkkhaki", "cornsilk4",
              "red", "brown", "orangered", "gold", "yellow",
              "forestgreen", "limegreen", "mediumblue", "royalblue", 
              "burlywood", "darkkhaki", "cornsilk4")	
PP <- 10

gap <- 0.1

line_seq_y <- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)



par(mar=c(6,7,3,1.5))
par(mgp=c(2.5, 1.25,0))


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,6*(1+gap)), ylim=c(0,1.002),
xaxs='i',
# yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="(C) High specificity (99%) target",
cex.lab=3, cex.axis=2, cex.main=2.5)



mtext(side = 2, line = 4.5, 
cex=2, 
text="sensitivity")



for(k in 1:6)
{
	polygon( x=gap*(k-1) + c(k-1,k,k,k-1), c(0,0,Sens_Combo[k,1]^PP,Sens_Combo[k,1]^PP), col=Sens_Combo_col[k], border=NA)

#	arrows( x0=gap*(k-1) + k - 0.5, x1=gap*(k-1) + k - 0.5, y0=Sens_Combo[k,2]^PP, y1=Sens_Combo[k,3]^PP, 
#              code=3, col="black", angle=90, lwd=2, length=0.15)  
}

arrows( x0=gap*(1-1) + 1 - 0.5, x1=gap*(1-1) + 1 - 0.5, y0=sens_xval_1ant[2]^PP, y1=sens_xval_1ant[3]^PP, 
        code=3, col="black", angle=90, lwd=2, length=0.15)  

arrows( x0=gap*(2-1) + 2 - 0.5, x1=gap*(2-1) + 2 - 0.5, y0=sens_xval_2ant[2]^PP, y1=sens_xval_2ant[3]^PP, 
        code=3, col="black", angle=90, lwd=2, length=0.15)  

arrows( x0=gap*(3-1) + 3 - 0.5, x1=gap*(3-1) + 3 - 0.5, y0=sens_xval_3ant[2]^PP, y1=sens_xval_3ant[3]^PP, 
        code=3, col="black", angle=90, lwd=2, length=0.15)  

arrows( x0=gap*(4-1) + 4 - 0.5, x1=gap*(4-1) + 4 - 0.5, y0=sens_xval_4ant[2]^PP, y1=sens_xval_4ant[3]^PP, 
        code=3, col="black", angle=90, lwd=2, length=0.15)  

arrows( x0=gap*(5-1) + 5 - 0.5, x1=gap*(5-1) + 5 - 0.5, y0=sens_xval_5ant[2]^PP, y1=sens_xval_5ant[3]^PP, 
        code=3, col="black", angle=90, lwd=2, length=0.15)  

arrows( x0=gap*(6-1) + 6 - 0.5, x1=gap*(6-1) + 6 - 0.5, y0=sens_xval_6ant[2]^PP, y1=sens_xval_6ant[3]^PP, 
        code=3, col="black", angle=90, lwd=2, length=0.15)  

points( x=gap*(1-1) + 1 - 0.5, y=sens_xval_1ant[1]^PP, pch=19, cex=2, col="black") 

points( x=gap*(2-1) + 2 - 0.5, y=sens_xval_2ant[1]^PP, pch=19, cex=2, col="black") 

points( x=gap*(3-1) + 3 - 0.5, y=sens_xval_3ant[1]^PP, pch=19, cex=2, col="black") 

points( x=gap*(4-1) + 4 - 0.5, y=sens_xval_4ant[1]^PP, pch=19, cex=2, col="black") 

points( x=gap*(5-1) + 5 - 0.5, y=sens_xval_5ant[1]^PP, pch=19, cex=2, col="black") 

points( x=gap*(6-1) + 6 - 0.5, y=sens_xval_6ant[1]^PP, pch=19, cex=2, col="black") 


for(i in 1:length(line_seq_y))
{
	points(x=c(0,100), y=rep(line_seq_y[i]^PP,2), type='l', col="grey", lty="dashed")
}



axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=2 ) 



axis(1, at=seq(from=0.5, by=1+gap, length=6), 
           label=c("Stri IgG", "+ RBDv2 IgG", "+ NPv1 IgG", "+ S2 IgG", "+ RBDv1 IgM",  "+ NPv1 IgM"), 
            las=2, cex.axis=1)




dev.off()






