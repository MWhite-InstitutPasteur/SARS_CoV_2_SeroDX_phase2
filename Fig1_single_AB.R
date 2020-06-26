library(binom)

AB_data <- read.csv("C:\\U\\CoronaVirus\\SeroSig_Phase2\\DATA\\SARS-CoV-2_serology.csv")


AB_data <- AB_data[-which(AB_data$days_post <= 10),]

N_positive <- length(which(AB_data$status == "positive"))
N_negative <- length(which(AB_data$status == "negative"))




####################################
####################################
##                                ##
##  ONE AT A TIME                 ##
##                                ##
####################################
####################################


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



################################################################################
################################################################################
##                                                                            ##
##   ####   ####  #####  #####  ##### ##     ####  ###### ####  ####  #   ##  ##
##  ##  ## ##  ## ##  ## ##  ## ##    ##    ##  ##   ##    ##  ##  ## ##  ##  ##
##  ##     ##  ## #####  #####  ####  ##    ######   ##    ##  ##  ## ### ##  ##
##  ##  ## ##  ## ## ##  ## ##  ##    ##    ##  ##   ##    ##  ##  ## ## ###  ##
##   ####   ####  ##  ## ##  ## ##### ##### ##  ##   ##   ####  ####  ##  ##  ##
##                                                                            ##
################################################################################
################################################################################

library(fields)

AB_cor <- matrix(NA, nrow=24, ncol=24)

cor_seq <- c(10:16, 43:47, 23:29, 56:60)


for(i in 1:24)
{
	for(j in 1:24)
	{
		AB_i <- AB_data[,cor_seq[i]]
		AB_j <- AB_data[,cor_seq[j]]

		index_ij <- which( is.na(AB_i)==FALSE & is.na(AB_j)==FALSE )

		AB_i <- AB_i[index_ij]
		AB_j <- AB_j[index_ij]

		AB_cor[i,j] = cor( AB_i, AB_j, method="spearman")
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

plot_names <- c("anti-Stri IgG", "anti-RBDv1 IgG", "anti-RBDv2 IgG",  "anti-S1 IgG", "anti-S2 IgG", 
                "anti-NPv1 IgG", "anti-NPv2 IgG", "anti-229E_NP IgG",
                "anti-NL63_NP IgG", "anti-fluA IgG", "anti-ade40 IgG", "anti-rub IgG",
                "anti-Stri IgM", "anti-RBDv1 IgM", "anti-RBDv2 IgM", "anti-S1 IgM", "anti-S2 IgM",  
                "anti-NPv1 IgM", "anti-NPv2 IgM", "anti-229E_NP IgM","anti-NL63_NP IgM", 
                "anti-fluA IgM", "anti-ade40 IgM", "anti-rub IgM")
     

plot_names_short <- c("Stri IgG",    "RBDv1 IgG", "RBDv2 IgG", "S1 IgG", "S2 IgG", 
                      "NPv1 IgG",    "NPv2 IgG",  "229E_NP IgG",
                      "NL63_NP IgG", "fluA IgG",  "ade40 IgG", "rub IgG",
                      "Stri IgM",    "RBDv1 IgM", "RBDv2 IgM", "S1 IgM", "S2 IgM", 
                      "NPv1 IgM",    "NPv2 IgM",  "229E_NP IgM","NL63_NP IgM", 
                      "fluA IgM",    "ade40 IgM", "rub IgM")
     

ant_cols <- c("red", "brown", "orangered", "gold", "yellow",
              "forestgreen", "limegreen", "mediumblue", "royalblue", 
              "burlywood", "darkkhaki", "cornsilk4",
              "red", "brown", "orangered", "gold", "yellow",
              "forestgreen", "limegreen", "mediumblue", "royalblue", 
              "burlywood", "darkkhaki", "cornsilk4")	 




tiff(file="Fig1_single_AB.tif", width=40, height=30, units="cm", res=500)

lay.mat <- rbind( c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ),              
                  c( 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7 ),
                  c( 8, 8, 9, 9,10,10,11,11,12,12,13,13 ),
                  c(14,14,14,14,14,14,14,14,14,14,14,14 ),
                  c(15,15,16,16,17,17,18,18,19,19,20,20 ),
                  c(21,21,22,22,23,23,24,24,25,25,26,26 ),
                  c(27,27,27,28,28,28,29,29,29,30,30,30 ),
                  c(27,27,27,28,28,28,29,29,29,31,31,31 ) )
layout(lay.mat, heights=c(1.5,8,8,1.5,8,8,14,2), widths=c(1,1,1,1,1,1))
layout.show(31)


############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))


plot.new()
title( "(A) IgG antibody responses", 
        cex.main=3.0, line=-2)




#####################################
##                                  
##  PANELS 1:7                     



par(mar=c(3,5,3,1))
par(mgp=c(2.5, 0.7, 0))

point.size = 0.75
lab.size   = 1.5
axis.size  = 1.1
main.size  = 2



plot_seq <- 10:16

line_seq_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01)


for( k in 1:7 )
{
	boxplot( AB_data[which(AB_data$site=="Bichat"),plot_seq[k]],
		   AB_data[which(AB_data$site=="Strasbourg"),plot_seq[k]],
		   AB_data[which(AB_data$site=="Cochin"),plot_seq[k]],
	 	   AB_data[which(AB_data$site=="EFS"),plot_seq[k]],
		   AB_data[which(AB_data$site=="PNC"),plot_seq[k]], 
		   AB_data[which(AB_data$site=="TRC"),plot_seq[k]], 
		   log="y",
		   pch=19, yaxt='n', xaxt='n',
		   ylim=c(1e-5, 0.03),
		   col=ant_cols[k],
		   ylab="", 
		   main=plot_names[k],
		   cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


	for(i in 1:length(line_seq_x))
	{
		points(x=c(0,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
	}

	mtext(side = 2, line = 3.5, 
	cex=0.8*axis.size, 
	text="     antibody dilution")



	axis(1, at=c(1, 3, 5), label=c("Bichat", "Cochin",  "PNC"), cex.axis=0.8*axis.size )

	axis(1, at=c(2, 4, 6), label=c( "Sburg", "TRC", "EFS"), cex.axis=0.8*axis.size )


	axis(2, at=c(0.00001, 0.0001, 0.001, 0.01), 
      	  label=c("0.00001", "0.0001", "0.001", "0.01"),
      	  las=2, cex.axis=axis.size )
}





#####################################
##                                  
##  PANELS 8:12                     



plot_seq <- 43:47

line_seq_x <- c(10, 30, 100, 300, 1000, 3000, 10000, 30000)


for( k in 1:5 )
{
	boxplot( AB_data[which(AB_data$site=="Bichat"),plot_seq[k]],
		   AB_data[which(AB_data$site=="Strasbourg"),plot_seq[k]],
		   AB_data[which(AB_data$site=="Cochin"),plot_seq[k]],
	 	   AB_data[which(AB_data$site=="EFS"),plot_seq[k]],
		   AB_data[which(AB_data$site=="PNC"),plot_seq[k]], 
		   AB_data[which(AB_data$site=="TRC"),plot_seq[k]], 
		   log="y",
		   pch=19, yaxt='n', xaxt='n',
		   ylim=c(10, 30000),
		   col=ant_cols[7+k],
		   ylab="", 
		   main=plot_names[7+k],
		   cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


	for(i in 1:length(line_seq_x))
	{
		points(x=c(0,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
	}

	mtext(side = 2, line = 3.5, 
	cex=0.8*axis.size, 
	text="     MFI")



	axis(1, at=c(1, 3, 5), label=c("Bichat", "Cochin",  "PNC"), cex.axis=0.8*axis.size )

	axis(1, at=c(2, 4, 6), label=c( "Sburg", "TRC", "EFS"), cex.axis=0.8*axis.size )


	axis(2, at=c(10, 30, 100, 300, 1000, 3000, 10000, 30000), 
      	  label=c("10", "", "100", "", "1000", "", "10000", ""),
      	  las=2, cex.axis=axis.size )
}


############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))


plot.new()
title( "(B) IgM antibody responses", 
        cex.main=3.0, line=-2)


#####################################
##                                  
##  PANELS 13:19                     


par(mar=c(3,5,3,1))
par(mgp=c(2.5, 0.7, 0))

plot_seq <- 23:29

line_seq_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01)


for( k in 1:7 )
{
	boxplot( AB_data[which(AB_data$site=="Bichat"),plot_seq[k]],
		   AB_data[which(AB_data$site=="Strasbourg"),plot_seq[k]],
		   AB_data[which(AB_data$site=="Cochin"),plot_seq[k]],
	 	   AB_data[which(AB_data$site=="EFS"),plot_seq[k]],
		   AB_data[which(AB_data$site=="PNC"),plot_seq[k]], 
		   AB_data[which(AB_data$site=="TRC"),plot_seq[k]], 
		   log="y",
		   pch=19, yaxt='n', xaxt='n',
		   ylim=c(1e-5, 0.03),
		   col=ant_cols[12+k],
		   ylab="", 
		   main=plot_names[12+k],
		   cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


	for(i in 1:length(line_seq_x))
	{
		points(x=c(0,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
	}

	mtext(side = 2, line = 3.5, 
	cex=0.8*axis.size, 
	text="     antibody dilution")



	axis(1, at=c(1, 3, 5), label=c("Bichat", "Cochin",  "PNC"), cex.axis=0.8*axis.size )

	axis(1, at=c(2, 4, 6), label=c( "Sburg", "TRC", "EFS"), cex.axis=0.8*axis.size )


	axis(2, at=c(0.00001, 0.0001, 0.001, 0.01), 
      	  label=c("0.00001", "0.0001", "0.001", "0.01"),
      	  las=2, cex.axis=axis.size )
}



#####################################
##                                  
##  PANELS 20:24                     

plot_seq <- 56:60

line_seq_x <- c(10, 30, 100, 300, 1000, 3000, 10000, 30000)


for( k in 1:5 )
{
	boxplot( AB_data[which(AB_data$site=="Bichat"),plot_seq[k]],
		   AB_data[which(AB_data$site=="Strasbourg"),plot_seq[k]],
		   AB_data[which(AB_data$site=="Cochin"),plot_seq[k]],
	 	   AB_data[which(AB_data$site=="EFS"),plot_seq[k]],
		   AB_data[which(AB_data$site=="PNC"),plot_seq[k]], 
		   AB_data[which(AB_data$site=="TRC"),plot_seq[k]], 
		   log="y",
		   pch=19, yaxt='n', xaxt='n',
		   ylim=c(10, 30000),
		   col=ant_cols[19+k],
		   ylab="", 
		   main=plot_names[19+k],
		   cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


	for(i in 1:length(line_seq_x))
	{
		points(x=c(0,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
	}

	mtext(side = 2, line = 3.5, 
	cex=0.8*axis.size, 
	text="     MFI")



	axis(1, at=c(1, 3, 5), label=c("Bichat", "Cochin",  "PNC"), cex.axis=0.8*axis.size )

	axis(1, at=c(2, 4, 6), label=c( "Sburg", "TRC", "EFS"), cex.axis=0.8*axis.size )


	axis(2, at=c(10, 30, 100, 300, 1000, 3000, 10000, 30000), 
      	  label=c("10", "", "100", "", "1000", "", "10000", ""),
      	  las=2, cex.axis=axis.size )
}





#####################################
#####################################
##                                  
##  PANEL 25                       
##  IgG ROC analysis


line_seq <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)

par(mar=c(5,7,2,1.0))
par(mgp=c(2.5, 1.25,0))


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="(C) IgG ROC analysis",
cex.lab=2, cex.axis=1.5, cex.main=1.75)

mtext(side = 1, line = 3.4, 
cex=2, 
text="1 - specificity")

mtext(side = 2, line = 5, 
cex=2, 
text="sensitivity")

for(i in 1:length(line_seq))
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


for(k in 1:12)
{
	points( x=1-spec_mat[,k], y=sens_mat[,k], 
    		  type='S', lwd=2, col=ant_cols[k] )
}




axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        cex.axis=2) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=2 ) 




#####################################
#####################################
##                                  
##  PANEL 26                       
##  IgM ROC analysis


line_seq <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)

par(mar=c(5,7,2,1.0))
par(mgp=c(2.5, 1.25,0))


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="(D) IgM ROC analysis",
cex.lab=2, cex.axis=1.5, cex.main=1.75)

mtext(side = 1, line = 3.4, 
cex=2, 
text="1 - specificity")

mtext(side = 2, line = 5, 
cex=2, 
text="sensitivity")

for(i in 1:length(line_seq))
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}




for(k in 1:12)
{
	points( x=1-spec_mat[,13+k], y=sens_mat[,13+k], 
    		  type='S', lwd=2, col=ant_cols[k] )
}



axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        cex.axis=2) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=2 ) 


###################################
## Panel 27: AUC


ant_cols <- c("red", "brown", "orangered", "gold", "yellow",
              "forestgreen", "limegreen", "mediumblue", "royalblue", 
              "burlywood", "darkkhaki", "cornsilk4",
              "red", "brown", "orangered", "gold", "yellow",
              "forestgreen", "limegreen", "mediumblue", "royalblue", 
              "burlywood", "darkkhaki", "cornsilk4")	

plot_names_short <- c("Stri IgG",    "RBDv1 IgG", "RBDv2 IgG", "S1 IgG", "S2 IgG", 
                      "NPv1 IgG",    "NPv2 IgG",  "229E_NP IgG",
                      "NL63_NP IgG", "fluA IgG",  "ade40 IgG", "rub IgG",
                      "Stri IgM",    "RBDv1 IgM", "RBDv2 IgM", "S1 IgM", "S2 IgM", 
                      "NPv1 IgM",    "NPv2 IgM",  "229E_NP IgM","NL63_NP IgM", 
                      "fluA IgM",    "ade40 IgM", "rub IgM")


PP <- 1

gap <- 0.1

line_seq_y <- c(0.0, 0.25, 0.5, 0.75,  1)



par(mar=c(5,7,2,1.0))
par(mgp=c(2.5, 1.25,0))


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,24*(1+gap)), ylim=c(0,1.002),
xaxs='i',
# yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="(E) Area under ROC curve",
cex.lab=2, cex.axis=1.5, cex.main=1.75)



mtext(side = 2, line = 4.5, 
cex=2, 
text="AUC")



for(k in 1:26)
{
	polygon( x=gap*(k-1) + c(k-1,k,k,k-1), c(0,0,AUC_one_ant[k]^PP,AUC_one_ant[k]^PP), col=ant_cols[k], border=NA)
}



for(i in 1:length(line_seq_y))
{
	points(x=c(0,100), y=rep(line_seq_y[i]^PP,2), type='l', col="grey", lty="dashed")
}



axis(2, at=c(0.0, 0.25, 0.5, 0.75, 1)^PP, 
        labels=c("0.0", "0.25", "0.5", "0.75", "1.0"), 
        las=2, cex.axis=2 ) 



axis(1, at=seq(from=0.5, by=1+gap, length=24), label=plot_names_short, las=2, cex.axis=0.5)


#####################################
#####################################
##                                  
##  PANEL 28                       
##  Correlation


par(mar=c(2,3,2,1))
par(mgp=c(1.5,0.75,0))


##N_cor_steps <- 100
##
##cor_cols <- rev(heat.colors(N_cor_steps))


cor_cols = c("springgreen4", "springgreen", "palegreen", "yellowgreen", "yellow", 
             "gold", "orange", "orangered", "firebrick1", "red3" )

N_cor_steps = length(cor_cols)

plot(x=100, y=100,
xlim=c(0,24), ylim=c(0,24),
xlab="", ylab="",
main="(F) Correlation between antibody titres",
xaxt='n', yaxt='n', xaxs='i', yaxs='i',
cex.lab=2, cex.axis=1.5, cex.main=1.75)


for(i in 1:24)
{
	for(j in 1:24)
	{
		polygon(x=c(i-1,i,i,i-1), y=c(j-1,j-1,j,j),
                    border=NA, col="blue" )

		polygon(x=c(i-1,i,i,i-1), y=c(j-1,j-1,j,j),
                    border=NA, col=cor_cols[ceiling(N_cor_steps*AB_cor[i,j])] )
	}
}


axis(1, at=seq(from=0.5, by=1, length=24), label=plot_names_short, las=2, cex.axis=0.5)
axis(2, at=seq(from=0.5, by=1, length=24), label=plot_names_short, las=2, cex.axis=0.5)



#############
## LEGEND

par(mar=c(2,1,1,1))


plot(x=100, y=100,
xlim=c(-10,100), ylim=c(0,1),
xlab="", ylab="",
main="",
xaxt='n', yaxt='n', xaxs='i', yaxs='i', bty='n')

polygon(y=c(0,1,1,0), x=c(-10,-10,0,0),
              border=NA, col="blue")	

for(i in 1:length(cor_cols))
{
	polygon(y=c(0,1,1,0), x=c(i-1,i-1,i,i)*100/N_cor_steps,
              border=NA, col=cor_cols[i])	
}


axis(1, at=100*c(-0.1,0,0.25,0.5,0.75,1), label=c("-100%", "0%", "25%", "50%", "75%", "100%"), cex.axis=1)





dev.off()




