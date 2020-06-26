library(binom)
library(MASS)
library(randomForest)
library(pROC)


N_tree <- 5000

N_rep <- 1000

N_SS_cut <- 5000

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

N_dil_cut <- 10000

dil_cut <- exp(seq(from=log(1e-5), to=log(0.3), length=N_dil_cut))


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



########################################
########################################
##                                    ##
##  #####   ####  #   ## ##### ##     ##
##  ##  ## ##  ## ##  ## ##    ##     ##
##  #####  ###### ### ## ####  ##     ##
##  ##     ##  ## ## ### ##    ##     ## 
##  ##     ##  ## ##  ## ##### #####  ##
##                                    ##
########################################
########################################

##################################################
## Cross-panel validation

#########################################
## 4 antigen algorithm: cross-panels

RF_Bichat_EFS <- RF_multi_ant_SS( AB_train = AB[which(AB_data$site %in% c("Cochin", "Strasbourg", "PNC", "TRC")),], 
                                  status_train = status[which(AB_data$site %in% c("Cochin", "Strasbourg", "PNC", "TRC"))], 
                                  AB_test = AB[which(AB_data$site %in% c("Bichat_early", "Bichat_late", "EFS")),], 
                                  status_test = status[which(AB_data$site %in% c("Bichat_early", "Bichat_late", "EFS"))], 
                                  antigen_list = c(1,2,3,4,5,6) )


RF_Bichat_PNC <- RF_multi_ant_SS( AB_train = AB[which(AB_data$site %in% c("Cochin", "Strasbourg", "EFS", "TRC")),], 
                                  status_train = status[which(AB_data$site %in% c("Cochin", "Strasbourg", "EFS", "TRC"))], 
                                  AB_test = AB[which(AB_data$site %in% c("Bichat_early", "Bichat_late", "PNC")),], 
                                  status_test = status[which(AB_data$site %in% c("Bichat_early", "Bichat_late", "PNC"))], 
                                  antigen_list = c(1,2,3,4,5,6) )


RF_Bichat_TRC <- RF_multi_ant_SS( AB_train = AB[which(AB_data$site %in% c("Cochin", "Strasbourg", "PNC", "EFS")),], 
                                  status_train = status[which(AB_data$site %in% c("Cochin", "Strasbourg", "PNC", "EFS"))], 
                                  AB_test = AB[which(AB_data$site %in% c("Bichat_early", "Bichat_late", "TRC")),], 
                                  status_test = status[which(AB_data$site %in% c("Bichat_early", "Bichat_late", "TRC"))], 
                                  antigen_list = c(1,2,3,4,5,6) )




RF_Stras_EFS <- RF_multi_ant_SS( AB_train = AB[which(AB_data$site %in% c("Cochin", "Bichat_early", "Bichat_late", "PNC", "TRC")),], 
                                 status_train = status[which(AB_data$site %in% c("Cochin", "Bichat_early", "Bichat_late", "PNC", "TRC"))], 
                                 AB_test = AB[which(AB_data$site %in% c("Strasbourg", "EFS")),], 
                                 status_test = status[which(AB_data$site %in% c("Strasbourg", "EFS"))], 
                                 antigen_list = c(1,2,3,4,5,6) )



RF_Stras_PNC <- RF_multi_ant_SS( AB_train = AB[which(AB_data$site %in% c("Cochin", "Bichat_early", "Bichat_late", "EFS", "TRC")),], 
                                 status_train = status[which(AB_data$site %in% c("Cochin", "Bichat_early", "Bichat_late", "EFS", "TRC"))], 
                                 AB_test = AB[which(AB_data$site %in% c("Strasbourg", "PNC")),], 
                                 status_test = status[which(AB_data$site %in% c("Strasbourg", "PNC"))], 
                                 antigen_list = c(1,2,3,4,5,6) )


RF_Stras_TRC <- RF_multi_ant_SS( AB_train = AB[which(AB_data$site %in% c("Cochin", "Bichat_early", "Bichat_late", "EFS", "PNC")),], 
                                 status_train = status[which(AB_data$site %in% c("Cochin", "Bichat_early", "Bichat_late", "EFS", "PNC"))], 
                                 AB_test = AB[which(AB_data$site %in% c("Strasbourg", "TRC")),], 
                                 status_test = status[which(AB_data$site %in% c("Strasbourg", "TRC"))], 
                                 antigen_list = c(1,2,3,4,5,6) )



RF_Cochin_EFS <- RF_multi_ant_SS( AB_train = AB[which(AB_data$site %in% c("Strasbourg", "Bichat_early", "Bichat_late", "PNC", "TRC")),], 
                                 status_train = status[which(AB_data$site %in% c("Strasbourg", "Bichat_early", "Bichat_late", "PNC", "TRC"))], 
                                 AB_test = AB[which(AB_data$site %in% c("Cochin", "EFS")),], 
                                 status_test = status[which(AB_data$site %in% c("Cochin", "EFS"))], 
                                 antigen_list = c(1,2,3,4,5,6) )



RF_Cochin_PNC <- RF_multi_ant_SS( AB_train = AB[which(AB_data$site %in% c("Strasbourg", "Bichat_early", "Bichat_late", "EFS", "TRC")),], 
                                 status_train = status[which(AB_data$site %in% c("Strasbourg", "Bichat_early", "Bichat_late", "EFS", "TRC"))], 
                                 AB_test = AB[which(AB_data$site %in% c("Cochin", "PNC")),], 
                                 status_test = status[which(AB_data$site %in% c("Cochin", "PNC"))], 
                                 antigen_list = c(1,2,3,4,5,6) )


RF_Cochin_TRC <- RF_multi_ant_SS( AB_train = AB[which(AB_data$site %in% c("Strasbourg", "Bichat_early", "Bichat_late", "EFS", "PNC")),], 
                                 status_train = status[which(AB_data$site %in% c("Strasbourg", "Bichat_early", "Bichat_late", "EFS", "PNC"))], 
                                 AB_test = AB[which(AB_data$site %in% c("Cochin", "TRC")),], 
                                 status_test = status[which(AB_data$site %in% c("Cochin", "TRC"))], 
                                 antigen_list = c(1,2,3,4,5,6) )






#########################################
## 1 antigen: cross-panels

dil_cut <- exp(seq(from=log(1e-5), to=log(0.3), length=N_dil_cut))


sens_1ant_panel <- matrix(NA, ncol=12, nrow=N_dil_cut)
spec_1ant_panel <- matrix(NA, ncol=12, nrow=N_dil_cut)

colnames(sens_1ant_panel) <- c( "Bichat_EFS", "Bichat_PNC", "Bichat_TRC",
                                "Stras_EFS", "Stras_PNC", "Stras_TRC",
                                "Cochin_EFS", "Cochin_PNC", "Cochin_TRC",
                                "Bichat_late_EFS", "Bichat_late_PNC", "Bichat_late_TRC" )
colnames(spec_1ant_panel) <- colnames(sens_1ant_panel)


for(i in 1:N_dil_cut)
{
	#####################
	## Bichat EFS

	Spike_pos <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("Bichat_early", "Bichat_late"))]
	Spike_neg <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("EFS"))]

	sens_1ant_panel[i,1] <- length(which( Spike_pos > dil_cut[i] ))/length(Spike_pos)
	spec_1ant_panel[i,1] <- length(which( Spike_neg < dil_cut[i] ))/length(Spike_neg)


	#####################
	## Bichat PNC

	Spike_pos <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("Bichat_early", "Bichat_late"))]
	Spike_neg <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("PNC"))]

	sens_1ant_panel[i,2] <- length(which( Spike_pos > dil_cut[i] ))/length(Spike_pos)
	spec_1ant_panel[i,2] <- length(which( Spike_neg < dil_cut[i] ))/length(Spike_neg)


	#####################
	## Bichat TRC

	Spike_pos <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("Bichat_early", "Bichat_late"))]
	Spike_neg <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("TRC"))]

	sens_1ant_panel[i,3] <- length(which( Spike_pos > dil_cut[i] ))/length(Spike_pos)
	spec_1ant_panel[i,3] <- length(which( Spike_neg < dil_cut[i] ))/length(Spike_neg)


	#####################
	## Stras EFS

	Spike_pos <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("Strasbourg"))]
	Spike_neg <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("EFS"))]

	sens_1ant_panel[i,4] <- length(which( Spike_pos > dil_cut[i] ))/length(Spike_pos)
	spec_1ant_panel[i,4] <- length(which( Spike_neg < dil_cut[i] ))/length(Spike_neg)


	#####################
	## Stras PNC

	Spike_pos <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("Strasbourg"))]
	Spike_neg <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("PNC"))]

	sens_1ant_panel[i,5] <- length(which( Spike_pos > dil_cut[i] ))/length(Spike_pos)
	spec_1ant_panel[i,5] <- length(which( Spike_neg < dil_cut[i] ))/length(Spike_neg)


	#####################
	## Stras TRC

	Spike_pos <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("Strasbourg"))]
	Spike_neg <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("TRC"))]

	sens_1ant_panel[i,6] <- length(which( Spike_pos > dil_cut[i] ))/length(Spike_pos)
	spec_1ant_panel[i,6] <- length(which( Spike_neg < dil_cut[i] ))/length(Spike_neg)


	#####################
	## Cochin EFS

	Spike_pos <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("Cochin"))]
	Spike_neg <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("EFS"))]

	sens_1ant_panel[i,7] <- length(which( Spike_pos > dil_cut[i] ))/length(Spike_pos)
	spec_1ant_panel[i,7] <- length(which( Spike_neg < dil_cut[i] ))/length(Spike_neg)


	#####################
	## Cochin PNC

	Spike_pos <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("Cochin"))]
	Spike_neg <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("PNC"))]

	sens_1ant_panel[i,8] <- length(which( Spike_pos > dil_cut[i] ))/length(Spike_pos)
	spec_1ant_panel[i,8] <- length(which( Spike_neg < dil_cut[i] ))/length(Spike_neg)


	#####################
	## Cochin TRC

	Spike_pos <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("Cochin"))]
	Spike_neg <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("PNC"))]

	sens_1ant_panel[i,9] <- length(which( Spike_pos > dil_cut[i] ))/length(Spike_pos)
	spec_1ant_panel[i,9] <- length(which( Spike_neg < dil_cut[i] ))/length(Spike_neg)


	#####################
	## Bichat_late EFS

	Spike_pos <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("Bichat_late"))]
	Spike_neg <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("EFS"))]

	sens_1ant_panel[i,10] <- length(which( Spike_pos > dil_cut[i] ))/length(Spike_pos)
	spec_1ant_panel[i,10] <- length(which( Spike_neg < dil_cut[i] ))/length(Spike_neg)


	#####################
	## Bichat PNC

	Spike_pos <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("Bichat_late"))]
	Spike_neg <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("PNC"))]

	sens_1ant_panel[i,11] <- length(which( Spike_pos > dil_cut[i] ))/length(Spike_pos)
	spec_1ant_panel[i,11] <- length(which( Spike_neg < dil_cut[i] ))/length(Spike_neg)


	#####################
	## Bichat TRC

	Spike_pos <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("Bichat_late"))]
	Spike_neg <- AB_data$Stri_IPP_IgG_dil[which(AB_data$site %in% c("TRC"))]

	sens_1ant_panel[i,12] <- length(which( Spike_pos > dil_cut[i] ))/length(Spike_pos)
	spec_1ant_panel[i,12] <- length(which( Spike_neg < dil_cut[i] ))/length(Spike_neg)
}

save.image("SupFig1_crosspanel_validation.RData")


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

tiff(file="SupFig1_crosspanel_validation.tif", width=40, height=30, units="cm", res=500)

lay.mat <- rbind( c( 1, 1, 7, 7, 7 ),
                  c( 2, 3, 8, 9,10 ),
                  c( 2, 3,11,12,13 ),
                  c( 4, 4,11,12,13 ),
                  c( 5, 6,11,12,13 ),
                  c( 5, 6,14,15,16 ) )
layout(lay.mat, heights=c(0.1,1,0.45,0.1,0.45,1), widths=c(3,3,2,2,2))
layout.show(16)


############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))


plot.new()
title( "(A) Binomial uncertainty", 
        cex.main=3.0, line=-2)



lab.size   = 2
axis.size  = 1.25
main.size  = 2

par(mar=c(5,7,2,1.5))
par(mgp=c(2.5, 1.25,0))

line_seq_x <- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)
line_seq_y <- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)

#####################################
#####################################
##                                 ## 
##  PANEL 1                        ##
##  Approach 1: frequentist        ##
##  Variation in specificity       ##
##                                 ##
#####################################
##################################### 

plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Variation in specificity",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3.5, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 4, 
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



points( x=1 - RF_multi_ant_ROC$specificities^PP, y=RF_multi_ant_ROC$sensitivities^PP, 
    	  type='S', lwd=2, col="aquamarine4" )


for(i in 1:(N_SS_multi_ant-1))
{
	polygon(x= c(1 - spec_multi_ant_binom[c(i,i+1),2]^PP, rev(1 - spec_multi_ant_binom[c(i,i+1),3]^PP)), 
		  y=c( RF_multi_ant_ROC$sensitivities[c(i,i+1)]^PP, rev(RF_multi_ant_ROC$sensitivities[c(i,i+1)]^PP) ),
	        col=rgb(69/256,139/256,116/256,0.2), border=NA)
}

points( x=1 - spec_1ant_binom[,1]^PP, y=sens_1ant_binom[,1]^PP, 
    	  type='S', lwd=2, col="orangered" )

for(i in 1:(N_dil_cut-1))
{
	polygon(x=c(1 - spec_1ant_binom[c(i,i+1),2]^PP, rev(1 - spec_1ant_binom[c(i,i+1),3]^PP)), 
		  y=c( sens_1ant_binom[c(i,i+1),1]^PP, rev(sens_1ant_binom[c(i,i+1),1]^PP) ),
	        col=rgb(190/256,90/256,90/256,0.2), border=NA)
}

axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 



#####################################
#####################################
##                                 ## 
##  PANEL 2                        ##
##  Approach 1: frequentist        ##
##  Variation in sensitivity       ##
##                                 ##
#####################################
##################################### 

plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Variation in sensitivity",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3.5, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 4, 
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



points( x=1 - RF_multi_ant_ROC$specificities^PP, y=RF_multi_ant_ROC$sensitivities^PP, 
    	  type='S', lwd=2, col="aquamarine4" )

for(i in 1:(N_SS_multi_ant-1))
{
	polygon(x=c(1 - RF_multi_ant_ROC$specificities[c(i,i+1)]^PP, rev(1 - RF_multi_ant_ROC$specificities[c(i,i+1)]^PP)), 
		  y=c( sens_multi_ant_binom[c(i,i+1),2]^PP, rev(sens_multi_ant_binom[c(i,i+1),3]^PP) ),
      	  col=rgb(69/256,139/256,116/256,0.2), border=NA)
}

points( x=1 - spec_1ant_binom[,1]^PP, y=sens_1ant_binom[,1]^PP, 
    	  type='S', lwd=2, col="orangered" )

for(i in 1:(N_dil_cut-1))
{
	polygon(x=c(1 - spec_1ant_binom[c(i,i+1),1]^PP, rev(1 - spec_1ant_binom[c(i,i+1),1]^PP)), 
		  y=c( sens_1ant_binom[c(i,i+1),2]^PP, rev(sens_1ant_binom[c(i,i+1),3]^PP) ),
	        col=rgb(190/256,90/256,90/256,0.2), border=NA)
}


axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 



############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))


plot.new()
title( "(B) Cross-validated uncertainty", 
        cex.main=3.0, line=-2)





par(mar=c(5,7,2,1.5))
par(mgp=c(2.5, 1.25,0))

line_seq_x <- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)
line_seq_y <- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)

#####################################
#####################################
##                                 ## 
##  PANEL 3                        ##
##  Approach 2: cross-validated    ##
##  Variation in specificity       ##
##                                 ##
#####################################
##################################### 

plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Variation in specificity",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3.5, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 4, 
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



points( x=1 - RF_multi_ant_xval_se_fix_sp_vary_quant[,1]^PP, y=SS_seq^PP, 
    	  type='S', lwd=2, col="aquamarine4" )


for(i in 1:(N_SS_cut-1))
{
	polygon(x=c(1 - RF_multi_ant_xval_se_fix_sp_vary_quant[c(i,i+1),2]^PP, rev(1 - RF_multi_ant_xval_se_fix_sp_vary_quant[c(i,i+1),3]^PP)), 
		  y=c( SS_seq[c(i,i+1)]^PP, rev(SS_seq[c(i,i+1)]^PP) ),
      	  col=rgb(69/256,139/256,116/256,0.2), border=NA)
}


points( x=1 - spec_1ant_binom[,1]^PP, y=sens_1ant_binom[,1]^PP, 
    	  type='S', lwd=2, col="orangered" )

for(i in 1:(N_dil_cut-1))
{
	polygon(x=c(1 - spec_1ant_xval_quant[c(i,i+1),2]^PP, rev(1 - spec_1ant_xval_quant[c(i,i+1),3]^PP)), 
		  y=c( sens_1ant_binom[c(i,i+1),1]^PP, rev(sens_1ant_binom[c(i,i+1),1]^PP) ),
	        col=rgb(190/256,90/256,90/256,0.2), border=NA)
}

axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 



#####################################
#####################################
##                                 ## 
##  PANEL 4                        ##
##  Approach 2: cross-validated    ##
##  Variation in specificity       ##
##                                 ##
#####################################
##################################### 

plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Variation in sensitivity",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3.5, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 4, 
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



points( x=1 - SS_seq^PP, y=RF_multi_ant_xval_se_vary_sp_fix_quant[,1]^PP, 
    	  type='S', lwd=2, col="aquamarine4" )


for(i in 1:(N_SS_cut-1))
{
	polygon(x=c(1 - SS_seq[c(i,i+1)]^PP, rev(1 - SS_seq[c(i,i+1)]^PP)), 
		  y=c( RF_multi_ant_xval_se_vary_sp_fix_quant[c(i,i+1),2]^PP, rev(RF_multi_ant_xval_se_vary_sp_fix_quant[c(i,i+1),3]^PP) ),
	        col=rgb(69/256,139/256,116/256,0.2), border=NA)
}

points( x=1 - spec_1ant_binom[,1]^PP, y=sens_1ant_binom[,1]^PP, 
    	  type='S', lwd=2, col="orangered" )

for(i in 1:(N_dil_cut-1))
{
	polygon(x=c(1 - spec_1ant_binom[c(i,i+1),1]^PP, rev(1 - spec_1ant_binom[c(i,i+1),1]^PP)), 
		  y=c( sens_1ant_xval_quant[c(i,i+1),2]^PP, rev(sens_1ant_xval_quant[c(i,i+1),3]^PP) ),
	        col=rgb(190/256,90/256,90/256,0.2), border=NA)
}

axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 




############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))


plot.new()
title( "(C) Cross-panel validation", 
        cex.main=3.0, line=-2)



par(mar=c(5,5.5,2,1))
par(mgp=c(2.5, 1.25,0))

lab.size   = 1.25
axis.size  = 1.25
main.size  = 1.5

#####################################
#####################################
##                                 ## 
##  PANEL 5                        ##
##  Positive: Bichat               ##
##  Negative: EFS                  ##
##                                 ##
#####################################
##################################### 


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Bichat & EFS",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 3.5, 
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


points( x=1 - RF_Bichat_EFS[,2]^PP, y=RF_Bichat_EFS[,1]^PP, 
    	  type='S', lwd=2, col="aquamarine4" )



points( x=1 - spec_1ant_panel[,1]^PP, y=sens_1ant_panel[,1]^PP, 
    	  type='S', lwd=2, col="red" )



axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 




#####################################
#####################################
##                                 ## 
##  PANEL 6                        ##
##  Positive: Bichat               ##
##  Negative: PNC                  ##
##                                 ##
#####################################
##################################### 


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Bichat & PNC",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 3.5, 
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


points( x=1 - RF_Bichat_PNC[,2]^PP, y=RF_Bichat_PNC[,1]^PP, 
    	  type='S', lwd=2, col="aquamarine4" )



points( x=1 - spec_1ant_panel[,2]^PP, y=sens_1ant_panel[,2]^PP, 
    	  type='S', lwd=2, col="red" )




axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 




#####################################
#####################################
##                                 ## 
##  PANEL 7                        ##
##  Positive: Bichat               ##
##  Negative: TRC                  ##
##                                 ##
#####################################
##################################### 


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Bichat & TRC",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 3.5, 
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


points( x=1 - RF_Bichat_TRC[,2]^PP, y=RF_Bichat_TRC[,1]^PP, 
    	  type='S', lwd=2, col="aquamarine4" )


points( x=1 - spec_1ant_panel[,3]^PP, y=sens_1ant_panel[,3]^PP, 
    	  type='S', lwd=2, col="red" )



axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 



#####################################
#####################################
##                                 ## 
##  PANEL 8                        ##
##  Positive: Strasbourg           ##
##  Negative: EFS                  ##
##                                 ##
#####################################
##################################### 


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Strasbourg & EFS",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 3.5, 
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


points( x=1 - RF_Stras_EFS[,2]^PP, y=RF_Stras_EFS[,1]^PP, 
    	  type='S', lwd=2, col="aquamarine4" )



points( x=1 - spec_1ant_panel[,4]^PP, y=sens_1ant_panel[,4]^PP, 
    	  type='S', lwd=2, col="red" )



axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 




#####################################
#####################################
##                                 ## 
##  PANEL 9                        ##
##  Positive: Strasbourg           ##
##  Negative: PNC                  ##
##                                 ##
#####################################
##################################### 


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Strasbourg & PNC",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 3.5, 
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


points( x=1 - RF_Stras_PNC[,2]^PP, y=RF_Stras_PNC[,1]^PP, 
    	  type='S', lwd=2, col="aquamarine4" )



points( x=1 - spec_1ant_panel[,5]^PP, y=sens_1ant_panel[,5]^PP, 
    	  type='S', lwd=2, col="red" )



axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 



#####################################
#####################################
##                                 ## 
##  PANEL 10                       ##
##  Positive: Strasbourg           ##
##  Negative: TRC                  ##
##                                 ##
#####################################
##################################### 


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Strasbourg & TRC",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 3.5, 
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


points( x=1 - RF_Stras_TRC[,2]^PP, y=RF_Stras_TRC[,1]^PP, 
    	  type='S', lwd=2, col="aquamarine4" )


points( x=1 - spec_1ant_panel[,6]^PP, y=sens_1ant_panel[,6]^PP, 
    	  type='S', lwd=2, col="red" )



axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 




#####################################
#####################################
##                                 ## 
##  PANEL 11                       ##
##  Positive: Cochin               ##
##  Negative: EFS                  ##
##                                 ##
#####################################
##################################### 


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Cochin & EFS",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 3.5, 
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


points( x=1 - RF_Cochin_EFS[,2]^PP, y=RF_Cochin_EFS[,1]^PP, 
    	  type='S', lwd=2, col="aquamarine4" )



points( x=1 - spec_1ant_panel[,7]^PP, y=sens_1ant_panel[,7]^PP, 
    	  type='S', lwd=2, col="red" )



axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 






#####################################
#####################################
##                                 ## 
##  PANEL 12                       ##
##  Positive: Cochin               ##
##  Negative: PNC                  ##
##                                 ##
#####################################
##################################### 


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Cochin & PNC",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 3.5, 
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


points( x=1 - RF_Cochin_PNC[,2]^PP, y=RF_Cochin_PNC[,1]^PP, 
    	  type='S', lwd=2, col="aquamarine4" )



points( x=1 - spec_1ant_panel[,8]^PP, y=sens_1ant_panel[,8]^PP, 
    	  type='S', lwd=2, col="red" )



axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 








#####################################
#####################################
##                                 ## 
##  PANEL 13                       ##
##  Positive: Cochin               ##
##  Negative: TRC                  ##
##                                 ##
#####################################
##################################### 


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="Cochin & TRC",
cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

mtext(side = 1, line = 3, 
cex=lab.size, 
text="1 - specificity")

mtext(side = 2, line = 3.5, 
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


points( x=1 - RF_Cochin_TRC[,2]^PP, y=RF_Cochin_TRC[,1]^PP, 
    	  type='S', lwd=2, col="aquamarine4" )


points( x=1 - spec_1ant_panel[,9]^PP, y=sens_1ant_panel[,9]^PP, 
    	  type='S', lwd=2, col="red" )



axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=axis.size) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=axis.size ) 




dev.off()




