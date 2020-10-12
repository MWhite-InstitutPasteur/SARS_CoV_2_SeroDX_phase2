# SARS_CoV_2_SeroDX_phase2

This folder contains code in R and C++ for replicating the analysis on Serological Signatures of SARS-CoV-2 infection as reported here.
https://www.medrxiv.org/content/10.1101/2020.05.07.20093963v2

SARS-CoV-2_serology.csv
A csv file containing measurements of antibody responses from our Luminex MAGPIX multiplex assay.

Fig1_single_AB.R
An R script for preliminaryof antibody responses, plotting ROC curves, calculating AUCs, and calculating correlations.

Fig2_multi_AB.R
Analysis of the assocation between multiple measurements of antibody response and SARS-CoV-2 infection status. A random forests classification algorithm (with cross-validation) is used.

Fig3_AB_kinetics.R
R script taking the posterior outputs of the C++ code for antibody kinetics and plotting antibody levels over time and sensitivity.

Fig4_kinetic_classifier.R
R script for model predicted sensitivity over time using a single antigen (spike) and multiplex combinations analysed with random forests algorithms. Depends on C++ output.

Fig5_sero_surveillance.R
R script for statistical assessment of the trade-off between sensitivity and specificity on serological surveys.

SupFig1_crosspanel_validation.R
R script for quantification of uncertainty of serological classification.

SARSCoV2_IgG_antigen_combination.xlsx
Excel spreadsheet detailing classification performance of all combinations of IgG responses to 7 SARS-CoV-2 antigens using a random forests classifier.

Source.cpp
C++ code for Bayesian statistical inference of anntibody kinetic models.

com.cpp
linpack.cpp
randlib.cpp
randlib.cpp
Randlib library for generation of random numbers for MCMC algorithm.

Chain_diagnostics.R
R script for plotting the outputted MCMC chains from C++. 

Stri_IPP_IgG_dil.txt
RBD_IPP_IgG_dil.txt
S1RBD_NA_IgG_dil.txt
S1_NA_IgG_dil.txt
S2_NA_IgG_dil.txt
NP_IPP_IgG_dil.txt
NP_NA_IgG_dil.txt
Input files for C++ statistical inference. Note that these files were genereated by prcoessing SARS-CoV-2_serology.csv.


