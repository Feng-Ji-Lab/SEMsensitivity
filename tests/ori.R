# Load packages
library(lavaan)
library(dplyr)
#library(tidyverse)
#library(gridExtra)
library(semfindr) #to install enter remotes::install_github("sfcheung/semfindr")
library(R.utils) #necessary for the negamax algorithm
#library(semPlot)
#library(mirt)
# library(simsem)



############################# USER INPUTS ######################

#Import data
df <- PoliticalDemocracy

###Build Model###
model <- '
  # measurement model
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  # regressions
  dem60 ~ ind60
  dem65 ~ ind60 + dem60
  # residual correlations
  y1 ~~ y5
  y2 ~~ y4 + y6
  y3 ~~ y7
  y4 ~~ y8
  y6 ~~ y8
  '
equalCons <- 0 #specify number of equality constraints in model (for calculations)
var_one = 'dem65' #first term
var_two = 'ind60' #second term
PAR <- c("dem65~ind60") #full relation
#################### SELECT INFLUENCE ANALYSIS PARAMETERS
###CHANGE ME - What parameter you are interested in###
test1 <- T #Exact method
test2 <- T #Approximate method
test3 <- T #Sign switch test
test4 <- T #Depth method
test5 <- T #Combined method
test6 <- T #Fit index tests
test7 <- T #Pruned Negamax method
test8 <- T #Depth method for switching sign
test9 <- T #Negamax method for switching sign
test10 <- T #Simple Depth for Fit Metrics
test11 <- T #anneal search method
test12 <- T #pso method
test13 <- T #brute force with cut method

##Change me - 'Hessian' or 'Covariance' Approximation Method
calcMeth <- 'Hessian'

###CHANGE ME - percentage of data points to drop
threshold <- 10


#################################################################################

#Fit SEM model
fit <- lavaan::sem(model, data=df)
summary(fit)

###Get Estimates of Parameters from SEM###
estimates <- parameterEstimates(fit)

###Determine The Value of The Parameter of Interest
conc <- data.frame(c(estimates$lhs), c(estimates$rhs), c(estimates$est), val = c(estimates$est))
int <- conc %>% filter_all(any_vars(. %in% c(var_one))) %>% filter_all(any_vars(. %in% c(var_two)))
par_value <- int$c.estimates.est. #this is the value of the parameter of interest

###Compute max number of points to be dropped###
max_final <- ceiling(threshold*nrow(df)/100) #perform rounding if necessary
N = nrow(df) #store number of observations in df for convenience

#Determine whether parameter is negative or positive in order
#to assess which direction to perturb it
if (par_value >= 0){
  signFactor <- T } else {signFactor <- F
}
# .SEMClusterEnv$getValueWithMemory = myMemorize(getValueInternal)

## BEGIN TESTS ##
###TEST 13 - brute force with cut method ###################################################
if (exists('test13') && test13 == T){
  Test13_result = Test13(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor)
  summary(Test13_result)
}

###TEST 12 - pso method ###################################################
if (exists('test12') && test12 == T){
  Test12_result = Test12(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor)
  summary(Test12_result)

}

###TEST 11 - anneal search method ###################################################
if (exists('test11') && test11 == T){
  Test11_result = Test11(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor)
  summary(Test11_result)
}
###TEST 1 - Naive Method with Exact Influence ###################################################
if (exists('test1') && test1 == T){
  Test1_result = Test1(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor)
  summary(Test1_result)
}


######## Test 2: Naive method with approximate influence #######

if(exists('test2') && test2 == T){
  Test2_result = Test2(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons=0, calcMeth = 'Hessian')
  summary(Test2_result)

}


###TEST 3 - Drop points until sign of parameter is switched used naive approximate method ########################

if(exists('test3') && test3 == T){
  Test3_result = Test3(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth)
  summary(Test3_result)

}

### TEST 4: Simple DEPTH METHOD ####

if(exists('test4') && test4 == T){
  Test4_result = Test4(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor)
  summary(Test4_result)

}


if(exists('test5') && test5 == T){
  Test5_result = Test5(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor)
  summary(Test5_result)

}


### Test 6 - Finding case deletions required to change fit metric###
## REQUIRES semfindr Package ##
if(exists('test6') && test6 == T){

  Test6a_result = Test6a(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth = "Hessian",measureTest="cfi", fitThreshold = 0.9, highGood = T)


  ##6b Using appx scores ##
    Test6b_result = Test6b(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth = "Hessian",measureTest="cfi", fitThreshold = 0.9, highGood = T)
}



################## Test 7 - Negamax Search Using Approx. Scores ######

if(exists('test7') && test7 == T){

  Test7_result = Test7(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth = "Hessian",adaptA = F, alpha= 0.25, maxTime = 300, pruneNum = 3)
  summary(Test7_result)
}





## Test 8: Use depth method to try to switch sign of parameter ####
if (exists('test8') && test8 == T){
  Test8_result = Test8(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth = "Hessian")
  summary(Test8_result)
}




################## Test 9 - Using Negamax to try to switch sign of parameter################

## NOTE: as of right now, this only uses approximate influence scores calculated via the Hessian (because of some singularity issues)
#There is also no option for automatic tuning
#And it requires test8 to be run to generate target data

if(exists('test9') && test9 == T){
  Test9_result = Test9(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth = "Hessian")
  summary(Test9_result)


}

###### COULD DEPTH METHOD USE EXACT COMPUTATION FOR LAST STEP IF REACHES A SINGULAR CALCULATION????
#TRY THIS OUT!!


###### DEPTH METHOD FOR FIT MEASURES
## Test 10: Use depth method to try to drop fit measure below threshold ####

if (exists('test10') && test10 == T){
  Test10_result = Test10(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth = "Hessian", depthmeasureTest = "cfi", depthfitThreshold = 0.9)
  summary(Test10_result)

}







