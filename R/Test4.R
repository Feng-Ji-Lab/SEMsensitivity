#' @title Simple depth method
#' @description This function determines a specific path in a Structural Equation Modeling (SEM) model value changing by removing samples iteratively, in which influences are determined by naive method and outputs relevant results. Use method 4 - Simple depth method.
#' @param df A data frame containing the dataset.
#' @param model A specified SEM model.
#' @param var_one The first variable of interest.
#' @param var_two The second variable of interest.
#' @param PAR The path of interest.
#' @param threshold The threshold for percentage of data dropped.
#' @param fit The SEM object.
#' @param estimates The estimates from the SEM model.
#' @param conc A data frame containing the parameter of interest.
#' @param int The value of the path of interest.
#' @param par_value The original value of the parameter of interest.
#' @param max_final The maximum number of influential data points to consider.
#' @param N The total number of data points.
#' @param signFactor A factor indicating the direction of parameter change (positive or negative).
#' @param equalCons The equality constraint used in the SEM model.
#' @param calcMeth The method used for approximation (default is 'Hessian').
#' @param ... Other arguments.
#' @importFrom lavaan sem
#' @import dplyr
#' @import semfindr
#' @importFrom stats vcov


#' @importFrom lavaan lavInspect
#' @importFrom lavaan parameterEstimates
#' @return A list of class \code{TestResult} containing:
#' \item{value}{The value of the parameter after dropping the influential points.}
#' \item{points}{The indices of the most influential data points.}
#' \item{methodname}{The name of the method used.}
#' \item{testindex}{The index of the test performed.}
#' @examples
#' \dontrun{
#' library(lavaan)
#' library(dplyr)
#' library(semfindr)
#' library(R.utils)
#' library(simsem)
#'
#' # Import data
#' df <- PoliticalDemocracy
#'
#' # Build Model
#' model <- '
#'   # measurement model
#'   ind60 =~ x1 + x2 + x3
#'   dem60 =~ y1 + y2 + y3 + y4
#'   dem65 =~ y5 + y6 + y7 + y8
#'   # regressions
#'   dem60 ~ ind60
#'   dem65 ~ ind60 + dem60
#'   # residual correlations
#'   y1 ~~ y5
#'   y2 ~~ y4 + y6
#'   y3 ~~ y7
#'   y4 ~~ y8
#'   y6 ~~ y8
#' '
#'
#' var_one <- 'dem65' # first term
#' var_two <- 'ind60' # second term
#' PAR <- c("dem65~ind60") # full relation
#' threshold <- 10
#'
#' # Fit SEM model
#' fit <- lavaan::sem(model, data = df)
#' summary(fit)
#'
#' # Get Estimates of Parameters from SEM
#' estimates <- parameterEstimates(fit)
#'
#' # Determine The Value of The Parameter of Interest
#' conc <- data.frame(lhs = estimates$lhs, rhs = estimates$rhs, est = estimates$est)
#' int <- conc %>% filter(lhs == var_one & rhs == var_two)
#' par_value <- int$est # this is the value of the parameter of interest
#'
#' # Compute max number of points to be dropped
#' max_final <- ceiling(threshold * nrow(df) / 100) # perform rounding if necessary
#' N <- nrow(df) # store number of observations in df for convenience
#'
#' # Determine whether parameter is negative or positive in order
#' # to assess which direction to perturb it
#' signFactor <- ifelse(par_value >= 0, TRUE, FALSE)
#'
#' Test4_result = Test4(df, model, var_one, var_two, PAR, threshold, fit, estimates,
#' conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth)
#' summary(Test4_result)
#' }
#' @export

Test4 = function(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons=0,calcMeth = "Hessian", ...){
  df_depth <- df #Store new copy of data frame (since we will modify this)
  df_depth$index = 1:nrow(df) #Create an index column to help keep track of which original points dropped (since data frame changes each loop)
  drops_index <- replicate(max_final, 0) #Store indices for drops
  tallyDropScore <- rep(0, nrow(df)) #Used to generate targets if using negamax search

  #Iteratively drop the best point, recalculating influence after each drop
  for (j in 1:max_final){

    model_j <- lavaan::sem(model, data=df_depth) #update model given any data removals

    #Calculate influence based on chosen approximation method
    if(calcMeth == 'Hessian'){
      j_sc <- lavaan::estfun.lavaan(model_j, scaling = T, ignore.constraints = T)
      j_asd <- lavInspect(model_j, "hessian")
      j_appr <- t(solve(j_asd[(equalCons + 1):dim(j_asd)[1], (equalCons + 1):dim(j_asd)[2]]) %*% t(j_sc))
      j_approx_infl <- -j_appr[,PAR]
    } else{
      sc2 <- lavaan::estfun.lavaan(model_j, scaling = F, ignore.constraints = T)
      j_approx_infl <- t(vcov(model_j)%*%t(sc2))
      j_approx_infl <- j_approx_infl[, PAR]
    }

    #Sort influence
    j_approx_sorted <- sort(j_approx_infl, decreasing = signFactor, index.return = TRUE)

    #Create target set for negamax search (not necessary if only using simple depth)
    tallyDropScore[j + 1] <- tallyDropScore[j] + j_approx_sorted$x[1]

    #Drop point
    depth_drop_actual <- j_approx_sorted$ix[1] #The single data point that will be dropped (most influential)
    drops_index[j] <- df_depth$index[depth_drop_actual] #Store the position of the drop (relative to the indexing of the original data frame)

    df_depth  <- df_depth[-c(depth_drop_actual),] #update df to remove point
  }

  #Fit model with full drops
  fit_depth <- lavaan::sem(model, data=df_depth)

  ###Get Estimates of Parameters from SEM###
  j_estimates <- parameterEstimates(fit_depth)

  ###Determine value of Parameter We Want####
  j_conc=data.frame(c(j_estimates$lhs), c(j_estimates$rhs), val = c(j_estimates$est))
  int_j <- j_conc %>% filter_all(any_vars(. %in% c(var_one))) %>% filter_all(any_vars(. %in% c(var_two)))
  depth_par_value = int_j$val


  resultList = list()
  class(resultList) <- "TestResult"
  resultList$value = depth_par_value
  resultList$points = drops_index
  resultList$methodname = "Simple DEPTH METHOD"
  resultList$testindex = 4

  resultList$PAR = PAR
  resultList$threshold = threshold
  resultList$N = N
  resultList$max_final = max_final
  resultList$par_value = par_value

  return(resultList)

}
