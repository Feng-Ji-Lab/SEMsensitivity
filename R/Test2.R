#' @title Naive Method with Approximate Influence
#' @description Remove a fixed percentage of samples (determined by the approximate influence) at a time and refit the model to observe the change in the parameter of interest. Use method 2 - Naive Method with Approximate Influence.
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
#' @param equalCons The equality constraint used in the SEM model (default is 0).
#' @param calcMeth The method used for approximation (default is 'Hessian').
#' @param ... Other arguments.
#' @importFrom lavaan sem
#' @import dplyr
#' @import semfindr


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
#' Test2_result <- Test2(df, model, var_one, var_two, PAR, threshold, fit, estimates,
#' conc, int, par_value, max_final, N, signFactor, equalCons, calcMeth)
#' summary(Test2_result)
#' }
#' @export

Test2 = function(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons=0, calcMeth = 'Hessian', ...){
  # #Perform calculations based on specified approximation method

  approx_infl <- getApproximationInfluence(df, calcMeth, N/(N-1), model, equalCons, PAR)

  #Sort influence and select most influential points to delete
  approx_sorted <- sort(approx_infl, decreasing = signFactor, index.return = TRUE)
  approx_drops <- approx_sorted$ix[1:max_final]

  ###Re-run LAVAAN with dropped points to test and determine new parameter value###
  new_par_value_appx = getValue(approx_drops,df, model,var_one, var_two)


  resultList = list()
  class(resultList) <- "TestResult"
  resultList$value = new_par_value_appx
  resultList$points = approx_drops
  resultList$methodname = "Naive method with approximate influence"
  resultList$testindex = 2

  resultList$PAR = PAR
  resultList$threshold = threshold
  resultList$N = N
  resultList$max_final = max_final
  resultList$par_value = par_value

  return(resultList)

}
