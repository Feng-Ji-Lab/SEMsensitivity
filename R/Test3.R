#' @title Specified Approximation Method
#' @description This function determines when a specific path in a Structural Equation Modeling (SEM) model changes sign by iteratively removing most influential data points (determined by naive method) and outputs relevant results. Use method 3 - Specified Approximation Method.
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
#' @importFrom lavaan lavInspect
#' @importFrom lavaan parameterEstimates
#' @importFrom stats vcov
#' @import dplyr
#' @import semfindr


#' @importFrom utils head
#' @return A list of class \code{TestResult3} containing:
#' \item{methodname}{The name of the method used.}
#' \item{testindex}{The index of the test performed.}
#' \item{max_drops}{The maximum number of data points dropped.}
#' \item{est_diff}{The expected change in parameter value.}
#' \item{act_diff}{The actual change in parameter value.}
#' \item{final_par_value}{The parameter value after dropping the influential points.}
#' \item{initial_par_value}{The original parameter value.}
#' \item{sign_switch_possible}{Logical indicating whether it was possible to change the sign of the parameter.}
#' \item{dropped_points}{The indices of the most influential data points.}
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
#' Test3_result = Test3(df, model, var_one, var_two, PAR, threshold, fit, estimates,
#' conc, int, par_value, max_final, N, signFactor,equalCons = 0, calcMeth)
#' summary(Test3_result)
#' }
#' @export

Test3 = function(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons=0,calcMeth = "Hessian", ...){
  if(calcMeth == 'Hessian'){
    ##Using Hessian
    nabla3 <- lavaan::estfun.lavaan(fit, scaling = T, ignore.constraints = T) #Scaled empirical estimating function
    hess3 <- lavInspect(fit, "hessian") #Hessian matrix
    appr3 <- t(solve(hess3[(equalCons + 1):dim(hess3)[1], (equalCons + 1):dim(hess3)[2]]) %*% t(nabla3)) #calculate influence matrix
    approx_infl3 <- -(N/(N-1))*appr3[,PAR]
  }
  else{
    #using covariance
    nabla4 <- lavaan::estfun.lavaan(fit, scaling = F, ignore.constraints = T) #Unscaled empirical estimating function
    approx_infl3 <- t(vcov(fit)%*%t(nabla4)) #Calculate influence matrix using variance-covariance matrix
    approx_infl3 <- approx_infl3[, PAR] #select appropriate data for given variable of interest
    approx_infl3<- (N/(N-1))*approx_infl3 #Adjustment factor
  }

  #Sort influence and select most influential points to delete
  approx_sorted_3 <- sort(approx_infl3, decreasing = signFactor, index.return = TRUE)

  ###For estimating change###
  tot = 0 #counting variable
  count = 0 #counting variable
  arg = 1 #truth checker

  #Iterate over most influential points, accumulating influence
  #this continues until predicted influence is large enough to switch sign
  for (i in 1:N) {
    tot = tot + approx_sorted_3$x[i]
    count = count + 1

    if(abs(tot) > abs(par_value)){break}

    if(sign(approx_sorted_3$x[i]) == sign(-par_value)){arg = 2} #If we have exhausted all influence values of one sign, quit
    if(arg == 2){break}
  }

  if (arg == 2) {
    sprintf("Failure, unable to switch sign.")
  }

  all_drops = head(approx_sorted_3$ix, count) #keep track of the points we dropped

  #calculate model after deletions
  second_df_drops <- df[-c(all_drops),]
  second_fit_drops = lavaan::sem(model, data=second_df_drops)
  test3Ints <- parameterEstimates(second_fit_drops)
  test3Int_frame =data.frame(c(test3Ints$lhs), c(test3Ints$rhs), val = c(test3Ints$est))
  test3Par <- test3Int_frame  %>% filter_all(any_vars(. %in% c(var_one))) %>% filter_all(any_vars(. %in% c(var_two)))
  test3ParValue = test3Par$val

  est_diff = -1*sum(approx_sorted_3$x[1:count]) #Estimated change in parameter (based on influence scores)
  act_diff = test3ParValue - par_value #Actual change in parameter
  print(sprintf("The maximum number of data points dropped was %d", count))
  print(sprintf("The expected change in parameter value was %.4f, while the actual change was %.4f", est_diff, act_diff))
  print(sprintf("The parameter at the biggest change is %.2f, while it was originally %.2f", test3ParValue, par_value))
  if (arg == 2) {
    print(sprintf("It was not possible to change the sign of the parameter"))
  }

  resultList <- list()
  class(resultList) <- "TestResult3"
  resultList$methodname <- "specified approximation method"
  resultList$testindex <- 3
  resultList$max_drops <- count
  resultList$est_diff <- est_diff
  resultList$act_diff <- act_diff
  resultList$final_par_value <- test3ParValue
  resultList$initial_par_value <- par_value
  resultList$sign_switch_possible <- ifelse(arg == 2, FALSE, TRUE)
  resultList$dropped_points <- all_drops

  resultList$PAR = PAR
  resultList$threshold = threshold
  resultList$N = N
  resultList$max_final = max_final
  resultList$par_value = par_value

  return(resultList)

}

#' @export
summary.TestResult3 <- function(object, ...) {
  if (!inherits(object, "TestResult3")) {
    stop("The object is not of class 'TestResult3'.")
  }

  cat("\nSummary of Specified Approximation Method: \n")
  cat("Method Name: ", object$methodname, "\n")
  # cat("Test Index: ", object$testindex, "\n")
  cat("Maximum Number of Data Points Dropped: ", object$max_drops, "\n")
  cat("Expected Change in Parameter Value: ", sprintf("%.4f", object$est_diff), "\n")
  cat("Actual Change in Parameter Value: ", sprintf("%.4f", object$act_diff), "\n")
  cat("Parameter Value at Biggest Change: ", sprintf("%.2f", object$final_par_value), "\n")
  cat("Original Parameter Value: ", sprintf("%.2f", object$initial_par_value), "\n")

  if (object$sign_switch_possible) {
    cat("It was possible to change the sign of the parameter.\n")
  } else {
    cat("It was not possible to change the sign of the parameter.\n")
  }

  cat("Dropped Points: ", paste(object$dropped_points, collapse = ", "), "\n")
}
summary <- function(object, ...) UseMethod("summary")
