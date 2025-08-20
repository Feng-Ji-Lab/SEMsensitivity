#' @title Use Depth Method to Try to Switch Sign of Parameter
#' @description This function uses a depth method to iteratively remove data points in order to switch the sign of a specific path in a Structural Equation Modeling (SEM) model.
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
#' @return A list of class \code{TestResult8} containing:
#' \item{deletedPoints}{The indices of the most influential data points.}
#' \item{initialValue}{The original value of the parameter.}
#' \item{finalValue}{The final value of the parameter.}
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
#'   Test8_result = Test8(df, model, var_one, var_two, PAR, threshold, fit, estimates,
#'   conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth = "Hessian")
#' summary(Test8_result)
#' }
#' @export

Test8 = function(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor, equalCons, calcMeth = "Hessian", ...){

  checkIndDep <- 10 #How often to update internal model of parameter change using an exact calculation

  df_u <- df #Copy of data frame to modify
  df_u$index <- 1:nrow(df) #Add an additional index to keep track of observations during deletions

  dfu_og <- df #Copy of data frame to serve as reference (with an index added)
  dfu_og$index <- 1:nrow(df)  #Add an additional index to keep track of observations during deletions

  drops8depth <- 0
  drops8index <- integer(nrow(df) - 1)
  totalDropScore <- rep(0, nrow(df)) #To keep track of accumulating influence
  tolerance <- 0.001 * par_value #Tolerance to help prevent singularities by avoiding going too close to zero

  #Continue to drop points and accumulate influence until zero is reached
  for (u in 1:(nrow(df) - 1)) {
    model_u <- lavaan::sem(model, data = df_u) #update model given data removals

    if (u %% checkIndDep == 0) { #every checkInd times through the loop, update the drop score for better accuracy
      ResetDep <- lavaan::sem(model, data = df_u)
      ResetDep <- parameterEstimates(ResetDep)
      ResetDep <- ResetDep[ResetDep$lhs == var_one & ResetDep$rhs == var_two, ]
      ResetDepValue <- ResetDep$est

      totalDropScore[u] <- par_value - ResetDepValue
    }

    #Calculate influence scores for current model
    if (calcMeth == 'Hessian') {
      u_sc <- lavaan::estfun.lavaan(model_u, scaling = TRUE, ignore.constraints = TRUE, remove.duplicated = TRUE)
      u_asd <- lavInspect(model_u, "hessian")
      u_appr <- t(solve(u_asd[(equalCons + 1):dim(u_asd)[1], (equalCons + 1):dim(u_asd)[2]]) %*% t(u_sc))
      infl_u <- -u_appr[, PAR]
    } else {
      scu <- lavaan::estfun.lavaan(model_u, scaling = FALSE, ignore.constraints = TRUE)
      layerscores_u <- t(vcov(model_u) %*% t(scu))
      infl_u <- layerscores_u[, PAR]
    }

    #Sort influence
    u_approx_sorted <- sort(infl_u, decreasing = signFactor, index.return = TRUE)
    udepth_drop_actual <- u_approx_sorted$ix[1] #The single datapoint that will be dropped
    totalDropScore[u + 1] <- totalDropScore[u] + u_approx_sorted$x[1] #add influence of this point to our tally
    drops8index[u] <- df_u$index[udepth_drop_actual] #Keep track of which data point we actually dropped (since size of df_u is changing)
    df_u <- df_u[-c(udepth_drop_actual), ] #update df to remove point

    #When we are sufficiently close to zero, break
    if (totalDropScore[u + 1] >= par_value - tolerance) {
      break
    }
  }

  if (u != nrow(df) - 1) {
    #Calculate results
    df_u <- df[-c(drops8index), ]
    fit_depthu <- lavaan::sem(model, data = df_u)

    ###Get Estimates of Parameters from SEM###
    u_estimates <- parameterEstimates(fit_depthu)

    ###Find value of parameter
    u_conc <- data.frame(c(u_estimates$lhs), c(u_estimates$rhs), val = c(u_estimates$est))
    int_u <- u_conc %>% filter_all(any_vars(. %in% c(var_one))) %>% filter_all(any_vars(. %in% c(var_two)))
    u_par_value <- int_u$val
  } else {
    u_par_value <- NA
  }

  resultList <- list()
  class(resultList) <- "TestResult8"
  resultList$deletedPoints <- drops8index[1:u]
  resultList$initialValue <- par_value
  resultList$finalValue <- u_par_value
  resultList$methodname <- "Use depth method to try to switch sign of parameter"
  resultList$testindex <- 8

  resultList$PAR = PAR
  resultList$threshold = threshold
  resultList$N = N
  resultList$max_final = max_final
  resultList$par_value = par_value

  return(resultList)
}


#' @export
# Summary function for TestResult8
summary.TestResult8 <- function(object, ...) {
  cat("Summary of Depth Method to Switch Sign of Parameter Results:\n")
  if (!is.null(object$finalValue)) {

    cat(sprintf("Method Name: %s \n",object$methodname))
    cat(sprintf("Path: %s \n",object$PAR))
    cat(sprintf("Original Value: %f \n",object$par_value))

    cat(sprintf("Original Number of Samples: %d \n",object$N))
    cat(sprintf("Drop Points Percentage: %d \n",object$threshold))
    cat(sprintf("Dropped Number of Samples: %d \n",object$max_final))
    cat(sprintf("New Parameter Value: %f \n",object$finalValue))

    cat("drop points list: \n")
    print(object$deletedPoints)

    cat(sprintf("Using the depth method, dropping %d data points, the parameter was switched from %.4f to %.4f\n", length(object$deletedPoints), object$initialValue, object$finalValue))
  } else {
    cat("With all data dropped, the parameter sign could not be switched.\n")
  }
}
summary <- function(object, ...) UseMethod("summary")
