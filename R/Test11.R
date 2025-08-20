#' @title Use Simulated Annealing Method to Try to Switch Sign of Parameter
#' @description This function uses the simulated annealing method to iteratively remove data points in order to switch the sign of a specific path in a Structural Equation Modeling (SEM) model.
#' @param df A data frame containing the dataset.
#' @param model A specified SEM model.
#' @param var_one The first variable of interest.
#' @param var_two The second variable of interest.
#' @param PAR The path of interest.
#' @param threshold The threshold for the percentage of data dropped.
#' @param fit The SEM object.
#' @param estimates The estimates from the SEM model.
#' @param conc A data frame containing the parameter of interest.
#' @param int The value of the path of interest.
#' @param par_value The original value of the parameter of interest.
#' @param max_final The maximum number of influential data points to consider.
#' @param N The total number of data points.
#' @param signFactor A factor indicating the direction of parameter change (positive or negative).
#' @param ... Other arguments.
#' @importFrom lavaan sem
#' @import dplyr
#' @import semfindr
#' @importFrom stats vcov
#' @importFrom lavaan lavInspect
#' @importFrom lavaan parameterEstimates
#' @importFrom stats runif
#' @return A list of class \code{TestResult11} containing:
#' \item{annealingDrops}{The indices of the most influential data points selected by the simulated annealing method.}
#' \item{initialValue}{The original value of the parameter.}
#' \item{finalValue}{The final value of the parameter after applying the simulated annealing method.}
#' \item{methodname}{The name of the method used.}
#' \item{testindex}{The index of the test performed.}
#' \item{PAR}{The path of interest that was evaluated.}
#' \item{threshold}{The threshold used for the percentage of dropped points.}
#' \item{N}{The total number of data points.}
#' \item{max_final}{The maximum number of points allowed to be dropped.}
#' \item{par_value}{The original value of the parameter.}
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
#' # Determine the value of the parameter of interest
#' conc <- data.frame(lhs = estimates$lhs, rhs = estimates$rhs, est = estimates$est)
#' int <- conc %>% filter(lhs == var_one & rhs == var_two)
#' par_value <- int$est # this is the value of the parameter of interest
#'
#' # Compute the max number of points to be dropped
#' max_final <- ceiling(threshold * nrow(df) / 100) # perform rounding if necessary
#' N <- nrow(df) # store the number of observations in df for convenience
#'
#' # Determine whether the parameter is negative or positive in order
#' # to assess which direction to perturb it
#' signFactor <- ifelse(par_value >= 0, TRUE, FALSE)
#'
#' Test11_result = Test11(df, model, var_one, var_two, PAR, threshold, fit, estimates,
#' conc, int, par_value, max_final, N, signFactor)
#' summary(Test11_result)
#' }
#' @export

Test11 = function(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor, ...) {
  simulated_annealing <- function(num_deletions, current_state = sample(1:N, num_deletions, replace = FALSE), initial_temp = 100, cooling_rate = 0.03, min_temp = 0.1, max_iter = 100) {
    num_points <- nrow(df)

    current_value <- getValue(current_state,df, model,var_one, var_two)
    best_state <- current_state
    best_value <- current_value

    temp <- initial_temp

    for (iter in 1:max_iter) {
      if (temp < min_temp) break
      new_state <- current_state
      add_index <- sample(setdiff(1:num_points, current_state), 1)
      remove_index <- sample(current_state, 1)
      new_state[new_state == remove_index] <- add_index

      new_value <- getValue(new_state,df, model,var_one, var_two)

      if (new_value < current_value || runif(1) < exp((current_value - new_value) / temp)) {
        current_state <- new_state
        current_value <- new_value
      }

      if (current_value < best_value) {
        best_state <- current_state
        best_value <- current_value
      }
      temp <- temp * (1 - cooling_rate)
    }

    return(list(best_state = best_state, best_value = best_value))
  }

  # Get initial index
  inflX <- replicate(N, 0)
  par_values <- replicate(N, 0)
  for (droppingP in 1:N) {
    par_values[droppingP] <- getValue(droppingP,df, model,var_one, var_two)
    inflX[droppingP] <- par_value - par_values[droppingP]
  }
  infl_sorted <- sort(inflX, decreasing = signFactor, index.return = TRUE)
  drops <- infl_sorted$ix[1:max_final]

  result <- simulated_annealing(threshold, drops)
  annealing_drops <- result$best_state
  annealing_value <- result$best_value

  resultList <- list()
  class(resultList) <- "TestResult11"
  resultList$annealingDrops <- annealing_drops
  resultList$initialValue <- par_value
  resultList$finalValue <- annealing_value
  resultList$methodname <- "Simulated Annealing Method to Switch Sign of Parameter"
  resultList$testindex <- 11
  resultList$PAR <- PAR
  resultList$threshold <- threshold
  resultList$N <- N
  resultList$max_final <- max_final
  resultList$par_value <- par_value

  return(resultList)
}

#' @export
summary.TestResult11 <- function(object, ...) {
  cat("Summary of Simulated Annealing Method to Switch Sign of Parameter Results:\n")
  if (!is.null(object$finalValue)) {

    cat(sprintf("Method Name: %s \n", object$methodname))
    cat(sprintf("Path: %s \n", object$PAR))
    cat(sprintf("Original Value: %f \n", object$initialValue))

    cat(sprintf("Original Number of Samples: %d \n", object$N))
    cat(sprintf("Drop Points Percentage: %d \n", object$threshold))
    cat(sprintf("Dropped Number of Samples: %d \n", length(object$annealingDrops)))
    cat(sprintf("New Parameter Value: %f \n", object$finalValue))

    cat("Simulated annealing drop points list: \n")
    print(object$annealingDrops)

    cat(sprintf("Using the simulated annealing method, dropping %d data points, the parameter was switched from %.4f to %.4f\n", length(object$annealingDrops), object$initialValue, object$finalValue))
  } else {
    cat("The simulated annealing method was unable to switch the parameter sign.\n")
  }
}
summary <- function(object, ...) UseMethod("summary")
