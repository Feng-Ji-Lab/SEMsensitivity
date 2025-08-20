#' @title Use Brute Search with Cut Method to Try to Switch Sign of Parameter
#' @description This function uses a brute-force search method with pruning (cutting) to iteratively search for combinations of data points that switch the sign of a specific path in a Structural Equation Modeling (SEM) model.
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
#' @return A list of class \code{TestResult13} containing:
#' \item{bruteSearchDrops}{The indices of the most influential data points selected by the brute search method.}
#' \item{initialValue}{The original value of the parameter.}
#' \item{finalValue}{The final value of the parameter after applying the brute search method.}
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
#' Test13_result = Test13(df, model, var_one, var_two, PAR, threshold, fit, estimates,
#' conc, int, par_value, max_final, N, signFactor)
#' summary(Test13_result)
#' }
#' @export

Test13 = function(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor, ...) {

  initialize_best_combinations <- function() {
    best_combinations <- list()
    best_values <- numeric()
    return(list(best_combinations = best_combinations, best_values = best_values))
  }

  update_best_combinations <- function(best_combinations, best_values, new_combination, new_value, max_best) {
    if (length(best_values) < max_best || new_value < max(best_values)) {
      if (length(best_values) >= max_best) {
        max_index <- which.max(best_values)
        best_combinations[[max_index]] <- new_combination
        best_values[max_index] <- new_value
      } else {
        best_combinations <- c(best_combinations, list(new_combination))
        best_values <- c(best_values, new_value)
      }
    }
    return(list(best_combinations = best_combinations, best_values = best_values))
  }

  generate_combinations <- function(n, k, max_best) {
    best <- initialize_best_combinations()

    for (i in 1:n) {
      value <- getValue(i,df, model,var_one, var_two)
      best <- update_best_combinations(best$best_combinations, best$best_values, c(i), value, max_best)
    }

    for (depth in 2:k) {
      new_best_combinations <- list()
      new_best_values <- numeric()
      for (combo in best$best_combinations) {
        for (i in setdiff(1:n, combo)) {
          new_combination <- c(combo, i)
          value <- getValue(new_combination,df, model,var_one, var_two)
          new_best <- update_best_combinations(new_best_combinations, new_best_values, new_combination, value, max_best)
          new_best_combinations <- new_best$best_combinations
          new_best_values <- new_best$best_values
        }
      }
      best$best_combinations <- new_best_combinations
      best$best_values <- new_best_values
    }

    return(best)
  }

  result <- generate_combinations(nrow(df), max_final, 5)
  brute_search_value <- min(result$best_values)
  brute_search_drops <- result$best_combinations[[which.min(result$best_values)]]

  print(sprintf("The brute search with cut method yields a new parameter value of %f", brute_search_value))

  resultList <- list()
  class(resultList) <- "TestResult13"
  resultList$bruteSearchDrops <- brute_search_drops
  resultList$initialValue <- par_value
  resultList$finalValue <- brute_search_value
  resultList$methodname <- "Brute Search with Cut Method to Switch Sign of Parameter"
  resultList$testindex <- 13
  resultList$PAR <- PAR
  resultList$threshold <- threshold
  resultList$N <- N
  resultList$max_final <- max_final
  resultList$par_value <- par_value

  return(resultList)
}


#' @export
summary.TestResult13 <- function(object, ...) {
  cat("Summary of Brute Search with Cut Method to Switch Sign of Parameter Results:\n")
  if (!is.null(object$finalValue)) {

    cat(sprintf("Method Name: %s \n", object$methodname))
    cat(sprintf("Path: %s \n", object$PAR))
    cat(sprintf("Original Value: %f \n", object$initialValue))

    cat(sprintf("Original Number of Samples: %d \n", object$N))
    cat(sprintf("Drop Points Percentage: %d \n", object$threshold))
    cat(sprintf("Dropped Number of Samples: %d \n", length(object$bruteSearchDrops)))
    cat(sprintf("New Parameter Value: %f \n", object$finalValue))

    cat("Brute search drop points list: \n")
    print(object$bruteSearchDrops)

    cat(sprintf("Using the brute search with cut method, dropping %d data points, the parameter was switched from %.4f to %.4f\n", length(object$bruteSearchDrops), object$initialValue, object$finalValue))
  } else {
    cat("The brute search with cut method was unable to switch the parameter sign.\n")
  }
}
summary <- function(object, ...) UseMethod("summary")
