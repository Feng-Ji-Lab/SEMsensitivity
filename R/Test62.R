#' @title Finding case deletions required to change fit metric using approximate method
#' @description Remove a fixed percentage of samples (determined by the approximate influence) at a time and refit the model to observe the change in the fit metric of interest. Use method 2 - Approximate Method.
#' @param df A data frame containing the dataset.
#' @param model A specified SEM model.
#' @param threshold The threshold for percentage of data dropped.
#' @param fit The SEM object.
#' @param max_final The maximum number of influential data points to consider.
#' @param N The total number of data points.
#' @param equalCons Logical; whether equality constraints exist in the model. The approximate method can only be run if no equality constraints are present (default is 0).
#' @param measureTest The fit measurement name. Can be "cfi", "chisq", "tli", "rmsea".
#' @param fitThreshold The threshold of the fit measurement to be a "good" model. For example, for CFI (measureTest = "cfi"), this threshold can be 0.9.
#' @param highGood A boolean argument stating if the fit measurement is higher the better. For CFI, this argument is TRUE.
#' @param ... Other arguments.
#' @importFrom lavaan sem
#' @import dplyr
#' @import semfindr
#' @return A list of class \code{TestResult62} containing:
#' \item{methodname}{The name of the method used.}
#' \item{testindex}{The index of the test performed.}
#' \item{original_fit_value}{The original value of the fit measurement.}
#' \item{final_fit_value}{The fit value after dropping the influential points using the approximate method.}
#' \item{num_drops}{The number of data points dropped.}
#' \item{threshold_crossed}{Logical indicating whether the threshold was crossed.}
#' \item{final_drops}{The indices of the most influential data points dropped using the approximate method.}
#' \item{measureTest}{The name of the fit measurement used.}
#' \item{fitThreshold}{The threshold of the fit measurement used.}
#' \item{highGood}{Logical indicating if a higher value is better for the fit measurement.}
#' \item{appx_threshold_tally}{The tally value required to cross the threshold using the approximate method.}
#' \item{model_appx_threshold_final}{The final fit value after dropping points sufficient to cross the threshold using the approximate method.}
#' \item{max_final}{The maximum number of influential data points to consider.}
#' \item{N}{The total number of data points.}
#' \item{equalCons}{Logical indicating whether equality constraints exist in the model.}
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
#' threshold <- 10
#'
#' # Fit SEM model
#' fit <- lavaan::sem(model, data = df)
#' summary(fit)
#'
#' # Compute max number of points to be dropped
#' max_final <- ceiling(threshold * nrow(df) / 100)
#' N <- nrow(df)
#'
#' Test62_result <- Test62(df, model, threshold, fit, max_final, N, equalCons = 0,
#'  measureTest = "cfi", fitThreshold = 0.9, highGood = T)
#' summary(Test62_result)
#' }
#' @export

Test62 <- function(df, model, threshold, fit, max_final, N, equalCons = 0, measureTest = "cfi", fitThreshold = 0.9, highGood = T, ...) {
  resultList <- list()
  class(resultList) <- "TestResult62"

  # Check if the approximate method can be run
  if (equalCons == 0) {
    # Get original fit value
    originalFit <- lavInspect(fit, "fit")
    originalFitValue <- originalFit[measureTest]

    # Calculate the difference between original fit and the threshold
    diffFit <- originalFitValue - fitThreshold

    # Get influence scores, sort and drop
    appxFitScore <- as.numeric(fit_measures_change_approx(fit, c(measureTest)))
    appxFitScoreSorted <- sort(appxFitScore, decreasing = highGood, index.return = TRUE)
    appxDropsFit <- appxFitScoreSorted$ix[1:max_final]

    # Estimate whether all points are needed to cross the threshold
    appxThresholdTally <- 0
    for (tally2 in 1:(nrow(df) - 1)) {
      appxThresholdTally <- appxThresholdTally + appxFitScoreSorted$x[tally2]
      if (appxThresholdTally >= diffFit) {
        break
      }
    }

    # Examine model after dropping points
    modelAppxDropFit <- lavaan::sem(model, df[-c(appxDropsFit), ])
    modelAppxDropFitValues <- lavInspect(modelAppxDropFit, "fit")
    modelAppxDropFitFinal <- modelAppxDropFitValues[measureTest]

    # Check if threshold was crossed
    threshold_crossed <- if (tally2 < (nrow(df) - 1)) {
      TRUE
    } else {
      FALSE
    }

    # If threshold was crossed, calculate the final fit after dropping required points
    if (threshold_crossed) {
      modelAppxThresholdFit <- lavaan::sem(model, df[-c(appxFitScoreSorted$ix[1:tally2]), ])
      modelAppxThresholdValues <- lavInspect(modelAppxThresholdFit, "fit")
      modelAppxThresholdFinal <- modelAppxThresholdValues[measureTest]
    } else {
      modelAppxThresholdFinal <- NA
      tally2 <- NA
    }

    # Add calculated values to the result list
    resultList = list()
    class(resultList) <- "TestResult62"
    resultList$methodname <- "Finding case deletions required to change fit metric using approximate method"
    resultList$testindex <- 62
    resultList$original_fit_value <- originalFitValue
    resultList$final_fit_value <- modelAppxDropFitFinal
    resultList$num_drops <- max_final
    resultList$threshold_crossed <- threshold_crossed
    resultList$final_drops <- appxDropsFit
    resultList$measureTest <- measureTest
    resultList$fitThreshold <- fitThreshold
    resultList$highGood <- highGood
    resultList$appx_threshold_tally <- tally2
    resultList$model_appx_threshold_final <- modelAppxThresholdFinal
    resultList$max_final <- max_final
    resultList$N <- N
    resultList$equalCons <- equalCons

  } else {
    # Cannot run approximate method due to equality constraints
    resultList$methodname <- "approximate method"
    resultList$testindex <- 62
    resultList$warning <- sprintf("Warning: can only run appx method for %s using semFindr if no equality constraints", measureTest)
    resultList$equalCons <- equalCons
  }

  return(resultList)
}

#' @export
summary.TestResult62 <- function(object, ...) {
  if (!inherits(object, "TestResult62")) {
    stop("The object is not of class 'TestResult62'.")
  }

  cat("\nSummary of Finding case deletions required to change fit metric using approximate method: \n")
  cat("Method Name: ", object$methodname, "\n")

  if (!is.null(object$warning)) {
    cat(object$warning, "\n")
  } else {
    cat("Original Fit Value: ", sprintf("%.4f", object$original_fit_value), "\n")
    cat("Final Fit Value After Dropping Points: ", sprintf("%.4f", object$final_fit_value), "\n")
    cat("Number of Data Points Dropped: ", object$num_drops, "\n")
    cat("Threshold Crossed: ", ifelse(object$threshold_crossed, "Yes", "No"), "\n")
    if (!is.na(object$appx_threshold_tally)) {
      cat("Number of Points Dropped to Cross Threshold: ", object$appx_threshold_tally, "\n")
      cat("Final Fit Value After Crossing Threshold: ", sprintf("%.4f", object$model_appx_threshold_final), "\n")
    } else {
      cat("Threshold Not Crossed.\n")
    }
    cat("Dropped Points: ", paste(object$final_drops, collapse = ", "), "\n")
    cat("Fit Measurement Used: ", object$measureTest, "\n")
    cat("Fit Threshold Used: ", sprintf("%.4f", object$fitThreshold), "\n")
    cat("Higher Value Better (for Fit Measurement): ", ifelse(object$highGood, "Yes", "No"), "\n")
    cat("Maximum Points Considered for Dropping: ", object$max_final, "\n")
    cat("Total Number of Data Points: ", object$N, "\n")
  }
}

summary <- function(object, ...) UseMethod("summary")
