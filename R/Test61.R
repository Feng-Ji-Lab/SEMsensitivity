#' @title Finding case deletions required to change fit metric using exact influences
#' @description Remove a fixed percentage of samples (determined by the exact influence) at a time and refit the model to observe the change in the path of interest. Use method 1 - Naive Method with Exact Influence.
#' @param df A data frame containing the dataset.
#' @param model A specified SEM model.
#' @param threshold The threshold for percentage of data dropped.
#' @param fit The SEM object.
#' @param max_final The maximum number of influential data points to consider.
#' @param N The total number of data points.
#' @param measureTest The fit measurement name. Can be "cfi", "chisq", "tli", "rmsea".
#' @param fitThreshold The threshold of the fit measurement to be a "good" model. For example, for CFI (measureTest = "cfi), this threshold can be 0.9.
#' @param highGood A boolean argument stating if the fit measurement is higher the better. For CFI, this argument is TRUE.
#' @param ... Other arguments.
#' @importFrom lavaan sem
#' @import dplyr
#' @import semfindr
#' @return A list of class \code{TestResult61} containing:
#' \item{methodname}{The name of the method used.}
#' \item{testindex}{The index of the test performed.}
#' \item{original_fit_value}{The original value of the fit measurement.}
#' \item{final_fit_value}{The fit value after dropping the influential points.}
#' \item{num_drops}{The number of data points dropped.}
#' \item{threshold_crossed}{Logical indicating whether the threshold was crossed.}
#' \item{final_drops}{The indices of the most influential data points dropped.}
#' \item{measureTest}{The name of the fit measurement used.}
#' \item{fitThreshold}{The threshold of the fit measurement used.}
#' \item{highGood}{Logical indicating if a higher value is better for the fit measurement.}
#' \item{exact_threshold_tally}{The tally value required to cross the threshold.}
#' \item{model_exact_threshold_final}{The final fit value after dropping points sufficient to cross the threshold.}
#' \item{max_final}{The maximum number of influential data points to consider.}
#' \item{N}{The total number of data points.}
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
#' Test61_result <- Test61(df, model, threshold, fit, max_final, N,
#'  measureTest = "cfi", fitThreshold = 0.9, highGood = T)
#' summary(Test61_result)
#' }
#' @export

Test61 <- function(df, model, threshold, fit, max_final, N, measureTest = "cfi", fitThreshold = 0.9, highGood = T, ...) {
  numDropsFit <- max_final

  # Determine original value of fit metric
  originalFit <- lavInspect(fit, "fit")
  originalFitValue <- originalFit[measureTest]

  # Determine difference between original value and threshold
  diffFit <- originalFitValue - fitThreshold

  # Get exact influence scores and sort them
  reruns <- lavaan_rerun(fit)
  exactFitScore <- as.numeric(fit_measures_change(reruns, fit_measures = measureTest))
  exactFitScoreSorted <- sort(exactFitScore, decreasing = highGood, index.return = TRUE)
  exactDropsFit <- exactFitScoreSorted$ix[1:numDropsFit]

  # Estimate whether all points are necessary to cross threshold
  exactThresholdTally <- 0
  for (tally in 1:(nrow(df) - 1)) {
    exactThresholdTally <- exactThresholdTally + exactFitScoreSorted$x[tally]
    if (exactThresholdTally >= diffFit) {
      break
    }
  }

  # Examine model after dropping points
  modelExactDropFit <- sem(model, df[-c(exactDropsFit), ])
  modelExactDropFitValues <- lavInspect(modelExactDropFit, "fit")
  modelExactDropFitFinal <- modelExactDropFitValues[measureTest]

  # Check if threshold was crossed
  threshold_crossed <- if (tally < (nrow(df) - 1)) {
    TRUE
  } else {
    FALSE
  }

  # If threshold was crossed, calculate the final fit after dropping required points
  if (threshold_crossed) {
    modelExactThresholdFit <- sem(model, df[-c(exactFitScoreSorted$ix[1:tally]), ])
    modelExactThresholdValues <- lavInspect(modelExactThresholdFit, "fit")
    modelExactTresholdFinal <- modelExactThresholdValues[measureTest]
  } else {
    modelExactTresholdFinal <- NA
    tally <- NA
  }

  resultList <- list()
  class(resultList) <- "TestResult61"
  resultList$methodname <- "naive method with exact influence"
  resultList$testindex <- 61
  resultList$original_fit_value <- originalFitValue
  resultList$final_fit_value <- modelExactDropFitFinal
  resultList$num_drops <- numDropsFit
  resultList$threshold_crossed <- threshold_crossed
  resultList$final_drops <- exactDropsFit
  resultList$measureTest <- measureTest
  resultList$fitThreshold <- fitThreshold
  resultList$highGood <- highGood
  resultList$exact_threshold_tally <- tally
  resultList$model_exact_threshold_final <- modelExactTresholdFinal
  resultList$max_final <- max_final
  resultList$N <- N

  return(resultList)
}

#' @export
summary.TestResult61 <- function(object, ...) {
  if (!inherits(object, "TestResult61")) {
    stop("The object is not of class 'TestResult61'.")
  }

  cat("\nSummary of Naive Method with Exact Influence: \n")
  cat("Method Name: ", object$methodname, "\n")
  cat("Original Fit Value: ", sprintf("%.4f", object$original_fit_value), "\n")
  cat("Final Fit Value After Dropping Points: ", sprintf("%.4f", object$final_fit_value), "\n")
  cat("Number of Data Points Dropped: ", object$num_drops, "\n")
  cat("Threshold Crossed: ", ifelse(object$threshold_crossed, "Yes", "No"), "\n")
  if (!is.na(object$exact_threshold_tally)) {
    cat("Number of Points Dropped to Cross Threshold: ", object$exact_threshold_tally, "\n")
    cat("Final Fit Value After Crossing Threshold: ", sprintf("%.4f", object$model_exact_threshold_final), "\n")
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

summary <- function(object, ...) UseMethod("summary")
