#' @title Use depth method to drop fit measure below threshold
#' @description This function uses the depth method to iteratively remove the most influential data points based on influence scores to reduce the fit measure below a specified threshold. The method is based on approximating the influence of each data point and is similar to TEST 8 but focused on fit metrics like CFI, RMSEA, etc.
#' @param df A data frame containing the dataset.
#' @param model A specified SEM model.
#' @param fit The SEM object after fitting the model.
#' @param max_final The maximum number of influential data points to consider for removal.
#' @param N The total number of data points.
#' @param measureTest The fit measurement name to be tested (e.g., "cfi", "tli", "rmsea"). Default is "cfi".
#' @param fitThreshold The threshold of the fit measurement. For example, for CFI (measureTest = "cfi"), the threshold could be 0.9. Default is 0.9.
#' @param highGood A boolean indicating whether higher values of the fit measure are better. For instance, for CFI, this should be set to TRUE. Default is TRUE.
#' @param ... Other arguments.
#' @return A list of class \code{TestResult10} containing:
#' \item{methodname}{The name of the method used.}
#' \item{testindex}{The index of the test performed.}
#' \item{original_fit_value}{The original value of the fit measurement before data points were removed.}
#' \item{final_fit_value}{The fit value after the influential data points were removed.}
#' \item{num_drops}{The number of data points dropped to achieve the desired fit value reduction.}
#' \item{depthdiffFit}{The difference between the original fit value and the threshold.}
#' \item{depthDropScore}{The cumulative drop in the fit value after removing points.}
#' \item{final_drops}{The indices of the most influential data points dropped.}
#' @examples
#' \dontrun{
#' library(lavaan)
#' library(dplyr)
#' library(semfindr)
#'
#' # Import data
#' df <- PoliticalDemocracy
#'
#' # Build SEM model
#' model <- '
#'   ind60 =~ x1 + x2 + x3
#'   dem60 =~ y1 + y2 + y3 + y4
#'   dem65 =~ y5 + y6 + y7 + y8
#'   dem60 ~ ind60
#'   dem65 ~ ind60 + dem60
#' '
#'
#' fit <- lavaan::sem(model, data = df)
#' max_final <- ceiling(10 * nrow(df) / 100)  # dropping 10% of the data points
#' N <- nrow(df)
#'
#' Test10_result <- Test10(df, model, fit, max_final, N,
#' measureTest = "cfi", fitThreshold = 0.9, highGood = TRUE)
#' summary(Test10_result)
#' }
#' @export

Test10 <- function(df, model, fit, max_final, N, measureTest = "cfi", fitThreshold = 0.9, highGood = TRUE, ...) {

  # Initialize necessary variables
  checkInd <- 10  # Interval for recalculating the drop score for accuracy
  depthnumDropsFit <- max_final  # Max number of data points to drop

  # Get original fit value
  depthoriginalFit <- lavInspect(fit, "fit")
  depthoriginalFitValue <- depthoriginalFit[measureTest]

  # Calculate the difference between the original fit and the threshold
  depthdiffFit <- depthoriginalFitValue - fitThreshold
  df_p <- df  # Create a copy of the dataset to manipulate
  depthDropScore <- 0  # Cumulative drop in fit score
  final_drops <- numeric(0)  # Initialize a vector to store dropped points

  for (p in 1:(nrow(df) - 1)) {
    model_p <- lavaan::sem(model, data = df_p)  # Update model with the current dataset after data removals

    if (p %% checkInd == 0) {  # Recalculate drop score every 'checkInd' iterations for accuracy
      fitReset <- lavaan::sem(model, data = df_p)
      fitResetValue <- lavInspect(fitReset, "fit")
      fitResetValue <- fitResetValue[measureTest]
      depthDropScore <- depthoriginalFitValue - fitResetValue
    }

    # Calculate influence scores using approximate influence method
    depthappxFitScore <- as.numeric(semfindr::fit_measures_change_approx(model_p, c(measureTest)))  # Get influence scores
    depthappxFitScoreSorted <- sort(depthappxFitScore, decreasing = highGood, index.return = TRUE)  # Sort influence scores

    # Remove the most influential data point
    df_p <- df_p[-c(depthappxFitScoreSorted$ix[1]),]  # Remove the most influential point from the dataset
    final_drops <- c(final_drops, depthappxFitScoreSorted$ix[1])  # Store dropped point indices
    depthDropScore <- depthDropScore + depthappxFitScoreSorted$x[1]  # Update the drop score

    # Stop if the cumulative drop score meets or exceeds the threshold difference
    if (depthDropScore >= depthdiffFit) {
      break
    }
  }

  # Final fit value after dropping the most influential points
  if (p < nrow(df) - 1) {
    finalfit <- lavaan::sem(model, data = df_p)
    depthfinalFit <- lavInspect(finalfit, "fit")
    depthfinalFitValue <- depthfinalFit[measureTest]
  }

  # Collect results into a list
  resultList <- list()
  class(resultList) <- "TestResult10"
  resultList$methodname <- "Use depth method to try to drop fit measure below threshold"
  resultList$testindex <- 10
  resultList$original_fit_value <- depthoriginalFitValue
  resultList$final_fit_value <- depthfinalFitValue
  resultList$num_drops <- p
  resultList$depthdiffFit <- depthdiffFit
  resultList$depthDropScore <- depthDropScore
  resultList$final_drops <- final_drops  # Store the indices of the most influential points dropped

  return(resultList)
}

#' @export
summary.TestResult10 <- function(object, ...) {
  if (!inherits(object, "TestResult10")) {
    stop("The object is not of class 'TestResult10'.")
  }

  cat("\nSummary of Depth Method to Drop Fit Measure Below Threshold: \n")
  cat("Method Name: ", object$methodname, "\n")
  cat("Original Fit Value: ", sprintf("%.4f", object$original_fit_value), "\n")
  cat("Final Fit Value After Dropping Points: ", sprintf("%.4f", object$final_fit_value), "\n")
  cat("Number of Data Points Dropped: ", object$num_drops, "\n")
  cat("Difference Between Original Fit and Threshold: ", sprintf("%.4f", object$depthdiffFit), "\n")
  cat("Cumulative Drop in Fit Value: ", sprintf("%.4f", object$depthDropScore), "\n")
  cat("Dropped Points: ", paste(object$final_drops, collapse = ", "), "\n")
}

# Create a generic summary function
summary <- function(object, ...) UseMethod("summary")
