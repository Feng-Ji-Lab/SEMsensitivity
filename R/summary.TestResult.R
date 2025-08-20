
#' @export
summary.TestResult <- function(object, ...) {
  ###Print Results###


  cat("Summary of SEM sensitivity analysis result: \n")

  cat(sprintf("Method Name: %s \n",object$methodname))
  cat(sprintf("Path: %s \n",object$PAR))
  cat(sprintf("Original Value: %f \n",object$par_value))

  cat(sprintf("Original Number of Samples: %d \n",object$N))
  cat(sprintf("Drop Points Percentage: %d \n",object$threshold))
  cat(sprintf("Dropped Number of Samples: %d \n",object$max_final))
  cat(sprintf("New Parameter Value: %f \n",object$value))

  cat("drop points list: \n")
  print(object$points)
}
summary <- function(object, ...) UseMethod("summary")
