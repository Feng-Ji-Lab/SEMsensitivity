#' @title Use Particle Swarm Optimization (PSO) to Try to Switch Sign of Parameter
#' @description This function uses the Particle Swarm Optimization (PSO) method to iteratively remove data points in order to switch the sign of a specific path in a Structural Equation Modeling (SEM) model.
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
#' @importFrom methods setRefClass
#' @return A list of class \code{TestResult12} containing:
#' \item{psoDrops}{The indices of the most influential data points selected by the PSO method.}
#' \item{initialValue}{The original value of the parameter.}
#' \item{finalValue}{The final value of the parameter after applying the PSO method.}
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
#' Test12_result = Test12(df, model, var_one, var_two, PAR, threshold, fit, estimates,
#' conc, int, par_value, max_final, N, signFactor)
#' summary(Test12_result)
#' }
#' @export

Test12 = function(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor, ...) {
  ParticleClass <- setRefClass(
    "ParticleClass",
    fields = list(
      position = "numeric",
      velocity = "numeric",
      best_position = "numeric",
      best_value = "numeric"
    ),
    where = globalenv(),
    methods = list(
      initialize = function(num_points, num_deletions, tposition = sample(1:num_points, num_deletions, replace = FALSE)) {
        .self$position <- tposition
        .self$velocity <- runif(num_deletions, -1, 1)
        .self$best_position <- .self$position
        .self$best_value <- Inf
      },
      update_velocity = function(global_best_position, w, c1, c2) {
        r1 <- runif(length(.self$position))
        r2 <- runif(length(.self$position))
        cognitive_component <- c1 * r1 * (.self$best_position - .self$position)
        social_component <- c2 * r2 * (global_best_position - .self$position)
        .self$velocity <- w * .self$velocity + cognitive_component + social_component
      },
      update_position = function(num_points) {
        .self$position <- round(.self$position + .self$velocity)
        .self$position[.self$position < 1] <- 1
        .self$position[.self$position > num_points] <- num_points
        .self$position <- unique(.self$position)
        if (length(.self$position) < max_final) {
          .self$position <- c(.self$position, sample(setdiff(1:num_points, .self$position), max_final - length(.self$position)))
        } else if (length(.self$position) > max_final) {
          .self$position <- sample(.self$position, max_final)
        }
      },
      evaluate = function(df, model, var_one, var_two) {
        current_value <- getValue(.self$position, df, model, var_one, var_two)
        if (current_value < .self$best_value) {
          .self$best_value <- current_value
          .self$best_position <- .self$position
        }
        return(current_value)
      }
    )
  )

  pso_minimize <- function(current_state = sample(1:num_points, num_deletions, replace = FALSE), num_deletions, num_particles = 10, max_iter = 30, w = 0.5, c1 = 1.5, c2 = 1.5) {
    num_points <- nrow(df)
    particles <- vector("list", num_particles)
    particles[[1]] <- ParticleClass$new(num_points, num_deletions, current_state)
    for (i in 2:num_particles) {
      particles[[i]] <- ParticleClass$new(num_points, num_deletions)
    }

    global_best_value <- Inf
    global_best_position <- NULL

    for (iter in 1:max_iter) {
      for (particle in particles) {
        current_value <- particle$evaluate(df, model, var_one, var_two)
        if (current_value < global_best_value) {
          global_best_value <- current_value
          global_best_position <- particle$position
        }
      }

      for (particle in particles) {
        particle$update_velocity(global_best_position, w, c1, c2)
        particle$update_position(num_points)
      }
    }

    return(list(best_position = global_best_position, best_value = global_best_value))
  }

  # Get initial index
  inflX <- replicate(N, 0)
  par_values <- replicate(N, 0)
  for (droppingP in 1:N) {
    par_values[droppingP] <- getValue(droppingP, df, model, var_one, var_two)
    inflX[droppingP] <- par_value - par_values[droppingP]
  }
  infl_sorted <- sort(inflX, decreasing = signFactor, index.return = TRUE)
  drops <- infl_sorted$ix[1:max_final]

  result <- pso_minimize(drops, num_deletions = max_final)
  pso_value <- result$best_value
  pso_drops <- result$best_position

  resultList <- list()
  class(resultList) <- "TestResult12"
  resultList$psoDrops <- pso_drops
  resultList$initialValue <- par_value
  resultList$finalValue <- pso_value
  resultList$methodname <- "Particle Swarm Optimization (PSO) to Switch Sign of Parameter"
  resultList$testindex <- 12
  resultList$PAR <- PAR
  resultList$threshold <- threshold
  resultList$N <- N
  resultList$max_final <- max_final
  resultList$par_value <- par_value

  return(resultList)
}



#' @export
summary.TestResult12 <- function(object, ...) {
  cat("Summary of Particle Swarm Optimization (PSO) Method to Switch Sign of Parameter Results:\n")
  if (!is.null(object$finalValue)) {

    cat(sprintf("Method Name: %s \n", object$methodname))
    cat(sprintf("Path: %s \n", object$PAR))
    cat(sprintf("Original Value: %f \n", object$initialValue))

    cat(sprintf("Original Number of Samples: %d \n", object$N))
    cat(sprintf("Drop Points Percentage: %d \n", object$threshold))
    cat(sprintf("Dropped Number of Samples: %d \n", length(object$psoDrops)))
    cat(sprintf("New Parameter Value: %f \n", object$finalValue))

    cat("PSO drop points list: \n")
    print(object$psoDrops)

    cat(sprintf("Using the PSO method, dropping %d data points, the parameter was switched from %.4f to %.4f\n", length(object$psoDrops), object$initialValue, object$finalValue))
  } else {
    cat("The PSO method was unable to switch the parameter sign.\n")
  }
}
summary <- function(object, ...) UseMethod("summary")
