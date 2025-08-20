#' @title Run Specified Test for SEM Path Sign or Fit Measure Change
#' @description This function runs a specified test by method name or method index to determine either (1) when a specific path in a Structural Equation Modeling (SEM) model changes sign (positive or negative) by iteratively removing data points, or (2) whether a fit measure (e.g., CFI, TLI) exceeds or falls below a specified threshold (e.g., 0.9). The function outputs relevant results for the specified test.
#' @param df A data frame containing the dataset.
#' @param model A specified SEM model.
#' @param var_one The first variable of interest.
#' @param var_two The second variable of interest.
#' @param PAR The path of interest.
#' @param threshold The threshold for percentage of data dropped.
#' @param fit The fit measure to assess the model.
#' @param estimates The estimates from the SEM model.
#' @param conc The convergence criterion for the SEM model.
#' @param int The interval for checking the parameter.
#' @param par_value The original value of the parameter of interest.
#' @param max_final The maximum number of influential data points to consider.
#' @param N The total number of data points.
#' @param signFactor A factor indicating the direction of parameter change (positive or negative).
#' @param equalCons The equality constraint used in the SEM model (default is 0).
#' @param calcMeth The method used for approximation (default is 'Hessian').
#' @param ratio The ratio used to determine the number of points to check using the exact method (default is 2) in combined method.
#' @param adaptA Logical indicating whether to use adaptive pruning (default is FALSE).
#' @param alpha Manual tuning parameter (default is 0.25).
#' @param maxTime Maximum time allowed for the search in seconds (default is 300).
#' @param pruneNum The number of branches to explore from each node (default is 3).
#' @param measureTest The fit measurement name to be tested (e.g., "cfi", "tli", "rmsea"). Default is "cfi".
#' @param fitThreshold The threshold of the fit measurement. For example, for CFI (measureTest = "cfi"), the threshold could be 0.9. Default is 0.9.
#' @param highGood A boolean indicating whether higher values of the fit measure are better. For instance, for CFI, this should be set to TRUE. Default is TRUE.
#' @param method The method name or index specifying which test to run. Default is "Naive Method with Approximate Influence".
#' @param ... Other arguments.
#' \describe{
#'   \item{1: }{\code{"Naive Method with Exact Influence"}}
#'   \item{2: }{\code{"Naive Method with Approximate Influence"}}
#'   \item{3: }{\code{"Specified Approximation Method"}}
#'   \item{4: }{\code{"Simple Depth Method"}}
#'   \item{5: }{\code{"Combined Method"}}
#'   \item{7: }{\code{"Negamax Search Algorithm Function"}}
#'   \item{8: }{\code{"Use Depth Method to Try to Switch Sign of Parameter"}}
#'   \item{9: }{\code{"Use Depth Method to Try to Switch Sign of Parameter with Negamax"}}
#'   \item{11: }{\code{"Simulated Annealing Method"}}
#'   \item{12: }{\code{"Particle Swarm Optimization"}}
#'   \item{13: }{\code{"Brute Search with Cut Method"}}
#'   \item{10: }{\code{"Use depth method to drop fit measure below threshold"}}
#'   \item{61: }{\code{"Finding case deletions required to change fit metric using exact influences"}}
#'   \item{62: }{\code{"Finding case deletions required to change fit metric using approximate method"}}
#' }
#' @return A list containing the results of the specified test.
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
#' fit <- sem(model, data = df)
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
#' # Run specified test by method name
#' result_by_name <- SEMsensitivity(df, model, var_one, var_two, PAR, threshold, fit,
#' estimates, conc, int, par_value, max_final, N, signFactor,
#' "Naive Method with Exact Influence")
#' summary(result_by_name)
#'
#' # Run specified test by method index
#' result_by_index <- SEMsensitivity(df, model, var_one, var_two, PAR, threshold, fit,
#' estimates, conc, int, par_value, max_final, N, signFactor, 1)
#' summary(result_by_index)
#' }
#' @export

SEMsensitivity <- function(df, model, var_one = NULL, var_two = NULL, PAR = NULL, threshold = 10, fit, estimates = NULL, conc = NULL, int = NULL, par_value = NULL, max_final, N, signFactor = NULL, equalCons = 0, calcMeth = NULL, ratio = NULL, adaptA = NULL, alpha = NULL,
                           maxTime = NULL, pruneNum = NULL,
                           measureTest = "cfi", fitThreshold = 0.9, highGood = TRUE,
                           method = 2, ...) {

  # Match the method to the corresponding test function
  method_functions <- list(
    "1" = Test1,
    "2" = Test2,
    "3" = Test3,
    "4" = Test4,
    "5" = Test5,
    "7" = Test7,
    "8" = Test8,
    "9" = Test9,
    "11" = Test11,
    "12" = Test12,
    "13" = Test13,
    "10" = Test10,
    "61" = Test61,
    "62" = Test62
  )

  method_names <- list(
    "1" = "Naive Method with Exact Influence",
    "2" = "Naive Method with Approximate Influence",
    "3" = "Specified Approximation Method",
    "4" = "Simple Depth Method",
    "5" = "Combined Method",
    "7" = "Negamax Search Algorithm Function",
    "8" = "Use Depth Method to Try to Switch Sign of Parameter",
    "9" = "Use Depth Method to Try to Switch Sign of Parameter with Negamax",
    "10" = "Use depth method to drop fit measure below threshold",
    "11" = "Simulated Annealing Method",
    "12" = "Particle Swarm Optimization",
    "13" = "Brute Search with Cut Method",
    "61" = "Finding case deletions required to change fit metric using exact influences",
    "62" = "Finding case deletions required to change fit metric using approximate method"
  )

  # Handle numeric and character inputs for method
  if (is.numeric(method)) {
    method <- as.character(method)
    if (!method %in% names(method_functions)) {
      stop(paste("Invalid method index. Please provide a valid method index. Available options are:", paste(names(method_functions), collapse = ", ")))
    }
    method_name <- method_names[[method]]
  } else if (is.character(method)) {
    # Clean up the string and match it to method names
    method_clean <- tolower(trimws(method))
    valid_methods <- sapply(method_names, tolower)
    if (!(method_clean %in% valid_methods)) {
      stop("Invalid method name. Please provide a valid method name or method index.")
    }
    method <- names(valid_methods[valid_methods == method_clean])
    method_name <- method_names[[method]]
  } else {
    stop("Invalid method specification. Provide method name or method index.")
  }

  # Call the corresponding function with the provided arguments
  result <- do.call(method_functions[[method]], list(
    df = df, model = model, var_one = var_one, var_two = var_two, PAR = PAR,
    threshold = threshold, fit = fit, estimates = estimates, conc = conc, int = int,
    par_value = par_value, max_final = max_final, N = N, signFactor = signFactor,
    equalCons = equalCons, calcMeth = calcMeth, ratio = ratio, adaptA = adaptA, alpha = alpha,
    maxTime = maxTime, pruneNum = pruneNum,
    measureTest = measureTest, fitThreshold = fitThreshold, highGood = highGood
  ))

  return(result)
}
