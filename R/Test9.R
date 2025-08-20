#' @title Use Depth Method to Try to Switch Sign of Parameter
#' @description This function uses a depth method combined with a Negamax search algorithm to iteratively remove data points in order to switch the sign of a specific path in a Structural Equation Modeling (SEM) model.
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
#' @importFrom lavaan lavInspect
#' @importFrom lavaan parameterEstimates
#' @importFrom R.utils insert
#' @import dplyr
#' @import semfindr
#' @importFrom stats vcov



#' @return A list of class \code{TestResult9} containing:
#' \item{deletedPoints}{The indices of the most influential data points.}
#' \item{initialValue}{The original value of the parameter.}
#' \item{finalValue}{The final value of the parameter.}
#' \item{methodname}{The name of the method used.}
#' \item{testindex}{The index of the test performed.}
#' \item{r_squared}{R-squared value for the approximation, if applicable.}
#' \item{predictedReduction}{Predicted reduction in parameter, if applicable.}
#' \item{message}{Summary message of the results.}
#' \item{failureMessage}{Failure message if the sign switch was unsuccessful.}
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
#' Test9_result <- Test9(df, model, var_one, var_two, PAR, threshold, fit, estimates,
#' conc, int, par_value, max_final, N, signFactor, equalCons, calcMeth = "Hessian")
#' summary(Test9_result)
#' }
#' @export

Test9 = function(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth = "Hessian", ...){

  r_squared = 1
  resultsBreak = NULL

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


  targets <- -totalDropScore[1:length(drops8index)+1]
  targets[length(drops8index)+1] <- targets[length(drops8index)]
  target <- -par_value

  df_neg <- df

  #### Negamax search algorithm
  negamaxSearch <- function(scores, currentSubset, startLayer, currentSum) {

    # Helper function to check if time limit has exceeded

    currentTime <- proc.time()

    if (startLayer == maxLayer + 1 || (currentTime[3]- startTime[3] >= maxTime)) {#If no more layers to check, quit
      return(list(score = currentSum, subset = currentSubset))
    }

    #If we have already found we can get to 0 in a certain number of layers, we no longer need to search
    #any paths longer than that, so we set the maxLayer accordingly
    if (currentSum <= target){
      maxLayer <- startLayer - 2
      return(list(score = currentSum, subset = currentSubset))
    }

    bestScore <- Inf
    bestSubset <- NULL

    if(length(currentSubset) == 0)
    {
      df_neg <- df
    } else {
      df_neg  <- df[-c(currentSubset),] #update df to remove points
    }

    model_k <- lavaan::sem(model, data=df_neg) #update model given data removals

    # calculate according the econometrica paper
    k_sc <- lavaan::estfun.lavaan(model_k, scaling = T, ignore.constraints = T)
    k_asd <- lavInspect(model_k, "hessian")

    #print(currentSubset)
    BREAKc <- FALSE
    tryCatch({
      k_appr <- t(solve(k_asd[(equalCons + 1):dim(k_asd)[1], (equalCons + 1):dim(k_asd)[2]]) %*% t(k_sc))
    }, error = function(e) {
      # Code to be executed if an error occurs
      print("An error occurred in Negamax:")
      print(conditionMessage(e))
      BREAKc <- TRUE
      resultsBreak <- list(currentSum, currentSubset)
      return(NULL)
    })

    if(BREAKc == FALSE){
      k_appr <- t(solve(k_asd[(equalCons + 1):dim(k_asd)[1], (equalCons + 1):dim(k_asd)[2]]) %*% t(k_sc))
      layerscores <- k_appr[,PAR]
    } else{layerscores <- rep(0, nrow(df))
    maxLayer <- 1}

    if(length(currentSubset) != 0)
    {
      if(nrow(df) %in% currentSubset){
        layerscores <- append(layerscores, 0, after = length(layerscores))
      }

      dropSubset <- sort(currentSubset) #this all tries to add the right zeros into the vector
      for(counter in 1:length(currentSubset)){
        dropSubset[counter] <- dropSubset[counter] - counter + 1
      }

      if(length(dropSubset) != 0)
      {
        if(nrow(df) %in% dropSubset){
          layerscores <- append(layerscores, 0, after = length(layerscores))
        }
        if(0 %in% dropSubset){ #if we need to insert an empty piece into the zeroth spot
          layerscores <- append(layerscores, 0, after = 0)
          dropSubset <- dropSubset[dropSubset != 0]
        }
      }
      layerscores <- insert(layerscores, dropSubset, vector("numeric", length = length(dropSubset)))
    }
    sortedScores <- sort(layerscores, index.return = TRUE)
    topPrune <- sortedScores$ix[1:pruneNum]
    for (k in topPrune){ #1:length(layerscores)) {
      # print(k)
      if(!BREAKc){
        # print("HELLO")
        if(layerscores[k] + currentSum <= targets[startLayer] + alphab[startLayer]) {
          updatedSubset <- c(currentSubset, k) #check current subset in further depth
          # print("HERE")
          if (!any(sapply(exploredSubsets, identical, sort(updatedSubset)))) {
            exploredSubsets <- append(exploredSubsets, list(sort(updatedSubset))) #global variable
            updatedSum <- currentSum + layerscores[k] #update sum
            updatedRemaining <- startLayer + 1 #update which layer we are on
            result <- negamaxSearch(layerscores, updatedSubset, updatedRemaining, updatedSum)
            value <- result$score
            if (value < bestScore) {
              bestScore <- value
              bestSubset <- result$subset
            }
          }
        }
      }
    }
    return(list(score = bestScore, subset = bestSubset))
  }

  # usage
  scores <- replicate(N, 0) #initialize
  currentSubset <- c() #initialize
  remainingElements <- 1:length(scores)
  currentSum <- 0 #start with no score
  startLayer <- 1 #start on layer 1
  alphab <- -0.25*targets #Tuning coefficient to prune highly inaccurate branches
  exploredSubsets <- list() #initialize empty subset
  pruneNum <- 3 #number of new nodes to consider along each branch of search tree
  maxLayer <- nrow(df) - 1 #max number of layers to search
  maxTime <- 300 #max number of seconds to allow search
  # Start time
  startTime <- proc.time()
  results <- negamaxSearch(scores, currentSubset, startLayer, currentSum)

  #re-run model with drops
  if(!is.null(results$subset)){
    fitNegamax = lavaan::sem(model, data=df[-c(results$subset),])
    estsNegamax <- parameterEstimates(fitNegamax)
    intNegamax <- data.frame(c(estsNegamax$lhs), c(estsNegamax$rhs), val = c(estsNegamax$est))
    intNegamax2 <- intNegamax %>% filter_all(any_vars(. %in% c(var_one))) %>% filter_all(any_vars(. %in% c(var_two)))
    parValueNegamax = intNegamax2$val

    #Print results
    print(sprintf("The Negamax method reduced the parameter value from %.4f to %.4f", par_value, parValueNegamax))
    print(sprintf("%d points were dropped", length(results$subset)))

    if(parValueNegamax > 0.01*par_value){
      print(sprintf("Failure to switch the sign of the parameter may be due to inaccurate approximations."))
      print(sprintf("R-squared for the approximation is %.4f", r_squared))
    }
  } else{print("Negamax failed. Check for singularity errors.")
    print(sprintf("Best solution dropped %d points resulting in predicted %.2f reduction in parameter",length(resultsBreak[[2]]),resultsBreak[1]))}

  resultList <- list()
  class(resultList) <- "TestResult9"
  resultList$deletedPoints <- if (!is.null(results$subset)) results$subset else resultsBreak[[2]]
  resultList$initialValue <- par_value
  resultList$finalValue <- if (!is.null(results$subset)) parValueNegamax else NA
  resultList$methodname <- "Use depth method to try to switch sign of parameter"
  resultList$testindex <- 9
  resultList$r_squared <- if (parValueNegamax > 0.01 * par_value) r_squared else NA
  resultList$predictedReduction <- if (is.null(results$subset)) resultsBreak[1] else NA

  if (!is.null(results$subset)) {
    resultList$message <- sprintf("The Negamax method reduced the parameter value from %.4f to %.4f", par_value, parValueNegamax)
    resultList$droppedPoints <- length(results$subset)
    if (parValueNegamax > 0.01 * par_value) {
      resultList$failureMessage <- "Failure to switch the sign of the parameter may be due to inaccurate approximations."
    }
  } else {
    resultList$message <- "Negamax failed. Check for singularity errors."
    resultList$droppedPoints <- length(resultsBreak[[2]])
  }

  return(resultList)

}

#' @export
# Summary function for TestResult9
summary.TestResult9 <- function(object, ...) {
  cat("Summary of Depth Method to Switch Sign of Parameter Results:\n")
  cat(object$message, "\n")
  cat(sprintf("%d points were dropped\n", object$droppedPoints))
  if (!is.null(object$finalValue)) {
    # cat(sprintf("R-squared for the approximation is %.4f\n", object$r_squared))
    if (!is.null(object$failureMessage)) {
      cat(object$failureMessage, "\n")
    }
  } else {
    cat(sprintf("Best solution dropped %d points resulting in predicted %.2f reduction in parameter\n", length(object$deletedPoints), object$predictedReduction))
  }
  cat("Deleted points indices: ", paste(object$deletedPoints, collapse = ", "), "\n")
}
summary <- function(object, ...) UseMethod("summary")
