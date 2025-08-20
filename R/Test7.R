#' @title Negamax Search Algorithm Function
#' @description Utilizes a Negamax search algorithm to iteratively drop data points and update the SEM, ultimately aiming to minimize the parameter value of interest. Use method 7 - Negamax Search Algorithm Function.
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
#' @param adaptA Logical indicating whether to use adaptive pruning (default is FALSE).
#' @param alpha Manual tuning parameter (default is 0.25).
#' @param maxTime Maximum time allowed for the search in seconds (default is 300).
#' @param pruneNum The number of branches to explore from each node (default is 3).
#' @param ... Other arguments.
#' @importFrom lavaan sem
#' @import dplyr
#' @importFrom stats vcov
#' @import semfindr



#' @importFrom R.utils insert
#' @importFrom lavaan lavInspect
#' @importFrom lavaan parameterEstimates
#' @return A list of class \code{TestResult7} containing:
#' \item{value}{The final value of the parameter after dropping the influential points.}
#' \item{methodname}{The name of the method used.}
#' \item{testindex}{The index of the test performed.}
#' \item{deletedPoints}{The indices of the most influential data points.}
#' \item{initialValue}{The original value of the parameter.}
#' \item{finalValue}{The final value of the parameter.}
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
#'   Test7_result = Test7(df, model, var_one, var_two, PAR, threshold, fit, estimates,
#'   conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth = "Hessian",
#'   adaptA = F, alpha= 0.25, maxTime = 300, pruneNum = 3)
#' summary(Test7_result)
#' }
#' @export

Test7 = function(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth = "Hessian",adaptA = F, alpha= 0.25, maxTime = 300, pruneNum = 3, ...){

  #### Negamax search algorithm function
  negamaxSearch <- function(currentSubset, layers, startLayer, currentSum, targets) {

    currentTime <- proc.time() #Get current time

    if (startLayer == layers + 1 || (currentTime[3]- startTime[3] >= maxTime)) {#If we are at the final layer (i.e. we have dropped max_final points), don't check deeper
      return(list(score = currentSum, subset = currentSubset))
    }

    #Depending on whether we have dropped any points along this branch, create corresponding data frame accounting for drops
    if(length(currentSubset) == 0){
      df_neg <- df

    } else {
      df_neg  <- df[-c(currentSubset),] #update df to remove point
    }

    if(length(exploredSubsets) == 0){
      #Define defaults to be replaced by the best stored values if none have been searched yet
      bestScore <- 0
      bestSubset <- NULL
    }

    #Update model given data removals
    model_k <- lavaan::sem(model, data=df_neg)

    #Calculate influence scores for current model
    if(calcMeth == 'Hessian'){
      k_sc <- lavaan::estfun.lavaan(model_k, scaling = T, ignore.constraints = T, remove.duplicated = TRUE)
      k_asd <- lavInspect(model_k, "hessian")
      k_appr <- t(solve(k_asd[(equalCons + 1):dim(k_asd)[1], (equalCons + 1):dim(k_asd)[2]]) %*% t(k_sc))
      layerscores <- -k_appr[,PAR]
    }
    else{
      sc9 <- lavaan::estfun.lavaan(model_k, scaling = F, ignore.constraints = T)
      layerscores <- t(vcov(model_k)%*%t(sc9))
      layerscores <- layerscores[, PAR]
    }

    #The layerscore vector also has to always be length N for other code to work, so we add zeros in the positions of cases that
    #have already been deleted. This code also helps with keeping track of deletion indices.
    if(length(currentSubset) != 0)
    { #If last point was deleted, add a zero at the end of the vector
      if(nrow(df) %in% currentSubset){
        layerscores <- append(layerscores, 0, after = length(layerscores))
      }

      dropSubset <- sort(currentSubset) #create an array of all the points being dropped, sorted

      #Iterate over all dropped points, keeping track of how many smaller points have been dropped to adjust indices
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

    #Sort scores
    sortedScores <- sort(layerscores, decreasing = signFactor, index.return = TRUE)

    #Pick the candidate nodes to evaluate
    topPrune <- sortedScores$ix[1:pruneNum]

    #Iterate over prospective drop points
    for (k in topPrune){
      #If auto-tuning is enabled, skip checking further depths if appears impossible to match target score
      #If pruning is enabled, determine whether to skip based on pruning tolerance
      skipLayer <- 0
      if(adaptA == T){
        if(((layers - startLayer + 2)*sortedScores$x[1] + currentSum) < abs(targets[max_final]) ){
          skipLayer <- 1
          alphab[startLayer+1] <- Inf
        }
      }else {
        if(layerscores[k] + currentSum < abs(targets[startLayer+1]) - alphab[startLayer+1]){
          skipLayer <- 1}}

      if(skipLayer == 1){ #FOR DEBUGGING
        print("Layer Skipped")
      }

      #If not skipping layer due to pruning
      if(skipLayer == 0){
        updatedSubset <- c(currentSubset, k) #update dropped subset to include current point, and check this further
        if (!any(sapply(exploredSubsets, identical, sort(updatedSubset)))) {#If current subset has not been explored before, search in deeper depth
          exploredSubsets <- append(exploredSubsets, list(sort(updatedSubset))) #global variable tracking searched subsets. Update to include current.
          updatedSum <- currentSum + layerscores[k] #update predicted influence sum
          updatedRemaining <- startLayer + 1 #update which layer we will be on when we start a new search
          #Perform search stemming from this node
          result <- negamaxSearch(updatedSubset, layers, updatedRemaining, updatedSum, targets)
          value <- result$score #record score when branch reaches terminus
          if (value > bestScore) {#If this is the best branch explored, save it
            bestScore <- value
            bestSubset <- result$subset
          }
        }
      }
    }
    return(list(score = bestScore, subset = bestSubset)) #when all branches are exhausted, exit
  }


  #### Decide pruning method
  # adaptA = F #Whether to use adaptive pruning (True or False)


  df_depth <- df #Store new copy of data frame (since we will modify this)
    df_depth$index = 1:nrow(df) #Create an index column to help keep track of which original points dropped (since data frame changes each loop)
    drops_index <- replicate(max_final, 0) #Store indices for drops
    tallyDropScore <- rep(0, nrow(df)) #Used to generate targets if using negamax search

    #Iteratively drop the best point, recalculating influence after each drop
    for (j in 1:max_final){

      model_j <- lavaan::sem(model, data=df_depth) #update model given any data removals

      #Calculate influence based on chosen approximation method
      if(calcMeth == 'Hessian'){
        j_sc <- lavaan::estfun.lavaan(model_j, scaling = T, ignore.constraints = T)
        j_asd <- lavInspect(model_j, "hessian")
        j_appr <- t(solve(j_asd[(equalCons + 1):dim(j_asd)[1], (equalCons + 1):dim(j_asd)[2]]) %*% t(j_sc))
        j_approx_infl <- -j_appr[,PAR]
      } else{
        sc2 <- lavaan::estfun.lavaan(model_j, scaling = F, ignore.constraints = T)
        j_approx_infl <- t(vcov(model_j)%*%t(sc2))
        j_approx_infl <- j_approx_infl[, PAR]
      }

      #Sort influence
      j_approx_sorted <- sort(j_approx_infl, decreasing = signFactor, index.return = TRUE)

      #Create target set for negamax search (not necessary if only using simple depth)
      tallyDropScore[j + 1] <- tallyDropScore[j] + j_approx_sorted$x[1]

      #Drop point
      depth_drop_actual <- j_approx_sorted$ix[1] #The single data point that will be dropped (most influential)
      drops_index[j] <- df_depth$index[depth_drop_actual] #Store the position of the drop (relative to the indexing of the original data frame)

      df_depth  <- df_depth[-c(depth_drop_actual),] #update df to remove point
    }


  #Regardless of pruning technique, need to have performed Simple Depth search to create target array
  targets <- -tallyDropScore[1:max_final]
  targets[max_final+1] <- 0

  #Manual tuning
  # alpha <- 0.25 #Manual tuning parameter
  alphab <- alpha*targets*-1 #Create tuned targets

  #number of layers to search (i.e. points to drop)
  num <- max_final

  #Max time to allow (in seconds)
  # maxTime <- 300

  #Store dataframe
  df_neg <- df

  #Decide number of branches to explore from each node
  # pruneNum <- 3

  #Initialize Best Score and subset
  bestScore <- 0
  bestSubset <- NULL

  # Run algorithm
  currentSubset <- c() #initialize the starting node for the search
  currentSum <- 0 #initialize score for starting node
  startLayer <- 1 #initialize starting layer
  exploredSubsets <- list() #initialize empty subset for tracking explored branches

  #Get start time
  startTime <- proc.time()

  #Run algorithm
  results <- negamaxSearch(currentSubset, num, startLayer, currentSum, targets)

  #Process results
  if(!is.null(results$subset)){
    #re-run model with drops
    fitNegamax = lavaan::sem(model, data=df[-c(results$subset),])
    estsNegamax <- parameterEstimates(fitNegamax)
    intNegamax <- data.frame(c(estsNegamax$lhs), c(estsNegamax$rhs), val = c(estsNegamax$est))
    intNegamax2 <- intNegamax %>% filter_all(any_vars(. %in% c(var_one))) %>% filter_all(any_vars(. %in% c(var_two)))
    parValueNegamax = intNegamax2$val

    #Print results
    print(sprintf("The Negamax method reduced the parameter value from %.4f to %.4f", par_value, parValueNegamax))
  }else{print("Negamax incomplete due to singularity.")}



  resultList <- list()
  class(resultList) <- "TestResult7"
  resultList$value <- parValueNegamax
  resultList$methodname <- "Negamax search algorithm function"
  resultList$testindex <- 7
  resultList$deletedPoints <- results$subset
  resultList$initialValue <- par_value
  resultList$finalValue <- parValueNegamax

  resultList$PAR = PAR
  resultList$threshold = threshold
  resultList$N = N
  resultList$max_final = max_final
  resultList$par_value = par_value

  return(resultList)

}
#' @export
# Summary function for TestResult7
summary.TestResult7 <- function(object, ...) {
  cat("Summary of Negamax Search Algorithm Results:\n")
  if (!is.null(object$value)) {

    cat(sprintf("Method Name: %s \n",object$methodname))
    cat(sprintf("Path: %s \n",object$PAR))
    cat(sprintf("Original Value: %f \n",object$par_value))

    cat(sprintf("Original Number of Samples: %d \n",object$N))
    cat(sprintf("Drop Points Percentage: %d \n",object$threshold))
    cat(sprintf("Dropped Number of Samples: %d \n",object$max_final))
    cat(sprintf("New Parameter Value: %f \n",object$parValueNegamax))

    cat("drop points list: \n")
    print(object$deletedPoints)

    cat(sprintf("The Negamax method reduced the parameter value from %.4f to %.4f\n", object$initialValue, object$finalValue))
  } else {
    cat("Negamax incomplete due to singularity.\n")
  }
}

summary <- function(object, ...) UseMethod("summary")

