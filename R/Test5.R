#' @title Combined Method
#' @description This function determines a specific path in a Structural Equation Modeling (SEM) model value changing by removing samples iteratively, in which influences are determined by both naive method and approximate method and outputs relevant results. Use method 5 - Combined Method.
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
#' @param ratio The ratio used to determine the number of points to check using the exact method (default is 2).
#' @param ... Other arguments.
#' @importFrom lavaan sem
#' @importFrom lavaan lavInspect
#' @importFrom lavaan parameterEstimates
#' @import dplyr
#' @import semfindr
#' @importFrom stats vcov



#' @return A list of class \code{TestResult} containing:
#' \item{value}{The value of the parameter after dropping the influential points.}
#' \item{points}{The indices of the most influential data points.}
#' \item{methodname}{The name of the method used.}
#' \item{testindex}{The index of the test performed.}
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
#' Test5_result = Test5(df, model, var_one, var_two, PAR, threshold, fit, estimates,
#' conc, int, par_value, max_final, N, signFactor, equalCons,calcMeth)
#' summary(Test5_result)
#' }
#' @export

Test5 = function(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons=0,calcMeth = "Hessian",ratio = 2, ...){

  ##### Test 5 - Combined Method ##########
  new_final <- ratio*max_final #how many points to check using exact method (we use a ratio of 2)

  #Get approximate influence scores
  if(calcMeth == 'Hessian'){
    sc5 <- lavaan::estfun.lavaan(fit, scaling = T, ignore.constraints = T)
    asd5 <- lavInspect(fit, "hessian")
    appr5 <- t(solve(asd5[(equalCons + 1):dim(asd5)[1], (equalCons + 1):dim(asd5)[2]]) %*% t(sc5))
    approx_infl_com <- -appr5[,PAR]
  }else{
    sc3 <- lavaan::estfun.lavaan(fit, scaling = F, ignore.constraints = T)
    approx_infl_com <- t(vcov(fit)%*%t(sc3))
    approx_infl_com <- approx_infl_com[, PAR]
  }

  #Sort influence scores
  approx_sorted_com <- sort(approx_infl_com, decreasing = signFactor, index.return = TRUE)
  approx_drops_com <- approx_sorted_com$ix[1:new_final]

  inflXb = replicate(new_final, 0) #Create empty array to store influence values
  par_valuesb <- replicate(N, 0) #Create empty array to store parameter values
  countcom <- 1 #counting variable

  #Iterate over prospective case deletions and check exact influence
  for (droppingQ in approx_drops_com){
    rerun <- lavaan::sem(model, data <- df[-c(droppingQ),]) #Rerun model with deletion
    estimates <- parameterEstimates(rerun) #Obtrain estimates

    ###Determine The Value of The Parameter of Interest
    conc <- data.frame(c(estimates$lhs), c(estimates$rhs), c(estimates$est), val = c(estimates$est))
    int <- conc %>% filter_all(any_vars(. %in% c(var_one))) %>% filter_all(any_vars(. %in% c(var_two)))
    par_valuesb[droppingQ] <- int$c.estimates.est.

    #Calculate influence and store
    inflXb[countcom] <- par_value - par_valuesb[droppingQ]
    countcom <- countcom + 1
  }

  #Re-sort influence based on exact influence scores
  infl_exa_com_sorted <- sort(inflXb, decreasing = signFactor, index.return = TRUE)

  #Determine which cases to delete (based on original ordering of observations)
  com_drops <- infl_exa_com_sorted$ix[1:max_final]
  com_drops <- approx_drops_com[com_drops]

  #Delete cases and re-test model
  df_drops_com <- df[-c(com_drops),]
  fit_drops_com = lavaan::sem(model, data=df_drops_com)

  #Determine new parameter value
  estCom<- parameterEstimates(fit_drops_com)
  comInt =data.frame(c(estCom$lhs), c(estCom$rhs), val = c(estCom$est))
  comInt2 <- comInt %>% filter_all(any_vars(. %in% c(var_one))) %>% filter_all(any_vars(. %in% c(var_two)))
  com_par_val = comInt2$val


  print(sprintf("The combined method yields a new parameter value of %f", com_par_val))

  resultList = list()
  class(resultList) <- "TestResult"
  resultList$value = com_par_val
  resultList$points = com_drops
  resultList$methodname = "Combined Method"
  resultList$testindex = 5

  resultList$PAR = PAR
  resultList$threshold = threshold
  resultList$N = N
  resultList$max_final = max_final
  resultList$par_value = par_value

  return(resultList)

}
