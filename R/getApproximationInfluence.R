#' @importFrom stats vcov
#' @noRd
#Perform calculations based on specified approximation method for specific data
getApproximationInfluence <- function(df, calcMeth, adjust, model, equalCons, PAR) {
  fit <- lavaan::sem(model, data=df)
  if(calcMeth == 'Hessian'){
    ##Using Hessian
    nabla <- lavaan::estfun.lavaan(fit, scaling = T, ignore.constraints = T) #Scaled empirical estimating function
    hess <- lavInspect(fit, "hessian") #Hessian matrix
    appr <- t(solve(hess[(equalCons + 1):dim(hess)[1], (equalCons + 1):dim(hess)[2]]) %*% t(nabla)) #calculate influence matrix
    approx_infl <- -1 * adjust*appr[,PAR]
  }else{
    #using covariance
    nabla2 <- lavaan::estfun.lavaan(fit, scaling = F, ignore.constraints = T) #Unscaled empirical estimating function
    approx_infl <- t(vcov(fit)%*%t(nabla2)) #Calculate influence matrix using variance-covariance matrix
    approx_infl <- approx_infl[, PAR] #select appropriate data for given variable of interest
    approx_infl<- adjust*approx_infl #Adjustment factor
  }
  return(approx_infl)
}
