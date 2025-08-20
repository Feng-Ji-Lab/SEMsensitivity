
#' @noRd
#Determine value of parameter for specific fit index
getValue <- function(drops_index,df, model,var_one, var_two) {
  options(warn=-1)
  df_drops <- df[-c(drops_index),]
  fit_drops = lavaan::sem(model, data=df_drops)
  new_estimates <- parameterEstimates(fit_drops)
  new_conc <- data.frame(c(new_estimates$lhs), c(new_estimates$rhs), val = c(new_estimates$est))
  new_int <- new_conc %>% filter_all(any_vars(. %in% c(var_one))) %>% filter_all(any_vars(. %in% c(var_two)))
  new_par_value = new_int$val
  return(new_par_value)
  # drops_index = sort(drops_index)
  #
  # return(.SEMClusterEnv$getValueWithMemory(drops_index,df, model,var_one, var_two))
}
