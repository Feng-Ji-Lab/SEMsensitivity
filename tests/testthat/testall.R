# Load packages
library(lavaan)
library(dplyr)
#library(tidyverse)
#library(gridExtra)
library(semfindr) #to install enter remotes::install_github("sfcheung/semfindr")
library(R.utils) #necessary for the negamax algorithm
#library(semPlot)
#library(mirt)
# library(simsem)

library(SEMsensitivity)


test_that('test-all', {
  ############################# USER INPUTS ######################

  #Import data
  df <- PoliticalDemocracy

  ###Build Model###
  model <- '
  # measurement model
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  # regressions
  dem60 ~ ind60
  dem65 ~ ind60 + dem60
  # residual correlations
  y1 ~~ y5
  y2 ~~ y4 + y6
  y3 ~~ y7
  y4 ~~ y8
  y6 ~~ y8
  '
  equalCons <- 0 #specify number of equality constraints in model (for calculations)
  var_one = 'dem65' #first term
  var_two = 'ind60' #second term
  PAR <- c("dem65~ind60") #full relation
  #################### SELECT INFLUENCE ANALYSIS PARAMETERS
  ###CHANGE ME - What parameter you are interested in###
  test1 <- T #Exact method
  test2 <- T #Approximate method
  test3 <- T #Sign switch test
  test4 <- T #Depth method
  test5 <- T #Combined method
  test6 <- T #Fit index tests
  test7 <- T #Pruned Negamax method
  test8 <- T #Depth method for switching sign
  test9 <- T #Negamax method for switching sign
  test10 <- F #Simple Depth for Fit Metrics
  test11 <- F #anneal search method
  test12 <- F #pso method
  test13 <- F #brute force with cut method

  ##Change me - 'Hessian' or 'Covariance' Approximation Method
  calcMeth <- 'Hessian'

  ###CHANGE ME - percentage of data points to drop
  threshold <- 10


  #################################################################################

  #Fit SEM model
  fit <- lavaan::sem(model, data=df)
  summary(fit)

  ###Get Estimates of Parameters from SEM###
  estimates <- parameterEstimates(fit)

  ###Determine The Value of The Parameter of Interest
  conc <- data.frame(c(estimates$lhs), c(estimates$rhs), c(estimates$est), val = c(estimates$est))
  int <- conc %>% filter_all(any_vars(. %in% c(var_one))) %>% filter_all(any_vars(. %in% c(var_two)))
  par_value <- int$c.estimates.est. #this is the value of the parameter of interest

  ###Compute max number of points to be dropped###
  max_final <- ceiling(threshold*nrow(df)/100) #perform rounding if necessary
  N = nrow(df) #store number of observations in df for convenience

  #Determine whether parameter is negative or positive in order
  #to assess which direction to perturb it
  if (par_value >= 0){
    signFactor <- T } else {signFactor <- F
    }
  # .SEMClusterEnv$getValueWithMemory = myMemorize(getValueInternal)

  ## BEGIN TESTS ##
  ###TEST 13 - brute force with cut method ###################################################
  if (exists('test13') && test13 == T){

    # 初始化最好的组合
    initialize_best_combinations <- function() {
      best_combinations <- list()
      best_values <- numeric()
      return(list(best_combinations = best_combinations, best_values = best_values))
    }

    # 更新最好的组合
    update_best_combinations <- function(best_combinations, best_values, new_combination, new_value, max_best) {
      if (length(best_values) < max_best || new_value < max(best_values)) {
        if (length(best_values) >= max_best) {
          max_index <- which.max(best_values)
          best_combinations[[max_index]] <- new_combination
          best_values[max_index] <- new_value
        } else {
          best_combinations <- c(best_combinations, list(new_combination))
          best_values <- c(best_values, new_value)
        }
      }
      return(list(best_combinations = best_combinations, best_values = best_values))
    }

    # 逐步生成组合并保留最好的组合
    generate_combinations <- function(n, k, max_best) {
      best <- initialize_best_combinations()

      # 初始化一位数组合
      for (i in 1:n) {
        value <- getValue(i)
        best <- update_best_combinations(best$best_combinations, best$best_values, c(i), value, max_best)
      }

      # 逐步增加组合的长度
      for (depth in 2:k) {
        new_best_combinations <- list()
        new_best_values <- numeric()
        for (combo in best$best_combinations) {
          for (i in setdiff(1:n, combo)) {
            new_combination <- c(combo, i)
            value <- getValue(new_combination)
            new_best <- update_best_combinations(new_best_combinations, new_best_values, new_combination, value, max_best)
            new_best_combinations <- new_best$best_combinations
            new_best_values <- new_best$best_values
          }
        }
        best$best_combinations <- new_best_combinations
        best$best_values <- new_best_values
      }

      return(best)
    }
    result <- generate_combinations(nrow(df), max_final, 5)
    brute_search_value = result$best_value
    brute_search_drops = result$best_position
    ###Print Results###
    print(sprintf("The brute search with cut method yields a new parameter value of %f", brute_search_value))
  }

  ###TEST 12 - pso method ###################################################
  if (exists('test12') && test12 == T){

    # 粒子类
    Particle <- setRefClass(
      "Particle",
      fields = list(
        position = "numeric",
        velocity = "numeric",
        best_position = "numeric",
        best_value = "numeric"
      ),
      methods = list(
        initialize = function(num_points, num_deletions, tposition = sample(1:num_points, num_deletions, replace = FALSE)) {
          position <<- tposition
          velocity <<- runif(num_deletions, -1, 1)
          best_position <<- position
          best_value <<- Inf
        },
        update_velocity = function(global_best_position, w, c1, c2) {
          r1 <- runif(length(position))
          r2 <- runif(length(position))
          cognitive_component <- c1 * r1 * (best_position - position)
          social_component <- c2 * r2 * (global_best_position - position)
          velocity <<- w * velocity + cognitive_component + social_component
        },
        update_position = function(num_points) {
          position <<- round(position + velocity)
          position[position < 1] <<- 1
          position[position > num_points] <<- num_points
          position <<- unique(position)
          if (length(position) < max_final) {
            position <<- c(position, sample(setdiff(1:num_points, position), max_final - length(position)))
          } else if (length(position) > max_final) {
            position <<- sample(position, max_final)
          }
        },
        evaluate = function(df, model, var_one, var_two) {
          current_value <- getValue(position)
          if (current_value < best_value) {
            best_value <<- current_value
            best_position <<- position
          }
          return(current_value)
        }
      )
    )

    # PSO算法
    pso_minimize <- function(current_state = sample(1:num_points, num_deletions, replace = FALSE), num_deletions, num_particles = 10, max_iter = 30, w = 0.5, c1 = 1.5, c2 = 1.5) {
      num_points <- nrow(df)
      particles <- vector("list", num_particles)
      particles[[1]] = Particle$new(num_points, num_deletions, current_state)
      for (i in 2:num_particles) {
        particles[[i]] <- Particle$new(num_points, num_deletions)
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

    # get init index
    inflX <- replicate(N, 0)
    par_values <- replicate(N, 0)
    for (droppingP in 1:N){
      par_values[droppingP] <- getValue(droppingP)
      inflX[droppingP] <- par_value -  par_values[droppingP]
    }
    infl_sorted <- sort(inflX, decreasing = signFactor, index.return = TRUE)
    drops <- infl_sorted$ix[1:max_final]

    result <- pso_minimize(drops, num_deletions = max_final)
    pso_value = result$best_value
    pso_drops = result$best_position
    ###Print Results###
    print(sprintf("The pso method yields a new parameter value of %f", pso_value))
  }

  ###TEST 11 - anneal search method ###################################################
  if (exists('test11') && test11 == T){
    # 模拟退火算法
    simulated_annealing <- function(num_deletions, current_state = sample(1:num_points, num_deletions, replace = FALSE), initial_temp = 100, cooling_rate = 0.03, min_temp = 0.1, max_iter = 100) {
      num_points <- nrow(df)

      # 初始化
      current_value <- getValue(current_state)
      best_state <- current_state
      best_value <- current_value

      temp <- initial_temp

      for (iter in 1:max_iter) {
        if (temp < min_temp) break
        # 生成邻居状态
        new_state <- current_state
        add_index <- sample(setdiff(1:num_points, current_state), 1)
        remove_index <- sample(current_state, 1)
        new_state[new_state == remove_index] <- add_index

        new_value <- getValue(new_state)

        # 接受准则
        if (new_value < current_value || runif(1) < exp((current_value - new_value) / temp)) {
          current_state <- new_state
          current_value <- new_value
        }

        # 更新最佳状态
        if (current_value < best_value) {
          best_state <- current_state
          best_value <- current_value
        }

        # 降低温度
        temp <- temp * (1 - cooling_rate)
      }

      return(list(best_state = best_state, best_value = best_value))
    }

    # get init index
    inflX <- replicate(N, 0)
    par_values <- replicate(N, 0)
    for (droppingP in 1:N){
      par_values[droppingP] <- getValue(droppingP)
      inflX[droppingP] <- par_value -  par_values[droppingP]
    }
    infl_sorted <- sort(inflX, decreasing = signFactor, index.return = TRUE)
    drops <- infl_sorted$ix[1:max_final]

    result <- simulated_annealing(threshold, drops)
    annealing_drops <- result$best_state
    annealing_value <- result$best_value
    print(sprintf("The simulated annealing method yields a new parameter value of %f", annealing_value))
  }
  ###TEST 1 - Naive Method with Exact Influence ###################################################
  if (exists('test1') && test1 == T){
    Test1_result = Test1(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor)
    summary(Test1_result)
  }


  ######## Test 2: Naive method with approximate influence #######

  if(exists('test2') && test2 == T){
    Test2_result = Test2(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons=0, calcMeth = 'Hessian')
    summary(Test2_result)

  }


  ###TEST 3 - Drop points until sign of parameter is switched used naive approximate method ########################

  if(exists('test3') && test3 == T){
    Test3_result = Test3(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth)
    summary(Test3_result)

  }

  ### TEST 4: Simple DEPTH METHOD ####

  if(exists('test4') && test4 == T){
    Test4_result = Test4(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor)
    summary(Test4_result)

  }


  if(exists('test5') && test5 == T){
    Test5_result = Test5(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor)
    summary(Test5_result)

  }


  ### Test 6 - Finding case deletions required to change fit metric###
  ## REQUIRES semfindr Package ##
  if(exists('test6') && test6 == T){

    Test61_result = Test61(df, model, threshold, fit, max_final, N, measureTest = "cfi", fitThreshold = 0.9, highGood = T)
    summary(Test61_result)


    ##6b Using appx scores ##
    Test62_result = Test62(df, model, threshold, fit, max_final, N, measureTest = "cfi", fitThreshold = 0.9, highGood = T)
    summary(Test62_result)
  }



  ################## Test 7 - Negamax Search Using Approx. Scores ######

  if(exists('test7') && test7 == T){

    Test7_result = Test7(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth = "Hessian",adaptA = F, alpha= 0.25, maxTime = 300, pruneNum = 3)
    summary(Test7_result)
  }





  ## Test 8: Use depth method to try to switch sign of parameter ####
  if (exists('test8') && test8 == T){
    Test8_result = Test8(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth = "Hessian")
    summary(Test8_result)
  }




  ################## Test 9 - Using Negamax to try to switch sign of parameter################

  ## NOTE: as of right now, this only uses approximate influence scores calculated via the Hessian (because of some singularity issues)
  #There is also no option for automatic tuning
  #And it requires test8 to be run to generate target data

  if(exists('test9') && test9 == T){
    Test9_result = Test9(df, model, var_one, var_two, PAR, threshold, fit, estimates, conc, int, par_value, max_final, N, signFactor,equalCons,calcMeth = "Hessian")
    summary(Test9_result)


  }

  ###### COULD DEPTH METHOD USE EXACT COMPUTATION FOR LAST STEP IF REACHES A SINGULAR CALCULATION????
  #TRY THIS OUT!!


  ###### DEPTH METHOD FOR FIT MEASURES
  ## Test 10: Use depth method to try to drop fit measure below threshold ####

  if (exists('test10') && test10 == T){
    Test10_result = Test10(df, model, fit, max_final, N, depthmeasureTest = "cfi", depthfitThreshold = 0.9)
    summary(Test10_result)

  }









})






