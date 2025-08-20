# SEMsensitivity

This package performs sensitivity analysis for Structural Equation Modeling (SEM). It determines which sample points need to be removed for the sign of a specific path in the SEM model to change, thus assessing the robustness of the model.

## Installation

The R package ``SEMsensitivity`` is publicly available on [Github](https://github.com/Feng-Ji-Lab/SEMsensitivity). It could be installed and run by using the following commands:

``` r
devtools::install_github("Feng-Ji-Lab/SEMsensitivity")
library(SEMsensitivity)
```

## Example

We take the dataset `PoliticalDemocracy` provided in `lavaan` as an example. First, define a SEM and fit it with lavaan. Then, indicate the path of interest and the two variables of interest. Then, run the `SEMsensitivity` function and see the results. The possible method names are shown below:

- Naive Method with Exact Influence
- Naive Method with Approximate Influence
- Specified Approximation Method
- Simple Depth Method
- Combined Method
- Negamax Search Algorithm Function
- Use Depth Method to Try to Switch Sign of Parameter
- Use Depth Method to Try to Switch Sign of Parameter with Negamax

Below is a sample code. 

First, install all needed packages and import dataset. 

``` r
library(lavaan)
library(dplyr)
library(semfindr)
library(R.utils)
library(SEMsensitivity)

df <- PoliticalDemocracy
```

State the model (same as the model in `lavaan` package).

```r
model <- '
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  dem60 ~ ind60
  dem65 ~ ind60 + dem60
  y1 ~~ y5
  y2 ~~ y4 + y6
  y3 ~~ y7
  y4 ~~ y8
  y6 ~~ y8
'
```

State the parameters and path of interest. Use the `lavaan` package to calculate the initial path parameters and extract them. 

``` r
var_one <- 'dem65'
var_two <- 'ind60'
PAR <- c("dem65~ind60")
threshold <- 10

fit <- sem(model, data = df)
summary(fit)

estimates <- parameterEstimates(fit)

conc <- data.frame(lhs = estimates$lhs, rhs = estimates$rhs, est = estimates$est)
int <- conc %>% filter(lhs == var_one & rhs == var_two)
par_value <- int$est

max_final <- ceiling(threshold * nrow(df) / 100)
N <- nrow(df)

signFactor <- ifelse(par_value >= 0, TRUE, FALSE)
```

Call the function `SEMsensitivity` to do sensitivity analysis. The method can be specified by either method name or index.

``` r
result_by_name <- SEMsensitivity(df, model, var_one, var_two, PAR, threshold, fit,
estimates, conc, int, par_value, max_final, N, signFactor,
"Naive Method with Exact Influence")
summary(result_by_name)

result_by_index <- SEMsensitivity(df, model, var_one, var_two, PAR, threshold, fit,
estimates, conc, int, par_value, max_final, N, signFactor, 1)
summary(result_by_index)
```

After running this, we can see the result as shown below. It shows the path of interest, original value, changed parameter value and the dropped samples indices list. 

``` 
Summary of SEM sensitivity analysis result: 
Method Name: Naive exact method 
Path: dem65~ind60 
Original Value: 0.572336 
Original Number of Samples: 75 
Drop Points Percentage: 10 
Dropped Number of Samples: 8 
New Parameter Value: 0.076319 
drop points list: 
[1] 30 38  2 47 36 14 53 26
```

We also provide an example of fit indices changing from a good-fit model to an inadequate model. 

```r
fit_result <- SEMsensitivity(df, model, fit = fit, max_final = max_final, N = N, measureTest = "cfi", fitThreshold = 0.9, highGood = TRUE, method = 10)
summary(fit_result)
```

The result is shown below:

```
Method Name:  Use depth method to try to drop fit measure below threshold 
Original Fit Value:  0.9954 
Final Fit Value After Dropping Points:  0.8823 
Number of Data Points Dropped:  19 
Difference Between Original Fit and Threshold:  0.0954 
Cumulative Drop in Fit Value:  0.0963 
Dropped Points:  74, 41, 12, 69, 63, 45, 66, 36, 18, 56, 53, 2, 59, 53, 4, 12, 53, 16, 5 
```

It provides the original fit value, and the changed fit value after dropping points, number of points dropped, and a list of indices of dropped points list. 

For more information, please refer to our [manual](SEMsensitivity_0.1.0.pdf).