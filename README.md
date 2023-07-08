
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Selective and robust external-control (srEC) integrative estimation

<!-- badges: start -->
<!-- badges: end -->

The srEC function is to develop an dynamic borrowing framework to
incorporate information from other external-control (EC) datasets with
the gold-standard randomized trials from [Gao et al.,
(2023)](https://arxiv.org/pdf/2306.16642.pdf). By adopting the
subject-level bias framework, the comparability of each EC subject is
assessed via a penalized estimation. The final integrative estimator
will only incorporate the external subjects with the estimated bias
being zero (i.e., a comparable subset).

## Installation with `devtools`:

``` r
devtools::install_github("Gaochenyin/SelectiveIntegrative")
```

## Main Paper:

C. Gao, S. Yang\*, M. Shan, W. Ye, I. Lipkovich, D. Faries (2023)
Integrating Randomized Placebo-Controlled Trial Data with External
Controls: A Semiparametric Approach with Selective Borrowing.
<https://arxiv.org/pdf/2306.16642>

## Example

``` r
library(SelectiveIntegrative)
set.seed(2333)
# generate data for the whole population
n_c.E <- 200; n_e.E <- 500
N <- n_c.E + n_e.E
X1 <- rnorm(N); X2 <- rnorm(N)
# omega: guage the unmeasured confounding for EC
omega <- 0.3
# generate the randomized trial population 
## generate the selection indicator for RT with expected sample size n_c.E
alpha0.opt <- uniroot(function(alpha0){
  mean(exp(alpha0 -2 * X1 - 2 * X2)/
  (1 + exp(alpha0 -2 * X1 - 2 * X2))) * N - n_c.E
}, interval = c(-50, 0))$root
eS <- exp(alpha0.opt -2 * X1 - 2 * X2)/
  (1 + exp(alpha0.opt -2 * X1 - 2 * X2))
delta <- rbinom(N, size = 1, prob = eS)
X.rt <- cbind(X1, X2)[delta == 1, ]
(n_c <- nrow(X.rt))
#> [1] 196
## generate the treatment assignment with marginal probability P.A
P.A <- 0.5
eta0.opt <- uniroot(function(eta0){
  mean(exp(eta0 - X.rt%*%c(1, 1))/
  (1 + exp(eta0 - X.rt%*%c(1, 1)))) - P.A
}, interval = c(-50, 0))$root
eA <- exp(eta0.opt - X.rt%*%c(1, 1))/
  (1 + exp(eta0.opt - X.rt%*%c(1, 1)))
A.rt <- rbinom(n_c, size = 1, prob = eA)
## generate the observed outcomes for RT
Y.rt <- as.vector(1 +  X.rt%*%c(1, 1) + A.rt * X.rt%*%c(.3, .3) + rnorm(n_c) + 
                    omega * rnorm(n_c)) # maintain a similar variation as the EC
data_rt <- list(X = X.rt, A = A.rt, Y = Y.rt)

# generate the external control population
X.ec <- cbind(X1, X2)[delta == 0, ]
(n_h <- nrow(X.ec))
#> [1] 504
A.ec <- 0
## generate the observed outcomes for EC (possibly confounded)
Y.ec <- as.vector(1 +  X.ec%*%c(1, 1) + omega * rnorm(n_h, mean = 1) + rnorm(n_h))
data_ec <- list(X = X.ec, A = A.ec, Y = Y.ec)
```

Now, we have generated the RT dataset `data_rt` and the EC dataset
`data_ec`. We are ready to implement our selective integrative
estimation by calling `srEC()`.

``` r
out <- srEC(data_rt = data_rt,
     data_ec = list(data_ec),
     method = 'gbm')
#> 载入需要的程辑包：ggplot2
#> 载入需要的程辑包：lattice
```

``` r
# AIPW
print(paste('AIPW: ', round(out$est$AIPW, 3), 
      ', S.E.: ', round(out$sd$AIPW, 3)))
#> [1] "AIPW:  -0.397 , S.E.:  2.158"
# ACW
print(paste('ACW: ', round(out$est$ACW, 3), 
      ', S.E.: ', round(out$sd$ACW, 3)))
#> [1] "ACW:  -0.831 , S.E.:  2.538"
# selective integrative estimation
print(paste('Our: ', round(out$est$ACW.final, 3), 
      ', S.E.: ', round(out$sd$ACW.final, 3)))
#> [1] "Our:  -0.377 , S.E.:  2.006"
```
