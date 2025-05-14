
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Selective and robust external-control (srEC) integrative estimation

<!-- badges: start -->
<!-- badges: end -->

The srEC package is to develop a data-adaptive borrowing framework to
incorporate external-controls (EC) with the trial data, which could
potentially facilitate the drug development. By adopting the
subject-level bias framework, the comparability of each EC subject is
assessed via a penalized estimation. The final integrative estimator
will only incorporate the ECs with the estimated bias being zero (i.e.,
a comparable subset).

## Installation with `devtools`:

``` r
devtools::install_github("Gaochenyin/SelectiveIntegrative")
```

## Example 1: data-adaptive borrowing for continuous outcomes

``` r
library(SelectiveIntegrative)
set.seed(2333)
# generate data for the whole population
n_c.E <- 3000; n_e.E <- 1000
N <- n_c.E + n_e.E
X1 <- rnorm(N); X2 <- rnorm(N)
# omega: guage the unmeasured confounding for EC
omega <- 3
# generate the randomized trial population 
## generate the selection indicator for RT with expected sample size n_c.E
alpha0.opt <- uniroot(function(alpha0){
  mean(exp(alpha0 -2 * X1 - 2 * X2)/
  (1 + exp(alpha0 -2 * X1 - 2 * X2))) * N - n_c.E
}, interval = c(-50, 10))$root
eS <- exp(alpha0.opt -2 * X1 - 2 * X2)/
  (1 + exp(alpha0.opt -2 * X1 - 2 * X2))
delta <- rbinom(N, size = 1, prob = eS)
X.rt <- cbind(1, X1, X2)[delta == 1, ]
(n_c <- nrow(X.rt))
#> [1] 2986
#> [1] 2986
## generate the treatment assignment with marginal probability P.A
P.A <- 2/3
eta0.opt <- uniroot(function(eta0){
  mean(exp(eta0 - X.rt%*%c(1, 1, 1))/
  (1 + exp(eta0 - X.rt%*%c(1, 1, 1)))) - P.A
}, interval = c(-50, 10))$root
eA <- exp(eta0.opt - X.rt%*%c(1, 1, 1))/
  (1 + exp(eta0.opt - X.rt%*%c(1, 1, 1)))
A.rt <- rbinom(n_c, size = 1, prob = eA)
## generate the observed outcomes for RT
Y.rt <- as.vector(1 +  X.rt%*%c(1, 1, 1) + A.rt * X.rt%*%c(0, .3, .3) + rnorm(n_c) + 
                    omega * rnorm(n_c)) # maintain a similar variation as the EC
data_rt <- list(X = X.rt, A = A.rt, Y = Y.rt)

# generate the external control population
X.ec <- cbind(1, X1, X2)[delta == 0, ]
(n_h <- nrow(X.ec))
#> [1] 1014
#> [1] 1014
A.ec <- 0
## generate the observed outcomes for EC (possibly confounded)
Y.ec <- as.vector(1 +  X.ec%*%c(1, 1, 1) + omega * rnorm(n_h, mean = 1) + rnorm(n_h))
data_ec <- list(X = X.ec, A = A.ec, Y = Y.ec)
```

Now, we have generated the RT dataset `data_rt` and the EC dataset
`data_ec`. We are ready to implement our selective integrative
estimation by calling `srEC()`.

``` r
out <- srEC(data_rt = data_rt,
     data_ec = list(data_ec),
     method = 'glm')
#> Loading required package: ggplot2
#> Loading required package: lattice
```

``` r
# AIPW
print(paste('AIPW: ', round(out$est$AIPW, 3), 
      ', S.E.: ', round(out$sd$AIPW, 3)))
#> [1] "AIPW:  -0.229 , S.E.:  0.161"
# ACW
print(paste('ACW: ', round(out$est$ACW, 3), 
      ', S.E.: ', round(out$sd$ACW, 3)))
#> [1] "ACW:  -0.52 , S.E.:  0.157"
# selective integrative estimation
print(paste('Our: ', round(out$est$ACW.final, 3), 
      ', S.E.: ', round(out$sd$ACW.final, 3)))
#> [1] "Our:  -0.235 , S.E.:  0.158"
```

## Example 2: data-adaptive borrowing for survival outcomes

``` r
data_rt <- data.frame(Y = pmin(event.timeR, cens.timeR),
                      D = as.numeric(event.timeR <= cens.timeR), 
                      X.R, A = A)
data_ec <- data.frame(Y = pmin(event.timeE, cens.timeE), 
                    D = as.numeric(event.timeE <= cens.timeE),
                    X.E)
head(data_rt); head(data_ec)
#>            Y D          X1          X2        X3 A
#> 1 0.21071985 0  0.04001744  0.31599294 0.3233869 1
#> 2 0.04700862 1 -1.12705449  0.23351994 0.3128803 0
#> 3 0.40752655 0  1.62962950 -1.16486405 0.8914390 1
#> 4 2.80010600 1 -2.06646239  0.28221082 2.3744909 1
#> 5 0.22524576 1 -1.28563288  0.83940768 0.8075524 1
#> 6 0.66835014 0  2.35225073  0.05560586 1.5876957 1
#>            Y D         X1          X2          X3
#> 1 0.34947475 1  1.6358954 -1.67439772 -1.04312329
#> 2 0.03582542 0 -1.9807189 -1.28329416  0.76766774
#> 3 0.24573391 1 -0.5391487 -0.61115066 -0.89693185
#> 4 0.07190597 1  1.1221320  0.06017221 -0.02517464
#> 5 1.01936194 0 -1.0958985  0.36200646 -0.40039740
#> 6 1.09258864 0  0.5181361  0.01479075 -1.27158364
```

Now, we have generated the RT dataset `data_rt` and the EC dataset
`data_ec` for the survival outcomes. We are ready to implement our
selective integrative estimation by calling `srEC_Surv()`.

``` r
out_Surv <- srEC_Surv(data_rt = data_rt,
                 data_ec = data_ec, 
                 k_grid = 1:4, 
                 X.pred.R = c('X1', 'X2', 'X3'),
                 X.pred.A = c('X1', 'X2', 'X3'),
                 X.pred.c = c('X1', 'X2', 'X3'),
                 X.pred.f = c('A', 'X1', 'X2', 'X3'))
```

``` r
# AIPW
out_Surv$aipcw.RCT
#> $est
#>      tau        S1           S0      RMST
#> [1,]   1 0.5426229  0.339123465 0.1017497
#> [2,]   2 0.3906198  0.022647286 0.3874857
#> [3,]   3 0.1626788  0.073588191 0.6160173
#> [4,]   4 0.1615468 -0.003297191 0.7429846
# ACW
out_Surv$aipcw
#> $est
#>      tau        S1         S0       RMST
#> [1,]   1 0.5426229 0.35309687 0.09476302
#> [2,]   2 0.3906198 0.10889752 0.33038719
#> [3,]   3 0.1626788 0.05829815 0.52343867
#> [4,]   4 0.1615468 0.04490456 0.63395013
# selective integrative estimation
out_Surv$aipcw_lasso
#> $est
#>      tau        S1         S0      RMST
#> [1,]   1 0.5426229 0.30719354 0.1177147
#> [2,]   2 0.3906198 0.07999399 0.3907423
#> [3,]   3 0.1626788 0.04636410 0.6042126
#> [4,]   4 0.1615468 0.03292925 0.7266787
```

## Referenced Papers:

1.  C. Gao, S. Yang\*, M. Shan, W. Ye, I. Lipkovich, D. Faries.
    Improving randomized controlled trial analysis with data-adaptive
    borrowing. *Biometrika* (2025).

2.  C. Gao, S. Yang\*, M. Shan, W. Ye, I. Lipkovich, D. Faries. Doubly
    protected estimation for survival outcomes utilizing external
    controls for randomized clinical trials. *International Conference
    on Machine Learning (ICML)* (2025).
