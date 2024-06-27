#' Selective and robust external-control borrowing strategy for evaluating treatment effect
#'
#' @description
#' `srEC()` is a penalized dynamic integrative framework to augment a randomized trial (RT)
#'  with a external control (EC) dataset, in which the subject-level compatibility of the EC
#'  is assessed by a well-crafted penalty (e.g., adaptive lasso penalty). The parameter of interest
#'  is the average treatment effect.
#'
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom caret train
#' @param data_rt A list contains X, Y and A for RT. The propensity scores P(A=1|X) can also be
#' contained as `prob_A` (or `NULL`).
#' @param data_ec A list contains X and Y for EC with treatment A=0.
#' @param rt.ctrl An object of [caret::trainControl()] for [caret::train()], which controls the
#' computational nuances for fitting the outcome model using the RT dataset.
#' @param hc.ctrl An object of [caret::trainControl()], which controls the
#' computational nuances for fitting the outcome model using the EC dataset.
#' @param ... Other options used to fit the predictive models. Passed on to [caret::train()].
#'
#' @returns A list with components:
#' * est: estimated average treatment effect by AIPW, ACW and the selective integrative estimator.
#' * sd: estimated standard errors for the aforementioned estimators.
#' * subset.idx: a subset of indices of the external controls which have been selected for the
#' final integrative estimation.
#' @export
srEC <- function(data_rt,
                 data_ec = NULL,
                 rt.ctrl = caret::trainControl(method = 'cv', number = 10),
                 hc.ctrl = caret::trainControl(method = 'cv', number = 10),
                 method = 'gbm', ...)
{
  # assign the variables for analyses
  X_c <- data_rt$X
  Y_c <- data_rt$Y; n_c <- length(Y_c)
  A_c <- data_rt$A
  prob_A <- data_rt$prob_A
  # data pre-processing
  if(is.data.frame(X_c)){
    # construct the basis function for prediction and calibration
    X_c <- model.matrix(~., X_c) # for calibration
  }else{X_c <- as.matrix(X_c)}
  Y0_c <- Y_c[which(A_c==0)]; X0_c <- X_c[which(A_c==0),, drop = F]
  Y1_c <- Y_c[which(A_c==1)]; X1_c <- X_c[which(A_c==1),, drop = F]
  # estimate the outcome model for RT
  fit.Y.RCT <- caret::train(Y ~ 0 + . + A * (.),
                     data = data.frame(Y = Y_c, A = A_c,
                                       X_c),
                     trControl = rt.ctrl,
                     method = method, ...)
  # estimate propensity score (if not provided)
  if(is.null(prob_A))
  {
    glm.ps <- glm(A_c ~ X_c, family = binomial())
    prob_A <- predict(glm.ps, type = 'response')
  }
  # doubly robust variance estimation
  mu1_AIPW_i <- {A_c * Y_c/prob_A +
      (prob_A - A_c)/prob_A * predict(fit.Y.RCT, newdata = data.frame(X_c, A = 1))}
  mu0_AIPW_i <- {(1-A_c) * Y_c/(1-prob_A) +
      (A_c - prob_A)/(1-prob_A) * predict(fit.Y.RCT, newdata = data.frame(X_c, A = 0))}
  tau_0_i <-  mu1_AIPW_i - mu0_AIPW_i
  tau_hat_AIPW <- sum(tau_0_i)/n_c

  # compute estimated variances for mu1_AIPW - mu0_ACW
  var_mu1_AIPW_hat <- sum({mu1_AIPW_i-sum(mu1_AIPW_i)/n_c}^2)*n_c^{-1}
  var_mu0_AIPW_hat <- sum({mu0_AIPW_i-sum(mu0_AIPW_i)/n_c}^2)*n_c^{-1}
  # compute estimated standard error for mu1_AIPW-mu0_AIPW
  sd_AIPW_hat <- sqrt(sum({tau_0_i - tau_hat_AIPW}^2)*n_c^{-1})
  sd_AIPW_hat


  # if no external controls are provided
  if(is.null(data_ec)){

    tau_hat_AIPW <- sum(tau_0_i)/n_c
    tau_upper_initial <- tau_hat_AIPW + 1.96 * sd_AIPW_hat/sqrt(n_c)
    tau_lower_initial <- tau_hat_AIPW - 1.96 * sd_AIPW_hat/sqrt(n_c)

    return(list(n_c = n_c,
                est = list(AIPW = tau_hat_AIPW,
                           ACW = tau_hat_AIPW,
                           ACW.lasso = tau_hat_AIPW),
                sd = list(AIPW = sd_AIPW_hat,
                          ACW = sd_AIPW_hat,
                          ACW.lasso = sd_AIPW_hat),
                CI = list(lower = tau_lower_initial,
                          upper = tau_upper_initial),
                subset.idx = integer(0)))

  }

  # detect the number of external controls
  K <- length(data_ec)

  # the function for the augmented calibration weighting estimator
  EST.ACW.FUN <- function(X_c, Y_c, A_c, prob_A,
                          data_ec,
                          bias_b_list = list(NULL),
                          ...){

    # obtain the ACW estimator for each external control dataset
    data_ec_i <- mapply(function(ec_i, bias_b){
      estimate_hc(X_c = X_c, X_h = ec_i$X,
                  Y_c = Y_c, Y_h = ec_i$Y,
                  A_c = A_c, prob_A = prob_A,
                  fit.Y.RCT = fit.Y.RCT,
                  hc.ctrl = hc.ctrl,
                  method = method, bias_b = bias_b, ...
      )
    }, ec_i = data_ec, bias_b = bias_b_list,
    SIMPLIFY = FALSE)

    # obtain the number of each external control dataset
    n_ei <- sapply(data_ec_i, function(ec_i)nrow(ec_i$X))

    # compute the estimated standard errors for mu1_ACW-mu0_ACW
    sd_ACW_hat_i <- sapply(data_ec_i, function(ec_i)
    {
      tau.ACW.lin <- c(mu1_AIPW_i,
                       rep(0, nrow(ec_i$X))) - ec_i$mu0_score # influence function
      sqrt(sum({tau.ACW.lin-
          c(rep(sum(tau.ACW.lin, na.rm = TRUE)/n_c, n_c),
            rep(0, nrow(ec_i$X)))}**2, na.rm = TRUE)*n_c^{-1})
    })
    sd_ACW_hat_i

    # construct consistent estimator for the subject-level bias b_k
    tau_ec_i <- lapply(data_ec_i, function(ec_i)c(mu1_AIPW_i,
                                                  rep(0, nrow(ec_i$X)))- ec_i$mu0_i)

    tau_ec_score <- lapply(data_ec_i, function(ec_i)c(mu1_AIPW_i,
                                                      rep(0, nrow(ec_i$X)))- ec_i$mu0_score)

    return(list(EST.ACW = data_ec_i, n_ei = n_ei,
                sd_ACW_hat_i = sd_ACW_hat_i,
                tau_ec_i = tau_ec_i,
                tau_ec_score = tau_ec_score))
  }

  EST.ACW.res <- EST.ACW.FUN(X_c = X_c, Y_c = Y_c,
                             A_c = A_c, prob_A = prob_A,
                             data_ec = data_ec)
  # compute estimate for tau by AIPW (RT-only)
  tau_hat_AIPW <- sum(tau_0_i)/n_c
  # compute estimate for tau by ACW (RT and EC)
  tau_hat_ACW <- sapply(EST.ACW.res$tau_ec_i,
                        function(tau_i)sum(tau_i, na.rm = TRUE)/n_c)

  # obtain the subject-level bias estimates
  bias_h <- lapply(EST.ACW.res$EST.ACW,
                   function(x)as.vector(x$bias_i))
  # compute the standard error of the bias based on IF
  bias_h.score <- lapply(EST.ACW.res$EST.ACW,
                         function(x)as.vector(x$bias_i.score))
  sd.bias_h <- lapply(EST.ACW.res$EST.ACW,
                      function(x)sqrt(mean({x$bias_i.score-
                          mean(x$bias_i.score)}**2)))
  bias_h <- mapply(function(x, y)x, x = bias_h)
  # subject-level selective borrowing framework
  Z <- factor(rep(1:sum(EST.ACW.res$n_ei)))
  Z.mat <- model.matrix(~Z+0)
  # obtain the OLS estimates for adaptive lasso
  gamma.hat <- unlist(lapply(EST.ACW.res$EST.ACW, function(x)x$gamma.hat))
  # tuning the parameters (omega and lambda)
  nu.vector <- c(1,2)#rev(seq(1, 2, length.out = 10)) with length 10 or 50
  cv.lasso.vec <- sapply(nu.vector, function(nu){

    # specify the lambda vector for cross-validation
    lambda.vector <- seq(nrow(Z.mat)**-(1+nu/5)/2,
                         nrow(Z.mat)**(-0.1),
                         length.out = 100)
    # construct the adaptive weights (subject-level)
    b_k.w <- c(1/abs(gamma.hat)**nu)

    # the penalty factors are re-scaled with sum equal to `nvars`
    b_k.w <- b_k.w/sum(b_k.w)*ncol(Z.mat)
    # replace any NA number
    b_k.w <- sapply(b_k.w, function(x)ifelse(is.na(x), Inf, x))

    # fit the adaptive lasso penalized estimates for the subject-level bias
    fitLasso <- glmnet(Z.mat, bias_h, family = 'gaussian',
                       alpha = 1,
                       penalty.factor=b_k.w,
                       intercept = FALSE,
                       standardize = FALSE,
                       standardize.response = FALSE,
                       lambda = lambda.vector)

    group_indicator <- rep(1:K, times = EST.ACW.res$n_ei)

    trade.off <- apply(predict(fitLasso, newx = Z.mat, type = 'coef'), 2,
                       function(x)
                       {
                         tau.ACW.lin <- unlist(bias_h)*
                           c(x[-1]==0)
                         bias.hi <- tapply(tau.ACW.lin, group_indicator,
                                           function(x)sum(x)/n_c) -0
                         sum(bias.hi**2)/n_c
                       }) + # bias
      apply(predict(fitLasso, newx = Z.mat, type = 'coef'), 2,
            function(x){
              tau.ACW.score.all <- unlist(lapply(bias_h.score, function(x)x[-(1:n_c)]))
              tau.ACW.lin <- tau.ACW.score.all*c(x[-1]==0)
              var.hi <- tapply(tau.ACW.lin, group_indicator,
                               function(x) sum(x-sum(x)/n_c)**2/n_c)

              tau.ACW.score.RCT <- lapply(bias_h.score, function(x)x[(1:n_c)])
              var.ci <- unlist(lapply(tau.ACW.score.RCT,
                                      function(x)sum(x-sum(x)/n_c)**2/n_c))

              mean(var.hi) + mean(var.ci)
            }) # variance
    min(trade.off)
  })
  # select the one that minimize the cv error
  nu.opt <- nu.vector[which.min(cv.lasso.vec)]
  b_k.OLS <- c(1/abs(gamma.hat)**nu.opt)
  b_k.w <- b_k.OLS/sum(b_k.OLS)*ncol(Z.mat)
  b_k.w <- sapply(b_k.w, function(x)ifelse(is.na(x), 1e10, x))

  # nu.opt/3 due to a smaller order of convergence rate of the machine learning models
  lambda.vector <- c(seq(nrow(Z.mat)**-(-.001+nu.opt/3)/2, # change to 3 for conservative convergence rate
                         nrow(Z.mat)**(-0.001),
                         length.out = 100))

  fitLasso <- glmnet(Z.mat, bias_h, family = 'gaussian',
                     alpha = 1,
                     penalty.factor=b_k.w,
                     intercept = FALSE,
                     standardize = FALSE,
                     standardize.response = FALSE,
                     lambda = lambda.vector)

  group_indicator <- rep(1:K, times = EST.ACW.res$n_ei)
  coef.mat <- predict(fitLasso, newx = Z.mat.tilde, type = 'coef')

  trade.off <- apply(coef.mat, 2,
                     function(x)
                     {
                       # bias
                       tau.ACW.lin <- unlist(bias_h)* c(x[-1]==0)
                       bias.hi <- tapply(tau.ACW.lin, group_indicator,
                                         function(x)sum(x)/(n_c))-#/(n_c+length(x))) -
                         0
                       bias.hi <- bias.hi**2
                       tau.ACW.score.all <- unlist(lapply(bias_h.score, function(x)x[-(1:n_c)]))
                       tau.ACW.lin <- tau.ACW.score.all*c(x[-1]==0)
                       tau.ACW.lin.list <- structure(split(tau.ACW.lin,
                                                           f = factor(rep(1:K,
                                                                          times = EST.ACW.res$n_ei))),
                                                     names = NULL)
                       tau.ACW.score.recont <- mapply(function(x,y) c(x[1:n_c],
                                                                      y),
                                                      x = bias_h.score, y = tau.ACW.lin.list,
                                                      SIMPLIFY = FALSE)
                       # variance
                       var.i <- unlist(mapply(function(x,y){
                         sum((x - mean(x))**2)/n_c
                       }, x = tau.ACW.score.recont,
                       y = EST.ACW.res$n_ei, SIMPLIFY = FALSE))
                       mean(var.i) + sum(bias.hi)
                     })

  lambda.opt <- lambda.vector[which.min(trade.off)]
  lambda.opt.idx <- which.min(trade.off)

  # obtain the penalized estimates for the subject-level bias
  b_k_lasso <- coef.mat[-1, lambda.opt.idx] # delete the intercept
  b_k_lasso
  # obtain the indices for borrowing
  hc.idx.lasso <- unname(which(b_k_lasso==0))
  sd_ACW_hat_i <- EST.ACW.res$sd_ACW_hat_i


  if(identical(integer(0), hc.idx.lasso)|identical(numeric(0), hc.idx.lasso)){
    # if no external controls are selected
    tau_hat_AIPW <- tau_final <- sum(tau_0_i)/n_c
    tau_hat_ACW <- sapply(EST.ACW.res$tau_hc_i,
                          function(tau_i)sum(tau_i)/n_c)
    sd_final <- sd_AIPW_hat
    Sigma2_tau <- sd_AIPW_hat**2
    tau_cont_upper <- tau_final +
      qnorm(1-0.025) * sd_final/sqrt(n_c)
    tau_cont_lower <- tau_final-
      qnorm(1-0.025) * sd_final/sqrt(n_c)
    return(list(n_c = n_c,
                est = list(AIPW = tau_hat_AIPW,
                           ACW = tau_hat_ACW,
                           ACW.lasso = tau_final),
                sd = list(AIPW = sd_AIPW_hat,
                          ACW = sd_ACW_hat_i,
                          ACW.lasso = sd_final),
                CI = list(lower = tau_cont_lower,
                          upper = tau_cont_upper),
                subset.idx = hc.idx.lasso))
  }else{
    # if some external controls are selected
    hc.val.lasso.list <- split(b_k_lasso, f = factor(rep(1:K,
                                                         times = EST.ACW.res$n_ei)))

    hc.val.lasso.list <- structure(hc.val.lasso.list, names = NULL)
    EST.ACW.res.lasso <- EST.ACW.FUN(X_c = X_c, Y_c = Y_c,
                                     A_c = A_c, prob_A = prob_A,
                                     data_ec = data_ec,
                                     bias_b_list = hc.val.lasso.list)
    # construct the eta.all for the post-lasso selection
    q_hat.all.lasso <- mapply(function(x, y)x$q_hat*(y==0), x = EST.ACW.res.lasso$EST.ACW,
                              y = hc.val.lasso.list, SIMPLIFY = FALSE)
    # q_hat.all.lasso <-  lapply(EST.ACW.res.lasso$EST.ACW,  function(x)x$q_hat)
    r_hat.all.lasso <-  lapply(EST.ACW.res.lasso$EST.ACW, function(x)x$r_X)

    tau.acw.lasso.list <- unlist(mapply(function(x, y)sum(x*c(rep(1, n_c), (y==0)),
                                                          na.rm = TRUE)/n_c,
                                           x = EST.ACW.res.lasso$tau_ec_i,
                                           y = hc.val.lasso.list, SIMPLIFY = FALSE))

    # }
    # acw.final for each external controls
    var.acw.lasso.list <- EST.ACW.res.lasso$sd_ACW_hat_i**2

    covar.acw.lasso.list <- sapply(EST.ACW.res.lasso$tau_ec_score, function(x)
    {
      sapply(EST.ACW.res.lasso$tau_ec_score, function(y)
      {
        sum((x[1:n_c] - sum(x[1:n_c], na.rm = TRUE)/n_c)*
              (y[1:n_c] - sum(y[1:n_c], na.rm = TRUE)/n_c),
            na.rm = TRUE)/n_c
      })
    })

    Sigma2.group <- matrix(covar.acw.lasso.list, ncol = K + 1,
                           nrow = K + 1)
    diag(Sigma2.group)[1] <- sd_AIPW_hat**2
    diag(Sigma2.group)[2:(K+1)] <- unlist(var.acw.lasso.list)
    Sigma2.group
  }

  # obtain the final selective integrative estimator for each external controls
  d_A_n <- c(rep(1, K+1)%*%MASS::ginv(Sigma2.group)%*%rep(1, K+1))^(-1) *
    MASS::ginv(Sigma2.group)%*%rep(1, K+1)

  tau_final <- c(tau_hat_AIPW, tau.acw.lasso.list)%*%d_A_n
  var.final <- c(d_A_n)%*%Sigma2.group%*%d_A_n
  sd_final <- sqrt(var.final)
  sd_final


  return(list(n_c = n_c,
              n_e = EST.ACW.res$n_ei,
              est = list(AIPW = tau_hat_AIPW,
                         ACW = tau_hat_ACW,
                         ACW.lasso = tau.acw.lasso.list,
                         ACW.final = tau_final
              ),
              sd = list(AIPW = sd_AIPW_hat,
                        ACW = sd_ACW_hat_i,
                        ACW.lasso = sqrt(var.acw.lasso.list),
                        ACW.final = sd_final
              ),
              subset.idx = hc.idx.lasso))
}
