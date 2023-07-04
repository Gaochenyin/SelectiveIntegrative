library(nleqslv)
library(Matrix)
library(pracma)
library(MatchIt)
library(MASS)
library(caret)
# library(gamlss)
#
remove_backticks = function(str) {
  r = sapply(str, function(x) {
    if(substring(x, 1, 1) == "`") x = substring(x, 2); l = nchar(x);
    if(substring(x, l, l) == "`") x = substring(x, 1, l - 1); x })
  r = as.character(r)
  r
}
cal_eqs <- function(par,
                    X_c, X_h)
{

  weight_cal <- as.numeric(exp(X_h%*%par))

  # dem <- rep(tapply(weight_cal, cluster, sum), times=ni)
  # special condition for handling cluster-specific nonignorable missingness
  # weight_cal <- weight_cal/dem*ni
  left_side <- apply(X_h*weight_cal, 2, sum)
  right_side <- apply(X_c, 2, sum)
  c(left_side-right_side)
}

cal_eqs.GMM <- function(par,
                    X_c, X_h)
{

  weight_cal <- as.numeric(exp(X_h%*%par))

  # dem <- rep(tapply(weight_cal, cluster, sum), times=ni)
  # special condition for handling cluster-specific nonignorable missingness
  # weight_cal <- weight_cal/dem*ni
  left_side <- apply(X_h*weight_cal, 2, sum)
  right_side <- apply(X_c, 2, sum)

  sum((left_side-right_side)**2)

  # c(left_side-right_side)
}

generate_hc <- function(X, Y0, n_h0,
                        confound_Y = 0,
                        delta_time = 0,
                        gamma_mean = NA)
{
  N <- nrow(X)
  # # w/o confounder
  # pi_delta_h_noC <- expit(c(X%*%c(-7.5, .1, -0.2)))
  # delta_h_noC <- rbinom(N, size = 1,
  #                   prob = pi_delta_h_noC)
  # n_h_noC <- sum(delta_h_noC)
  # n_h_noC
  # w confounder
  alpha.opt <- uniroot(function(alpha)
    {sum(expit(c(X%*%c(alpha, .1, -0.2)))) - n_h0},
          c(-50, 50))$root
  pi_delta_h <- expit(c(X%*%c(alpha.opt, .1, -0.2)))
  delta_h <- rbinom(N, size = 1,
                    prob = pi_delta_h)
  n_h <- sum(delta_h)
  # generate concurrency bias
  T_Y <- sample(0:2, n_h, replace = T,
                prob = c(1/3, 1/3, 1/3))

  Y0_h <- Y0[delta_h==1] +
    delta_time * T_Y +
    confound_Y * rnorm(n_h, mean = 1);

  # generate measurement error
  gamma_X <- if(is.na(gamma_mean)){0}else{
    rnorm(n_h, gamma_mean)
  }
  X0_h <- cbind(1, X[delta_h==1,-1]  + gamma_X) # intercept has no measurement error
  return(list(X0 = X0_h,
              Y0 = Y0_h))

}
estimate_hc <- function(X_h, Y_h,
                        X_c, Y_c, A_c, prob_A,
                        fit.Y.RCT, hc.ctrl,
                        bias_b = NULL, method = 'rpart', tuneGrid = NULL,
                        ...)
{
  # r_X: the ratio of var(Y(0)|X,delta_c=1)/var(Y(0)|X,delta_h=1)
  # X_c <- as.matrix(X_c); X_h <- as.matrix(X_h)
  # data pre-possessing
  p <- ncol(X_h);
  if(is.data.frame(X_h)){
    X_h <- model.matrix(~.,#+ I(AGE_standardized**2) + I(BMI_standardized**2),
                        X_h)
  }else{X_h <- as.matrix(X_h)}
  n_c <- nrow(X_c); n_h <- nrow(X_h)
  # the default value for bias
  if(is.null(bias_b)){bias_b <- rep(0, n_h)}
  n_c0 <- sum(A_c==0)
  mu_i <- as.numeric(n_c + n_h)
  # estimate outcome model (combined)
  X0_ch <- rbind(X_c[A_c==0,], X_h)
  X0_ch <- as.matrix(X0_ch)
  Y0_ch <- c(Y_c[A_c==0], Y_h)
  # X0_ch <- as.matrix(rbind(X_c, X_h))
  # Y0_ch <- c(Y_c - X_c%*%psi_A, Y_h)
  library(MASS)
  # combined modeling
  # lm.fit(x = X0_ch, y = Y0_ch)$coefficients
  # beta_hat <- ginv(t(X0_ch)%*%X0_ch)%*%t(X0_ch)%*%Y0_ch
  # combined weighted modeling
  # beta_hat_0.h <- beta_hat_0.c <- nleqslv(x = c(beta_hat),
  #         fn = function(x)
  #         {
  #           apply(c(Y0_ch - X0_ch%*%x)*
  #                   c(1/(1 - prob_A[A_c==0]), rep(1, n_h))*X0_ch,
  #                 2, mean)
  #         })$x

  # if concurrent control is small, might have problem
  # separated modeling
  # beta_hat_0.h <- nleqslv(x = c(beta_hat),
  #         fn = function(x)
  #         {
  #           apply(c(Y_h - X_h%*%x)*X_h,
  #                 2, mean)
  #         })$x
  #
  # beta_hat_0.c <- nleqslv(x = c(beta_hat),
  #         fn = function(x)
  #         {
  #           apply(c(Y_c[A_c==0] - X_c[A_c==0,]%*%x)*
  #                   c(1/(1 - prob_A[A_c==0]))*X_c[A_c==0,],
  #                 2, mean)
  #         })$x
  # only the RCT modelling
  # beta_hat_0.c <- beta_A0_c
  # beta_hat_0.h <- ginv(t(X_h)%*%X_h)%*%t(X_h)%*%Y_h

  # use machine learning model to fitting mu_1 and mu_0
  fit.Y.HC <- train(Y ~ .+0, data = data.frame(Y = Y_h, X_h),
                     trControl = hc.ctrl, ...)


  # estimate odds-ratio
  library(BB)
  if(ncol(X_c)<n_h)
  {
    lambda_hat <- dfsane(par = rep(0, dim(X_c)[2]), # initial points
                         fn=cal_eqs,
                         X_c = X_c, X_h = X_h,#[bias_b==0,],
                         control = list(trace=T,
                                        NM=T,
                                        BFGS=F,
                                        tol=1.e-8,
                                        maxit = 500))$par

    # cal_eqs(lambda_hat,
    #         X_c, X_h[bias_b==0,])
  }else{
    # GMM calibration
    lambda_hat <-  optim(par = rep(0, dim(X_c)[2]), # initial points
          fn=cal_eqs.GMM,
          X_c = X_c, X_h = X_h,#[bias_b==0,],
          control = list(trace=T,
                         abstol=1.e-8,
                         maxit = 500),
          method = 'BFGS')$par
    # cal_eqs.GMM(lambda_hat,
    #         X_c, X_h[bias_b==0,])
    # lambda_hat <- dfsane(par = rep(0, dim(X_c)[2]), # initial points
    #                      fn=cal_eqs,
    #                      X_c = X_c, X_h = X_h,
    #                      control = list(trace=T,
    #                                     NM=T,
    #                                     BFGS=F,
    #                                     tol=1.e-8,
    #                                     maxit = 500))$par
  }


  q_hat <- as.numeric(exp(X_h%*%lambda_hat))
  q_hat.c <- as.numeric(exp(X_c%*%lambda_hat))


  # # -----------------
  # # estimate the mu_acw
  # ## the first n_c are the concurrent controls
  # ## the last n_h are the historical controls
  # mu_i[1:n_c] <- (X_c%*%beta_hat_0.c) #* n_c^{-1}
  # mu_i[(n_c+1):(n_c+n_h)] <- q_hat * (Y_h - predict(fit.Y.RCT, data.frame(X_h)))#X_h%*%beta_hat_0.c) #* n_c^{-1}
  # mu_acw <- sum(mu_i) * n_c**(-1)
  # # mu0_acw <- (sum(q_hat))^{-1}*sum(q_hat*X0_h%*%beta_hat) +
  # #   mean(Y0_h-X0_h%*%beta_hat)
  # # compute the score for estimating beta
  # # dot.tau.ACW.beta0 <- apply(X_c, 2, sum)/n_c -
  # #   apply(q_hat * X_h, 2, sum)/n_c
  # dot.tau.ACW.beta0 <- apply(X_c, 2, sum)/n_c
  # S.ACW.beta0.RCT <- c(Y_c-A_c*X_c%*%psi_A-predict(fit.Y.RCT, data.frame(X_c)))*X_c
  # dot.S.ACW.beta0 <- -  -(t(X_c)%*%X_c)/n_c
  #
  # # S.ACW.beta0.RCT <- c((1-A_c)/(1-prob_A)*(Y_c-X_c%*%beta_hat_0.c))*X_c
  # # dot.S.ACW.beta0 <- -  -(t((1-A_c)/(1-prob_A)*X_c)%*%X_c)/n_c
  #
  # # compute the score for estimating q
  # dot.tau.ACW.q <- apply(q_hat*c(Y_h - predict(fit.Y.RCT, data.frame(X_h)))*X_h,
  #                        2, sum)/n_c
  S.ACW.q <- rbind(-X_c, q_hat*X_h)
  dot.S.ACW.q <- t(q_hat*X_h)%*%X_h/n_c
  #
  # dot.tau.ACW.beta0.RWD <- -apply(q_hat * X_h, 2, sum)/n_c
  # # compute the adjusted influence function
  # tau.ACW.c <- c(X_c%*%beta_hat_0.c -
  #                  S.ACW.beta0.RCT%*%t(dot.tau.ACW.beta0%*%ginv(dot.S.ACW.beta0)),
  #                rep(0, n_h))
  # # tau.ACW.h <- c(rep(0, n_c),
  # #                 q_hat * (Y_h - X_h%*%beta_hat_0.h))-
  # #   S.ACW.q%*%t(dot.tau.ACW.q%*%ginv(dot.S.ACW.q))
  #
  # tau.ACW.h <- c(-S.ACW.beta0.RCT%*%t(dot.tau.ACW.beta0.RWD%*%ginv(dot.S.ACW.beta0)),
  #                 q_hat * (Y_h - predict(fit.Y.RCT, data.frame(X_h))))-
  #   S.ACW.q%*%t(dot.tau.ACW.q%*%ginv(dot.S.ACW.q))
  #
  # var_hat_c <- sum({tau.ACW.c-sum(tau.ACW.c)/n_c}**2)
  # var_hat_h <- sum({tau.ACW.h-sum(tau.ACW.h)/n_c}**2)
  #
  # tau_ACWi <- tau.ACW.c + tau.ACW.h
  # var_hat <- var_hat_c + var_hat_h
  # # var_hat_c <- sum((X_c%*%beta_hat_0.c - mu_acw)^2)
  # # var_hat_h <- sum(q_hat^2 * (Y_h - X_h%*%beta_hat_0.h)^2)
  # # var_hat <- var_hat_h + var_hat_c
  # # --------------------

  # a new set of semi-parameteric estimator

  # fit gamlss for conditional variance
  # gam.RCT.A0 <- gamlss(Y~.+0, data = data.frame(Y = Y_c[A==0],
  #                                               X_c[A==0]))
  # gam.HC <- gamlss(Y~.+0, data = data.frame(Y = Y_h,
  #                                            X_h))
  # r_X.c <- predict(gam.RCT.A0, 'sigma', newdata = data.frame(X_c[A==0]))**2/
  #   predict(gam.HC, 'sigma', newdata = data.frame(X_c[A==0]))**2
  #
  # r_X.h <- predict(gam.RCT.A0, 'sigma', newdata = data.frame(X_h))**2/
  #   predict(gam.HC, 'sigma', newdata = data.frame(X_h))**2

  r_X <- r_X.c <- r_X.h <- mean({predict(fit.Y.RCT, data.frame(X_c[A_c==0,], A = 0)) - mean(Y_c[A_c==0])}**2)/
    tryCatch(mean({predict(fit.Y.HC, data.frame(X_h[bias_b==0, ])) - mean(Y_h[bias_b==0])}**2),
             error = function(e)1e6) # if no selected, down weight the HC

  # prob_b_zero_c <- prob_b_zero_h <- sum(bias_b==0)/n_h

  # fit a glm for b==0 against X_h
  fit.b0.hc <- glm(b0 ~ . +0, data = data.frame(b0 = (bias_b == 0),
                                                X_h),
                   family = 'binomial')
  prob_b_zero_c <- predict(fit.b0.hc, newdata = data.frame(X_c),
                           type = 'response')
  prob_b_zero_h <- predict(fit.b0.hc, newdata = data.frame(X_h),
                           type = 'response')
  # fit a glm for A against X_c
  fit.ps <- glm(A~.+0, data = data.frame(A = A_c,
                                       X_c),
                family = 'binomial')
  prob_A.h <- predict(fit.ps, newdata = data.frame(X_h),
                      type = 'response')


  # construct the influence function
  dot.mu0_c.q <- apply(c({q_hat.c * r_X * prob_b_zero_c}*
                           (1-A_c)*(Y_c - predict(fit.Y.RCT, data.frame(X_c, A = 0)))/
                           {{r_X * prob_b_zero_c + (1-prob_A)*q_hat.c}**2})*X_c, 2, sum)/n_c

  dot.mu0_h.q <- apply(c({q_hat * prob_b_zero_h * r_X**2}*
                           prob_b_zero_h * (Y_h - predict(fit.Y.RCT, data.frame(X_h, A = 0)))/
                           {{r_X * prob_b_zero_h + (1-prob_A.h)*q_hat}**2})*X_h, 2, sum)/n_c
  mu0_c <- c(predict(fit.Y.RCT, data.frame(X_c, A = 0)) +
               q_hat.c/{prob_b_zero_c * r_X + (1-prob_A)*q_hat.c}*(1-A_c)*(Y_c - predict(fit.Y.RCT, data.frame(X_c, A = 0))),
             rep(0, n_h))
  mu0_h <- c(rep(0, n_c),
             q_hat*r_X/{prob_b_zero_h * r_X + (1-prob_A.h)*q_hat} *(bias_b==0)*
               (Y_h - predict(fit.Y.RCT, data.frame(X_h, A = 0))))
  mu0_i <- mu0_c + mu0_h
  mu0_score <-  mu0_i -  S.ACW.q%*%t(dot.mu0_c.q%*%ginv(dot.S.ACW.q)) -
    S.ACW.q%*%t(dot.mu0_h.q%*%ginv(dot.S.ACW.q))

  # mu0_score <- mu0_score * c(rep(1, n_c), bias_b==0)
  # sum(mu1_AIPW_i)/n_c - sum(mu0_i)/n_c
  # dot.tau.ACW.q <- apply(q_hat*c(Y_h - X_h%*%beta_hat_0.h)*X_h,
  #                        2, sum)/n_c

  mu0_ACW.hat <- (sum(mu0_i))/n_c
  # fit a model for the bias
  # fit.bias.lm <- lm(bias ~ .+0,
  #                   data = data.frame(bias = Y_h - predict(fit.Y.RCT.A0, data.frame(X_h)),
  #                                     X_h))
  # b.X <- fit.bias.lm$fitted.values
  bias_i <- q_hat * c(Y_h - predict(fit.Y.RCT, data.frame(X_h, A = 0)))*(n_c+n_h)/n_c #-
    # {q_hat*(n_c+n_h)/n_c - 1}*b.X

   # compute the score for estimating q
  dot.bias.ACW.q <- apply(q_hat*c(Y_h - predict(fit.Y.RCT, data.frame(X_h, A = 0)))*X_h,
                          2, sum)*(n_c+n_h)/n_c**2
  bias_i.score <- c(rep(0, n_c), bias_i) -
    S.ACW.q%*%t(dot.bias.ACW.q%*%ginv(dot.S.ACW.q))
  gamma.hat <- predict(fit.Y.HC, data.frame(X_h)) -
    predict(fit.Y.RCT, data.frame(X_h, A = 0))

  # compute the score for the bias.dr
  # S.ACW.q%*%ginv(dot.S.ACW.q)%*%t(dot.tau.ACW.q) # score for calibration weights
  #
  #
  # dot.bias.q <- apply(q_hat*c(Y_h - X_h%*%beta_hat_0.h - b.X)*X_h,
  #                     2, sum)/n_c
  #
  # dot.bias.beta0 <- -apply(q_hat * X_h, 2, sum)/n_c
  #
  # dot.bias.gamma <- -{apply(q_hat * X_h*(n_c+n_h)/n_c, 2, sum)-
  #     apply(rbind(X_c, X_h), 2, sum)}/(n_c+n_h)
  #
  #
  # dot.S.bias.beta0 <- -  -((1-A_c)*t(X_c)%*%X_c)/n_c
  # dot.S.bias.gamma <- -(t(X_h)%*%X_h)/n_c
  #
  # bias.ipw.score <- bias.ipw - score.q[-(1:n_c)]


  list(n_c = n_c, n_h = n_h,
       X_h = X_h, Y_h = Y_h,
       q_hat = q_hat,
       # var_hat = var_hat * n_c^{-2},
       # var_hat_h = var_hat_h * n_c^{-2},
       # var_hat_c =var_hat_c * n_c^{-2},
       # beta_hat_0.c = beta_hat_0.c,
       # beta_hat_0.h = beta_hat_0.h,
       # tau_ACW.RCT = tau.ACW.c,
       # tau_ACW.RWD = tau.ACW.h,
       # tau_ACWi = tau_ACWi,
       mu0_i = mu0_i, mu0_score = mu0_score,
       bias_i = bias_i,
       bias_i.score = bias_i.score,
       gamma.hat = gamma.hat, r_X = r_X)
}


estimate_hc_learner <- function(X_h, Y_h,
                                X_c, Y_c, A_c, prob_A)
{
  # X_c <- as.matrix(X_c); X_h <- as.matrix(X_h)
  # data pre-possessing
  p <- ncol(X_h)
  if(is.data.frame(X_h)){
    X_h <- model.matrix(~.,#+ I(AGE_standardized**2) + I(BMI_standardized**2),
                        X_h)
  }else{X_h <- as.matrix(X_h)}
  n_c <- nrow(X_c); n_h <- nrow(X_h)
  n_c0 <- sum(A_c==0)
  mu_i <- as.numeric(n_c + n_h)
  # estimate outcome model (combined)
  X0_ch <- rbind(X_c[A_c==0,], X_h)
  X0_ch <- as.matrix(X0_ch)
  Y0_ch <- c(Y_c[A_c==0], Y_h)

  Y.fit.rf.combined <- randomForest::randomForest(x = X_c[A_c==0,],
                                                  y = Y_c[A_c==0],
                                                  ntree = 500)

  library(MASS)
  # estimate odds-ratio
  library(BB)
  if(ncol(X_c)<nrow(X_h))
  {
    lambda_hat <- dfsane(par = rep(0, dim(X_c)[2]), # initial points
                         fn=cal_eqs,
                         X_c = X_c, X_h = X_h,
                         control = list(trace=T,
                                        NM=T,
                                        BFGS=F,
                                        tol=1.e-8,
                                        maxit = 500))$par
  }else{
    # GMM calibration
    lambda_hat <-  optim(par = rep(0, dim(X_c)[2]), # initial points
                         fn=cal_eqs.GMM,
                         X_c = X_c, X_h = X_h,
                         control = list(trace=T,
                                        abstol=1.e-8,
                                        maxit = 500))$par

    # lambda_hat <- dfsane(par = rep(0, dim(X_c)[2]), # initial points
    #                      fn=cal_eqs,
    #                      X_c = X_c, X_h = X_h,
    #                      control = list(trace=T,
    #                                     NM=T,
    #                                     BFGS=F,
    #                                     tol=1.e-8,
    #                                     maxit = 500))$par
  }


  q_hat <- as.numeric(exp(X_h%*%lambda_hat))
  # estimate the mu_acw
  ## the first n_c are the concurrent controls
  ## the last n_h are the historical controls

  mu_i[1:n_c] <- predict(Y.fit.rf.combined,
                         newdata = X_c) #* n_c^{-1}
  mu_i[(n_c+1):(n_c+n_h)] <- q_hat * (Y_h -
                                        predict(Y.fit.rf.combined,
                                                newdata = X_h)) #* n_c^{-1}
  mu_acw <- sum(mu_i) * n_c**(-1)

  # compute the adjusted influence function
  tau.ACW.c <- c(predict(Y.fit.rf.combined, newdata = X_c),
                 rep(0, n_h))
  # tau.ACW.h <- c(rep(0, n_c),
  #                 q_hat * (Y_h - X_h%*%beta_hat_0.h))-
  #   S.ACW.q%*%t(dot.tau.ACW.q%*%ginv(dot.S.ACW.q))

  tau.ACW.h <- c(rep(0, n_c),
                 q_hat * (Y_h - predict(Y.fit.rf.combined, newdata = X_h)))

  var_hat_c <- sum({tau.ACW.c-sum(tau.ACW.c)/n_c}**2)
  var_hat_h <- sum({tau.ACW.h-sum(tau.ACW.h)/n_c}**2)

  tau_ACWi <- tau.ACW.c + tau.ACW.h
  var_hat <- var_hat_c + var_hat_h
  # var_hat_c <- sum((X_c%*%beta_hat_0.c - mu_acw)^2)
  # var_hat_h <- sum(q_hat^2 * (Y_h - X_h%*%beta_hat_0.h)^2)
  # var_hat <- var_hat_h + var_hat_c




  list(n_c = n_c, n_h = n_h,
       X_h = X_h, Y_h = Y_h,
       q_hat = q_hat, mu_i = mu_i,
       var_hat = var_hat * n_c^{-2},
       var_hat_h = var_hat_h * n_c^{-2},
       var_hat_c =var_hat_c * n_c^{-2},
       tau_ACW.RCT = tau.ACW.c,
       tau_ACW.RWD = tau.ACW.h,
       tau_ACWi = tau_ACWi)
}

# sub-function
AIPW <- function(X_c, Y_c, A_c, prob_A = NULL)
{
  # data pre-processing
  X_c <- as.matrix(X_c)
  n_c <- length(Y_c)
  # K <- length(data_hc)
  Y0_c <- Y_c[which(A_c==0)]; X0_c <- X_c[which(A_c==0),]
  Y1_c <- Y_c[which(A_c==1)]; X1_c <- X_c[which(A_c==1),]

  # estimate outcome model (concurrent only)
  beta_hat_A0_c <- ginv(t(X0_c)%*%X0_c)%*%t(X0_c)%*%Y0_c
  beta_hat_A1_c <- ginv(t(X1_c)%*%X1_c)%*%t(X1_c)%*%Y1_c



  # estimate propensity score (if not provided)
  if(is.null(prob_A))
  {
    glm.ps <- glm(A_c ~ X_c, family = binomial())
    prob_A <- predict(glm.ps, type = 'response')
  }

  # # improve efficiency
  # library(nleqslv)
  # nleqslv(x = c(beta_hat_A1_c),
  #         fn = function(x)
  #         {
  #          apply(c(A_c/prob_A* (1-prob_A)/prob_A *
  #                    (Y_c - X_c%*%x))*X_c,
  #                2, mean)
  #         })
  #
  #
  # nleqslv(x = c(beta_hat_A0_c),
  #         fn = function(x)
  #         {
  #           apply(c((1-A_c)/(1-prob_A) *
  #                     (prob_A)/(1-prob_A) *
  #                     (Y_c - X_c%*%x))*X_c,
  #                 2, mean)
  #         })

  # doubly robust variance estimation
  mu1_AIPW_i <- {A_c * Y_c/prob_A +
      (prob_A - A_c)/prob_A*(X_c%*%beta_hat_A1_c)}#/n_c
  mu0_AIPW_i <- {(1-A_c) * Y_c/(1-prob_A) +
      (A_c - prob_A)/(1-prob_A) * (X_c%*%beta_hat_A0_c)}#/n_c
  tau_0_i <-  mu1_AIPW_i - mu0_AIPW_i
  tau_AIPW <- sum(tau_0_i)/n_c
  sqrt(sum({tau_0_i - tau_AIPW}^2)*n_c^{-1})
  # compute \partial \tau/\partial \beta
  dot.tau.beta1 <- apply((1-A_c/prob_A)*X_c, 2, sum)/n_c
  dot.tau.beta0 <- -apply((1-(1-A_c)/(1-prob_A))*X_c, 2, sum)/n_c

  dot.S.beta0 <- -t((1-A_c)*X_c)%*%((1-A_c)*X_c)/n_c
  dot.S.beta1 <- -t((A_c)*X_c)%*%((A_c)*X_c)/n_c

  S.beta1 <- c(A_c*(Y_c-X_c%*%beta_hat_A1_c))*X_c
  S.beta0 <- c((1-A_c)*(Y_c-X_c%*%beta_hat_A0_c))*X_c

  tau_AIPW.beta1 <- S.beta1%*%t(dot.tau.beta1%*%ginv(dot.S.beta1))
  tau_AIPW.beta0 <- S.beta0%*%t(dot.tau.beta0%*%ginv(dot.S.beta0))
  tau_0_i <- tau_0_i - tau_AIPW.beta1 - tau_AIPW.beta0
  # sum(tau_0_i)/n_c
  # compute s.d. for mu1_AIPW-mu0_AIPW
  sd_AIPW_hat <- sqrt(sum({tau_0_i - tau_AIPW}^2)*n_c^{-1})

  # compute s.d. for mu1_AIPW - mu0_ACW
  sum({mu1_AIPW_i-sum(mu1_AIPW_i)/n_c}^2)*n_c^{-1}
  sum({mu0_AIPW_i-sum(mu0_AIPW_i)/n_c}^2)*n_c^{-1}

  return(list(AIPW = tau_AIPW,
              sd = sd_AIPW_hat/sqrt(n_c)))
}



# main function (super learner)
AIPW_ACW_lasso_SuperLearner <- function(X_c, Y_c, A_c, prob_A = NULL,
                                        data_hc,
                                        lambda = 'min',
                                        strata = 1, lasso = 'group',
                                        match.adjust = FALSE,
                                        refit = TRUE)
{
  # data pre-processing
  if(is.data.frame(X_c)){
    # construct the basis function for prediction and calibration
    X_c <- model.matrix(~.,# + I(AGE_standardized**2) + I(BMI_standardized**2),
                        X_c) # for calibration
  }else{X_c <- as.matrix(X_c)}
  n_c <- length(Y_c)
  K <- length(data_hc)
  Y0_c <- Y_c[which(A_c==0)]; X0_c <- X_c[which(A_c==0),, drop = F]
  Y1_c <- Y_c[which(A_c==1)]; X1_c <- X_c[which(A_c==1),, drop = F]


  # estimate outcome model (concurrent only)
  RCT.dat <- as.data.frame(cbind(Y_c, X_c, group = A_c))

  ## constant treatment effects
  # ## joint estimation
  # fit.RCT <- lm(Y_c~.+0, data = RCT.dat)
  # # fit.RCT <- lm.fit(x = cbind(X_c,
  # #                              group = A_c), y = Y_c)
  # beta_hat <- fit.RCT$coefficients
  # beta_hat_A1_c <- beta_hat_A1_c <- beta_hat[names(beta_hat)!='group']
  # beta_hat_A1_c['`(Intercept)`'] <- beta_hat_A0_c['`(Intercept)`'] + beta_hat['group']
  # psi_A <- beta_hat_A1_c - beta_hat_A0_c

  ## hetero treatment effects
  # X.design <- model.matrix(~ 0 + . + group +
  #                            group:(AGE_standardized + BMI_standardized),
  #                          as.data.frame(cbind(X_c, group = A_c)))
  #
  # nleqslv(x = rep(0, ncol(X.design)),
  #                 fn = function(x)
  #                 {
  #                   X.design.A_c <- cbind(X.design[,-grep("group$",  colnames(X.design))],
  #                                         A_c*X.design[,grep("group$",  colnames(X.design))])
  #
  #                   res.Y <- Y_c - X.design.A_c%*%x
  #
  #                   apply(c(res.Y)*X.design.A_c, 2, sum)
  #                 })$x
  Y.fit.rf <- randomForest::randomForest(x = cbind(X_c, group = A_c),
                             y = Y_c,
                             ntree = 500)

  # fit.RCT <- lm(Y_c ~ 0 + . + group +
  #                 group:(.-`(Intercept)`),
  #               data = RCT.dat)
  # beta_hat <- fit.RCT$coefficients
  # beta_hat_A1_c <- beta_hat_A0_c <- beta_hat[-grep("group$",  names(beta_hat))]
  #
  # # interaction term and intercept
  # interaction.idx <- grep(":group$",  names(beta_hat))
  # interaction.name <- sub("(*):group", "\\1", names(beta_hat))[interaction.idx]
  # beta_hat_A1_c[interaction.name] <- beta_hat_A0_c[interaction.name] +
  #   beta_hat[interaction.idx]
  # beta_hat_A1_c['`(Intercept)`'] <- beta_hat_A0_c['`(Intercept)`'] + beta_hat['group']
  #
  #
  # psi_A <- beta_hat_A1_c - beta_hat_A0_c
  # estimate propensity score (if not provided)
  if(is.null(prob_A))
  {
    glm.ps <- glm(A_c ~ X_c_A0, family = binomial())
    prob_A <- predict(glm.ps, type = 'response')
  }

  ## super-learner
  # library(AIPW)
  # AIPW.learner <- AIPW::AIPW$new(Y = Y_c,
  #                A = A_c,
  #                W = X_c[,-1],
  #                Q.SL.library="SL.mean",g.SL.library="SL.mean",
  #                k_split=2,
  #                save.sl.fit = TRUE)
  # AIPW.learner <- aipw_wrapper(Y = Y_c,
  #              A = A_c,
  #              W = X_c[,-1],
  #              Q.SL.library="SL.mean",g.SL.library="SL.mean",
  #              k_split=1)
  # mu1_AIPW_i <- AIPW.learner$obs_est$aipw_eif1
  # mu0_AIPW_i <- AIPW.learner$obs_est$aipw_eif0


  # # step 1: initial estimate of psi0
  # psi0 <- nleqslv(x = rep(0, ncol(X_c)),
  #         fn = function(x)
  #         {
  #           apply(X_c*
  #           c(A_c - prob_A)*
  #             c(Y_c - A_c*X_c%*%x), 2, mean)
  #         })$x
  # # step 2: obtain the estimator for the pseudo outcome Y(0)
  # Y0_c.pseudo <- Y_c - A_c*X_c%*%psi0
  # beta_hat_A0_c <- ginv(t(X_c)%*%X_c)%*%t(X_c)%*%Y0_c.pseudo
  # # step 3: re-estimate the psi
  # psi <- nleqslv(x = psi0,
  #         fn = function(x)
  #         {
  #           apply(X_c*
  #                   c(A_c - prob_A)*
  #                   c(Y_c - A_c*X_c%*%x- X_c%*%beta_hat_A0_c), 2, mean)
  #         })$x
  # # step 4: revise the estimator for the benchmark
  # Y0_c.pseudo <- Y_c - A_c*X_c%*%psi
  # beta_hat_A0_c <- ginv(t(X_c)%*%X_c)%*%t(X_c)%*%Y0_c.pseudo
  # beta_hat_A1_c <- beta_hat_A0_c + psi


  # beta_hat_A0_c <- ginv(t(X0_c)%*%X0_c)%*%t(X0_c)%*%Y0_c
  # beta_hat_A0_c <- nleqslv(x = c(beta_hat_A0_c),
  #         fn = function(x)
  #         {
  #           apply(c(Y0_c - X0_c%*%x)*
  #                   c(1/(1 - prob_A[A_c==0]))*X0_c,
  #                 2, mean)
  #         })$x

  # beta_hat_A1_c <- ginv(t(X1_c)%*%X1_c)%*%t(X1_c)%*%Y1_c
  # beta_hat_A1_c <- nleqslv(x = c(beta_hat_A1_c),
  #         fn = function(x)
  #         {
  #           apply(c(Y1_c - X1_c%*%x)*
  #                   c(1/prob_A[A_c==1])*X1_c,
  #                 2, mean)
  #         })$x
  # beta_hat_A0_c <- c(1, rep(.5, 4), rep(-.5, 4), rep(1, 4))
  # beta_hat_A1_c <- c(1, rep(.5, 4), rep(-.5, 4), rep(1, 4))

  # doubly robust variance estimation
  mu1_AIPW_i <- predict(Y.fit.rf,
                        newdata = cbind(X_c, group = 1))
  mu0_AIPW_i <- predict(Y.fit.rf,
                        newdata = cbind(X_c, group = 0))
  tau_0_i <-  mu1_AIPW_i - mu0_AIPW_i
  tau_hat_AIPW <- sum(tau_0_i)/n_c
  # compute s.d. for mu1_AIPW-mu0_AIPW
  sd_AIPW_hat <- sqrt(sum({tau_0_i - tau_hat_AIPW}^2)*n_c^{-1})
  sd_AIPW_hat


  # if no historical controls are matched
  if(nrow(data_hc[[1]]$X0)==0){

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



  # strata <- 5
  # stratify the historical controls based on the fitted residual
  data_hc_stratified <- lapply(data_hc, function(hc_i){
    X_h <- hc_i$X0
    Y_h <- hc_i$Y0
    if(is.data.frame(X_h)){
      X_h <- model.matrix(~.,# + I(AGE_standardized**2) + I(BMI_standardized**2),
                          X_h)
    }else{X_h <- as.matrix(X_h)}
    residual.RW <- Y_h - predict(Y.fit.rf,
                                 newdata = cbind(X_h, group = 0))
    # strata <- length(residual.RW)
    # strata <- 5
    ## stratify based on quantile
    # strata.idx <- split(order(abs(residual.RW)),
    #                     cut(seq_along(residual.RW), strata, labels = FALSE))
    #
    # strata.idx <- structure(strata.idx, names = NULL)
    #
    # lapply(strata.idx, function(idx)
    # {
    #   list(X0 = X_h[idx,],
    #        Y0 = Y_h[idx])
    # })
    ## stratify based on K-means
    # plot(Y_h, col = strata.idx)
    # plot(Y_h, col = data.list$RWD$T_i+1)
    strata.idx <- kmeans(residual.RW, centers = strata)$cluster
    lapply(1:strata, function(s)
    {
      idx <- which(strata.idx == s)
      list(X0 = X_h[idx,],
           Y0 = Y_h[idx])
    })
  })
  # unwrap the stratified
  data_hc <- do.call(c, data_hc_stratified)
  # compute s.d. for mu1_AIPW - mu0_ACW
  var_mu1_AIPW_hat <- sum({mu1_AIPW_i-sum(mu1_AIPW_i)/n_c}^2)*n_c^{-1}
  var_mu0_AIPW_hat <- sum({mu0_AIPW_i-sum(mu0_AIPW_i)/n_c}^2)*n_c^{-1}


  EST.ACW.FUN <- function(X_c, Y_c, A_c, prob_A,
                          data_hc){
    # linear expansion of augmented calibration estimator
    data_hc_i <- lapply(data_hc, function(hc_i)
    {
      estimate_hc_learner(X_c = X_c, X_h = hc_i$X0,
                          Y_c = Y_c, Y_h = hc_i$Y0,
                          A_c = A_c, prob_A = prob_A)
    })

    # cov_AIPW_ACW_i <- lapply(data_hc_i, function(hc_i)sum((X_c%*%beta_hat_A1_c -
    #                                                          sum(mu1_AIPW_i)/n_c)*
    #                                                         (X_c%*%hc_i$beta_hat_0.c -
    #                                                            sum(hc_i$mu_i)/n_c))*n_c^{-1})


    sd_ACW_hat_i <-  sapply(data_hc_i, function(hc_i)
    {
      tau.ACW.lin <- c(mu1_AIPW_i,
                       rep(0, hc_i$n_h))- hc_i$tau_ACWi
      sqrt(sum({tau.ACW.lin-sum(tau.ACW.lin)/n_c}**2)*n_c^{-1})
    })

    # sd_ACW_hat_i <- mapply(function(hc_i, cov){
    #   sqrt(hc_i$var_hat*n_c + var_mu1_AIPW_hat-2*cov)
    # }, hc_i = data_hc_i, cov = cov_AIPW_ACW_i)

    # lapply(data_hc_i, function(x)x$var_hat_h*n_c)
    # lapply(data_hc_i, function(x)x$var_hat_c*n_c)
    # lapply(data_hc_i, function(x)x$var_hat*n_c)
    # construct consistent estimator of b_k
    tau_hc_i <- lapply(data_hc_i, function(hc_i)c(mu1_AIPW_i,
                                                  rep(0, hc_i$n_h))- hc_i$mu_i)

    # tau_hc_i <- lapply(data_hc_i, function(hc_i)c(mu1_AIPW_i,
    #                                               rep(0, hc_i$n_h))- hc_i$tau_ACWi)

    # b_hat <- lapply(tau_hc_i,
    #                 function(tau_i)sum(tau_i)-2)%>%unlist()
    # b_hat <- lapply(tau_hc_i,
    #                 function(tau_i)sum(tau_i)/n_c-sum(tau_0_i)/n_c)%>%unlist()

    # construct consistent variance estimator
    # # One attempt: (tau_i-tau)^2 (Lunceford and Davidian (2004))
    # tau_hc_i_s1 <- lapply(tau_hc_i,
    #                      function(tau_i)tau_i/sd(tau_i))
    # another variance estimation by M-estimation
    # construct the co-variance matrix for tau_hat
    ## var(tau_0): V_c^[0]: sd_AIPW_hat**2
    V_c.0 <- sd_AIPW_hat**2
    ## var(tau_k): V_h^[k]: hc_i$var_hat*n_c, V_c^[k]: var_mu1_AIPW_hat-2*cov
    V_h.k <- sapply(data_hc_i, function(hc_i)hc_i$var_hat_h*n_c)

    # V_c.k <- mapply(function(hc_i, cov) hc_i$var_hat_c*n_c +
    #                   var_mu1_AIPW_hat-2*cov,
    #                 hc_i = data_hc_i,
    #                 cov = cov_AIPW_ACW_i)
    V_c.k <- sapply(data_hc_i, function(x){
      tau.ACW.RCT <- c(mu1_AIPW_i,
                       rep(0, x$n_h))-
        x$tau_ACW.RCT
      # tau.ACW.RWD <- x$tau_ACW.RWD

      sum({tau.ACW.RCT-sum(tau.ACW.RCT)/n_c}^2)*n_c^{-1}
      # sum(tau.ACW.RCT**2)*n_c^{-1}#+
      # sum({tau.ACW.RWD-sum(tau.ACW.RWD)/n_c}^2)*n_c^{-1}
    })

    V_c.kk.dot <- sapply(data_hc_i, function(x)
    {
      sapply(data_hc_i, function(y)
      {
        ## for dataset x
        # compute the adjusted influence function
        tau.ACW.RCT.x <- c(mu1_AIPW_i,
                           rep(0, x$n_h)) - x$tau_ACW.RCT
        # tau.ACW.RWD.x <- -S.ACW.beta0.RWD%*%t(dot.tau.ACW.beta0%*%ginv(dot.S.ACW.beta0))

        ## for dataset y
        # compute the adjusted influence function
        tau.ACW.RCT.y <- c(mu1_AIPW_i,
                           rep(0, y$n_h))- y$tau_ACW.RCT
        # tau.ACW.RWD.y <- -S.ACW.beta0.RWD%*%t(dot.tau.ACW.beta0%*%ginv(dot.S.ACW.beta0))

        if(identical(x,y)){
          sum({tau.ACW.RCT.x-sum(tau.ACW.RCT.x)/n_c}*
                {tau.ACW.RCT.y-sum(tau.ACW.RCT.y)/n_c})*n_c^{-1}
        }else{
          sum({tau.ACW.RCT.x[1:n_c]-sum(tau.ACW.RCT.x)/n_c}*
                {tau.ACW.RCT.y[1:n_c]-sum(tau.ACW.RCT.y)/n_c})*n_c^{-1}
        }



        # sum(tau.ACW.RCT.x*tau.ACW.RCT.y)*n_c^{-1}
        # sum({tau.ACW.RWD.x-sum(tau.ACW.RWD.x)/n_c}*
        #       {tau.ACW.RWD.y-sum(tau.ACW.RWD.y)/n_c})*n_c^{-1}

      })
    })

    ## cov(tau_0, tau_k): V_c^[0,k]
    V_c.0k <- sapply(data_hc_i, function(x){

      # compute the adjusted influence function
      tau.ACW.RCT.0 <- tau_0_i
      tau.ACW.RCT.x <- mu1_AIPW_i - x$tau_ACW.RCT[1:n_c]
      # tau.ACW.RWD <- -S.ACW.beta0.RWD%*%t(dot.tau.ACW.beta0%*%ginv(dot.S.ACW.beta0))


      sum({tau.ACW.RCT.0-sum(tau.ACW.RCT.0)/n_c}*
            {tau.ACW.RCT.x-sum(tau.ACW.RCT.x)/n_c})*n_c^{-1}
      sum(tau.ACW.RCT.0*tau.ACW.RCT.x)*n_c^{-1}#+
      # sum({tau.ACW.RWD-sum(tau.ACW.RWD)/n_c}^2)*n_c^{-1}

    })

    V_h.0k.list <- sapply(data_hc_i, function(x){
      V_c.vec <- -c(mu1_AIPW_i -
                      x$tau_ACW.RCT[1:n_c])*(x$tau_ACW.RWD[1:n_c]/n_c)
      matrix(rep(V_c.vec, x$n_h),
             ncol = x$n_h)
    }, simplify = FALSE)

    V_h.0k <- do.call(cbind, V_h.0k.list)
    # diag(V_c.kk.dot) <- V_c.k


    # ## cov(tau_0, tau_k): V_c^[0,k]
    # V_c.0k <- sapply(data_hc_i, function(x){
    #   var_mu1_AIPW_hat +
    #     sum((X_c%*%x$beta_hat_0.c - sum(x$mu_i)/n_c)*
    #           (X_c%*%beta_hat_A0_c - sum(mu0_AIPW_i)/n_c))/n_c -
    #     sum((X_c%*%x$beta_hat_0.c - sum(x$mu_i)/n_c)*
    #           (X_c%*%beta_hat_A1_c - sum(mu1_AIPW_i)/n_c))/n_c -
    #     sum((X_c%*%beta_hat_A0_c - sum(mu0_AIPW_i)/n_c)*
    #           (X_c%*%beta_hat_A1_c - sum(mu1_AIPW_i)/n_c))/n_c
    # })

    # ## cov(tau_k, tau_k'): V_c^[k,k']
    # V_c.kk.dot <- sapply(data_hc_i, function(x)
    # {
    #   sapply(data_hc_i, function(y)
    #   {
    #     var_mu1_AIPW_hat +
    #       sum((X_c%*%x$beta_hat_0.c - sum(x$mu_i)/n_c)*
    #             (X_c%*%y$beta_hat_0.c - sum(y$mu_i)/n_c))/n_c-
    #       sum((X_c%*%x$beta_hat_0.c - sum(x$mu_i)/n_c)*
    #             (X_c%*%beta_hat_A1_c - sum(mu1_AIPW_i)/n_c))/n_c-
    #       sum((X_c%*%beta_hat_A1_c - sum(mu1_AIPW_i)/n_c)*
    #             (X_c%*%y$beta_hat_0.c - sum(y$mu_i)/n_c))/n_c
    #   })
    # })
    n_hi <- sapply(data_hc_i, function(hc_i)hc_i$n_h)
    # # approximation
    # sd_ACW_hat_i_s1 <- mapply(function(V_c, V_h, n_h){
    #   var_hi <- c(rep(V_c, n_c),
    #               rep(V_h, n_h))
    #   sqrt(var_hi)}, V_c = V_c.k, V_h = V_h.k, n_h = n_hi)
    # tau_hc_i_s1 <- mapply(function(x,y)x/y,
    #                       x = tau_hc_i, y = sd_ACW_hat_i_s1)
    # # direct standardization
    # tau_hc_i_s2 <- mapply(function(x,y)x/y,
    #                       x = tau_hc_i, y = sd_ACW_hat_i)

    # need to use Sigma_tau^{-1/2} for standardization
    # use rearranged tau for adaptive lasso
    Sigma2_tau_c <- cbind(c(unname(V_c.0), unname(V_c.0k)),
                          rbind(unname(V_c.0k),
                                unname(V_c.kk.dot)))

    # right.upper.corner <- kronecker(rep(1,K*strata+1), V_h.0k)
    # left.bottom.corner <- t(right.upper.corner)
    Sigma_Y <- as.matrix(bdiag(kronecker(Sigma2_tau_c,
                                         eye(n_c)),
                               diag(rep(V_h.k, times = n_hi))))
    # Sigma_Y <- rbind(cbind(kronecker(Sigma2_tau_c,
    #                                  eye(n_c)), right.upper.corner),
    #                  cbind(left.bottom.corner, diag(rep(V_h.k, times = n_hi))))

    library(expm)
    # expm::sqrtm(Sigma2_tau_c)
    # pracma::sqrtm(Sigma2_tau_c)$Binv
    # approximated inverse
    Sigma2_tau_c_1.2_inv <- expm::sqrtm(Sigma2_tau_c)#$Binv
    Sigma2_V_c_1.2 <- kronecker(Sigma2_tau_c_1.2_inv,
                                eye(n_c))
    Sigma2_V_h_1.2 <- diag(rep(V_h.k**(-1/2), times = n_hi))
    Sigma_Y_1.2_inv <- as.matrix(bdiag(Sigma2_V_c_1.2,
                                       Sigma2_V_h_1.2))
    # rearrange the response
    tau.hat <- c(tau_0_i,
                 lapply(tau_hc_i, function(x)x[1:n_c])%>%unlist(),
                 lapply(tau_hc_i, function(x)x[-(1:n_c)])%>%unlist())

    # tapply(tau.hat, c(rep(0, n_c),
    #                   rep(1:K, each = n_c),
    #                   rep((K+1):(2*K), times = n_hi)), sum)/n_c

    tau.tilde <- Sigma_Y_1.2_inv%*%tau.hat


    return(list(EST.ACW = data_hc_i, n_hi = n_hi,
                sd_ACW_hat_i = sd_ACW_hat_i,
                Sigma2_tau_c = Sigma2_tau_c, Sigma_Y = Sigma_Y,
                V_h.k = V_h.k,
                Sigma_Y_1.2_inv = Sigma_Y_1.2_inv,
                tau.hat = tau.hat, tau_hc_i = tau_hc_i,
                tau.tilde = tau.tilde))
  }

  EST.ACW.res <- EST.ACW.FUN(X_c = X_c, Y_c = Y_c,
                             A_c = A_c, prob_A = prob_A,
                             data_hc = data_hc)

  ###################################
  # sd: sd_AIPW_hat
  tau_hat_AIPW <- sum(tau_0_i)/n_c
  # sd: sd_ACW_hat_i
  tau_hat_ACW <- sapply(EST.ACW.res$tau_hc_i,
                        function(tau_i)sum(tau_i)/n_c)



  # Sigma_Y <- as.matrix(bdiag(kronecker(Sigma2_tau_c,
  #                                      eye(n_c)),
  #                            diag(rep(V_h.k, times = n_hi))))

  if(lasso == 'group')
  {
    # group-level integration
    Z <- factor(c(rep(0, n_c*(K*strata+1)),
                  rep(1:(K*strata), times = EST.ACW.res$n_hi)))
    Z.mat <- model.matrix(~Z+0)
    Z.mat.tilde <- EST.ACW.res$Sigma_Y_1.2_inv%*%Z.mat

    ## OLS for gamma for adaptive lasso
    gamma.hat <- ginv(t(Z.mat.tilde)%*%Z.mat.tilde)%*%
      t(Z.mat.tilde)%*%EST.ACW.res$tau.tilde
    # tuning parameter
    omega <- 1

    # construct the adaptive weights (group-level)
    b_k.w <- c(0, 1/abs(gamma.hat[-1])**omega)
    # the penalty factors are rescaled to sum to nvars
    b_k.w <- b_k.w/sum(b_k.w)*ncol(Z.mat.tilde)
    b_k.w <- sapply(b_k.w, function(x)ifelse(is.na(x), Inf, x))
    # adaptive Lasso
    library(glmnet)
    cv.lasso <- cv.glmnet(Z.mat.tilde, EST.ACW.res$tau.tilde, family = 'gaussian',
                          alpha = 1,
                          penalty.factor=b_k.w,
                          intercept = FALSE,
                          standardize = FALSE,
                          standardize.response = FALSE,#,
                          # lambda = c(seq(1,
                          #                1e-5/sqrt(nrow(Z.mat.tilde)),
                          #                length.out = 500))
                          # lambda = c(seq(1/sqrt(nrow(Z.mat.tilde)),
                          #                1.1e-2/sqrt(nrow(Z.mat.tilde)),
                          #                length.out = 500),
                          #            seq(1e-8/sqrt(nrow(Z.mat.tilde)),
                          #                1e-10/sqrt(nrow(Z.mat.tilde)),
                          #                length.out = 100))
                          lambda = seq(nrow(Z.mat.tilde)**-(1+omega)/2,
                                       nrow(Z.mat.tilde)**(-1/2),
                                       length.out = 100)
    )
    if(lambda == 'min')
    {
      lambda.opt <- cv.lasso$lambda.min
    }else{
      lambda.opt <- cv.lasso$lambda.1se
    }
    b_k_lasso <- coef(cv.lasso, s='lambda.min')[-1]
  }else{
    if(lasso == 'individual'){
      ## individual-level integration
      Z <- factor(c(rep(0, n_c*(K*strata+1)),
                    rep(1:sum(EST.ACW.res$n_hi))))

      # Z <- factor(c(rep(0, n_c),
      #               rep(1:(K*strata), each = n_c),
      #               rep(((K*strata)+1):((K*strata)+sum(n_hi)))))
      Z.mat <- model.matrix(~Z+0)
      Z.mat.tilde <- EST.ACW.res$Sigma_Y_1.2_inv%*%Z.mat

      ## OLS for gamma for adaptive lasso
      gamma.hat <- ginv(t(Z.mat.tilde)%*%Z.mat.tilde)%*%
        t(Z.mat.tilde)%*%EST.ACW.res$tau.tilde
      # adaptive Lasso
      library(glmnet)
      # set.seed(1234)
      # lambda/sqrt{n}=0, lambda = infty
      # tapply(tau.hat, c(rep(1:(1+K*strata), each = n_c),
      #                   rep(((1+K*strata)+1):((1+2*K*strata)), times = n_hi)), mean)
      # tuning parameter
      # ------------------------------------
      ## two-dimensional cross validation (ada-lasso)
      omega.vector <- c(1,2)#rev(seq(1, 2, length.out = 10)) # choose for 10 or 50

      # tuning the parameter
      cv.lasso.vec <- sapply(omega.vector, function(omega){

        lambda.vector <- seq(nrow(Z.mat.tilde)**-(1+omega)/2,
                             nrow(Z.mat.tilde)**(-1/2),
                             length.out = 100)
        # construct the adaptive weights (subject-level)
        b_k.w <- c(0,
                   1/abs(gamma.hat[-1])**omega)
        # b_k.w <- c(rep(0, K*strata+1),
        #            1/abs(gamma.hat[-(1:(K*strata+1))])**omega)

        # the penalty factors are rescaled to sum to nvars
        b_k.w <- b_k.w/sum(b_k.w)*ncol(Z.mat.tilde)
        b_k.w <- sapply(b_k.w, function(x)ifelse(is.na(x), Inf, x))


        fitLasso <- glmnet(Z.mat.tilde, EST.ACW.res$tau.tilde, family = 'gaussian',
                           alpha = 1,
                           penalty.factor=b_k.w,
                           intercept = FALSE,
                           standardize = FALSE,
                           standardize.response = FALSE,
                           lambda = lambda.vector)

        group_indicator <- c(rep(1:(K*strata), each = n_c),
                             rep(1:(K*strata), times = EST.ACW.res$n_hi))

        trade.off <- apply(predict(fitLasso, newx = Z.mat.tilde, type = 'coef'), 2,
                           function(x)
                           {
                             tau.ACW.lin <- EST.ACW.res$tau.hat[-c(1:n_c)]*
                               c(rep(x[2]!=0, (K*strata)*n_c), x[-(1:2)]!=0)
                             bias.hi <- tapply(tau.ACW.lin, group_indicator,
                                               function(x)sum(x)/n_c) -
                               tau_hat_AIPW
                             mean(bias.hi**2)
                           })+ # bias
          apply(predict(fitLasso, newx = Z.mat.tilde, type = 'coef'), 2,
                function(x){
                  # variance should use the score function
                  tau.ACW.score.all <- lapply(EST.ACW.res$EST.ACW, function(x)x$tau_ACWi)%>%unlist()
                  tau.ACW.lin <- tau.ACW.score.all*#EST.ACW.res$tau.hat[-c(1:n_c)]*
                    c(rep(x[2]!=0, (K*strata)*n_c), x[-(1:2)]!=0)
                  var.hi <- tapply(tau.ACW.lin, group_indicator,
                                   function(x) sum(x-sum(x)/n_c)**2/n_c**2)
                  mean(var.hi)
                }) # variance
        min(trade.off)
        # lambda.opt <- lambda.vector[which.min(trade.off)]


        # cv.lasso <- cv.glmnet(Z.mat.tilde, EST.ACW.res$tau.tilde, family = 'gaussian',
        #                       alpha = 1,
        #                       penalty.factor=b_k.w,
        #                       intercept = FALSE,
        #                       standardize = FALSE,
        #                       standardize.response = FALSE,#,
        #                       # lambda = c(seq(1,
        #                       #                1e-5/sqrt(nrow(Z.mat.tilde)),
        #                       #                length.out = 500))
        #                       # lambda = c(seq(1/sqrt(nrow(Z.mat.tilde)),
        #                       #                1.1e-2/sqrt(nrow(Z.mat.tilde)),
        #                       #                length.out = 500),
        #                       #            seq(1e-8/sqrt(nrow(Z.mat.tilde)),
        #                       #                1e-10/sqrt(nrow(Z.mat.tilde)),
        #                       #                length.out = 100))
        #                       lambda = seq(nrow(Z.mat.tilde)**-(1+omega)/2,
        #                                    nrow(Z.mat.tilde)**(-1/2),
        #                                    length.out = 100),
        #                       nfolds = 10)
        # ifelse(lambda == 'min', cv.lasso$cvm[cv.lasso$index['min',]],
        #        cv.lasso$cvm[cv.lasso$index['1se',]])
      })
      # ------------------------------------
      # select the one that minimize the cv error
      omega.opt <- omega.vector[which.min(cv.lasso.vec)]
      # construct the adaptive weights (subject-level)
      # b_k.w <- c(rep(0, K*strata+1),
      #            1/abs(gamma.hat[-(1:(K*strata+1))])**omega.opt)
      b_k.OLS <- c(0,
                   1/abs(gamma.hat[-1])**omega.opt)
      # the penalty factors are rescaled to sum to nvars (internally)
      b_k.w <- b_k.OLS/sum(b_k.OLS)*ncol(Z.mat.tilde)
      b_k.w <- sapply(b_k.w, function(x)ifelse(is.na(x), Inf, x))

      # # cross-validation tuning against AIPW (custom metric)
      library(caret)
      # custom metric
      # custom_CV <- function(data, lev = NULL, model = NULL)
      # {
      #   # customize metric based on rowIndex
      #   # xi_0 <- sapply(data$rowIndex, function(x)
      #   # {
      #   #   if(x<=((K+1)*n_c)){tau_hat_AIPW}else{
      #   #     b_k.OLS[x-(K+1)*n_c]}
      #   # })
      #   out <- mean((data$pred - data$obs)**2)
      #
      #   stats <- defaultSummary(data, lev = lev, model = model)
      #   c(custom = out, stats)
      # }
      #
      # fitControl <- trainControl(## 10-fold CV
      #   method = "cv",
      #   number = 2,
      #   summaryFunction = custom_CV)
      # fitGrid <- expand.grid(lambda = rev(seq(nrow(Z.mat.tilde)**-(1+omega.opt)/2,
      #                                     nrow(Z.mat.tilde)**(-1/2),
      #                                     length.out = 100)),
      #                        alpha = 1)
      #
      # fitGlmnet <- train(x = Z.mat.tilde,
      #       y = as.numeric(EST.ACW.res$tau.tilde),
      #       method = 'glmnet',
      #       trControl = fitControl,
      #       tuneGrid = fitGrid,
      #       # for glmnet
      #       penalty.factor=b_k.w,
      #       intercept = FALSE,
      #       standardize = FALSE,
      #       standardize.response = FALSE)
      # predict(fitGlmnet$finalModel,
      #         newx = Z.mat.tilde, type = 'response', s= 100)
      # predict(fitGlmnet$finalModel, type = 'coef', s= 1e-100)

      # library(recipes)
      # recipe(x = Z.mat.tilde)%>%
      #   add_role(as.numeric(EST.ACW.res$tau.tilde),
      #            new_role = 'outcome')
      lambda.vector <- seq(nrow(Z.mat.tilde)**-(1+omega.opt)/2,
                           nrow(Z.mat.tilde)**(-1/2),
                           length.out = 100)
      fitLasso <- glmnet(Z.mat.tilde, EST.ACW.res$tau.tilde, family = 'gaussian',
                         alpha = 1,
                         penalty.factor=b_k.w,
                         intercept = FALSE,
                         standardize = FALSE,
                         standardize.response = FALSE,
                         lambda = lambda.vector)

      group_indicator <- c(rep(1:(K*strata), each = n_c),
                           rep(1:(K*strata), times = EST.ACW.res$n_hi))

      trade.off <- apply(predict(fitLasso, newx = Z.mat.tilde, type = 'coef'), 2,
                         function(x)
                         {
                           tau.ACW.lin <- EST.ACW.res$tau.hat[-c(1:n_c)]*
                             c(rep(x[2]!=0, (K*strata)*n_c), x[-(1:2)]!=0)
                           bias.hi <- tapply(tau.ACW.lin, group_indicator,
                                             function(x)sum(x)/n_c) -
                             tau_hat_AIPW
                           mean(bias.hi**2)
                         })+ # bias
        apply(predict(fitLasso, newx = Z.mat.tilde, type = 'coef'), 2,
              function(x){
                # variance should use the score function
                tau.ACW.score.all <- lapply(EST.ACW.res$EST.ACW, function(x)x$tau_ACWi)%>%unlist()
                tau.ACW.lin <- tau.ACW.score.all*#EST.ACW.res$tau.hat[-c(1:n_c)]*
                  c(rep(x[2]!=0, (K*strata)*n_c), x[-(1:2)]!=0)
                var.hi <- tapply(tau.ACW.lin, group_indicator,
                                 function(x) sum(x-sum(x)/n_c)**2/n_c**2)
                mean(var.hi)
              }) # variance

      lambda.opt <- lambda.vector[which.min(trade.off)]

      # cv.lasso <- cv.glmnet(Z.mat.tilde, EST.ACW.res$tau.tilde, family = 'gaussian',
      #                       alpha = 1,
      #                       penalty.factor=b_k.w,
      #                       intercept = FALSE,
      #                       standardize = FALSE,
      #                       standardize.response = FALSE,#,
      #                       # lambda = c(seq(1,
      #                       #                1e-5/sqrt(nrow(Z.mat.tilde)),
      #                       #                length.out = 500))
      #                       # lambda = c(seq(1/sqrt(nrow(Z.mat.tilde)),
      #                       #                1.1e-2/sqrt(nrow(Z.mat.tilde)),
      #                       #                length.out = 500),
      #                       #            seq(1e-8/sqrt(nrow(Z.mat.tilde)),
      #                       #                1e-10/sqrt(nrow(Z.mat.tilde)),
      #                       #                length.out = 100))
      #                       lambda = seq(nrow(Z.mat.tilde)**-(1+omega.opt)/2,
      #                                    nrow(Z.mat.tilde)**(-1/2),
      #                                    length.out = 100),
      #                       nfolds = 10)
      #
      # if(lambda == 'min')
      # {
      #   lambda.opt <- cv.lasso$lambda.min
      # }else{
      #   lambda.opt <- cv.lasso$lambda.1se
      # }

      b_k_lasso <- coef(fitLasso, s= lambda.opt)[-1]
      b_k_lasso
    }
  }

  b_k_lasso
  # subtract the idx for the benchmarks
  hc.idx.lasso <- unname(which(b_k_lasso==0)-1)
  sd_ACW_hat_i <- EST.ACW.res$sd_ACW_hat_i
  if(identical(integer(0), hc.idx.lasso)|identical(numeric(0), hc.idx.lasso)){
    # select null set
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
    if(lasso == 'individual'){
      hc.val.lasso.list <- split(b_k_lasso[-1], f = factor(rep(1:(K*strata),
                                                               times = EST.ACW.res$n_hi)))

      hc.val.lasso.list <- structure(hc.val.lasso.list, names = NULL)
      if(refit){


        X0.lasso <- do.call(rbind, mapply(function(hc_i, hc.lasso.val){
          hc.lasso.idx <-  which(hc.lasso.val == 0)
          hc_i$X0[hc.lasso.idx,]
        }, hc_i = data_hc, hc.lasso.val = hc.val.lasso.list,
        SIMPLIFY = FALSE))

        Y0.lasso <- mapply(function(hc_i, hc.lasso.val){
          hc.lasso.idx <-  which(hc.lasso.val == 0)
          hc_i$Y0[hc.lasso.idx]
        }, hc_i = data_hc, hc.lasso.val = hc.val.lasso.list,
        SIMPLIFY = FALSE)%>%unlist()

        # data_hc.lasso <- mapply(function(hc_i, hc.lasso.val){
        #   hc.lasso.idx <-  which(hc.lasso.val == 0)
        #   list(X0 = hc_i$X0[hc.lasso.idx,],
        #        Y0 = hc_i$Y0[hc.lasso.idx])
        # }, hc_i = data_hc, hc.lasso.val = hc.val.lasso.list,
        # SIMPLIFY = FALSE)


        data_hc.lasso <- list(list(X0 = X0.lasso, Y0 = Y0.lasso))

        EST.ACW.res.lasso <- EST.ACW.FUN(X_c = X_c, Y_c = Y_c,
                                         A_c = A_c, prob_A = prob_A,
                                         data_hc = data_hc.lasso)

        # construct the eta.all for the post-lasso selection
        q_hat.all <-  lapply(EST.ACW.res$EST.ACW, function(x)x$q_hat)
        q_hat.all.star <-  lapply(EST.ACW.res.lasso$EST.ACW, function(x)x$q_hat)


        # construct the eta for historical controls
        ## model-assisted calibration
        # hc.val.lasso.list
        # mapply(data_hc, function(hc_i, hc.val.lasso){
        #   X_h <- hc_i$X0
        #   if(is.data.frame(X_h)){
        #     X_h <- model.matrix(~., X_h)
        #   }else{X_h <- as.matrix(X_h)}
        #
        #   Y_h0.pred <- X_h%*%beta_hat_A0_c + hc.val.lasso/sqrt(n_c)
        #   Y_c0.pred <- X_c%*%beta_hat_A0_c
        #
        #   lambda.hat.refit <- dfsane(par = rep(0, 2), # initial points
        #          fn=cal_eqs,
        #          X_c = cbind(1, Y_c0.pred), X_h = cbind(1, Y_h0.pred),
        #          control = list(trace=T,
        #                         NM=T,
        #                         BFGS=F,
        #                         tol=1.e-8,
        #                         maxit = 500))$par
        #   q_hat.refit <- as.numeric(exp(cbind(1, Y_h0.pred)%*%lambda.hat.refit))
        #
        # }, hc_i = data_hc, hc.val.lasso = hc.val.lasso.list)

        eta_h.all <- rep(0, times = sum(EST.ACW.res$n_hi))
        # eta_h.all[hc.idx.lasso] <- q_hat.all.star%>%unlist()
        eta_h.all[hc.idx.lasso] <- (q_hat.all%>%unlist())[hc.idx.lasso]
        eta_h.all <- eta_h.all/q_hat.all%>%unlist()/n_c

        # construct the eta for the RCT
        ## selected components
        Sigma2_V_c_1 <- kronecker(ginv(EST.ACW.res$Sigma2_tau_c),
                                  eye(n_c))
        ## subject-level estimator
        eta_c.all <- sum(apply(Sigma2_V_c_1, 2, sum))^{-1}*
          apply(Sigma2_V_c_1, 2, sum)
        # eta_c.all <- c(rep(0, n_c), rep(1/n_c, n_c*(K*strata)))
        eta.ind.all <- c(eta_c.all, eta_h.all)

        # update the Sigma_Y
        EST.ACW.res$Sigma_Y <- as.matrix(bdiag(kronecker(EST.ACW.res$Sigma2_tau_c,
                                                         eye(n_c)),
                                               diag(rep(EST.ACW.res.lasso$V_h.k,
                                                        times = EST.ACW.res$n_hi))))

        # approximated inverse
        Sigma2_tau_c_1.2_inv <- expm::sqrtm(EST.ACW.res$Sigma2_tau_c)#$Binv
        Sigma2_V_c_1.2 <- kronecker(Sigma2_tau_c_1.2_inv,
                                    eye(n_c))
        Sigma2_V_h_1.2 <- diag(rep((EST.ACW.res.lasso$V_h.k)**(-1/2),
                                   times = EST.ACW.res$n_hi))
        EST.ACW.res$Sigma_Y_1.2_inv <- as.matrix(bdiag(Sigma2_V_c_1.2,
                                                       Sigma2_V_h_1.2))

        EST.ACW.res$tau.tilde <- EST.ACW.res$Sigma_Y_1.2_inv%*%EST.ACW.res$tau.hat

        # # updata Z.matrix
        # Z <- factor(c(rep(0, n_c*(1+1)),
        #               rep(1:length(hc.idx.lasso))))
        #
        # # Z <- factor(c(rep(0, n_c),
        # #               rep(1:(K*strata), each = n_c),
        # #               rep(((K*strata)+1):((K*strata)+sum(n_hi)))))
        # Z.mat <- model.matrix(~Z+0)
        # Z.mat.tilde <- EST.ACW.res.lasso$Sigma_Y_1.2_inv%*%Z.mat




        # # compute the final estimator
        # tau_final <- EST.ACW.res$tau.hat%*%eta.ind.all
        # var_final <- t(eta.ind.all)%*%EST.ACW.res.lasso$Sigma_Y%*%(eta.ind.all)*n_c#*nrow(Z.mat.tilde)
        # sd_final <- sqrt(var_final)
        #
        #
        #
        # # subsitute sqrt(eta%*%Sigma_Y%*%eta) with sd_final/sqrt(n_c)
        # tau_upper <- tau_final + 1.96 * sd_final/sqrt(n_c)
        # tau_lower <- tau_final - 1.96 * sd_final/sqrt(n_c)
      }else{ # not refit
        # data_hc.lasso <- data_hc
        # EST.ACW.res.lasso <- EST.ACW.res

        # construct the eta for AIPW
        eta.AIPW <- c(rep(1/n_c, n_c), # AIPW
                      rep(0, n_c*(K*strata)), # ACW main part (0)
                      rep(0, length(b_k_lasso[-1]))) # ACW augmented part (0)

        # construct the eta for ACW
        eta.ACW.list <- lapply(1:(K*strata), function(idx){
          # RCT group
          eta.ACW_c <- c(rep(0, n_c*idx),
                         rep(1/n_c, n_c),
                         rep(0, n_c*(K*strata-idx)))
          # RWD group
          eta.ACW_h <- numeric(length(b_k_lasso[-1]))
          n_hi <- EST.ACW.res$n_hi
          start.idx <- sum(n_hi[0:(idx-1)])
          eta.ACW_h[start.idx +
                      which(hc.val.lasso.list[[idx]]==0)] <- 1/n_c
          c(eta.ACW_c, eta.ACW_h)
        })

        # AIPW var
        var.tau.AIPW <- eta.AIPW%*%EST.ACW.res$Sigma_Y%*%eta.AIPW*n_c
        # AIPW-ACW var
        var.tau.AIPW_ACW <- sapply(eta.ACW.list, function(eta.ACW){
          eta.AIPW%*%EST.ACW.res$Sigma_Y%*%eta.ACW*n_c
        })
        # ACW-ACW var
        var.tau.ACW_ACW <- sapply(eta.ACW.list, function(eta.ACW.x){
          sapply(eta.ACW.list, function(eta.ACW.y){
            eta.ACW.x%*%EST.ACW.res$Sigma_Y%*%eta.ACW.y*n_c})
        })
        # construct the matrix
        Sigma2_tau <- cbind(c(var.tau.AIPW, var.tau.AIPW_ACW),
                            rbind(t(var.tau.AIPW_ACW), var.tau.ACW_ACW))
        d_A_n <- c(rep(1, 1+K*strata)%*%ginv(Sigma2_tau)%*%rep(1, 1+K*strata))^{-1} *
          ginv(Sigma2_tau)%*%rep(1, 1+K*strata)

        # construct the eta for the combined estimator
        eta.AIPW.all <- eta.AIPW*d_A_n[1]

        eta.ACW.all <- do.call(rbind, mapply(function(weight, eta.ACW){
          weight*eta.ACW
        }, weight = d_A_n[-1],
        eta.ACW = eta.ACW.list,
        SIMPLIFY = FALSE))%>%apply(2, sum)
        # eta.ind.all.var <- sum(apply(Sigma2_V_c_1, 2, sum))^{-1}*
        #   c(apply(Sigma2_V_c_1, 2, sum),
        #     sapply(hc.idx.lasso.ind, function(x)ifelse(x==0, 1, 0)))

        eta.ind.all <- eta.AIPW.all + eta.ACW.all


      }
    }
    if(lasso == 'group')
    {

      # ------------------------------
      # group-level estimator
      n.lasso <- length(hc.idx.lasso) + 1 # selected HC plus RCT
      # compute variance between tau_k and tau_k'
      Sigma2_tau <- sapply(EST.ACW.res.lasso$EST.ACW[hc.idx.lasso], function(x)
      {
        sapply(EST.ACW.res.lasso$EST.ACW[hc.idx.lasso], function(y)
        {
          cov <- var_mu1_AIPW_hat +
            sum((X_c%*%x$beta_hat_0.c - sum(x$mu_i)/n_c)*
                  (X_c%*%y$beta_hat_0.c - sum(y$mu_i)/n_c))/n_c-
            sum((X_c%*%x$beta_hat_0.c - sum(x$mu_i)/n_c)*
                  (X_c%*%beta_hat_A1_c - sum(mu1_AIPW_i)/n_c))/n_c-
            sum((X_c%*%beta_hat_A1_c - sum(mu1_AIPW_i)/n_c)*
                  (X_c%*%y$beta_hat_0.c - sum(y$mu_i)/n_c))/n_c
          # sqrt(cov)
        })
      })
      Sigma2_tau <- as.matrix(Sigma2_tau)
      diag(Sigma2_tau) <- sd_ACW_hat_i[hc.idx.lasso]**2
      # compute variance between tau_k and tau_0
      sigma2_tau_k_0 <- sapply(data_hc_i[hc.idx.lasso], function(x){
        cov <- var_mu1_AIPW_hat +
          sum((X_c%*%x$beta_hat_0.c - sum(x$mu_i)/n_c)*
                (X_c%*%beta_hat_A0_c - sum(mu0_AIPW_i)/n_c))/n_c -
          sum((X_c%*%x$beta_hat_0.c - sum(x$mu_i)/n_c)*
                (X_c%*%beta_hat_A1_c - sum(mu1_AIPW_i)/n_c))/n_c -
          sum((X_c%*%beta_hat_A0_c - sum(mu0_AIPW_i)/n_c)*
                (X_c%*%beta_hat_A1_c - sum(mu1_AIPW_i)/n_c))/n_c
        # sqrt(cov)
      })
      # concatenate the variance togather
      Sigma2_tau <- cbind(c(sd_AIPW_hat**2, unname(sigma2_tau_k_0)),
                          rbind(unname(sigma2_tau_k_0), Sigma2_tau))
      # Sigma2_tau <- Sigma_tau**2#/n_c
      # compute the weights
      d_A_n <- c(rep(1, n.lasso)%*%ginv(Sigma2_tau)%*%rep(1, n.lasso))^{-1} *
        ginv(Sigma2_tau)%*%rep(1, n.lasso)
      # final estimator
      tau_hat_selected <- c(tau_hat_AIPW,
                            tau_hat_ACW[hc.idx.lasso])
      tau_final <- tau_hat_selected%*%d_A_n
      var_final <- t(d_A_n)%*%Sigma2_tau%*%(d_A_n)
      sd_final <- sqrt(var_final)
      # construct eta
      eta_1_c <- lapply(1:(K*strata), function(k)
      {
        if(k%in%hc.idx.lasso){
          # map the k-th control to the post-lasso selection
          idx <- which(hc.idx.lasso==k)
          rep(d_A_n[idx+1]/(n_c), n_c)
        }else{
          rep(0, n_c)
        }})
      eta_1_h <- lapply(1:(K*strata), function(k)
      {
        if(k%in%hc.idx.lasso){
          # map the k-th control to the post-lasso selection
          idx <- which(hc.idx.lasso==k)
          rep(d_A_n[idx+1]/(n_c), n_hi[k])
        }else{
          rep(0, n_hi[k])
        }
      })
      eta.ind.all <- c(rep(d_A_n[1]/n_c, n_c),
                       eta_1_c%>%unlist(),
                       eta_1_h%>%unlist())
      eta.ind.all%*%tau.hat # equal to tau_hat_selected%*%d_A_n
      # hc.idx.deleted <- setdiff(1:(K*strata), hc.idx.lasso)
      # -----------------------------------------------------
    }
  }




  # extract the results A^c
  # AIPW and ACW
  # compute each tau_hat
  # sd: sd_AIPW_hat
  # tau_hat_AIPW <- sum(tau_0_i)/n_c
  # # sd: sd_ACW_hat_i
  # tau_hat_ACW <- sapply(tau_hc_i,
  #                       function(tau_i)sum(tau_i)/n_c)

  # compute the final estimator
  tau_final <- EST.ACW.res$tau.hat%*%eta.ind.all
  # var_ACW_lasso <- (sum(apply(kronecker(ginv(Sigma2_tau_c),
  #                                       eye(n_c)), 2, sum))/n_c)**{-1}

  # revised variance estimation
  # tau.hat.RWD <- lapply(data_hc_i, function(x)x$X_h%*%beta_hat_A0_c)%>%unlist()
  #
  # var_final <- (sum(apply(kronecker(ginv(Sigma2_tau_c),
  #                      eye(n_c)), 2, sum))/n_c+
  #     sum(tau.hat.RWD[hc.idx.lasso]**2)/length(hc.idx.lasso)/n_c)**{-1}

  # var_final <- t(eta.ind.all)%*%Sigma_Y%*%(eta.ind.all)*n_c#*nrow(Z.mat.tilde)

  var_final <- t(eta.ind.all)%*%EST.ACW.res$Sigma_Y%*%(eta.ind.all)*n_c#*nrow(Z.mat.tilde)

  sd_final <- sqrt(var_final)

  # use conditional distribution to revise the variance
  A_n <- which(b_k_lasso!=0);A_n.c <- which(b_k_lasso==0)
  # compute the sign matrix
  if(length(A_n)>=5) # if the non-zero parameters are too many (at most 2**10 candidates)
  {
    ## conditional on the sign (subject-level)
    sgn.mat <- matrix(abs(b_k_lasso[A_n])/b_k_lasso[A_n], nrow = 1)}else{
      ## construct the sign matrix (group-level)
      sgn.mat <- as.matrix(do.call(expand.grid,
                                   rep(list(c(-1,1)), length(A_n))))
    }

  L_U.sign <- function(sgn){

    # sgn <- abs(b_k_lasso[A_n])/b_k_lasso[A_n] # uncomment it for less conservative CI

    sgn.diag <- if(length(sgn)==1){sgn
    }else{diag(sgn)}
    # segement the standardized Z matrix by A_n and A_n^c
    Z.mat.tilde.A_n <- Z.mat.tilde[,A_n,drop=F];
    Z.mat.tilde.A_n.c <- Z.mat.tilde[,-A_n,drop=F]
    # segement the weight
    w.A_n <- b_k.w[A_n]; w.A_n.c <- b_k.w[-A_n]
    w.A_n.diag <- if(length(w.A_n)==1){w.A_n}else{diag(w.A_n)}
    w.A_n.diag.inv <- if(length(w.A_n)==1){w.A_n^{-1}}else{diag(w.A_n^{-1})}
    w.A_n.c.diag <- if(length(w.A_n.c)==1){w.A_n.c}else{diag(w.A_n.c)}
    w.A_n.c.diag.inv <- if(length(w.A_n.c)==1){w.A_n.c^{-1}}else{diag(w.A_n.c^{-1})}
    # lambda:
    lambda_n <- lambda.opt * nrow(Z.mat.tilde)
    # construct the tau.tilde = Y
    tau.tilde <- EST.ACW.res$tau.tilde

    # construct the polyhedra constraint
    C_1 <- -sgn.diag%*%ginv(t(Z.mat.tilde.A_n)%*%Z.mat.tilde.A_n)%*%
      t(Z.mat.tilde.A_n)
    h_1 <- -lambda_n*sgn.diag%*%
      ginv(t(Z.mat.tilde.A_n)%*%Z.mat.tilde.A_n)%*%w.A_n.diag%*%sgn
    #%*%tau.tilde
    Z.mat.tilde.A_n.2.inv <- ginv(t(Z.mat.tilde.A_n)%*%Z.mat.tilde.A_n)
    C_01 <- lambda_n^{-1}*w.A_n.c.diag.inv%*%
      (t(Z.mat.tilde.A_n.c) -
         t(Z.mat.tilde.A_n.c)%*%Z.mat.tilde.A_n%*%
         Z.mat.tilde.A_n.2.inv%*%
         t(Z.mat.tilde.A_n))
    C_02 <- lambda_n^{-1}*w.A_n.c.diag.inv%*%
      (-t(Z.mat.tilde.A_n.c) +
         t(Z.mat.tilde.A_n.c)%*%Z.mat.tilde.A_n%*%
         Z.mat.tilde.A_n.2.inv%*%
         t(Z.mat.tilde.A_n))
    h_01 <- 1-w.A_n.c.diag.inv%*%t(Z.mat.tilde.A_n.c)%*%Z.mat.tilde.A_n%*%
      Z.mat.tilde.A_n.2.inv%*%w.A_n.diag%*%sgn
    h_02 <- 1+w.A_n.c.diag.inv%*%t(Z.mat.tilde.A_n.c)%*%Z.mat.tilde.A_n%*%
      Z.mat.tilde.A_n.2.inv%*%w.A_n.diag%*%sgn
    C_s <- rbind(C_01, C_02,
                 C_1)%*%EST.ACW.res$Sigma_Y_1.2_inv
    h_s <- rbind(h_01, h_02,
                 h_1)




    # ---------------------------
    # sanity check for KKT condition
    t(Z.mat.tilde.A_n)%*%
      (tau.tilde - Z.mat.tilde.A_n%*%b_k_lasso[A_n]) # equal to 0
    lambda_n <- lambda.opt * nrow(Z.mat.tilde)
    lambda_n * w.A_n * sign(b_k_lasso[A_n])

    # should be bounded by lambda_n * w.A_n.c
    t(Z.mat.tilde.A_n.c)%*%
      (tau.tilde - Z.mat.tilde.A_n%*%b_k_lasso[A_n])
    lambda_n * w.A_n.c
    # ---------------------------

    f_eta <- (c(eta.ind.all%*%EST.ACW.res$Sigma_Y%*%eta.ind.all)**{-1})*
      EST.ACW.res$Sigma_Y%*%eta.ind.all
    r <- EST.ACW.res$tau.hat-f_eta%*%(eta.ind.all%*%EST.ACW.res$tau.hat)

    den <- C_s%*%f_eta
    num <- h_s - C_s%*%r

    idx.neg <- which(den<0); idx.pos <- which(den>0)
    L_s <- max(num[idx.neg]/den[idx.neg])
    U_s <- min(num[idx.pos]/den[idx.pos])

    # if(L_s>U_s)
    # {
    #   temp <- L_s
    #   L_s <- U_s
    #   U_s <- temp
    # }
    # if(L_s>U_s) L_s <- U_s <- NA
    return(c(L_s = L_s,
             U_s = U_s))
  }

  L_U.mat <- apply(sgn.mat, 1, L_U.sign)

  L_s <- min(L_U.mat['L_s',], na.rm = T)
  U_s <- max(L_U.mat['U_s',], na.rm = T)

  # genetic algorithm for maximize (or minimize) lower bound and upper bound (to avoid computational efficiency)
  library(genalg)
  # U_s.best <- -rbga.bin(size = length(sgn),
  #          evalFunc =  function(sgn) -L_U.sign(sgn)['U_s'],
  #          popSize=200, iters=200, verbose = TRUE)$best
  # L_s.best <- rbga.bin(size = length(sgn),
  #                 evalFunc =  function(sgn) L_U.sign(sgn)['L_s'],
  #                 popSize=200, iters=200, verbose = TRUE)$best

  library(truncnorm)
  # subsitute sqrt(eta%*%Sigma_Y%*%eta) with sd_final/sqrt(n_c)
  tau_upper_initial <- tau_final + 1.96 * sd_final/sqrt(n_c)
  tau_lower_initial <- tau_final - 1.96 * sd_final/sqrt(n_c)

  # tau_cont_upper <- tau_upper_initial
  # tau_cont_lower <- tau_lower_initial

  # if the value is on the boundary, skip
  if(abs(eta.ind.all%*%EST.ACW.res$tau.hat - L_s)<1e-5 |
     abs(eta.ind.all%*%EST.ACW.res$tau.hat - U_s)<1e-5) # on the boundary
  {
    tau_cont_upper <- tau_upper_initial
    tau_cont_lower <- tau_lower_initial
  }else{
    # tau_cont_upper <- tau_upper_initial
    # tau_cont_lower <- tau_lower_initial
    # truncated mean function
    truncnorm.mean <- function(L, U, mean, sd)
    {
      truncnorm::etruncnorm(a = L, b = U,
                            mean = mean, sd = sd)*
        (pnorm(U, mean = mean, sd = sd) -
           pnorm(L, mean = mean, sd = sd))
    }

    truncnorm.dev <- function(x,
                              L, U, mean, sd)
    {
      denom <- (pnorm(U, mean = mean, sd = sd) -
                  pnorm(L, mean = mean, sd = sd))**2

      num <- 1/sd**2 *{truncnorm.mean(L = L, U = x,
                                      mean = mean, sd = sd)-
          mean * (pnorm(x, mean = mean, sd = sd) -
                    pnorm(L, mean = mean, sd = sd))} *
        (pnorm(U, mean = mean, sd = sd) -
           pnorm(L, mean = mean, sd = sd)) -
        1/sd**2 *{truncnorm.mean(L = L, U = U,
                                 mean = mean, sd = sd)-
            mean * (pnorm(U, mean = mean, sd = sd) -
                      pnorm(L, mean = mean, sd = sd))} *
        (pnorm(x, mean = mean, sd = sd) -
           pnorm(L, mean = mean, sd = sd))
      num/denom
    }


    tau_cont_upper <- nleqslv(x = c(tau_upper_initial),
                              fn=function(x)
                              {
                                ptruncnorm(eta.ind.all%*%EST.ACW.res$tau.hat,
                                           a=L_s, b=U_s,
                                           mean = x, sd = sd_final/sqrt(n_c))-
                                  0.05/2
                              },
                              # jac = function(x)
                              # {
                              #   truncnorm.dev(x = eta.ind.all%*%EST.ACW.res$tau.hat,
                              #                 L = L_s, U = U_s,
                              #                 mean = x,
                              #                 sd = sd_final/sqrt(n_c))
                              # },
                              method = 'Newton',
                              control = list(allowSingular = TRUE))$x

    tau_cont_lower <- nleqslv(x = c(tau_lower_initial),
                              fn=function(x)
                              {
                                ptruncnorm(eta.ind.all%*%EST.ACW.res$tau.hat,
                                           a=L_s, b=U_s,
                                           mean = x, sd = sd_final/sqrt(n_c))-
                                  (1-0.05/2)
                              },
                              # jac = function(x)
                              # {
                              #   truncnorm.dev(x = eta.ind.all%*%EST.ACW.res$tau.hat,
                              #                 L = L_s, U = U_s,
                              #                 mean = x,
                              #                 sd = sd_final/sqrt(n_c))
                              # },
                              method = 'Newton',
                              control = list(allowSingular = TRUE))$x



    # tau_cont_upper <- uniroot(function(x)
    # {
    #   (ptruncnorm(eta%*%tau.hat,
    #              a=L_s, b=U_s,
    #              mean = x, sd = sqrt(eta%*%Sigma_Y%*%eta))-
    #     0.05/2)
    # }, interval = c(L_s-0.1*abs(tau_lower_initial),
    #                 U_s+0.1*abs(tau_upper_initial)))$root
    #
    # tau_cont_lower <- uniroot(function(x)
    # {
    #   ptruncnorm(eta%*%tau.hat,
    #              a=L_s, b=U_s,
    #              mean = x, sd = sqrt(eta%*%Sigma_Y%*%eta))-
    #     (1-0.05/2)
    # }, interval = c(L_s-0.1*abs(tau_lower_initial),
    #                 U_s+0.1*abs(tau_upper_initial)))$root

  }
  return(list(n_c = n_c,
              est = list(AIPW = tau_hat_AIPW,
                         ACW = tau_hat_ACW,
                         ACW.lasso = tau_final),
              sd = list(AIPW = sd_AIPW_hat,
                        ACW = sd_ACW_hat_i,
                        ACW.lasso = sd_final),
              CI = list(lower = tau_cont_lower,
                        upper = tau_cont_upper),
              subset.idx = hc.idx.lasso,
              model = EST.ACW.res))
}
DATA_standardizer <- function(data)
{

  ## continuous
  data$AGEYR <- (data$AGEYR-mean(data$AGEYR, na.rm = TRUE))/sd(data$AGEYR, na.rm = TRUE) # age
  data$HGTCMBL <- (data$HGTCMBL-mean(data$HGTCMBL, na.rm = TRUE))/sd(data$HGTCMBL, na.rm = TRUE) # height
  # data$WGTKGBL <-  (data$WGTKGBL-mean(data$WGTKGBL, na.rm = TRUE))/sd(data$WGTKGBL, na.rm = TRUE) # weight
  data$BMIBL <- (data$BMIBL-mean(data$BMIBL, na.rm = TRUE))/sd(data$BMIBL, na.rm = TRUE) # BMI
  # data$EDCCNT_standarized <- (data$EDCCNT-mean(data$EDCCNT, na.rm = TRUE))/sd(data$EDCCNT, na.rm = TRUE) # years of education
  # data$DURDISONSYR <- (data$DURDISONSYR-mean(data$DURDISONSYR, na.rm = TRUE))/sd(data$DURDISONSYR, na.rm = TRUE) # time since onset of AD symptoms
  data$DURDIAGYR <- (data$DURDIAGYR-mean(data$DURDIAGYR, na.rm = TRUE))/sd(data$DURDIAGYR, na.rm = TRUE) # time since diagnosis of AD
  data$mmse_base <- (data$mmse_base-mean(data$mmse_base, na.rm = TRUE))/sd(data$mmse_base, na.rm = TRUE) # mini mental State examination

  ## categorical
  data$SEXSNM <- factor(data$SEXSNM, levels = c('M', 'F')) # gender, Male is 0
  # data$RACECAT <- as.factor(data$RACECAT) # race
  # data$TBBL <- as.factor(data$TBBL) # tobacco use
  data$ALCHLBL_C <- as.factor(data$ALCHLBL_C) # alcohol use
  # data$WRKSTAT1LNM <- as.factor(data$WRKSTAT1LNM) # work status
  # data$HISSTSBL <- as.factor(data$HISSTSBL) # modified Hachinski Ischemia Scale
  # data$GDSSEVBL <- as.factor(data$GDSSEVBL) # geriatric depression score
  # data$admed <- as.factor(data$admed) # AChEI or memantime use # 0, 1 and APOEGN
  data$hypertension <- as.factor(data$hypertension)
  data$diabetes <- as.factor(data$diabetes)
  data$treat <- as.factor(data$treat)

  return(data[, c(cov.names, y.names, 'TRT')])
}

# a real data application
RDA.LZAN <- function(seed, n_c_selected, y.names,
                     method = 'lm', ...){

  # week 80, adas_cog13 as the primary outcome
  outcomes_partial <- subset(outcomes, subset = viswk == 80,
                             select = c('USUBJID', y.names))
  data.combined <- left_join(baselines, outcomes_partial,
                             by = 'USUBJID')
  # separate the data
  data.RCT <- subset(data.combined, subset = cohort == 'LZAN')
  data.HC <- subset(data.combined, subset = cohort == 'Geras_EU')
  # selected name of baseline covariates
  cov.names <- c( # continuous
    'AGEYR', 'HGTCMBL',
    'BMIBL', 'EDCCNT',
    'DURDIAGYR',
    'mmse_base',
    # 'WGTKGBL_standarized',
    # categorical
    'SEXSNM',  'ALCHLBL_C', 'hypertension', 'diabetes', 'treat'
  ) #'RACECAT' #'WRKSTAT1LNM', 'HISSTSBL', 'GDSSEVBL' 'TBBL', 'DURDISONSYR_standarized', 'admed'


  # standardize the baseline covariates
  data.RCT.standardized <- DATA_standardizer(data = data.RCT)
  data.HC.standardized <- DATA_standardizer(data = data.HC)
  data.RCT.standardized$TRT <- sapply(data.RCT.standardized$TRT,
                                      function(x)ifelse(x=='Placebo',0,1))%>%as.vector()

  data.RCT.standardized$treat <- factor(data.RCT.standardized$treat,
                                        labels = c('Mild AD', 'Moderate'))
  # select the mild one
  data.RCT.standardized <-
    data.RCT.standardized[data.RCT.standardized$treat=='Mild AD',]
                          # !colnames(data.HC.standardized)%in%c('treat')]
  # recode the treatment and alchol use for HC
  data.HC.standardized$TRT <- 0
  data.HC.standardized$ALCHLBL_C <- sapply(data.HC.standardized$ALCHLBL_C,
                                           function(x)ifelse(x == 'Not at all', 'No', 'Yes'))%>%
    as.factor()
  # select the mild one
  data.HC.standardized <-
    data.HC.standardized[data.HC.standardized$treat=='Mild',]
                         # !colnames(data.HC.standardized)%in%c('treat')]


  # delete the row with missingness
  data.RCT.standardized <- na.omit(data.RCT.standardized) #370
  data.HC.standardized <- na.omit(data.HC.standardized) # 459



  # cov.names <- c( # continuous
  #   'AGEYR', 'HGTCMBL',
  #   'BMIBL', 'EDCCNT',
  #   'DURDIAGYR',
  #   'mmse_base',
  #   # 'WGTKGBL_standarized',
  #   # categorical
  #   'SEXSNM',  'ALCHLBL_C', 'hypertension', 'diabetes', 'treat',
  #   'TRT'
  # )
  # data.RCT.outcomes <- data.RCT[data.RCT$treat=='Mild AD', cov.names]
  # data.HC.outcomes <- data.HC[data.HC$treat=='Mild', cov.names]
  #
  # data.RCT.outcomes <- na.omit(data.RCT.outcomes) #370
  # data.HC.outcomes <- na.omit(data.HC.outcomes) # 459
  #
  # main <- function(arr){
  #  c(mean(arr),
  #    sd(arr))
  # }
  # apply(data.HC.outcomes[, 1:6],
  #       2, main)

  # # preliminary analysis
  # data.combined.standardized <-
  #   rbind(cbind(origin = 0, data.HC.standardized),
  #         cbind(origin = 1, data.RCT.standardized))
  # # focus on mild AD
  # data.combined.standardized.mild <- data.combined.standardized[data.combined.standardized$treat == 'Mild',
  # !colnames(data.combined.standardized)%in%c('treat')]
  # lm(iadrs_13~.,
  #    data = data.combined.standardized.mild)
  # for adas_cog 13:
  ## HC: 0
  ## RCT(0): 0 + 2.2658
  ## RCT(1): 0 + 2.2658 -2.4890 = -0.2232

  # select samples from the RCT, control group
  set.seed(seed)
  idx.1 <- which(data.RCT.standardized$TRT == 1)
  idx.0 <- which(data.RCT.standardized$TRT == 0)

  # for adas_cog13
  # control (129) and treated (120)

  # cat(length(idx.0))
  idx.0.selected <- sample(idx.0, n_c_selected, replace = F)
  selected_idx <- c(idx.1, idx.0.selected)

  # selected name of baseline covariates
  cov.names.selected <- c( # continuous
    'AGEYR', 'HGTCMBL',
    'BMIBL', 'EDCCNT',
    'DURDIAGYR',
    'mmse_base',
    # 'WGTKGBL_standarized',
    # categorical
    'SEXSNM',  'ALCHLBL_C', 'hypertension', 'diabetes'#, 'treat'
  ) #'RACECAT' #'WRKSTAT1LNM', 'HISSTSBL', 'GDSSEVBL' 'TBBL', 'DURDISONSYR_standarized', 'admed'

  # organize the dataset for subsequent analysis
  data.RCT.selected <- DATA_standardizer(data = data.RCT.standardized[selected_idx,])
  data.HC.selected <- DATA_standardizer(data = data.HC.standardized)

  LZAN_list <- list(X = subset(data.RCT.selected,
                               select = cov.names.selected),
                    Y = eval(parse(text=paste0('data.RCT.selected$', y.names))),
                    A = data.RCT.selected$TRT)
  n_t <- sum(LZAN_list$A); n_c <- sum(1-LZAN_list$A)
  n_RCT <- n_t + n_c

  # construct HC
  data.EU <- list(X0 = subset(data.HC.selected,
                              select = cov.names.selected),
                  Y0 = eval(parse(text=paste0('data.HC.selected$', y.names))))
  # data.EU$X0$treat
  n_h <- length(data.EU$Y0)
  # just one historical controls
  data_hc.init <- list(data.EU)
  # no PS matching
  data_hc <- data_hc.init
  # add PS matching
  # add the 1:1 optimal pair propensity matching
  # data_hc <- lapply(data_hc.init, function(hc_i){
  #   XY_h <- cbind(hc_i$X0, Y = hc_i$Y0)
  #   dat_X_ch <- data.frame(rbind(cbind(delta_c = 1,
  #                                      LZAN_list$X, Y = LZAN_list$Y),
  #                                cbind(delta_c = 0, XY_h)))
  #
  #   matched.opt <- matchit(delta_c~.+0-Y, data=dat_X_ch, method='optimal')
  #   X_h_post <- subset(dat_X_ch[matched.opt$match.matrix,],
  #                      select = c(-Y,-delta_c))
  #   Y_h_post <- dat_X_ch[matched.opt$match.matrix,'Y']
  #   n_c1 <- sum(LZAN_list$A); n_c0 <- sum(1-LZAN_list$A)
  #   # randomly selected to balance RCT
  #   delta_n_c <- n_c1 - n_c0
  #   if(delta_n_c<0){return(list(X0 = X_h_post,
  #                               Y0 = Y_h_post))}else{
  #                                 select.idx <- sample(1:(n_c1+ n_c0), delta_n_c)
  #                                 return(list(X0 = X_h_post[select.idx,],
  #                                             Y0 = Y_h_post[select.idx]))}
  #
  # })

  # begin our analysis
  # set.seed(123456)
  # RDA.results <- AIPW_ACW_lasso_SuperLearner(X_c = LZAN_list$X,
  #                Y_c = LZAN_list$Y,
  #                A_c = LZAN_list$A,
  #                prob_A = rep(n_t/(n_t+n_c),
  #                             times = n_t+n_c), # randomly assign
  #                data_hc = data_hc,
  #                lambda = 'min', strata = 1,
  #                lasso = 'individual', refit = FALSE)
  RDA.results <- AIPW_ACW_lasso(X_c = LZAN_list$X,
                                Y_c = LZAN_list$Y,
                                A_c = LZAN_list$A,
                                prob_A = rep(n_t/(n_t+n_c),
                                             times = n_t+n_c), # randomly assign
                                data_hc = data_hc,
                                lambda = 'min', strata = 1,
                                lasso = 'individual', refit = FALSE,
                                method = method, ...)
  # Bayesian estimation
  data.list <- list(RCT = list(A = LZAN_list$A,
                               X = LZAN_list$X,
                               Y = LZAN_list$Y),
                    RWD = list(X = data.EU$X0,
                               Y = data.EU$Y0))

  source('Bayesian//utils_Bayes.R')
  # Bayesian analysis subject to some unforseen error
  PPP.res <- tryCatch(PPP_estimator(data.list, lcnt0 = c( 'AGEYR_standarized',
                                                          'HGTCMBL_standarized',
                                                          'BMIBL_standarized',
                                                          'EDCCNT_standarized',
                                                          'DURDIAGYR_standarized',
                                                          'mmse_base_standarized'
                                                          # 'WGTKGBL_standarized'
  ),
  lord0 = c('SEXSNM',
            'ALCHLBL_C',
            'hypertension',
            'diabetes')),
  error = function(e)list(est = RDA.results$est$AIPW,
                          SE = RDA.results$sd$AIPW/sqrt(RDA.results$n_c)))


  list(
    results = RDA.results,
    n_h = n_h,
    y.names = y.names,
    subset.idx = RDA.results$subset.idx,
    PPP = c(est = PPP.res$est,
            SE = PPP.res$SE))

}


# A main simulation
main <- function(seed,
                 n_c1 = 300, n_c0 = 100, n_h0 = 3e3,
                 tau = 2,
                 confound_Y1 = 0, confound_Y2 = 0,
                 delta_time1 = 0, delta_time2 = 0,
                 gamma_mean1 = NA, gamma_mean2 = NA,
                 lambda = 'min',
                 method = 'separate', strata = 1, lasso = 'group')
{
  set.seed(seed)
  # finite population
  N <- 1e6
  X1 <- rnorm(N); X2 <- runif(N)
  X <- cbind('(Intercept)' = 1, X1, X2)
  Y0 <- X1 + X2 + 1 + rnorm(N, sd = 1)
  # generate concurrent control
  RCT.idx <- sample(1:N, n_c0 + n_c1)
  n_c <- n_c0 + n_c1
  Y_c <- Y0[RCT.idx]; X_c <- X[RCT.idx, ]

  # constant propensity score
  prob_A <- rep(n_c1/n_c, n_c)

  # # change to non-constant propensity score
  # alpha.opt <- nleqslv(x = 1,
  #                      fn=function(x)
  #                      {
  #                        p <- expit(X_c%*%c(x, 1, 1))
  #                        n_c1 - n_c * mean(p)
  #                      })$x
  # prob_A <- c(expit(X_c%*%c(alpha.opt, 1, 1)))

  # # ---------------- old setup
  # # nc <- 200; nh <- 1000
  # delta_c <- rbinom(N, size = 1,
  #                   prob = expit(-8+.3*X1 + .4*X2))
  # n_c <- sum(delta_c)
  # Y_c <- Y0[delta_c==1]; X_c <- X[delta_c==1,]
  # # assign treatment in the concurrent trial
  # prob_A <- rep(.5, n_c) # constant trt assignment
  # # -----------------

  # tau <- 2     # constant trt effects
  A_c <- rbinom(n_c, 1, prob_A)
  Y0_c <- Y_c[A_c==0]; X0_c <- X_c[A_c==0,]
  Y1_c <- Y_c[A_c==1] +
    X1[A_c=1] + X2[A_c=1] + tau; X1_c <- X_c[A_c==1,]
  Y_c[A_c==1] <- Y1_c; Y_c[A_c==0] <- Y0_c
  # generate historical controls
  data_hc_1 <- replicate(1, generate_hc(X = X, Y0 = Y0,
                                        n_h0 = n_h0,
                                        confound_Y = confound_Y1,
                                        delta_time = delta_time1,
                                        gamma_mean = gamma_mean1),
                         simplify = F)

  data_hc_2 <- replicate(1, generate_hc(X = X, Y0 = Y0,
                                        n_h0 = n_h0,
                                        confound_Y = confound_Y2,
                                        delta_time = delta_time2,
                                        gamma_mean = gamma_mean2),
                         simplify = F)
  # concurrent trial: X_c, Y_c and A_c
  if(method == 'joint'){
    data_hc_joint <- list(X0 = do.call(rbind, lapply(c(data_hc_1, data_hc_2),
                                                     function(dat){dat$X0})),
                          Y0 = unlist(lapply(c(data_hc_1, data_hc_2),
                                             function(dat){dat$Y0})))

    # historical controls: X0 and Y0
    data_hc <- list(data_hc_joint)
    # AIPW_ACW_lasso(X_c = X_c, Y_c = Y_c,
    #                A_c = A_c, prob_A = prob_A,
    #                data_hc = list(data_hc_joint),
    #                lambda = lambda)
  }
  if(method == 'separate'){
    data_hc <- c(data_hc_1, data_hc_2)
    # # historical controls: X0 and Y0
    # AIPW_ACW_lasso(X_c = X_c, Y_c = Y_c,
    #                A_c = A_c, prob_A = prob_A,
    #                data_hc = c(data_hc_1, data_hc_2),
    #                lambda = lambda)
  }
  if(method == 'nearest')
  {
    library(bayesmeta)
    kl.val <- lapply(c(data_hc_1, data_hc_2), function(dat){
      col.mean1 <- colMeans(dat$X0[,-1])
      col.mean2 <- colMeans(X_c[,-1])
      cov1 <- cov(dat$X0[,-1]);cov2 <- cov(X_c[,-1])
      kldiv(col.mean1, col.mean2,
            cov1, cov2)
      })
    select.idx <- which.min(kl.val%>%unlist())

    # with selected historical controls: X0 and Y0
    data_hc <- c(data_hc_1, data_hc_2)[select.idx]

    # AIPW_ACW_lasso(X_c = X_c, Y_c = Y_c,
    #                A_c = A_c, prob_A = prob_A,
    #                data_hc = c(data_hc_1, data_hc_2)[select.idx],
    #                lambda = lambda)
  }
  if(method == 'optimal')
  {
    data_hc <- lapply(c(data_hc_1, data_hc_2), function(data_hc){
      XY_h <- cbind(data_hc$X0, Y = data_hc$Y0)
      dat_X_ch <- data.frame(rbind(cbind(delta_c = 1, X_c, Y = Y_c),
                                   cbind(delta_c = 0, XY_h)))

      matched.opt <- matchit(delta_c~.+0-Y, data=dat_X_ch, method='optimal')
      X_h_post <- subset(dat_X_ch[matched.opt$match.matrix,],
                         select = c(-Y,-delta_c))
      X_h_post <- as.matrix(X_h_post)
      Y_h_post <- dat_X_ch[matched.opt$match.matrix,'Y']
      # randomly selected to balance RCT
      delta_n_c <- n_c1 - n_c0
      select.idx <- sample(1:(n_c1+n_c0), delta_n_c)

      list(X0 = X_h_post[select.idx,],
           Y0 = Y_h_post[select.idx])
    })

    # # with selected historical controls: X0 and Y0
    # AIPW_ACW_lasso(X_c = X_c, Y_c = Y_c,
    #                A_c = A_c, prob_A = prob_A,
    #                data_hc = data_hc,
    #                lambda = lambda)
  }
  AIPW_ACW_lasso(X_c = X_c, Y_c = Y_c,
                 A_c = A_c, prob_A = prob_A,
                 data_hc = data_hc,
                 lambda = lambda, strata = strata, lasso = lasso)
}

main_real <- function(seed,
                      N_TD = 200, N_CD = 100,
                      N_K = 3e3,
                      omega = .3, # confounder
                      beta_T = .3, # treatment effect
                      delta1 = .3, # concurrent bias
                      gamma_mu = .3 # measurement error
)
{
  M_normal <- round(calc_theory(Dist = 'Gaussian',
                                params = c(0,1)), 8)
  M <- replicate(12, M_normal)
  ncont <- ncol(M)

  # create a correlation matrix
  set.seed(seed)
  Rey <- diag(1, nrow = (ncont))
  for (i in 1:nrow(Rey)) {
    for (j in 1:ncol(Rey)) {
      if (i > j) Rey[i, j] <- (1/2)**(i-j)
      Rey[j, i] <- Rey[i, j]
    }
  }
  min(eigen(Rey, symmetric = TRUE)$values)

  # simulate variables with the error loop for the finite population
  N <- 1e5
  X.info <- rcorrvar(n = N, k_cont = ncont, method = "Polynomial",
                     means =  M[1, ], vars =  (M[2, ])^2,
                     skews = M[3, ], skurts = M[4, ],
                     fifths = M[5, ], sixths = M[6, ],
                     rho = Rey, seed = seed, errorloop = TRUE)
  X <- X.info$continuous_variables
  X <- as.matrix(X)
  # add intercept
  X <- cbind(1, X)
  dim(X)
  # slope for each covariate is selected empirically
  beta <- c(1, rep(.5, 4), rep(-.5, 4), rep(1, 4))
  # inclusion and exclusion criteria
  ## random sample for now
  ## RCT
  # N_TD <- 200; N_CD <- 100
  # omega <- .3; beta_T <- .3
  X_D.idx <- sample(1:N, size = N_TD + N_CD)
  X_D <- X[X_D.idx, ]
  # assign treatment (N_TD v.s. N_CD)
  alpha.opt <- nleqslv(x = 1,
                       fn=function(x)
                       {
                         p <- expit(X_D%*%c(x, rep(1,4), rep(.5,4), rep(-.5, 4)))
                         N_TD - (N_TD + N_CD) * mean(p)
                       }, method = 'Newton')$x
  prob_A <- expit(X_D%*%c(alpha.opt, rep(1,4), rep(.5,4), rep(-.5, 4)))
  mean(prob_A)
  W_D <- rbinom(N_TD + N_CD, size = 1, prob = prob_A)
  # split the covariate
  X_TD.idx <- which(W_D == 1)
  X_TD <- X_D[X_TD.idx, ]
  X_CD.idx <- which(W_D == 0)
  X_CD <- X_D[X_CD.idx, ]
  ## generate the outcomes
  Y_TD <- as.matrix(X_TD)%*%beta + beta_T +
    omega * rnorm(nrow(X_TD)) + rnorm(nrow(X_TD))
  Y_CD <- as.matrix(X_CD)%*%beta +
    omega * rnorm(nrow(X_CD)) + rnorm(nrow(X_CD))
  Y_D <- numeric(N_TD + N_CD)
  Y_D[X_TD.idx] <- Y_TD;Y_D[X_CD.idx] <- Y_CD;


  # external control
  # N_K <- 3e3
  # delta1 <- 0.3 # concurrency bias
  # concurrecy lack indicator
  T_K <- sample(0:2, N_K, replace = T,
                prob = c(1/3, 1/3, 1/3))
  gamma_K <- if(is.null(gamma_mu)|is.na(gamma_mu)){
    0}else{rnorm(N_K, mean = gamma_mu)} # measurement error
  X_K.idx <- sample(1:N, size = N_K)
  X_K <- X[X_K.idx, ]
  X_K <- X_K + gamma_K # add measurement error
  X_lin <- as.matrix(X_K)%*%beta + omega * rnorm(N_K, mean= 0.3)
  # outcome validity deviations:
  # (X_lin) to exp(X_lin)
  Y_K <- c((X_lin) + delta1 * T_K + rnorm(N_K))
  # randomly split the historical controls into K folds
  library(caret)
  library(dplyr)
  folds.idx <- createFolds(1:N_K, k = K, list = F)
  data_hc <- lapply(1:K, function(k)
  {
    idx <- which(folds.idx==k)
    list(X0 = as.matrix(X_K[idx,]), Y0 = Y_K[idx])
  })

  # AIPW(X_c = X_D, Y_c = Y_D, A_c = W_D,
  #      prob_A = c(prob_A))

  AIPW_ACW_lasso(X_c = X_D, Y_c = Y_D, A_c = W_D,
                 prob_A = c(prob_A),
                 data_hc = data_hc)
}


# evaluate real data simulation
sim.eval.real <- function(N_TD, N_CD,
                          N_H, beta_T, delta1, omega, gamma_mu,
                          error_Y = 0,
                          boot = FALSE,
                          suffix = 'lm'
                          )
{
  if(boot){
    file.name <- paste0('sim_real_N_CD_', N_CD,
                        '_beta_T_', beta_T,
                        '_delta1_', delta1,
                        '_omega_', omega,
                        '_gamma_mu_', gamma_mu,
                        # '_error_Y_', error_Y,
                        '_rep1000_strata1_individual_norefit_HTE_boot.RData')
  }else{

    if(suffix == 'lm'){
      if(omega==0){prop <- NULL}
      file.name <- paste0('sim_real_N_CD_', N_CD,
                          '_beta_T_', beta_T,
                          '_delta1_', delta1,
                          '_omega_', omega, prop,
                          '_gamma_mu_', gamma_mu,
                          '_error_Y_', error_Y,
                          # '_individual_norefit_HTE_newEIF_nonlinear_logit.RData' # incorrect mu
                          # '_individual_norefit_HTE_newEIF_notlogit.RData' # incorrect q
                          '_rep1000_individual_norefit_HTE_newEIF_logit.RData' # all correct
                          )
    }else{
      if(omega==0){prop <- NULL}
      file.name <- paste0('sim_real_N_CD_', N_CD,
                          '_beta_T_', beta_T,
                          '_delta1_', delta1,
                          '_omega_', omega, prop,
                          '_gamma_mu_', gamma_mu,
                          '_error_Y_', error_Y,
                          # '_individual_norefit_HTE_newEIF_nonlinear_logit_gbm.RData' # incorrect mu
                          # '_individual_norefit_HTE_newEIF_notlogit_gbm.RData' # incorrect q
                          '_rep1000_individual_norefit_HTE_newEIF_logit_gbm.RData' # all correct
                          )
    }

  }


  # sim_real_N_CD_100_beta_T_0_delta1_0_omega_0_gamma_mu_sd_error_Y_0_rep1000_strata1_individual_norefit_HTE_logit
  # file.name <- paste0('sim_real_N_TD_', N_TD,
  #                     '_N_CD_', N_CD,
  #                     '_N_H_', N_H,
  #                     '_beta_T_', beta_T,
  #                     '_delta1_', delta1,
  #                     '_omega_',omega,
  #                     '_gamma_mu_', gamma_mu,
  #                     '_rep1000_separate.RData')
  load(file = file.name)
  # load(file.path('RData_real',file.name))
  # load(file.path('RData_sim',file.name))
  # load(file.path('RData_sim_new',file.name))
  # combined subset
  # lapply(res.real_dat, function(x)tryCatch(x$subset.idx,
  #                                          error = function(e)NA))
  prob.comb <- mean(lapply(res.real_dat, function(x)tryCatch(length(x$subset.idx)/(200-N_CD),
                                           error = function(e)NA))%>%unlist(), na.rm = TRUE)

  # prob.comb <- 1 - mean(lapply(res.real_dat, function(x)tryCatch(identical(x$subset.idx, numeric(0)),
  #                                          error = function(e)NA))%>%unlist(), na.rm = TRUE)
  # point estimator
  tau.AIPW <- lapply(res.real_dat, function(x)tryCatch(x$est$AIPW,
                                                       error = function(e)NA))%>%unlist()
  tau.AIPW.sd <- lapply(res.real_dat, function(x)tryCatch(x$sd$AIPW,
                                                          error = function(e)NA))%>%unlist()
  tau.ACW <- lapply(res.real_dat, function(x)tryCatch(x$est$ACW,
                                                      error = function(e)NA))
  tau.ACW.mat <- do.call(rbind, tau.ACW)
  tau.ACW.sd <- lapply(res.real_dat, function(x)tryCatch(x$sd$ACW,
                                                         error = function(e)NA))
  tau.ACW.sd.mat <- do.call(rbind, tau.ACW.sd)
  tau.ACW.lasso <- lapply(res.real_dat, function(x)tryCatch(x$est$ACW.lasso,
                                                            error = function(e)NA))%>%unlist()
  tau.ACW.lasso.sd <- lapply(res.real_dat, function(x)tryCatch(x$sd$ACW.lasso,
                                           error = function(e)NA))%>%unlist()

  # add estimate for tau.final


  tau.ACW.final <- lapply(res.real_dat, function(x)tryCatch(as.vector(unlist(x$est)['ACW.final']),
                                                            error = function(e)NA))%>%unlist()
  tau.ACW.final.sd <- lapply(res.real_dat, function(x)tryCatch(as.vector(unlist(x$sd)['ACW.final']),
                                                               error =function(e)NA))%>%unlist()


  n_c <- lapply(res.real_dat, function(x)tryCatch(x$n_c,
                                                  error = function(e)NA))%>%unlist()
  # varince estimation
  var.est.AIPW <- lapply(res.real_dat, function(x)tryCatch(x$sd$AIPW**2,
                                                           error = function(e)NA))%>%unlist()
  var.est.ACW <- do.call(rbind, lapply(res.real_dat, function(x)tryCatch(x$sd$ACW**2,
                                                                         error = function(e)NA)))
  var.est.ACW.lasso <- lapply(res.real_dat, function(x)tryCatch(x$sd$ACW.lasso**2,
                                                                error = function(e)NA))%>%unlist()
  var.est.final <- tau.ACW.final.sd**2
  # construct the CI
  ## AIPW
  tau.AIPW.u <- tau.AIPW + qnorm(1-0.025) * tau.AIPW.sd/sqrt(n_c)
  tau.AIPW.l <- tau.AIPW - qnorm(1-0.025) * tau.AIPW.sd/sqrt(n_c)
  # ACW
  tau.ACW.mat.u <- tau.ACW.mat + qnorm(1-0.025) * tau.ACW.sd.mat/sqrt(n_c)
  tau.ACW.mat.l <- tau.ACW.mat - qnorm(1-0.025) * tau.ACW.sd.mat/sqrt(n_c)

  # ACW.lasso unconditional
  tau.ACW.lasso.u <- tau.ACW.lasso + qnorm(1-0.025) * tau.ACW.lasso.sd/sqrt(n_c)
  tau.ACW.lasso.l <- tau.ACW.lasso - qnorm(1-0.025) * tau.ACW.lasso.sd/sqrt(n_c)

  # ACW.lasso conditional
  tau.ACW.lasso.u.cont <- lapply(res.real_dat, function(x)tryCatch(x$CI$upper,
                                                                   error = function(e)NA))%>%unlist()
  tau.ACW.lasso.l.cont <- lapply(res.real_dat, function(x)tryCatch(x$CI$lower,
                                                                   error = function(e)NA))%>%unlist()
  # ACW.final unconditional
  tau.ACW.final.u <- tau.ACW.final + qnorm(1-0.025) * tau.ACW.final.sd/sqrt(n_c)
  tau.ACW.final.l <- tau.ACW.final - qnorm(1-0.025) * tau.ACW.final.sd/sqrt(n_c)

  # apply(tau.ACW.mat.u - tau.ACW.mat.l, 2, mean)
  # mean(tau.ACW.lasso.u - tau.ACW.lasso.l)
  # mean(tau.ACW.lasso.u.cont - tau.ACW.lasso.l.cont)

  return(list(n_c = n_c, est = list(AIPW = tau.AIPW,
                                    ACW = tau.ACW.mat,
                                    ACW.lasso = tau.ACW.lasso,
                                    ACW.final = tau.ACW.final),
              var.est = list(AIPW = var.est.AIPW,
                             ACW = var.est.ACW,
                             ACW.lasso = var.est.ACW.lasso,
                             ACW.final = var.est.final),
              CI = list(AIPW = cbind(lower = tau.AIPW.l,
                                     upper = tau.AIPW.u),
                        ACW = list(lower = tau.ACW.mat.l,
                                   upper = tau.ACW.mat.u),
                        ACW.lasso = cbind(lower = tau.ACW.lasso.l,
                                          upper = tau.ACW.lasso.u),
                        ACW.final = cbind(lower = tau.ACW.final.l,
                                          upper = tau.ACW.final.u)
                        # ACW.lasso.cont = cbind(lower = tau.ACW.lasso.l.cont,
                        #                   upper = tau.ACW.lasso.u.cont)
              ),
              prob = prob.comb))
}

# plot the dataframe
plot.df <- function(dat, y.label)
{
  # "#F8766D" "#A3A500" "#00BF7D" "#00B0F6" "#E76BF3"
  library(RColorBrewer)
  library("ggsci")

  K <- length(grep('ACW', row.names(dat))) - 1 # subtract ACW.lasso

  df <- data.frame(cbind(N_CD = N_CD.vec, t(dat)))
  df <- reshape2::melt(df, id.vars = 'N_CD',
                       variable.name = 'name')
  # if(K==1)
  # {
  #
  #   df$name <- factor(df$name, levels = c('AIPW', paste0('ACW',1:K), 'ACW.lasso',
  #                                         'PPP'),
  #                     labels =  c('AIPW.RCT', 'ACW', 'ACW.lasso',
  #                                 'PPP'))
  # }else{
  #   df$name <- factor(df$name, levels = c('AIPW', paste0('ACW',1:K), 'ACW.lasso',
  #                                         'PPP'))
  # }
  # df <- subset(df, subset = name%in%c('AIPW', 'PPP', 'ACW.lasso'))

  library(ggplot2)
  ggplot(df, aes(x = N_CD, y = value,
                 color = name))+
    geom_point(size = 3)+
    geom_line(size = 1)+theme_light()+
    guides(color=guide_legend(title=""),
           alpha='none')+
    scale_color_discrete(label = c(TeX('$\\hat{\\tau}_{aipw}$'),
                                   TeX('$\\hat{\\tau}_{acw}$'),
                                   TeX('$\\hat{\\tau}_{acw}^{lasso}$'),
                                   TeX('$\\hat{\\tau}_{acw,gbm}^{lasso}$'),
                                   TeX('$\\hat{\\tau}_{ppp}$')))+
    theme(legend.text = element_text(size = 15),
          legend.key.size = unit(1, 'cm'))+
    xlab(TeX('$n_{c0}$')) + ylab(y.label)#+
    # scale_color_manual(values = c("#F8766D",
    #                               "#A3A500",
    #                               "#00BF7D",
    #                               "#00B0F6"))




}

# handy functions
expit <- function(x){exp(x)/{1+exp(x)}}


# generate data based on real-data moments and correlation (or sampling score)
getdata <- function(n_c = 100, # RCT control
                    n_t = 200, # RCT treat
                    n_k = 3e3, # RWE data
                    corY_Xn = 0, bias_h = 0.5,# unmeasured confounder (correlation and biasedness)
                    eff_size = 0, # treatment effect
                    time_bias = 0, # time inconcurrency,
                    error_X = 0, # measurement error
                    error_Y = 0, # outcome validity
                    M = 1, linear = TRUE)
  {
  if(linear){
    file.name <- paste0("data_nc.",n_c,"_nt.",n_t,"_nk.",n_k,
                        "_effsize_", eff_size,
                        "_timebias_",time_bias,
                        "_bXn_", corY_Xn, "_bias_h_", bias_h,
                        "_error_X_", error_X,
                        '_corr',
                        # "_error_Y_", error_Y,
                        "_rep_", M,
                        ".RData")
    folder.name <- paste0("data_nc.",n_c,"_nt.",n_t,"_nk.",n_k,
                          "_effsize_", eff_size,
                          "_timebias_",time_bias,
                          "_bXn_", corY_Xn, "_bias_h_", bias_h,
                          "_error_X_", error_X, '_corr'
                          # "_error_Y_", error_Y
    )
  }else{
    file.name <- paste0("data_nc.",n_c,"_nt.",n_t,"_nk.",n_k,
                        "_effsize_", eff_size,
                        "_timebias_",time_bias,
                        "_bXn_", corY_Xn, "_bias_h_", bias_h,
                        "_error_X_", error_X,
                        "_error_Y_", error_Y,
                        '_HTE',
                        "_rep_", M,
                        ".RData")
    folder.name <- paste0("data_nc.",n_c,"_nt.",n_t,"_nk.",n_k,
                          "_effsize_", eff_size,
                          "_timebias_",time_bias,
                          "_bXn_", corY_Xn, "_bias_h_", bias_h,
                          "_error_X_", error_X,
                          "_error_Y_", error_Y,
                          '_HTE'
    )
  }


  # create the folder if not exists
  ifelse(!dir.exists(file.path('RWD_sim', folder.name)),
         dir.create(file.path('RWD_sim', folder.name)),
         FALSE)
  # if the data has been generated, load it and return
  if(file.exists(file.path('RWD_sim',
                           folder.name,
                           file.name))){
    load(file.path('RWD_sim',
                   folder.name,
                   file.name))}else{
    corXjXn <- 0 # U was independent with all covariates (or correlated with 0.3)
    corY_Ybase <- 0.5 # baseline outcome

    # read in the data set
    d <- haven::read_sas('..//real_data//hybrid_monarch_nsai_spotlight.sas7bdat')


    # Xn for RCT N(0,1) (unobserved confounder in RCT)
    d$Xn <- rnorm(nrow(d),0,1)
    # out_RCTcome
    Y <- 'PFS_standardized'
    # Categorical Covariates, note: "group" is removed as a covariate as it is the treatment assignment/data source indicator
    lord0 <- c(NULL,"race","ECOGBL","NATDISI","TGRADE","num_sites_mets","IPDDISTG","PGRSTAT")
    # continuous Covariates, note: "PFS1" is removed as a covariate because it represents the out_RCTcome
    lcnt0 <- c(NULL,"AGE_standardized","BMI_standardized","Xn","PFS_standardized")
    ###Process categorical covariates
    d$race <- ifelse(d$race=='white','WHITE','NONWHITE')
    d$ECOGBL <- ifelse(d$ECOGBL=="",NA,ifelse(d$ECOGBL==">1","1",d$ECOGBL))
    d$TGRADE <- ifelse(d$TGRADE=='Unknown',NA,d$TGRADE)
    d$num_sites_mets <- ifelse(d$num_sites_mets=="",NA,d$num_sites_mets)
    d$NATDISI <- ifelse(d$NATDISI=='VISCERAL','VISCERAL','BONE ONLY')
    d$IPDDISTG <- ifelse(d$IPDDISTG%in%c("missing","Not documented","0"), NA, d$IPDDISTG)
    d$PGRSTAT <- ifelse(d$PGRSTAT%in%c("","other"),"negative",d$PGRSTAT)
    ###Create standardized version of continuous covariates, standardized to combined population
    d$AGE_standardized=(d$AGE-mean(d$AGE,na.rm=TRUE))/sd(d$AGE,na.rm=TRUE)
    d$BMI_standardized=(d$BMI-mean(d$BMI,na.rm=TRUE))/sd(d$BMI,na.rm=TRUE)
    d$PFS_standardized=(d$PFS1-mean(d$PFS1,na.rm=TRUE))/sd(d$PFS1,na.rm=TRUE)

    ### Extract joint distribution of X1,...Xn from RCT
    ch <- 'Monarch3'
    lord <- lord0
    lcnt <- lcnt0[!(lcnt0%in%c("TimeEff","SettingEff"))]

    tm <- d %>% filter(cohort==ch)
    tm <- tm %>% mutate_if(is.character, as.factor)
    X <- tm[,c(lord,lcnt,Y)]
    for(col in lord) X[[col]] <- as.numeric(X[[col]])
    # obtain the correlation matrix for X
    rho <- cor(X[,c(lord,lcnt)], method="pearson",use="pairwise.complete.obs")

    # frml <- paste('~',paste(c(lord,lcnt[lcnt!='PFS1']),collapse='+'))
    # previous_na_action=options('na.action')
    # options(na.action='na.pass')
    # Xtgt <- model.matrix(as.formula(frml),tm)[,-1]
    # options(na.action=previous_na_action$na.action)

    # Adjust the correlation matrix to reflect the relationship between Xn and X1,...,X_n-1
    # The target correlation matrix must be ordered ordinal, continuous, Poisson, Negative Binomial
    Xn.idx <- which(names(X)=="Xn")
    Xj.idx <- which(names(X)=="BMI_standardized")
    # Assume that the Xn is at most correlated with one of the covariates
    rho[, Xn.idx] <- rho[Xn.idx, ] <- 0; rho[Xn.idx, Xn.idx] <- 1
    rho[Xn.idx, Xj.idx] <- rho[Xj.idx, Xn.idx] <- corXjXn
    # Assume the baseline is not correlated with any covariate
    rho['PFS_standardized', ] <- rho[,'PFS_standardized'] <- 0
    rho['PFS_standardized', 'PFS_standardized'] <- 1

    ##### Simulate X1,...,Xn for RCT using SimCorrVar function (Plasmode)
    method = "Polynomial"
    # ordinal marginal distribution and support
    marginal <- lapply(lord, function(lord.name)
    {
      tab <- as.numeric(table(X[,lord.name]))
      tab <- cumsum(tab)/sum(tab)
      tab[-length(tab)]
    })
    support <- lapply(lord, function(lord.name)
    {
      sort(unique(X[[lord.name]]))
    })

    # continuous moments
    cmoms_RCT <- do.call(rbind, lapply(lcnt, function(lcnt.name)
    {
      calc_moments(na.omit(X[[lcnt.name]]))
    }))

    # # correction terms
    # ## first attempt optimization
    # lSix0 <- lapply(1:nrow(cmoms_RCT), function(i){
    #   fc <- find_constants(
    #     method = method,
    #     skews = cmoms_RCT[i,3], skurts = cmoms_RCT[i,4],
    #     fifths = cmoms_RCT[i,5], sixths = cmoms_RCT[i,6],
    #     Six = seq(-10, 100, length.out = 10), seed = M)
    #   ifelse(length(fc$SixCorr1)==0, 0, fc$SixCorr1)
    # })
    # ## revised optimization
    # lSix <- lapply(1:nrow(cmoms_RCT), function(i){
    #   fc <- find_constants(
    #     method = method,
    #     skews = cmoms_RCT[i,3], skurts = cmoms_RCT[i,4],
    #     fifths = cmoms_RCT[i,5], sixths = cmoms_RCT[i,6],
    #     Six = seq(lSix0[[i]]-1, lSix0[[i]]+1,by = .1))
    #   ifelse(length(fc$SixCorr1)==0, 0, fc$SixCorr1)
    # })

    sim <- rcorrvar(n = n_t+n_c, k_cont = length(lcnt),k_cat = length(lord),  # change n
                    method = method,
                    means =cmoms_RCT[,1],
                    vars = cmoms_RCT[,2]^2,
                    skews = cmoms_RCT[,3],
                    skurts = cmoms_RCT[,4],
                    fifths = cmoms_RCT[,5],
                    sixths = cmoms_RCT[,6],
                    Six = rep(list(seq(0, 100, length.out = 100)),
                              length(lcnt)+length(lord)),
                    marginal=marginal, support=support,
                    Sigma = NULL,
                    rho = rho, cstart = NULL, seed = M, errorloop = FALSE,
                    epsilon = 0.001, maxit = 1000, extra_correct = TRUE)
    out_RCT <- cbind(sim$ordinal_variables,
                     sim$continuous_variables)
    colnames(out_RCT) <- c(lord,lcnt)
    # Check marginal distribution of categorical covaraites and compare to observed distribution
    # for (v in lord0) {print(v);print(prop.table(table(tm[[v]])))}
    # for (v in lord0) {print(v);print(prop.table(table(out_RCT[[v]])))}

    # Check marginal distribution of continuous covariates and compare to observed distribution
    # for (v in lcnt0) {print(v);print(c(summary(tm[[v]]), sd(tm[[v]],na.rm=TRUE)))}
    # for (v in lcnt0) {print(v);print(c(summary(out_RCT[[v]]), sd(out_RCT[[v]],na.rm=TRUE)))}

    #Check correlation structure of all covariates
    # library(corrplot)
    # corrplot(cor(out_RCT))
    # corrplot(rho)
    # rename it as baseline
    colnames(out_RCT)[11] <- "PFS_baseline"
    # re-map the ordinal covariates back
    for(col in lord) {
      out_RCT[[col]]=factor(out_RCT[[col]],labels=levels(tm[[col]]))
    }

    ## Randomly assign treatment for RCT on a 2:1 ratio (not involve the covariate)
    out_RCT$group <- rbinom(n_c + n_t, size = 1,
                            prob = n_t/(n_t+n_c))

    ## Generate out_RCTcome as a linear combination of covariates
    #summary(lm(PFS_standardized~race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+AGE_standardized+BMI_standardized,data=tm,na.action=na.exclude))

    # drop group

    if(linear){
      ## fit the model empirically
      fit.RCT <- lm(PFS_standardized~race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+
                      AGE_standardized+BMI_standardized,data=tm,na.action=na.exclude)
      betas <- fit.RCT$coefficients
      cond_var <- summary(fit.RCT)$sigma
      marginal_sd <- summary(lm(PFS_standardized~group, data=tm, na.action=na.exclude))$sigma

      ## construct the outcome
      b_Trt <- eff_size; b_base <- corY_Ybase; b_Xn <- corY_Xn
      X.mat <- model.matrix(~race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+
                              AGE_standardized+BMI_standardized, data=out_RCT)

      # Conditional outcome generation
      out_RCT$PFS_standardized <- rnorm(nrow(out_RCT),
                                        X.mat%*%betas+
                                          out_RCT$group*b_Trt+
                                          out_RCT$PFS_baseline*b_base+
                                          out_RCT$Xn*b_Xn, sd=marginal_sd)
      coef.b <- lm(PFS_standardized~race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+
                     AGE_standardized+BMI_standardized+PFS_baseline + group,
                   data=out_RCT,na.action=na.exclude)$coefficients
      coef.b.trt <- coef.b['group']
      mean.diff <- (mean(out_RCT$PFS_standardized[out_RCT$group==1])-
                      mean(out_RCT$PFS_standardized[out_RCT$group==0]))
    }else{

      # fit the model empirically by nonlinear model (e.g., random-forest or spline model)
      # library(randomForest)
      library(gam)
      fit_gam <- gam(PFS_standardized~race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+
                      AGE_standardized+BMI_standardized,#+
                      # I(AGE_standardized^2) + I(BMI_standardized^2) +
                      #  I(AGE_standardized^3) + I(BMI_standardized^3)+
                      #  exp(AGE_standardized) + exp(BMI_standardized) +
                      #  sin(AGE_standardized) + sin(BMI_standardized) +
                      #  cos(AGE_standardized) + cos(BMI_standardized)
                     data=tm,na.action=na.exclude)
      marginal_sd <- summary(lm(PFS_standardized~group, data=tm, na.action=na.exclude))$sigma

      ## construct the outcome
      b_Trt <- eff_size; b_base <- corY_Ybase; b_Xn <- corY_Xn
      # conditional outcome generation
      ## Y(0)
      mu0.X <- predict(fit_gam,
              newdata = out_RCT[,c(lord, 'AGE_standardized', 'BMI_standardized')])
      ## HTE
      tau.X <- eff_size +
        c(as.matrix(out_RCT[, c('AGE_standardized', 'BMI_standardized')])%*%c(1,1)) - # add the HTE
        c(cmoms_RCT[which(lcnt%in%c('AGE_standardized', 'BMI_standardized')),1]%*%c(1,1)) # adjust the mean

      ## generate the outcome for RWD
      out_RCT$PFS_standardized <- rnorm(nrow(out_RCT),
                                        mu0.X + out_RCT$group*tau.X +
                                          out_RCT$PFS_baseline*b_base+
                                          out_RCT$Xn*b_Xn,
                                        sd=marginal_sd)
    }



    # summmary.fit <- data.frame(coef.b.trt=unname(coef.b.trt),
    #                            mean.diff = mean.diff)
    # write.csv(summmary.fit,paste("~/real_data/data/Null/summary",M, ".csv",sep = ""),row.names = FALSE)

    #summary(lm(PFS_standardized~PFS_baseline, data=out_RCT, na.action=na.exclude))
    #cor(out_RCT$PFS_standardized,out_RCT$PFS_baseline)
    #cor(out_RCT$PFS_baseline,out_RCT$PFS_standardized)
    #cor(out_RCT$Xn,out_RCT$PFS_standardized)

    ####################################################################################################################
    ### using calibration method to generate X1,...Xn from RWD
    # fit P(delta_c = 1 | X)

    ### Extract joint distribution of X1,...Xn from RWD
    ch <- 'Spotlight'
    lord <- lord0
    lcnt <- lcnt0[!(lcnt0%in%c("TimeEff","SettingEff"))]

    # Xn for RWD N(0.3,1)
    # change Xn for RWD N(0.5,1). 04/18/22
    d$Xn <- rnorm(nrow(d), bias_h, 1)
    tm <- d %>% filter(cohort==ch)
    tm <- tm %>% mutate_if(is.character, as.factor)
    X <- tm[,c(lord,lcnt,Y)]

    #X <- tm %>% dplyr::select( all_of(lord),all_of(lcnt),all_of(Y) )

    for(col in lord) X[[col]] <- as.numeric(X[[col]])
    rho <- cor(X[,c(lord,lcnt)], method="pearson",use="pairwise.complete.obs")

    #Adjust the correlation matrix to reflect the relationship between Xn and X1,...,X_n-1
    #The target correlation matrix must be ordered ordinal, continuous, Poisson, Negative Binomial
    Xn.idx <- which(names(X)=="Xn")
    Xj.idx <- which(names(X)=="BMI_standardized")
    # Assume that the Xn is at most correlated with one of the covariates
    rho[, Xn.idx] <- rho[Xn.idx, ] <- 0; rho[Xn.idx, Xn.idx] <- 1
    rho[Xn.idx, Xj.idx] <- rho[Xj.idx, Xn.idx] <- corXjXn
    # Assume the baseline is not correlated with any covariate
    rho['PFS_standardized', ] <- rho[,'PFS_standardized'] <- 0
    rho['PFS_standardized', 'PFS_standardized'] <- 1
    rho

    ##### Simulate X1,...,Xn for RWD using SimCorrVar function (Plasmode)
    # ordinal marginal distribution and support
    method = "Polynomial"
    # ordinal marginal distribution and support
    marginal2 <- lapply(lord, function(lord.name)
    {
      tab <- as.numeric(table(X[,lord.name]))
      tab <- cumsum(tab)/sum(tab)
      tab[-length(tab)]
    })

    # continuous moments
    cmoms_RWD <- do.call(rbind, lapply(lcnt, function(lcnt.name)
    {
      calc_moments(na.omit(X[[lcnt.name]]))
    }))

    # correction terms
    # ## first attempt optimization
    # lSix0 <- lapply(1:nrow(cmoms_RWD), function(i){
    #   fc <- find_constants(
    #     method = method,
    #     skews = cmoms_RWD[i,3], skurts = cmoms_RWD[i,4],
    #     fifths = cmoms_RWD[i,5], sixths = cmoms_RWD[i,6],
    #     Six = seq(-10, 100, length.out = 10), seed = M)
    #   ifelse(length(fc$SixCorr1)==0, 0, fc$SixCorr1)
    # })
    # ## revised optimization
    # lSix <- lapply(1:nrow(cmoms_RWD), function(i){
    #   fc <- find_constants(
    #     method = method,
    #     skews = cmoms_RWD[i,3], skurts = cmoms_RWD[i,4],
    #     fifths = cmoms_RWD[i,5], sixths = cmoms_RWD[i,6],
    #     Six = seq(lSix0[[i]]-1, lSix0[[i]]+1,by = .1))
    #   ifelse(length(fc$SixCorr1)==0, 0, fc$SixCorr1)
    # })

    # need to change the mean of RWE?

    sim <- rcorrvar(n = n_k, k_cont = length(lcnt),k_cat = length(lord),
                    method = method,
                    means =cmoms_RWD[,1],
                    vars = cmoms_RWD[,2]^2,
                    skews = cmoms_RWD[,3],
                    skurts = cmoms_RWD[,4],
                    fifths = cmoms_RWD[,5],
                    sixths = cmoms_RWD[,6],
                    Six = rep(list(seq(0, 100, length.out = 100)),
                              length(lcnt)+length(lord)),
                    marginal=marginal2, support=support,
                    Sigma = NULL,
                    rho = rho, cstart = NULL, seed = M, errorloop = FALSE,
                    epsilon = 0.001, maxit = 1000, extra_correct = TRUE)


    out_RWD <- cbind(sim$ordinal_variables,
                     sim$continuous_variables)
    colnames(out_RWD) <- c(lord,lcnt); colnames(out_RWD)[11] <- "PFS_baseline"


    # Check marginal distribution of categorical covaraites and compare to observed distribution
    # for (v in lord0) {print(v);print(prop.table(table(tm[[v]])))}
    # for (v in lord0) {print(v);print(prop.table(table(out_RWD[[v]])))}

    #Check marginal distribution of continuous covariates and compare to observed distribution
    # for (v in lcnt0) {print(v);print(c(summary(tm[[v]]), sd(tm[[v]],na.rm=TRUE)))}
    # for (v in lcnt0) {print(v);print(c(summary(out_RWD[[v]]), sd(out_RWD[[v]],na.rm=TRUE)))}
    # mean(out_RWD$PFS_baseline)
    #
    #Check correlation structure of all covariates
    # dat=out_RWD %>% mutate_if(is.factor, as.numeric)
    # print(cor(dat))
    # print(rho)

    for(col in lord) {
      out_RWD[[col]] <- factor(out_RWD[[col]],labels=levels(tm[[col]]))
    }

    X <- model.matrix(~race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+
                        AGE_standardized+BMI_standardized,data=out_RWD)

    # Generate variable for time inconcurrency
    if(time_bias=='unif')
    {
      out_RWD$time <- runif(dim(out_RWD)[1])
      ## Generate out_RWDcome as a linear combination of covariates
      time_bias.vec <- 0.3*out_RWD$time

    }else{
      out_RWD$time <- sample(c(0,1,2), size=dim(out_RWD)[1], replace=TRUE)
      ## Generate out_RWDcome as a linear combination of covariates
      time_bias.vec <-  time_bias*out_RWD$time
    }

    # linear combination or not
    if(linear){

      out_RWD$PFS_standardized <- rnorm(nrow(out_RWD), X%*%betas +
                                          b_base*out_RWD$PFS_baseline +
                                          b_Xn*out_RWD$Xn +
                                          time_bias.vec, sd=marginal_sd)
    }else{

      mu0.X.RWD <- predict(fit_gam,
                           newdata = out_RWD[,c(lord, 'AGE_standardized', 'BMI_standardized')])
      out_RWD$PFS_standardized <- rnorm(nrow(out_RWD), mu0.X.RWD +
                                          b_base*out_RWD$PFS_baseline +
                                          b_Xn*out_RWD$Xn +
                                          time_bias.vec, sd=marginal_sd)
    }






    #cor(out_RWD$PFS_baseline,out_RWD$PFS_standardized)

    #Add in measurement error
    if(error_X=='sd'){
      out_RWD$PFS_baseline=out_RWD$PFS_baseline + rnorm(nrow(out_RWD),0,sd(out_RWD$PFS_baseline))
    }else if(error_X=='bias'){
      out_RWD$PFS_baseline=out_RWD$PFS_baseline + rnorm(nrow(out_RWD),0.3*sd(out_RWD$PFS_baseline),sd(out_RWD$PFS_baseline))
    }else if(error_X=='-sd'){
      out_RWD$PFS_baseline=out_RWD$PFS_baseline - rnorm(nrow(out_RWD),0,sd(out_RWD$PFS_baseline))
    }else if(error_X=='-bias'){
      out_RWD$PFS_baseline=out_RWD$PFS_baseline - rnorm(nrow(out_RWD),0.3*sd(out_RWD$PFS_baseline),sd(out_RWD$PFS_baseline))
    }

    if(error_Y=='sd'){
      out_RWD$PFS_standardized=out_RWD$PFS_standardized + rnorm(nrow(out_RWD),0,sd(out_RWD$PFS_standardized))
    }else if(error_Y=='bias'){
      out_RWD$PFS_standardized=out_RWD$PFS_standardized + rnorm(nrow(out_RWD),0.3*sd(out_RWD$PFS_standardized),sd(out_RWD$PFS_standardized))
    }else if(error_Y=='-sd'){
      out_RWD$PFS_standardized=out_RWD$PFS_standardized - rnorm(nrow(out_RWD),0,sd(out_RWD$PFS_standardized))
    }else if(error_Y=='-bias'){
      out_RWD$PFS_standardized=out_RWD$PFS_standardized - rnorm(nrow(out_RWD),0.3*sd(out_RWD$PFS_standardized),sd(out_RWD$PFS_standardized))
    }else if(error_Y=='exp'){
      out_RWD$PFS_standardized=rnorm(dim(out_RWD)[1],
                                     exp(mu0.X.RWD +
                                           b_base*out_RWD$PFS_baseline +
                                           b_Xn*out_RWD$Xn) +
                                       time_bias.vec, sd=marginal_sd)
    }

    X_c <- subset(out_RCT, select = -c(group, PFS_standardized, Xn));
    Y_c <- out_RCT$PFS_standardized; A_c <- out_RCT$group

    X_h <- subset(out_RWD, select = -c(PFS_standardized, time, Xn))
    Y_h <- out_RWD$PFS_standardized

    data.list <- list(RCT = list(X = X_c, Y = Y_c, A = A_c),
                      RWD = list(X = X_h, Y = Y_h, T_i = out_RWD$time))
    save(data.list,
         file = file.path('RWD_sim',
                          folder.name,
                          file.name))

  }

  return(data.list)


}


getdata_logit <- function(n_c = 100, # RCT control
                          n_t = 200, # RCT treat
                          n_k = 3e3, # RWE data
                          corY_Xn = 0, bias_h = 0.5, prop = NULL,# unmeasured confounder (correlation and biasedness)
                          eff_size = 0, # treatment effect
                          time_bias = 0, # time inconcurrency,
                          error_X = 0, # measurement error
                          error_Y = 0, # outcome validity
                          M = 1, linear = TRUE)
{
  if(linear){
    file.name <- paste0("data_nc.",n_c,"_nt.",n_t,"_nk.",n_k,
                        "_effsize_", eff_size,
                        "_timebias_",time_bias,
                        "_bXn_", paste(corY_Xn, collapse = '_'), "_bias_h_", bias_h, prop,
                        "_error_X_", error_X,
                        '_corr',
                        # "_error_Y_", error_Y,
                        "_rep_", M,
                        ".RData")
    folder.name <- paste0("data_nc.",n_c,"_nt.",n_t,"_nk.",n_k,
                          "_effsize_", eff_size,
                          "_timebias_",time_bias,
                          "_bXn_", paste(corY_Xn, collapse = '_'), "_bias_h_", bias_h, prop,
                          "_error_X_", error_X, '_corr'
                          # "_error_Y_", error_Y
    )
  }else{
    file.name <- paste0("data_nc.",n_c,"_nt.",n_t,"_nk.",n_k,
                        "_effsize_", eff_size,
                        "_timebias_",time_bias,
                        "_bXn_", paste(corY_Xn, collapse = '_'), "_bias_h_", bias_h, prop,
                        "_error_X_", error_X,
                        "_error_Y_", error_Y,
                        '_HTE',
                        "_rep_", M,
                        ".RData")
    folder.name <- paste0("data_nc.",n_c,"_nt.",n_t,"_nk.",n_k,
                          "_effsize_", eff_size,
                          "_timebias_",time_bias,
                          "_bXn_", paste(corY_Xn, collapse = '_'), "_bias_h_", bias_h, prop,
                          "_error_X_", error_X,
                          "_error_Y_", error_Y,
                          '_HTE_logit' # both correct
                          # '_nonlinear_logit' # incorrect mu
                          # '_nonlinear_cal' # incorrect q
    )
  }


  # create the folder if not exists
  ifelse(!dir.exists(file.path('RWD_sim', folder.name)),
         dir.create(file.path('RWD_sim', folder.name)),
         FALSE)
  # if the data has been generated, load it and return
  if(file.exists(file.path('RWD_sim',
                           folder.name,
                           file.name))){
    load(file.path('RWD_sim',
                   folder.name,
                   file.name))}else{
                     set.seed(M)
                     corXjXn <- 0 # U was independent with all covariates (or correlated with 0.3)
                     corY_Ybase <- 0.5 # baseline outcome

                     # read in the data set
                     d <- haven::read_sas('..//real_data//hybrid_monarch_nsai_spotlight.sas7bdat')


                     # Xn for RCT N(0,1) (unobserved confounder in RCT)
                     d$Xn <- rnorm(nrow(d),0,1)
                     # out_RCTcome
                     Y <- 'PFS_standardized'
                     # Categorical Covariates, note: "group" is removed as a covariate as it is the treatment assignment/data source indicator
                     lord0 <- c(NULL,"race","ECOGBL","NATDISI","TGRADE","num_sites_mets","IPDDISTG","PGRSTAT")
                     # continuous Covariates, note: "PFS1" is removed as a covariate because it represents the out_RCTcome
                     lcnt0 <- c(NULL,"AGE_standardized","BMI_standardized","Xn","PFS_standardized")
                     ###Process categorical covariates
                     d$race <- ifelse(d$race=='white','WHITE','NONWHITE')
                     d$ECOGBL <- ifelse(d$ECOGBL=="",NA,ifelse(d$ECOGBL==">1","1",d$ECOGBL))
                     d$TGRADE <- ifelse(d$TGRADE=='Unknown',NA,d$TGRADE)
                     d$num_sites_mets <- ifelse(d$num_sites_mets=="",NA,d$num_sites_mets)
                     d$NATDISI <- ifelse(d$NATDISI=='VISCERAL','VISCERAL','BONE ONLY')
                     d$IPDDISTG <- ifelse(d$IPDDISTG%in%c("missing","Not documented","0"), NA, d$IPDDISTG)
                     d$PGRSTAT <- ifelse(d$PGRSTAT%in%c("","other"),"negative",d$PGRSTAT)
                     ###Create standardized version of continuous covariates, standardized to combined population
                     d$AGE_standardized=(d$AGE-mean(d$AGE,na.rm=TRUE))/sd(d$AGE,na.rm=TRUE)
                     d$BMI_standardized=(d$BMI-mean(d$BMI,na.rm=TRUE))/sd(d$BMI,na.rm=TRUE)
                     d$PFS_standardized=(d$PFS1-mean(d$PFS1,na.rm=TRUE))/sd(d$PFS1,na.rm=TRUE)

                     ### Extract joint distribution of X1,...Xn from the real data
                     lord <- lord0
                     lcnt <- lcnt0[!(lcnt0%in%c("TimeEff","SettingEff"))]

                     tm <- d #%>% filter(cohort==ch)
                     tm <- tm %>% mutate_if(is.character, as.factor)
                     X <- tm[,c(lord,lcnt,Y)]
                     for(col in lord) X[[col]] <- as.numeric(X[[col]])
                     # obtain the correlation matrix for X
                     rho <- cor(X[,c(lord,lcnt)], method="pearson",use="pairwise.complete.obs")

                     # frml <- paste('~',paste(c(lord,lcnt[lcnt!='PFS1']),collapse='+'))
                     # previous_na_action=options('na.action')
                     # options(na.action='na.pass')
                     # Xtgt <- model.matrix(as.formula(frml),tm)[,-1]
                     # options(na.action=previous_na_action$na.action)

                     # Adjust the correlation matrix to reflect the relationship between Xn and X1,...,X_n-1
                     # The target correlation matrix must be ordered ordinal, continuous, Poisson, Negative Binomial
                     Xn.idx <- which(names(X)=="Xn")
                     Xj.idx <- which(names(X)=="BMI_standardized")
                     # Assume that the Xn is at most correlated with one of the covariates
                     rho[, Xn.idx] <- rho[Xn.idx, ] <- 0; rho[Xn.idx, Xn.idx] <- 1
                     rho[Xn.idx, Xj.idx] <- rho[Xj.idx, Xn.idx] <- corXjXn
                     # Assume the baseline is not correlated with any covariate
                     rho['PFS_standardized', ] <- rho[,'PFS_standardized'] <- 0
                     rho['PFS_standardized', 'PFS_standardized'] <- 1

                     ##### Simulate X1,...,Xn for RCT using SimCorrVar function (Plasmode)
                     method = "Polynomial"
                     # ordinal marginal distribution and support
                     marginal <- lapply(lord, function(lord.name)
                     {
                       tab <- as.numeric(table(X[,lord.name]))
                       tab <- cumsum(tab)/sum(tab)
                       tab[-length(tab)]
                     })
                     support <- lapply(lord, function(lord.name)
                     {
                       sort(unique(X[[lord.name]]))
                     })

                     # continuous moments
                     cmoms_RCT <- do.call(rbind, lapply(lcnt, function(lcnt.name)
                     {
                       calc_moments(na.omit(X[[lcnt.name]]))
                     }))

                     # # correction terms
                     # ## first attempt optimization
                     # lSix0 <- lapply(1:nrow(cmoms_RCT), function(i){
                     #   fc <- find_constants(
                     #     method = method,
                     #     skews = cmoms_RCT[i,3], skurts = cmoms_RCT[i,4],
                     #     fifths = cmoms_RCT[i,5], sixths = cmoms_RCT[i,6],
                     #     Six = seq(-10, 100, length.out = 10), seed = M)
                     #   ifelse(length(fc$SixCorr1)==0, 0, fc$SixCorr1)
                     # })
                     # ## revised optimization
                     # lSix <- lapply(1:nrow(cmoms_RCT), function(i){
                     #   fc <- find_constants(
                     #     method = method,
                     #     skews = cmoms_RCT[i,3], skurts = cmoms_RCT[i,4],
                     #     fifths = cmoms_RCT[i,5], sixths = cmoms_RCT[i,6],
                     #     Six = seq(lSix0[[i]]-1, lSix0[[i]]+1,by = .1))
                     #   ifelse(length(fc$SixCorr1)==0, 0, fc$SixCorr1)
                     # })
                     N <- n_c+n_t+n_k
                     sim <- rcorrvar(n = N, k_cont = length(lcnt),k_cat = length(lord),  # change n
                                     method = method,
                                     means =cmoms_RCT[,1],
                                     vars = cmoms_RCT[,2]^2,
                                     skews = cmoms_RCT[,3],
                                     skurts = cmoms_RCT[,4],
                                     fifths = cmoms_RCT[,5],
                                     sixths = cmoms_RCT[,6],
                                     Six = rep(list(seq(0, 100, length.out = 100)),
                                               length(lcnt)+length(lord)),
                                     marginal=marginal, support=support,
                                     Sigma = NULL,
                                     rho = rho, cstart = NULL, seed = M, errorloop = FALSE,
                                     epsilon = 0.001, maxit = 1000, extra_correct = TRUE)
                     out_population <- cbind(sim$ordinal_variables,
                                             sim$continuous_variables)
                     colnames(out_population) <- c(lord,lcnt)
                     # Check marginal distribution of categorical covaraites and compare to observed distribution
                     # for (v in lord0) {print(v);print(prop.table(table(tm[[v]])))}
                     # for (v in lord0) {print(v);print(prop.table(table(out_RCT[[v]])))}

                     # Check marginal distribution of continuous covariates and compare to observed distribution
                     # for (v in lcnt0) {print(v);print(c(summary(tm[[v]]), sd(tm[[v]],na.rm=TRUE)))}
                     # for (v in lcnt0) {print(v);print(c(summary(out_RCT[[v]]), sd(out_RCT[[v]],na.rm=TRUE)))}

                     #Check correlation structure of all covariates
                     # library(corrplot)
                     # corrplot(cor(out_RCT))
                     # corrplot(rho)
                     # rename it as baseline
                     colnames(out_population)[which(colnames(out_population)=='PFS_standardized')] <- "PFS_baseline"
                     # renmar some covariate
                     b_Trt <- eff_size; b_base <- corY_Ybase; b_Xn <- corY_Xn
                     # re-map the ordinal covariates back
                     for(col in lord) {
                       out_population[[col]]=factor(out_population[[col]],
                                             labels=levels(tm[[col]]))
                     }
                     ## fit propensity score for selection
                     X.mat <- model.matrix(~race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+
                                             AGE_standardized+BMI_standardized, data=out_population)
                     X.mat.logit <- model.matrix(~race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+
                                                   AGE_standardized+BMI_standardized,#+
                                                   # AGE_standardized:BMI_standardized +
                                                   # I(AGE_standardized^2) + I(BMI_standardized^2) +
                                                   # I(AGE_standardized^3) + I(BMI_standardized^3),
                                                 data=out_population)
                     ## empirical fit beta
                     fit.delta_c.coef <- glm(relevel(cohort, ref = 'Spotlight')~
                                               race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+
                                               AGE_standardized+BMI_standardized,#+
                                             #   AGE_standardized:BMI_standardized +
                                             # I(AGE_standardized^2) + I(BMI_standardized^2) +
                                             #  I(AGE_standardized^3) + I(BMI_standardized^3),
                                             data=tm,
                          na.action=na.exclude,
                         family = 'binomial')$coefficients
                     fit.delta_c.beta <- fit.delta_c.coef[names(fit.delta_c.coef)!='(Intercept)']
                     ## numerical method for intercept
                     beta_0.c <- uniroot(function(beta_0)
                       {
                       prob.delta_c <- exp(X.mat.logit%*%c(beta_0, fit.delta_c.beta))/
                         (1+exp(X.mat.logit%*%c(beta_0, fit.delta_c.beta)))
                       sum(prob.delta_c)-(n_t+n_c)
                     }, c(-50, 10))$root

                     ## prob for selection indicator delta_c (add unmeasured confounder)
                     prob.delta_c <- exp(X.mat.logit%*%c(beta_0.c, fit.delta_c.beta) +
                                           b_Xn*out_population$Xn)/
                       (1+exp(X.mat.logit%*%c(beta_0.c, fit.delta_c.beta) +
                                b_Xn*out_population$Xn))
                     ## generate for RCT population
                     delta_c.vec <- rbinom(n = N, size = 1, prob = prob.delta_c)
                     out_RCT <- out_population[which(delta_c.vec==1),]
                     ## Randomly assign treatment for RCT on a 2:1 ratio (not involve the covariate)
                     n_RCT <- nrow(out_RCT)
                     out_RCT$group <- rbinom(n_RCT, size = 1,
                                             prob = n_t/(n_t+n_c))

                     ## Generate out_RCTcome as a linear combination of covariates
                     #summary(lm(PFS_standardized~race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+AGE_standardized+BMI_standardized,data=tm,na.action=na.exclude))

                     # drop group

                     if(linear){
                       ## fit the model empirically
                       fit.RCT <- lm(PFS_standardized~race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+
                                       AGE_standardized+BMI_standardized,data=tm,na.action=na.exclude,
                                     subset = (cohort == 'Monarch3'))
                       betas <- fit.RCT$coefficients
                       # cond_var <- summary(fit.RCT)$sigma
                       marginal_sd <- summary(lm(PFS_standardized~group, data=tm, na.action=na.exclude,
                                                 subset = (cohort == 'Monarch3')))$sigma

                       ## construct the outcome
                       b_Trt <- eff_size; b_base <- corY_Ybase; b_Xn <- corY_Xn
                       X.mat <- model.matrix(~race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+
                                               AGE_standardized+BMI_standardized, data=out_RCT)

                       # Conditional outcome generation
                       out_RCT$PFS_standardized <- rnorm(nrow(out_RCT),
                                                         X.mat%*%betas+
                                                           out_RCT$group*b_Trt+
                                                           out_RCT$PFS_baseline*b_base+
                                                           out_RCT$Xn*b_Xn, sd=marginal_sd)
                       coef.b <- lm(PFS_standardized~race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+
                                      AGE_standardized+BMI_standardized+PFS_baseline + group,
                                    data=out_RCT,na.action=na.exclude)$coefficients
                       coef.b.trt <- coef.b['group']
                       mean.diff <- (mean(out_RCT$PFS_standardized[out_RCT$group==1])-
                                       mean(out_RCT$PFS_standardized[out_RCT$group==0]))
                     }else{

                       # fit the model empirically by nonlinear model (e.g., random-forest or spline model)
                       # library(randomForest)
                       library(gam)
                       fit_gam <- gam(PFS_standardized~race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+
                                        AGE_standardized+BMI_standardized,#+
                                      #   AGE_standardized:BMI_standardized +
                                      # I(AGE_standardized^2) + I(BMI_standardized^2) +
                                      #  I(AGE_standardized^3) + I(BMI_standardized^3),
                                      data = tm, na.action=na.exclude,
                                      subset = (cohort == 'Monarch3'))
                       marginal_sd <- summary(lm(PFS_standardized~group, data=tm,
                                                 subset = (cohort == 'Monarch3'), na.action=na.exclude))$sigma

                       # conditional outcome generation
                       ## Y(0)
                       mu0.X <- predict(fit_gam,
                                        newdata = out_RCT[,c(lord,
                                                             'AGE_standardized', 'BMI_standardized')])
                       ## HTE
                       tau.X <- eff_size +
                         c(as.matrix(out_RCT[, c('AGE_standardized', 'BMI_standardized')])%*%c(1,1)) - # add the HTE
                         c(apply(out_RCT[, c('AGE_standardized', 'BMI_standardized')], 2, mean)%*%c(1,1))# adjust the mean

                       ## generate the outcome for RWD
                       out_RCT$PFS_standardized <- rnorm(nrow(out_RCT),
                                                         mu0.X + out_RCT$group*tau.X +
                                                           out_RCT$PFS_baseline*b_base+
                                                           out_RCT$Xn*b_Xn,
                                                         sd=marginal_sd)
                     }



                     # summmary.fit <- data.frame(coef.b.trt=unname(coef.b.trt),
                     #                            mean.diff = mean.diff)
                     # write.csv(summmary.fit,paste("~/real_data/data/Null/summary",M, ".csv",sep = ""),row.names = FALSE)

                     #summary(lm(PFS_standardized~PFS_baseline, data=out_RCT, na.action=na.exclude))
                     #cor(out_RCT$PFS_standardized,out_RCT$PFS_baseline)
                     #cor(out_RCT$PFS_baseline,out_RCT$PFS_standardized)
                     #cor(out_RCT$Xn,out_RCT$PFS_standardized)

                     ####################################################################################################################

                     ## generate covariate information for RWD
                     # ## empirical fit beta
                     # fit.delta_h.coef <- glm(relevel(cohort, ref = 'Monarch3')~
                     #                           race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+
                     #                           AGE_standardized+BMI_standardized,data=tm,na.action=na.exclude,
                     #                         family = 'binomial')$coefficients
                     # fit.delta_h.beta <- fit.delta_h.coef[names(fit.delta_h.coef)!='(Intercept)']
                     # ## numerical method for intercept
                     # beta_0.h <- uniroot(function(beta_0)
                     # {
                     #   prob.delta_h <- exp(X.mat%*%c(beta_0, fit.delta_h.beta))/
                     #     (1+exp(X.mat%*%c(beta_0, fit.delta_h.beta)))
                     #   sum(prob.delta_h)-(n_k)
                     # }, c(-50, 10))$root
                     #
                     # ## prob for selection indicator delta_c
                     # prob.delta_h <- exp(X.mat%*%c(beta_0.h, fit.delta_h.beta))/
                     #   (1+exp(X.mat%*%c(beta_0.h, fit.delta_h.beta)))
                     ## generate for RCT population
                     out_RWD <- out_population[which(delta_c.vec==0),]
                     # out_RWD <- out_population[which(rbinom(n = N, size = 1, prob = prob.delta_h)==1),]

                     # colnames(out_RWD) <- c(lord,lcnt); colnames(out_RWD)[11] <- "PFS_baseline"


                     # Check marginal distribution of categorical covaraites and compare to observed distribution
                     # for (v in lord0) {print(v);print(prop.table(table(tm[[v]])))}
                     # for (v in lord0) {print(v);print(prop.table(table(out_RWD[[v]])))}

                     #Check marginal distribution of continuous covariates and compare to observed distribution
                     # for (v in lcnt0) {print(v);print(c(summary(tm[[v]]), sd(tm[[v]],na.rm=TRUE)))}
                     # for (v in lcnt0) {print(v);print(c(summary(out_RWD[[v]]), sd(out_RWD[[v]],na.rm=TRUE)))}
                     # mean(out_RWD$PFS_baseline)
                     #
                     #Check correlation structure of all covariates
                     # dat=out_RWD %>% mutate_if(is.factor, as.numeric)
                     # print(cor(dat))
                     # print(rho)


                     # add unmeasured confounder
                     out_RWD$Xn <- out_RWD$Xn + bias_h
                     if(!is.null(prop)){
                       # randomly delete Xn to zero
                       zero.idx <- sample(length(out_RWD$Xn),
                                          length(out_RWD$Xn)*prop)
                       out_RWD$Xn[zero.idx] <- 0
                     }

                     for(col in lord) {
                       out_RWD[[col]] <- factor(out_RWD[[col]],labels=levels(tm[[col]]))
                     }

                     X <- model.matrix(~race+ECOGBL+NATDISI+TGRADE+num_sites_mets+IPDDISTG+PGRSTAT+
                                         AGE_standardized+BMI_standardized,data=out_RWD)

                     # Generate variable for time inconcurrency
                     if(time_bias=='unif')
                     {
                       out_RWD$time <- runif(dim(out_RWD)[1])
                       ## Generate out_RWDcome as a linear combination of covariates
                       time_bias.vec <- 0.3*out_RWD$time

                     }else{
                       out_RWD$time <- sample(c(0,1,2), size=dim(out_RWD)[1], replace=TRUE)
                       ## Generate out_RWDcome as a linear combination of covariates
                       time_bias.vec <-  time_bias*out_RWD$time
                     }

                     # linear combination or not
                     if(linear){

                       out_RWD$PFS_standardized <- rnorm(nrow(out_RWD), X%*%betas +
                                                           b_base*out_RWD$PFS_baseline +
                                                           b_Xn*out_RWD$Xn +
                                                           time_bias.vec, sd=marginal_sd)
                     }else{

                       mu0.X.RWD <- predict(fit_gam,
                                            newdata = out_RWD[,c(lord, 'AGE_standardized', 'BMI_standardized')])
                       out_RWD$PFS_standardized <- rnorm(nrow(out_RWD), mu0.X.RWD +
                                                           b_base*out_RWD$PFS_baseline +
                                                           b_Xn*out_RWD$Xn +
                                                           time_bias.vec, sd=marginal_sd)
                     }






                     #cor(out_RWD$PFS_baseline,out_RWD$PFS_standardized)

                     #Add in measurement error
                     if(error_X=='sd'){
                       out_RWD$PFS_baseline=out_RWD$PFS_baseline + rnorm(nrow(out_RWD),0,sd(out_RWD$PFS_baseline))
                     }else if(error_X=='bias'){
                       out_RWD$PFS_baseline=out_RWD$PFS_baseline + rnorm(nrow(out_RWD),0.3*sd(out_RWD$PFS_baseline),sd(out_RWD$PFS_baseline))
                     }else if(error_X=='-sd'){
                       out_RWD$PFS_baseline=out_RWD$PFS_baseline - rnorm(nrow(out_RWD),0,sd(out_RWD$PFS_baseline))
                     }else if(error_X=='-bias'){
                       out_RWD$PFS_baseline=out_RWD$PFS_baseline - rnorm(nrow(out_RWD),0.3*sd(out_RWD$PFS_baseline),sd(out_RWD$PFS_baseline))
                     }

                     if(error_Y=='sd'){
                       out_RWD$PFS_standardized=out_RWD$PFS_standardized + rnorm(nrow(out_RWD),0,sd(out_RWD$PFS_standardized))
                     }else if(error_Y=='bias'){
                       out_RWD$PFS_standardized=out_RWD$PFS_standardized + rnorm(nrow(out_RWD),0.3*sd(out_RWD$PFS_standardized),sd(out_RWD$PFS_standardized))
                     }else if(error_Y=='-sd'){
                       out_RWD$PFS_standardized=out_RWD$PFS_standardized - rnorm(nrow(out_RWD),0,sd(out_RWD$PFS_standardized))
                     }else if(error_Y=='-bias'){
                       out_RWD$PFS_standardized=out_RWD$PFS_standardized - rnorm(nrow(out_RWD),0.3*sd(out_RWD$PFS_standardized),sd(out_RWD$PFS_standardized))
                     }else if(error_Y=='exp'){
                       out_RWD$PFS_standardized=rnorm(dim(out_RWD)[1],
                                                      exp(mu0.X.RWD +
                                                            b_base*out_RWD$PFS_baseline +
                                                            b_Xn*out_RWD$Xn) +
                                                        time_bias.vec, sd=marginal_sd)
                     }

                     X_c <- subset(out_RCT, select = -c(group, PFS_standardized, Xn));
                     Y_c <- out_RCT$PFS_standardized; A_c <- out_RCT$group

                     X_h <- subset(out_RWD, select = -c(PFS_standardized, time, Xn))
                     Y_h <- out_RWD$PFS_standardized

                     data.list <- list(RCT = list(X = X_c, Y = Y_c, A = A_c),
                                       RWD = list(X = X_h, Y = Y_h, T_i = out_RWD$time))
                     save(data.list,
                          file = file.path('RWD_sim',
                                           folder.name,
                                           file.name)
                          )

                   }

  return(data.list)


}


get_summary <- function(res)
{
  tau.AIPW <- res$est$AIPW
  tau.ACW.mat <-  res$est$ACW #do.call(rbind, res$est$ACW)
  tau.ACW.lasso <- res$est$ACW.lasso
  tau.ACW.final <- res$est$ACW.final

  # gbm can have extreme prediction error
  tau.ACW.lasso[abs(tau.ACW.lasso)>3] <- NA
  tau.ACW.final[abs(tau.ACW.final)>3] <- NA

  res$var.est$ACW.final[abs(res$var.est$ACW.final)>30] <- NA
  # point estimator
  ## bias
  bias.tau <- list(AIPW = mean(tau.AIPW, na.rm = T)-0,
                   ACW = apply(tau.ACW.mat, 2, mean, na.rm = T) - 0,
                   ACW.lasso = mean(tau.ACW.lasso, na.rm = T)-0,
                   ACW.final =mean(tau.ACW.final, na.rm = T)-0)
  ## var
  var.tau <- list(AIPW = var(tau.AIPW, na.rm = T),
                  ACW = apply(tau.ACW.mat, 2, var, na.rm = T),
                  ACW.lasso = var(tau.ACW.lasso, na.rm = T),
                  ACW.final = var(tau.ACW.final, na.rm = T))
  ## mse
  mse.tau <- list(AIPW = mean((tau.AIPW-0)**2, na.rm = T),
                  ACW = apply((tau.ACW.mat-0)**2, 2, mean, na.rm = T),
                  ACW.lasso = mean((tau.ACW.lasso-0)**2, na.rm = T),
                  ACW.final = mean((tau.ACW.final-0)**2, na.rm = T))
  ## type I error
  type.1.error <- list(AIPW = mean(res$CI$AIPW[,'lower']>0|
                                     res$CI$AIPW[,'upper']<0, na.rm = T),
                       ACW = apply(res$CI$ACW$lower>0|res$CI$ACW$upper<0,
                                   2, mean, na.rm = T),
                       ACW.lasso = mean(res$CI$ACW.lasso[,'lower']>0|
                                          res$CI$ACW.lasso[,'upper']<0, na.rm = T),
                       ACW.final = mean(res$CI$ACW.final[,'lower']>0|
                                          res$CI$ACW.final[,'upper']<0, na.rm = T)
                       # ACW.lasso.conservative = mean(res$CI$ACW.lasso.cont[,'lower']>0|
                       #                         res$CI$ACW.lasso.cont[,'upper']<0, na.rm = T)
  )
  ## CP1
  CP1 <- list(AIPW = mean(res$CI$AIPW[,'lower']<0&
                            res$CI$AIPW[,'upper']>0, na.rm = T),
              ACW = apply(res$CI$ACW$lower<0&res$CI$ACW$upper>0,
                          2, mean, na.rm = T),
              ACW.lasso = mean(res$CI$ACW.lasso[,'lower']<0&
                                 res$CI$ACW.lasso[,'upper']>0, na.rm = T),
              ACW.final = mean(res$CI$ACW.final[,'lower']<0&
                                 res$CI$ACW.final[,'upper']>0, na.rm = T)#,
              # ACW.lasso.conservative = mean(res$CI$ACW.lasso.cont[,'lower']<0&
              #                         res$CI$ACW.lasso.cont[,'upper']>0, na.rm = T)
  )
  ## var.est
  var.est <- list(AIPW = mean(res$var.est$AIPW/res$n_c, na.rm = TRUE),
                  ACW = apply(as.matrix(res$var.est$ACW)/res$n_c, 2, mean, na.rm = TRUE),
                  ACW.lasso = mean(res$var.est$ACW.lasso/res$n_c, na.rm = TRUE),
                  ACW.final = mean(res$var.est$ACW.final/res$n_c, na.rm = TRUE))

  CI.len <- list(AIPW = mean(res$CI$AIPW[,'upper'] -
                               res$CI$AIPW[,'lower'],
                             na.rm = T),
                 ACW = apply(res$CI$ACW$upper -
                               res$CI$ACW$lower,
                             2, mean, na.rm = T),
                 ACW.lasso = mean(res$CI$ACW.lasso[,'upper'] -
                                    res$CI$ACW.lasso[,'lower'], na.rm = T),
                 ACW.final = mean(res$CI$ACW.final[,'upper'] -
                                    res$CI$ACW.final[,'lower'], na.rm = T)#,
                 # ACW.lasso = mean(res$CI$ACW.lasso.cont[,'upper'] -
                 #                    res$CI$ACW.lasso.cont[,'lower'],
                 #                  na.rm = T)
  )

  ## prob
  prob.comb <- res$prob
  return(list(bias = bias.tau,
              var = var.tau, var.est = var.est,
              mse = mse.tau,
              CP = CP1,
              typeIerror = type.1.error,
              prob.comb = prob.comb))
}

# control sample is 200
get_benchmark <- function(delta1 = 0,
                          omega = 0,
                          gamma_mu = 0,
                          error_Y = 0)
{
  res.beta0.c200 <- sim.eval.real(N_CD = 2e2,
                                  beta_T = 0,
                                  delta1 = delta1,
                                  omega = omega,
                                  gamma_mu = gamma_mu,
                                  error_Y = error_Y)
  summary.beta0.c200 <- get_summary(res.beta0.c200)
  # load another RData
  res.beta1.c200 <- sim.eval.real(N_CD = 2e2,
                                  beta_T = 0.3,
                                  delta1 = delta1,
                                  omega = omega,
                                  gamma_mu = gamma_mu,
                                  error_Y = error_Y)
  ## CI len
  ## power
  power <- list(AIPW = 1 - mean(res.beta1.c200$CI$AIPW[,'lower']<0&
                                  res.beta1.c200$CI$AIPW[,'upper']>0, na.rm = T),
                ACW = 1 - apply(res.beta1.c200$CI$ACW$lower<0&res.beta1.c200$CI$ACW$upper>0,
                                2, mean, na.rm = T),
                # ACW.lasso = 1 - mean(res$CI$ACW.lasso[,'lower']<0&
                #                        res$CI$ACW.lasso[,'upper']>0, na.rm = T),
                ACW.lasso = 1 - mean(res.beta1.c200$CI$ACW.final[,'lower']<0&
                                       res.beta1.c200$CI$ACW.final[,'upper']>0, na.rm = T))

  CP2 <- list(AIPW = mean(res.beta1.c200$CI$AIPW[,'lower']<0.3&
                                  res.beta1.c200$CI$AIPW[,'upper']>0.3, na.rm = T),
                ACW = apply(res.beta1.c200$CI$ACW$lower<0.3&res.beta1.c200$CI$ACW$upper>0.3,
                                2, mean, na.rm = T),
                # ACW.lasso = 1 - mean(res$CI$ACW.lasso[,'lower']<0&
                #                        res$CI$ACW.lasso[,'upper']>0, na.rm = T),
                ACW.lasso = mean(res.beta1.c200$CI$ACW.final[,'lower']<0.3&
                                       res.beta1.c200$CI$ACW.final[,'upper']>0.3, na.rm = T))


  lapply(list(c(summary.beta0.c200, power = list(power),
                CP2 = list(CP2))),
         function(dat){c(bias = dat$bias$AIPW,
                         var = dat$var$AIPW,
                         mse = dat$mse$AIPW,
                         typeI = dat$typeIerror$AIPW,
                         power = dat$power$AIPW,
                         CP2 = dat$CP2$AIPW
         )})%>%unlist()
}
