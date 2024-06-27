remove_backticks = function(str) {
  r = sapply(str, function(x) {
    if(substring(x, 1, 1) == "`") x = substring(x, 2); l = nchar(x);
    if(substring(x, l, l) == "`") x = substring(x, 1, l - 1); x })
  r = as.character(r)
  r
}

# calibration equation
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
  sum(c(left_side-right_side)**2)
}

# estimate the outcome modeling for the external control dataset
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

  # use machine learning model to fitting mu_1 and mu_0
  fit.Y.HC <- train(Y ~ .+0, data = data.frame(Y = Y_h, X_h),
                     trControl = hc.ctrl)


  # estimate odds-ratio
  if(ncol(X_c)<n_h)
  {

    lambda_hat <- tryCatch(BB::dfsane(par = rep(0, dim(X_c)[2]), # initial points
                        fn=cal_eqs,
                        X_c = X_c, X_h = X_h,#[bias_b==0,],
                        control = list(trace=T,
                                       NM=T,
                                       BFGS=F,
                                       tol=1.e-8,
                                       maxit = 500))$par,
             warning = function(w){optim(par = rep(0, dim(X_c)[2]), # initial points
                                         fn=cal_eqs.GMM,
                                         X_c = X_c, X_h = X_h,#[bias_b==0,],
                                         control = list(trace=T,
                                                        abstol=1.e-8,
                                                        maxit = 500),
                                         method = 'BFGS')$par})

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

  S.ACW.q <- rbind(-X_c, q_hat*X_h)
  dot.S.ACW.q <- t(q_hat*X_h)%*%X_h/n_c

  r_X <- mean((predict(fit.Y.RCT, data.frame(X_c[A_c==0,], A = 0)) - mean(Y_c[A_c==0]))**2)/
    tryCatch(mean((predict(fit.Y.HC, data.frame(X_h[bias_b==0, ])) - mean(Y_h[bias_b==0]))**2),
             error = function(e)1e6) # if no selected, down weight the HC with extreme values

  # prob_b_zero_c <- prob_b_zero_h <- sum(bias_b==0)/n_h

  # fit a glm for b==0 against X_h
  # fit.b0.hc <- glm(b0 ~ . +0, data = data.frame(b0 = (bias_b == 0),
  #                                               X_h),
  #                  family = 'binomial')
  # prob_b_zero_c <- predict(fit.b0.hc, newdata = data.frame(X_c),
  #                          type = 'response')
  # prob_b_zero_h <- predict(fit.b0.hc, newdata = data.frame(X_h),
  #                          type = 'response')
  b_c.idx <- 0
  b_h.idx <- bias_b

  # fit a glm for A against X_c
  fit.ps <- glm(A~.+0, data = data.frame(A = A_c,
                                       X_c),
                family = 'binomial')
  prob_A.h <- predict(fit.ps, newdata = data.frame(X_h),
                      type = 'response')


  # construct the influence function
  dot.mu0_c.q <- apply(c((q_hat.c * r_X * (b_c.idx == 0))*
                           (1-A_c)*(Y_c - predict(fit.Y.RCT, data.frame(X_c, A = 0)))/
                           ((r_X * (b_c.idx == 0) + (1-prob_A)*q_hat.c)**2))*X_c, 2, sum,
                       na.rm = TRUE)/n_c

  dot.mu0_h.q <- apply(c((q_hat * (b_h.idx == 0) * r_X**2)*
                           (b_h.idx == 0) * (Y_h - predict(fit.Y.RCT, data.frame(X_h, A = 0)))/
                           ((r_X * (b_h.idx == 0) + (1-prob_A.h)*q_hat)**2))*X_h, 2, sum,
                       na.rm = TRUE)/n_c
  mu0_c <- c(predict(fit.Y.RCT, data.frame(X_c, A = 0)) +
               q_hat.c/((b_c.idx == 0) * r_X + (1-prob_A)*q_hat.c)*
               (1-A_c)*(Y_c - predict(fit.Y.RCT, data.frame(X_c, A = 0))),
             rep(0, n_h))
  mu0_h <- c(rep(0, n_c),
             q_hat*r_X/((b_h.idx == 0) * r_X + (1-prob_A.h)*q_hat) *(b_h.idx == 0)*
               (Y_h - predict(fit.Y.RCT, data.frame(X_h, A = 0))))
  mu0_i <- mu0_c + mu0_h
  mu0_score <-  mu0_i -  S.ACW.q%*%t(dot.mu0_c.q%*%MASS::ginv(dot.S.ACW.q)) -
    S.ACW.q%*%t(dot.mu0_h.q%*%MASS::ginv(dot.S.ACW.q))

  mu0_ACW.hat <- (sum(mu0_i))/n_c

  # construct the pseudo-observations for the bias parameters
  bias_i <- Y_h - predict(fit.Y.RCT, data.frame(X_h, A = 0))

  # bias_i <- q_hat * c(Y_h - predict(fit.Y.RCT, data.frame(X_h, A = 0)))*(n_c+n_h)/n_c #-
  #   # {q_hat*(n_c+n_h)/n_c - 1}*b.X

  # compute the score for estimating q
  bias_i.score <- c(rep(0, n_c), bias_i)

  # dot.bias.ACW.q <- apply(q_hat*c(Y_h - predict(fit.Y.RCT, data.frame(X_h, A = 0)))*X_h,
  #                         2, sum)*(n_c+n_h)/n_c**2
  # bias_i.score <- c(rep(0, n_c), bias_i) -
  #   S.ACW.q%*%t(dot.bias.ACW.q%*%MASS::ginv(dot.S.ACW.q))

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

# handy functions
expit <- function(x){exp(x)/{1+exp(x)}}

