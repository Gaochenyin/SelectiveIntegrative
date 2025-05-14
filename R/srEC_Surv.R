#' Selective and robust external-control borrowing strategy for evaluating treatment effect
#'
#' @description
#' `srEC_Surv()` is a penalized dynamic integrative framework for survival outcomes that augments
#'  a randomized trial (RT) with a external control (EC) dataset. The subject-level compatibility of the EC
#'  is assessed by a well-crafted penalty (e.g., adaptive lasso penalty). The parameter of interest
#'  is the average treatment effect in terms of the restricted mean survival time (RMST).
#'
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom survival coxph Surv
#' @param data_rt A list contains X, Y, D and A for RT. The propensity scores P(A=1|X) can also be
#' contained as `prob_A` (or `NULL`).
#' @param data_ec A list contains X, Y, and D for EC with treatment A=0.
#' @param k_grid A vector of time points to evaluate the survival curves.
#' @param X.pred.R A list of covariate name to model the propensity for R
#' @param X.pred.A A list of covariate names to model the propensity for A
#' @param X.pred.c A list of covariate names to model the survival curve of the event
#' @param X.pred.f A list of covariate names to model the survival curve of the censoring
#'
#' @returns A list with components:
#' * aipcw.RCT: RCT-only estimate
#' * aipcw: RCT and full EC estimate
#' * aipcw_lasso: RCT with selective EC estimate
#' * subset.idx: a subset of indices of the external controls which have been selected for the
#' final integrative estimation.
#' @export
srEC_Surv <- function(data_rt,
                      data_ec,
                      k_grid,
                      X.pred.R = c('X1', 'X2', 'X3'),
                      X.pred.A = c('X1', 'X2', 'X3'),
                      X.pred.c = c('X1', 'X2', 'X3'),
                      X.pred.f = c('A', 'X1', 'X2', 'X3')
                      ){

  # model fitting
  dat.combine <- rbind(data.frame(R = 1, data_rt),
                       data.frame(R = 0, data_ec, A = 0))

  N_1 <- sum(dat.combine$R == 1 & dat.combine$A == 1)
  N_c <- sum(dat.combine$R == 1)

  ## propensity score for R
  fit.R <- glm(as.formula(paste('R~',
                                paste(X.pred.R, collapse = '+'))),
               data = dat.combine,
               family = 'binomial')
  prob.R <- fit.R$fitted.values
  q.R <- prob.R / (1 - prob.R)

  ## propensity score for A
  ps.A <- data_rt$prob_A
  ## RCT
  ## propensity score
  if(is.null(ps.A)){
    fit.ps <- glm(as.formula(paste('A~', paste(X.pred.A, collapse = '+'))),
                  data = data_rt,
                  family = 'binomial')
    ps.A <-  fit.ps$fitted.values
  }

  # Cox PH survival modeling
  fit.event <- survival::coxph(as.formula(paste('Surv(Y,D == 1)~R+(',
                                      paste0(X.pred.f,
                                             collapse = '+'),
                                      ')')),
                     data = dat.combine)
  fit.cens <- survival::coxph(as.formula(paste('Surv(Y,D == 0)~R+(',
                                     paste0(X.pred.c,
                                            collapse = '+'),
                                     ')')),
                    data = dat.combine)

  ## AIPCW for each time point
  time_grid <- sort(dat.combine$Y)

  # estimation  for the survival curves S1
  aipcw.S1.est <- function(dat,
                           fit.surv.obj,
                           fit.cens.obj,
                           q.R, ps.A,
                           time_grid, # estimation
                           k_grid # evaluation
  ){

    event.base <- survival::basehaz(fit.surv.obj, centered = FALSE)
    cens.base <- survival::basehaz(fit.cens.obj, centered = FALSE)

    risk.event <- predict(fit.surv.obj,
                          newdata = data.frame(R = 1,
                                               dat[, X.pred.f]),
                          type = 'risk')

    risk.event.A1 <- predict(fit.surv.obj,
                             newdata = data.frame(R = 1,
                                                  A = 1,
                                                  dat[, setdiff(X.pred.f, 'A')]),
                             type = 'risk')

    risk.cens <- predict(fit.cens.obj,
                         newdata = data.frame(R = 1,
                                              dat[, X.pred.c]),
                         type = 'risk')

    # basic cumulative hazard and hazard
    event.cumhaz <- event.base$hazard[findInterval(time_grid,
                                                   event.base$time)]
    cens.cumhaz <- cens.base$hazard[findInterval(time_grid,
                                                 cens.base$time)]

    event.hazard.hat <- outer(risk.event,
                              c(event.cumhaz[1], diff(event.cumhaz)/diff(time_grid)))
    cens.hazard.hat <- outer(risk.cens,
                             c(cens.cumhaz[1], diff(cens.cumhaz)/diff(time_grid)))

    S.event.pred <- exp(-outer(risk.event,
                               event.cumhaz, "*"))
    S.cens.pred <- exp(-outer(risk.cens,
                              cens.cumhaz, "*"))

    # evaluate S1 and S1c on the grid
    S1.pred <- exp(-outer(risk.event.A1,
                          event.base$hazard[findInterval(k_grid,
                                                         event.base$time)]))
    Sc1.pred <- exp(-outer(risk.cens,
                           cens.base$hazard[findInterval(k_grid,
                                                         cens.base$time)]))
    # trim off small values
    threshold <- 1e-8
    S1.pred[S1.pred<threshold] <- threshold
    Sc1.pred[Sc1.pred<threshold] <- threshold
    S.event.pred[S.event.pred<threshold] <- threshold
    S.cens.pred[S.cens.pred<threshold] <- threshold

    # construct the IPCW
    S1.hat.IPCW <- sapply(1:length(k_grid), function(tk_idx){
      t_k <- k_grid[tk_idx]
      (dat$A/ps.A/Sc1.pred[, tk_idx] * (dat$Y > t_k) +
          (1 - dat$A/ps.A) * S1.pred[, tk_idx])
    })


    # construct the AIPW
    ## S1
    S1.hat.EE.mat <- sapply(1:length(k_grid), function(tk_idx){
      t_k <- k_grid[tk_idx]
      part1_sum <- rep(0, length(ps.A))

      for (j in 1:length(time_grid)) {
        t_j <- time_grid[j]
        h <- -as.numeric(dat$A/ps.A/S.cens.pred[, j] *
                           S1.pred[, tk_idx]/
                           S.event.pred[, j]) * (t_j <= t_k)
        part1 <- h * (
          (dat$Y == t_j) * (dat$D == 1) -
            (dat$Y >= t_j) * event.hazard.hat[, j]
        )
        # interval length
        if(j==1){
          len <- t_j
        }else{
          len <- t_j - time_grid[j-1]
        }
        part1_sum <- part1_sum + part1 * len
      }

      (dat$A/ps.A/Sc1.pred[, tk_idx] * (dat$Y > t_k) +
          (1 - dat$A/ps.A) * S1.pred[, tk_idx]
        # + part1_sum
      )
    })
    return(list(OI = S1.pred,
                IPCW.RCT = S1.hat.IPCW,
                AIPW.RCT = S1.hat.EE.mat))
  }


  # estimation for the survival curves S0
  aipcw.S0.est <- function(dat.RCT, dat.RWD,
                           fit.surv.obj,
                           fit.cens.obj,
                           ps.A,
                           time_grid, # estimation for integral
                           k_grid, X.pred.R = c('X1', 'X2', 'X3'), # evaluation
                           prob.R = NULL, # refit after seleciton
                           select_subset = NULL,
                           model_list = NULL
  ){
    dat.combine <- rbind(data.frame(R = 1, dat.RCT),
                         data.frame(R = 0, dat.RWD, A = 0))

    if(is.null(select_subset)){
      select_subset <- (dat.combine$R == 0) # include all ECs
      prob_subset <- 1
    }else{
      # fit a model for prob_subset
      fit.b <- glm(b~., data = data.frame(b = select_subset[dat.combine$R==0],
                                          dat.combine[dat.combine$R==0, X.pred.R]),
                   family = 'binomial')
      prob_subset <- predict(fit.b,
                             newdata = dat.combine,
                             type = 'response',
                             allow.new.levels=TRUE)
      # incorporate that in the calibration weight
      # prob_subset <- 1

      # prob_subset <- sum(select_subset)/nrow(dat.RWD)
    }
    # PS
    ps.A <- c(ps.A, rep(0, nrow(dat.RWD)))

    # CW
    if(is.null(prob.R)){
      fit.R <- glm(as.formula(paste('R~',
                                    paste(X.pred.R, collapse = '+'))),
                   data = dat.combine,
                   subset = dat.combine$R|select_subset,
                   family = 'binomial')
      prob.R <- predict(fit.R, newdata = dat.combine,
                        type = 'response')

      # cap the value of probR
      prob.R[prob.R>.99] <- .99
      prob.R[prob.R<.01] <- .01
    }

    q.R <- prob.R/(1 - prob.R)

    # modeling
    if(is.null(model_list)){
      # covariates with all zeros
      # baseline hazard functions
      event.base.R1 <- event.base.R0 <-
        survival::basehaz(fit.surv.obj, centered = FALSE)

      cens.base.R1 <- cens.base.R0 <-
        survival::basehaz(fit.cens.obj, centered = FALSE)

      model_list$risk.event.R1 <- predict(fit.surv.obj,
                                          newdata = data.frame(R = 1,
                                                               dat.combine[, X.pred.f]),
                                          type = 'risk')
      model_list$risk.event.R0 <- predict(fit.surv.obj,
                                          newdata = data.frame(R = 0,
                                                               dat.combine[, X.pred.f]),
                                          type = 'risk')

      risk.event.R1A0 <- predict(fit.surv.obj,
                                 newdata = data.frame(R = 1,
                                                      A = 0,
                                                      dat.combine[, setdiff(X.pred.f, 'A')]),
                                 type = 'risk')

      risk.event.R0A0 <- predict(fit.surv.obj,
                                 newdata = data.frame(R = 0,
                                                      A = 0,
                                                      dat.combine[, setdiff(X.pred.f, 'A')]),
                                 type = 'risk')



      risk.cens <- predict(fit.cens.obj,
                           newdata = dat.combine[, c('R', X.pred.c)],
                           type = 'risk')

      risk.cens.R1 <- predict(fit.cens.obj,
                              newdata = data.frame(R = 1,
                                                   dat.combine[, setdiff(X.pred.c, 'A')]),
                              type = 'risk')

      risk.cens.R0 <- predict(fit.cens.obj,
                              newdata = data.frame(R = 0,
                                                   dat.combine[, setdiff(X.pred.c, 'A')]),
                              type = 'risk')
      # basic cumulative hazard and hazard
      event.cumhaz.R0 <- event.base.R0$hazard[findInterval(time_grid,
                                                           event.base.R0$time)]
      event.cumhaz.R1 <- event.base.R1$hazard[findInterval(time_grid,
                                                           event.base.R1$time)]

      cens.cumhaz.R0 <- cens.base.R0$hazard[findInterval(time_grid,
                                                         cens.base.R0$time)]
      cens.cumhaz.R1 <- cens.base.R1$hazard[findInterval(time_grid,
                                                         cens.base.R1$time)]

      model_list$eventR1.hazard.hat <- outer(risk.event.R1A0,
                                             c(event.cumhaz.R1[1], diff(event.cumhaz.R1)))
      model_list$eventR0.hazard.hat <- outer(risk.event.R0A0,
                                             c(event.cumhaz.R0[1], diff(event.cumhaz.R0)))

      # predicted probabilities on the time grid
      model_list$SR1.event.pred <- exp(-outer(model_list$risk.event.R1,
                                              event.cumhaz.R1, "*"))
      model_list$SR0.event.pred <- exp(-outer(model_list$risk.event.R0,
                                              event.cumhaz.R0, "*"))

      model_list$SR1.cens.pred <- exp(-outer(risk.cens.R1,
                                             cens.cumhaz.R1, "*"))
      model_list$SR0.cens.pred <- exp(-outer(risk.cens.R0,
                                             cens.cumhaz.R0, "*"))

      # evaluate S1 and S1c on the grid
      model_list$S0R0.pred <- exp(-outer(risk.event.R0A0,
                                         event.base.R0$hazard[findInterval(k_grid,
                                                                           event.base.R0$time)]))
      model_list$S0R1.pred <- exp(-outer(risk.event.R1A0,
                                         event.base.R1$hazard[findInterval(k_grid,
                                                                           event.base.R1$time)]))

      model_list$ScR0.pred <- exp(-outer(risk.cens.R0,
                                         cens.base.R0$hazard[findInterval(k_grid,
                                                                          cens.base.R0$time)]))
      model_list$ScR1.pred <- exp(-outer(risk.cens.R1,
                                         cens.base.R1$hazard[findInterval(k_grid,
                                                                          cens.base.R1$time)]))

    }



    # trim off small values
    threshold <- 1e-8
    model_list$S0R0.pred[model_list$S0R0.pred<threshold] <- threshold
    model_list$S0R1.pred[model_list$S0R1.pred<threshold] <- threshold
    model_list$ScR0.pred[model_list$ScR0.pred<threshold] <- threshold
    model_list$ScR1.pred[model_list$ScR1.pred<threshold] <- threshold

    model_list$SR1.event.pred[model_list$SR1.event.pred<threshold] <- threshold
    model_list$SR1.cens.pred[model_list$SR1.cens.pred<threshold] <- threshold
    model_list$SR0.event.pred[model_list$SR0.event.pred<threshold] <- threshold
    model_list$SR0.cens.pred[model_list$SR0.cens.pred<threshold] <- threshold

    # construct the RCT-only IPCW
    S0.hat.IPCW <- sapply(1:length(k_grid), function(tk_idx){
      t_k <- k_grid[tk_idx]
      (dat.combine$R * (1 - dat.combine$A)/(1 - ps.A)/model_list$ScR1.pred[, tk_idx] *
          (dat.combine$Y > t_k) +
          dat.combine$R *
          {dat.combine$A - ps.A}/{1 - ps.A} *
          model_list$S0R1.pred[, tk_idx])
    })


    # construct the RCT-only AIPW
    S0.hat.EE.mat.RCT <- sapply(1:length(k_grid), function(tk_idx){
      t_k <- k_grid[tk_idx]
      part1_sum <- rep(0, length(ps.A))

      for (j in 1:length(time_grid)) {
        t_j <- time_grid[j]

        h.m1 <- -as.numeric(dat.combine$R*(1-dat.combine$A)/model_list$SR1.cens.pred[, j]/
                              (1 - ps.A) *
                              model_list$S0R1.pred[, tk_idx]/
                              model_list$SR1.event.pred[, j]) * (t_j <= t_k)

        part1 <- h.m1 * (
          (dat.combine$Y == t_j) * (dat.combine$D == 1) -
            (dat.combine$Y >= t_j) * model_list$eventR1.hazard.hat[, j]
        )

        # interval length
        if(j==1){
          len <- t_j
        }else{
          len <- t_j - time_grid[j-1]
        }
        part1_sum <- part1_sum + part1 * len
      }

      (dat.combine$R * (1 - dat.combine$A)/(1 - ps.A)/model_list$ScR1.pred[, tk_idx] *
          (dat.combine$Y > t_k) +
          dat.combine$R *
          {dat.combine$A - ps.A}/{1 - ps.A} *
          model_list$S0R1.pred[, tk_idx] +
          part1_sum)
    })
    # apply(S0.hat.EE.mat.RCT, 2, sum)/N.R

    # construct the EC-only AIPW
    S0.hat.EE.mat.EC <- sapply(1:length(k_grid), function(tk_idx){
      t_k <- k_grid[tk_idx]
      part2_sum <- rep(0, length(ps.A))

      for (j in 1:length(time_grid)) {
        t_j <- time_grid[j]

        h.m2 <- -as.numeric((1 - dat.combine$R)*q.R/model_list$SR0.cens.pred[, j] *
                              model_list$S0R0.pred[, tk_idx]/
                              model_list$SR0.event.pred[, j]) * (t_j <= t_k)

        part2 <- h.m2 * (
          (dat.combine$Y == t_j) * (dat.combine$D == 1) -
            (dat.combine$Y >= t_j) * model_list$eventR0.hazard.hat[, j]
        )
        # interval length
        if(j==1){
          len <- t_j
        }else{
          len <- t_j - time_grid[j-1]
        }
        part2_sum <- part2_sum + part2 * len
      }

      (1 - dat.combine$R) * q.R/model_list$ScR0.pred[, tk_idx]*
        (dat.combine$Y > t_k) +

        dat.combine$R * model_list$S0R1.pred[, tk_idx] -
        (1 - dat.combine$R) * q.R *
        model_list$S0R0.pred[, tk_idx] +

        part2_sum
    })

    # apply(S0.hat.EE.mat.EC, 2, sum)/N.R

    # construct the AIPW
    r.X <-  {model_list$S0R1.pred * (1 - model_list$S0R1.pred)}/
      {model_list$S0R0.pred * (1 - model_list$S0R0.pred)}
    r.X[is.na(r.X)] <- 1

    S0.hat.EE.mat <- sapply(1:length(k_grid), function(tk_idx){
      t_k <- k_grid[tk_idx]
      part1_sum <- rep(0, length(ps.A))
      part2_sum <- rep(0, length(ps.A))

      D.X <- r.X[, tk_idx] * prob_subset +
        (1-ps.A) * q.R

      for (j in 1:length(time_grid)) {
        t_j <- time_grid[j]

        h.m1 <- -as.numeric(dat.combine$R*(1-dat.combine$A)/
                              model_list$SR1.cens.pred[, j] *
                              model_list$S0R1.pred[, tk_idx]/
                              model_list$SR1.event.pred[, j]) * (t_j <= t_k) *
          q.R/D.X

        h.m2 <- -as.numeric((1 - dat.combine$R) * select_subset * q.R/
                              model_list$SR0.cens.pred[, j] *
                              model_list$S0R0.pred[, tk_idx]/
                              model_list$SR0.event.pred[, j]) * (t_j <= t_k)*
          r.X[, tk_idx]/D.X

        part1 <- h.m1 * (
          (dat.combine$Y == t_j) * (dat.combine$D == 1) -
            (dat.combine$Y >= t_j) * model_list$eventR1.hazard.hat[, j]
        )

        part2 <- h.m2 * (
          (dat.combine$Y == t_j) * (dat.combine$D == 1) -
            (dat.combine$Y >= t_j) * model_list$eventR0.hazard.hat[, j]
        )
        # interval length
        if(j==1){
          len <- t_j
        }else{
          len <- t_j - time_grid[j-1]
        }
        part1_sum <- part1_sum + part1 * len
        part2_sum <- part2_sum + part2 * len
      }

      (dat.combine$R * (1 - dat.combine$A) *
          (dat.combine$Y > t_k) * q.R/D.X/model_list$ScR1.pred[, tk_idx] +

          (1 - dat.combine$R) * select_subset * q.R *
          (dat.combine$Y > t_k) * r.X[, tk_idx]/D.X/model_list$ScR0.pred[, tk_idx] +

          {
            dat.combine$R - dat.combine$R * (1 - dat.combine$A) * q.R/D.X -
              (1  - dat.combine$R) * select_subset * q.R * r.X[, tk_idx]/D.X
          } * model_list$S0R1.pred[, tk_idx] +



          part1_sum + part2_sum
      )

    })


    kappa.R1 <- sapply(1:length(k_grid), function(tk_idx){
      t_k <- k_grid[tk_idx]
      part1_sum <- rep(0, length(ps.A))

      for (j in 1:length(time_grid)) {
        t_j <- time_grid[j]

        h.m1 <- -as.numeric(dat.combine$R*(1-dat.combine$A)/model_list$SR1.cens.pred[, j]/
                              (1 - ps.A)/prob.R *
                              model_list$S0R1.pred[, tk_idx]/
                              model_list$SR1.event.pred[, j]) * (t_j <= t_k)

        part1 <- h.m1 * (
          (dat.combine$Y == t_j) * (dat.combine$D == 1) -
            (dat.combine$Y >= t_j) * model_list$eventR1.hazard.hat[, j]
        )

        # interval length
        if(j==1){
          len <- t_j
        }else{
          len <- t_j - time_grid[j-1]
        }
        part1_sum <- part1_sum + part1 * len
      }

      (model_list$S0R1.pred[, tk_idx] +
          dat.combine$R * (1 - dat.combine$A)/(1 - ps.A)/prob.R *
          ((dat.combine$Y > t_k)/model_list$ScR1.pred[, tk_idx] -
             model_list$S0R1.pred[, tk_idx]) +
          part1_sum)
    })

    kappa.R0 <- sapply(1:length(k_grid), function(tk_idx){
      t_k <- k_grid[tk_idx]
      part2_sum <- rep(0, length(ps.A))

      for (j in 1:length(time_grid)) {
        t_j <- time_grid[j]

        h.m2 <- -as.numeric((1 - dat.combine$R)/(1-prob.R)/model_list$SR0.cens.pred[, j] *
                              model_list$S0R0.pred[, tk_idx]/
                              model_list$SR0.event.pred[, j]) * (t_j <= t_k)

        part2 <- h.m2 * (
          (dat.combine$Y == t_j) * (dat.combine$D == 1) -
            (dat.combine$Y >= t_j) * model_list$eventR0.hazard.hat[, j]
        )
        # interval length
        if(j==1){
          len <- t_j
        }else{
          len <- t_j - time_grid[j-1]
        }
        part2_sum <- part2_sum + part2 * len
      }

      model_list$S0R0.pred[, tk_idx] +
        (1 - dat.combine$R)/(1-prob.R) *
        ((dat.combine$Y > t_k)/model_list$ScR0.pred[, tk_idx] -
           model_list$S0R0.pred[, tk_idx]) +
        part2_sum
    })

    return(list(OI = model_list$S0R1.pred,# * dat.combine$R,
                OI.e = model_list$S0R0.pred,# * (1-dat.combine$R),
                IPCW.RCT = S0.hat.IPCW,
                AIPW.RCT = S0.hat.EE.mat.RCT,
                AIPW.EC = S0.hat.EE.mat.EC,
                ACW = S0.hat.EE.mat,
                ScR0.hat = model_list$ScR0.pred,
                bias.dr = kappa.R1 - kappa.R0,
                model_list = model_list))
  }

  S1.EE.hat.list <- aipcw.S1.est(dat = data_rt,
                                 fit.surv.obj = fit.event,
                                 fit.cens.obj = fit.cens,
                                 q.R = q.R, ps.A = ps.A,
                                 time_grid = time_grid,
                                 k_grid = k_grid)
  S1.EE.hat <- S1.EE.hat.list$AIPW.RCT
  S1.IPCW <- S1.EE.hat.list$IPCW.RCT
  S1.OI <- S1.EE.hat.list$OI

  S0.EE.hat.list <- aipcw.S0.est(dat.RCT = data_rt,
                                 dat.RWD = data_ec,
                                 fit.surv.obj = fit.event,
                                 fit.cens.obj = fit.cens,
                                 prob.R = prob.R, ps.A = ps.A,
                                 time_grid = time_grid,
                                 k_grid = k_grid)

  S0.EE.RCT.hat <- S0.EE.hat.list$AIPW.RCT
  S0.EE.EC.hat <- S0.EE.hat.list$AIPW.EC
  S0.EE.hat <- S0.EE.hat.list$ACW
  S0.OI <- S0.EE.hat.list$OI
  S0.IPCW <- S0.EE.hat.list$IPCW.RCT
  S0.OI.e <- S0.EE.hat.list$OI.e
  # selective borrowing
  bias.hat <- (t(apply(S0.OI, 1,
                       function(x)RMST.surv(x, k_grid)[, 'MST'])) -
                 t(apply(S0.OI.e, 1,
                         function(x)RMST.surv(x, k_grid)[, 'MST'])))
  ## should be conducted subject-wisely, not on time
  bias.hat_MST <- bias.hat[, length(k_grid)]


  # xi.hat
  ## DR-Learner
  xi.hat <- t(apply(S0.EE.hat.list$bias.dr, 1,
                    function(x)RMST.surv(x, k_grid)[, 'MST']))
  xi.hat_MST <- xi.hat[, length(k_grid)] # based on the integral from 0 to k_grid


  # compute the variance of xi_hat
  Z <- factor(
    c(rep(0, N_c),
      1:sum(dat.combine$R == 0))
  )
  Z.mat <- model.matrix(~Z+0)

  nu <- 1
  # constuct the adaptive penalty
  b_k.w_MST <- 1/abs(bias.hat_MST[dat.combine$R==0])**nu
  # the penalty factors are re-scaled with sum equal to `nvars`
  b_k.w_MST <- b_k.w_MST/sum(b_k.w_MST) * (ncol(Z.mat) - 1)
  # replace any NA number
  b_k.w_MST[is.na(b_k.w_MST)] <- Inf
  # add 0 for RCT (no penalty)
  b_k.w_MST <- c(0, b_k.w_MST)

  # adaptive lasso penalty
  fit.penalty <- glmnet(Z.mat,
                        xi.hat_MST,
                        family = 'gaussian',
                        alpha = 1,
                        penalty.factor = b_k.w_MST,
                        intercept = FALSE,
                        standardize = FALSE,
                        standardize.response = FALSE,
                        nlambda = 10,
                        lambda.min.ratio = 0.005
  )

  lambda.vector <- fit.penalty$lambda
  # SCAD penalty
  trade.off <- apply(predict(fit.penalty, newx = Z.mat, type = 'coef'), 2,
                     function(x)
                     {
                       select_subset <- c(rep(0, N_c),
                                          x[-c(1:2)] == 0)

                       xi.lin <- aipcw.S0.est(dat.RCT = data_rt,
                                              dat.RWD = data_ec,
                                              fit.surv.obj = fit.event,
                                              fit.cens.obj = fit.cens,
                                              prob.R = prob.R, # refit for calibration
                                              ps.A = ps.A,
                                              time_grid = time_grid,
                                              k_grid = k_grid,
                                              select_subset = select_subset,
                                              X.pred.R = X.pred.R,
                                              model_list = S0.EE.hat.list$model_list
                       )$ACW

                       bias.ec <- RMST.surv(apply(xi.lin, 2, sum)/N_c, k_grid)[length(k_grid), 'MST'] -
                         RMST.surv(apply(S0.EE.RCT.hat, 2, sum)/N_c, k_grid)[length(k_grid), 'MST']

                       var.ec <- var(apply(xi.lin, 1, function(xi.row){
                         RMST.surv(xi.row, k_grid)[length(k_grid), 'MST']
                       }))/N_c

                       sum(bias.ec**2) + var.ec
                       # var.ec
                     })
  lambda.opt.idx <- which.min(trade.off)
  lambda.opt <- lambda.vector[lambda.opt.idx]

  # pick the coefficient matrix
  coef.mat <- predict(fit.penalty, newx = Z.mat, type = 'coef')[, lambda.opt.idx]

  # select the comparable subset
  b_k_lasso <- coef.mat[-c(1:2)] # delete the intercept and the main study
  select_subset <- c(rep(0, N_c),
                     b_k_lasso == 0)
  # q.R_lasso <- q.R; q.R_lasso[b_k_lasso!=0] <- 0
  # data_ec_lasso <- data_ec; data_ec_lasso[b_k_lasso!=0, ] <- 0
  S0.EE.hat.list_lasso <- list(subset = as.vector(which(b_k_lasso == 0)),
                               S0.hat = aipcw.S0.est(dat.RCT = data_rt,
                                                     dat.RWD = data_ec,
                                                     fit.surv.obj = fit.event,
                                                     fit.cens.obj = fit.cens,
                                                     prob.R = prob.R,
                                                     ps.A = ps.A,
                                                     time_grid = time_grid,
                                                     k_grid = k_grid,
                                                     select_subset = select_subset,
                                                     X.pred.R = X.pred.R,
                                                     model_list = S0.EE.hat.list$model_list)$ACW)

  S0.EE.hat_lasso <- S0.EE.hat.list_lasso$S0.hat
  # combined estimation
  ## set R = 1 for the comparable external subset
  dat.update <- dat.combine
  dat.update[c(rep(0, N_c), b_k_lasso) == 0, 'R'] <- 1

  fit.event.update <- survival::coxph(as.formula(paste('Surv(Y,D == 1)~R*(',
                                             paste0(X.pred.f,
                                                    collapse = '+'),
                                             ')')),
                            data = dat.update)

  aipcw.RMST.est <- function(dat.RCT, dat.RWD,
                             S1.hat.X, S0.hat.X){

    N.R <- nrow(dat.RCT)
    N.E <- nrow(dat.RWD)

    S1.hat.X <- rbind(S1.hat.X,
                      matrix(rep(0, ncol(S1.hat.X) * N.E),
                             nrow = N.E))
    # S1.hat.EE.mat[abs(S1.hat.EE.mat)>5] <- 5
    # S0.hat.EE.mat[abs(S0.hat.EE.mat)>5] <- 5
    # average over the data
    if(is.null(dim(S1.hat.X))&
       is.null(dim(S0.hat.X))){
      S1.hat.EE <- S1.hat.X
      S0.hat.EE <- S0.hat.X

    }else{
      S1.hat.EE <- apply(S1.hat.X,
                         2, sum)/N.R
      S0.hat.EE <- apply(S0.hat.X,
                         2, sum)/N.R
    }


    # return a list
    return(list(est = RMST(c(S1.hat.EE),
                            c(S0.hat.EE),
                            k_grid = c(k_grid))))




  }



  aipcw.RCT.est <- aipcw.RMST.est(dat.RCT = data_rt,
                                  dat.RWD = data_ec,
                                  S1.hat.X = S1.EE.hat,
                                  S0.hat.X = S0.EE.RCT.hat)

  aipcw.est <- aipcw.RMST.est(dat.RCT = data_rt,
                              dat.RWD = data_ec,
                              S1.hat.X = S1.EE.hat,
                              S0.hat.X = S0.EE.hat)
  aipcw.est_lasso <- aipcw.RMST.est(dat.RCT = data_rt,
                                    dat.RWD = data_ec,
                                    S1.hat.X = S1.EE.hat,
                                    S0.hat.X = S0.EE.hat_lasso)

  return(list(aipcw.RCT = aipcw.RCT.est,
              aipcw = aipcw.est,
              aipcw_lasso = aipcw.est_lasso,
              subset.idx = S0.EE.hat.list_lasso$subset))
}
