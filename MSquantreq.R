library(gamlss)
library(quantreg)
library(ald)

#'@title Fit Markov-Switching with EM-(linear-)Quantile-Regression
#'@description This function estimates iteratively the underlying Markov-Chain
#'  via Expectation Maximization wrapping the weighted Quantile Regressions on
#'  each state. Each state is estimated seperately at quantile tau, weighting
#'  each observation by the probability to pertain to the state under
#'  evaluation.
#'@license This algorithm is a refactored version of code that was provided by Timo Adams.
#' only slight alterations such as comments, plots and formula interaction have been added.
#'@remark The extension to non-linear quantile regression fails, as weighting is
#'  not implemented in the RQSS package correctly.
#'@param N number of states
#'@param stat assumption of stationary Markov chain (has implication for initial
#'  state probability delta)
#'@param data dataframe containing y and x which is refered to in formula
#'@param formula list of length N, containing the state depending formula
#' which is used for rq fit
#'@param tau quantile that is estimated
#'@param conv_tol numerical. Specify the conversion precision in llh-llh_old.
#'@param conv_print logical. While iterating print the current llh and conv.crit.
#'
#'@return mod: list of N models, that can be used with
#'  fitted(mod[i]) to predict y at the tau quantile with
#'  their most likely state value, if the probs are close to 0 or 1
#'
#'@export
MS.fitquantreg = function(data,
                          formula = list(y ~ x1, y ~ x1),
                          N = 2,
                          stat = FALSE,
                          tau = 0.5,
                          max_iter = 250,
                          conv_tol = 1e-08,
                          conv_print = TRUE,
                          plotit = FALSE) {

  mcall = match.call()

  # user specified dependent
  y = data[,as.character(formula[[1]][[2]])]

  # (initial parameters) --------------------------------------------
  # FIXME: initial parameter guesses scalable to multiple statess N

  delta = NULL # rstateparam(N)$delta
  gamma = rstateparam(N)$gamma # matrix(c(0.95, 0.05, 0.05, 0.95), ncol = 2)
  mod = fv_mu = fv_sigma = vector("list")
  old_llh = 0

  for (j in 1:N) {   # adjusted to accomodate n >= 3
    fv_mu[[j]] = seq(
      from = seq(-1, 1, length = N)[j], #  c(-0.25, 0.75)[j],
      to = seq(0, 2, length= N)[j], # c(0, 1.5)[j],
      length = length(y)
    )
    fv_sigma[[j]] = rep(0.1, N)[j] # c(0.1, 0.1)[j]

  }
  while (TRUE) {
    for (i in 1:max_iter) {
      # (a) Initalize/Update parameters --------------------------------

      delta_next = delta
      gamma_next = gamma
      fv_mu_next = fv_mu
      fv_sigma_next = fv_sigma

      # stationarity assumption on MC
      if (is.null(delta)) {
        delta = solve(t(diag(N) - gamma + 1), rep(1, N))
      }

      # (b) Denisty value of y under current parameters ----------------
      lalpha = lbeta = matrix(NA, N, length(y))
      state_density = matrix(NA, nrow = length(y), ncol =  N)
      for (j in 1:N) {
        for (k in 1:length(y)) { # can be apply
          state_density[k, j] = dALD(y[k],
                                     mu = fv_mu[[j]][k],
                                     sigma = fv_sigma[[j]],
                                     p = tau)
        }
      }
      state_density = ifelse(!is.na(state_density), state_density, 1)

      # (c.1) forward algorithm ----------------------------------------
      # alpha_1 = P(State = i)*P(y) (ZML p38+p49)
      foo = delta * state_density[1,]

      # Avoid underflow LI by scaling alpha vector by means of ZML p48/p49.
      sumfoo = sum(foo)
      lscale = log(sumfoo)
      foo = foo / sumfoo
      lalpha[, 1] = log(foo) + lscale

      # P(y_t| State = j, gamma, density estimate)
      for (j in 2:length(y)) {
        foo = foo %*% gamma * state_density[j,]

        # scaling & approx LI
        sumfoo = sum(foo)
        lscale = lscale + log(sumfoo)
        foo = foo / sumfoo
        lalpha[, j] = log(foo) + lscale
      }

      # (c.2) backward algorithm (see p. 67 in ZML) --------------------
      foo = rep(1 / N, N)
      lbeta[, length(y)] = rep(0, N)
      lscale = log(N)
      foo = foo / sum(foo)

      for (j in (length(y) - 1):1) {
        foo = gamma %*% (state_density[j + 1,] * foo)
        lbeta[, j] = log(foo) + lscale
        sumfoo = sum(foo)
        foo = foo / sumfoo
        lscale = lscale + log(sumfoo)
      }

      # (c.3) loglikelihood -------------------------------------------

      # This is a workaround based on zucchini. Allows short coding as alpha is scaled.
      # originally sum of alphas would have resulted in logLI
      llh = max(lalpha[, length(y)]) +
        log(sum(exp(lalpha[, length(y)] - max(lalpha[, length(y)]))))


      # (d) Stateindicator --------------------------------------------
      # E(weights) = E(P(State_t = i| y, x)) in log for all i
      weights = exp(lalpha + lbeta - llh)
      # originally in code :
      # weights = matrix(NA, nrow =  N, ncol =  length(y))
      # for (j in 1:length(y)) {
      #   weights[, j] = exp(lalpha[, j] + lbeta[, j] - llh)
      # }

      # look at the weight's components
      # print(qplot(y = lalpha[1,]))
      # print(qplot(y = lbeta[1,]))
      # print(qplot(y = lalpha[1,] + lbeta[1,]- llh))
      # print(-llh)

      # (e.1) Gamma  (transitional prob. estimates) ----------------------------
      lstate_density = log(state_density)

      # how 'likely' is an observation in the respective state?
      # print(qplot(y = state_density[,1]))
      # print(qplot(y = state_density[,2]))

      # p72 ZML second equation 2. uses indicator (v_ij) of last state implicitly
      for (j in 1:N) {
        for (k in 1:N) {
          gamma_next[j, k] = gamma[j, k] *
            sum(exp(lalpha[j, 1:(length(y) - 1)] +
                      lstate_density[2:length(y), k] +
                      lbeta[k, 2:length(y)] -
                      llh))
        }
      }
      gamma_next = gamma_next / apply(gamma_next, 1, sum)

      # (e.2) Delta (initial state distrib) ------------------------------------
      if (stat == TRUE) { # stationary assumption MC
        delta_next = solve(t(diag(N) - gamma_next + 1), rep(1, N))
      } else {
        # current (ith step) initial state probability
        # accord to paper eq 6 & eq 8 after transformation. sum(u_i) = llh
        delta_next = exp(lalpha[, 1] + lbeta[, 1] - llh)
        delta_next = delta_next / sum(delta_next)
      }

      # (e.3) Model fit --------------------------------------------------------
      for (j in 1:N) {
        data$w = weights[j, ]

        if(plotit==TRUE){
          # plot the current realization of weights: indicates algorithm's guess
          print(qplot(y = weights[j,], color = weights[j,]))
        }

        # (RQ-Version)
        mod[[j]] = quantreg::rq(
          formula = formula[[j]],
          tau = tau,
          weights = w,
          data = data
        )

        # #(RQSS-Version) # weighted estimation is not working! even with RQSSkoenker
        # # requires source('RQSS_koenker'). The weight estimation fails considerably.
        # mod[[j]] = RQSS(formula = formula,
        #                 tau = tau,
        #                 weights = weights[j,],
        #                 data = data ,
        #                 lambda = 10 # als global argument und nicht term specific übergeben?
        # )

        # MLEs for parameter mu of the ALD
        fv_mu_next[[j]] = fitted(mod[[j]])

        # MLEs for parameter sigma of the ALD
        diff = y - fv_mu_next[[j]]
        fv_sigma_next[[j]] = (sum(weights[j, which(diff > 0)] * tau * diff[diff > 0]) -
                                sum(weights[j, which(diff < 0)] * (1 - tau) * diff[diff < 0])) / length(y)
        # for stability reasons
        if (fv_sigma_next[[j]] < 0.00001) {
          fv_sigma_next[[j]] = 0.00001
        }
      }

      # look at the actual fit
      if(plotit == TRUE){
        print(ggplot(aes(x1), data=data)+
                        geom_point(aes(x1,y))+
                        geom_line(aes(x1,y= fitted(mod[[1]])), color= 'lightblue')+
                        geom_line(aes(x1,y= fitted(mod[[2]])), color= 'steelblue')
                )
      }


      # (f.1) Iteration/Convergence counter --------------------------
      conv_crit = abs(llh - old_llh)
      if (conv_print == TRUE) { # not a perfect fit (otherwise it would have converged perfectly)
        print(paste("--- ite. =",i,
                    "- con.-cri. =", round(conv_crit, 6),
                    "- log.-llh. =", round(max(llh), 6),"---"
        ))
      }

      # Convergence was met or max_iter reached?
      if (conv_crit < conv_tol | i == max_iter) {
        if (i == max_iter) {
          print(paste("--- no convergence within", max_iter, "ite. ---"))
        } else{
          print(paste("--- convergence after", i, "ite. ---"))
        }

        return(list(
          delta = delta_next,
          gamma = gamma_next,
          mod = mod,
          llh = llh,
          state_probs = t(weights), # TODO return transpose and change dependencies in fitEVAL
          mcall = mcall,
          tau = tau
        ))

      }

      # (f.2) Update initial parameters with current estimates -----------------
      delta = delta_next
      gamma = gamma_next
      fv_mu = fv_mu_next
      fv_sigma = fv_sigma_next
      old_llh = llh

    }
  }
}


