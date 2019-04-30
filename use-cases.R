rm(list = ls())
source('./MSquantreq.R')
source('./Data.R')
source('./FitEval.R')


# (multistate clean) -----------------------------------------------------------
# note that stateformula also specifies number of states
state_formula = list()
state_formula[[1]] = list(mu = ~ 0.5 * x1 - 0.5, sigma = ~ 0.1)
state_formula[[2]] = list(mu = ~ x1 + 0.5, sigma = ~ 0.1)
state_formula[[3]] = list(mu = ~ 0.4*x1 + 3, sigma = ~ 0.1)

set.seed(3)
fitevals = fiteval(
  data = data.frame(x1 = runif(200, min = -2, max = 2)),
  state_formula,
  tau = c(0.3, 0.5),
  n_obs = c(100, 200),
  repetitions = 10,
  delta = NULL,
  gamma = NULL,
  fixstatechain = FALSE
)
#saveRDS(fitevals, "fiteval_3states.rds")
#fitevals = readRDS("fiteval_3states.rds")
#saveRDS(fitevals, "fiteval_3states_pred_config3.rds")
#fitevals = readRDS("fiteval_3states_pred_config3.rds")

names(fitevals)
names(fitevals$MSs)
names(fitevals$fits)
names(fitevals$QR_evals)
names(fitevals$EM_MC_evals)

# 'interface' for all following plots
n_obs = 200
tau = 0.3
repetition = 9

MSobs = paste('n_', n_obs, sep = '')
fitconfig = paste('n_', n_obs, '_tau_', tau, sep = '')
reps = paste('repetiton_', repetition, sep = '')

# (msquantreg realizations) --------------
# note if QR_eval(pred_config=2), (==rowwise which.max state_prob is model's state
# prediction), the lines spike. The reason is that the state is misclassified -
# as consequence, the bias is calculated with the model's prediction (with faulty state)
# against true quantile (true state). However, the lines connect by true state.
# More or less, the same applies to (pred_config=3), where the prediction is
# the weighted average (according to states' prob) of all three model predictions

# bias = E(T)- g(theta)
# my measure: sum[w(distance(x1_1, x1_2))*(T-g(theta))],
# where w is dnorm normalized for w to sum to 1.
# is red line in x1-biasvec plot

# MSE =  E[(T- g(theta))²]
# my measure: sum[w(distance(x1_1, x1_2))*(T-g(theta))²]

# note penalty lambda is gam default and not cross validated yet!

# CONSIDER DIFFERENT PRED_CONFIGS!
plotbiasrealizations(
  stackedbiasframe=fitevals$QR_evals[[fitconfig]][['stackedbiasframe']])

plotMCrealizations(
  stackedgammaforb=fitevals$EM_MC_evals[[fitconfig]][['stackedgammaforb']],
  stackeddeltaeucl=fitevals$EM_MC_evals[[fitconfig]][['stackeddeltaeucl']])



# (single realizations) -------------
MS.plot1realization(
  MS= fitevals$MSs[[MSobs]][[repetition]],
  fit=fitevals$fits[[fitconfig]][[repetition]],
  fiteval=fitevals$QR_evals[[fitconfig]][[repetition]],
  stateind = 1)

MS.plot1realization(
  MS= fitevals$MSs[[MSobs]][[repetition]],
  fit=fitevals$fits[[fitconfig]][[repetition]],
  fiteval=fitevals$QR_evals[[fitconfig]][[repetition]],
  stateind = 2)

MS.plot1realization(
  MS= fitevals$MSs[[MSobs]][[repetition]],
  fit=fitevals$fits[[fitconfig]][[repetition]],
  fiteval=fitevals$QR_evals[[fitconfig]][[repetition]],
  stateind = 3)

# ----------
fit.plot1resid(
  #TODO shape: model_prediction of state?
  MS= fitevals$MSs[[MSobs]][[reps]],
  fit=fitevals$fits[[fitconfig]][[reps]],
  fiteval=fitevals$QR_evals[[fitconfig]][[reps]],
  stateind = 1)



# (ADAM example) ---------------------------------------------------------------
state_formula = list()
state_formula[[1]] = list(mu = ~ 0.5 * x1 - 0.5, sigma = ~ 0.1)
state_formula[[2]] = list(mu = ~ x1 + 0.5, sigma = ~ 0.1)

n = 200
set.seed(3)
fitevals = fiteval(
  data = data.frame(x1 = runif(n, min = -2, max = 2)),
  state_formula,
  tau = c(0.25, 0.5, 0.75),
  n_obs = c(2000),
  repetitions = 10,
  delta = rep(0.5, 2),
  gamma = matrix(c(0.95, 0.05, 0.05, 0.95), ncol = 2),
  fixstatechain = FALSE
)

#saveRDS(fitevals, "fiteval_3ADAM.rds")
#fitevals = readRDS("fiteval_3ADAM.rds")

n_obs = 2000
tau = 0.25
repetition = 3

MSobs = paste('n_', n_obs, sep = '')
fitconfig = paste('n_', n_obs, '_tau_', tau, sep = '')
reps = paste('repetiton_', repetition, sep = '')

ggplot(aes(x=x1, y=biasvec), data = fitevals$QR_evals[[fitconfig]][['stackedbiasframe']])+
  facet_grid(~state)+
  geom_point(aes(color=realization, alpha = 0.01))+
  geom_line(aes(color=realization))+
  stat_smooth(formula = y~x, method = loess)+
  geom_hline(yintercept = 0)

plotMCrealizations(
  stackedgammaforb=fitevals$EM_MC_evals[[fitconfig]][['stackedgammaforb']],
  stackeddeltaeucl=fitevals$EM_MC_evals[[fitconfig]][['stackeddeltaeucl']])

MS.plot1realization(
  MS= fitevals$MSs[[MSobs]][[repetition]],
  fit=fitevals$fits[[fitconfig]][[repetition]],
  fiteval=fitevals$QR_evals[[fitconfig]][[repetition]],
  stateind = 1)




# (serious example) ------------------------------------------------------------
# algorithm might fail in some repetition due to singularity in Msquanreg's qr
n =  2000
fitevals = fiteval(
  data = data.frame(x1 = runif(n, min = -2, max = 2)),
  state_formula,
  tau = c(0.5),
  n_obs = c(200),
  repetitions = 100,
  delta = NULL,
  gamma = NULL,
  fixstatechain = FALSE
)


# (non linear/ sparse data generator) ------------------------------------------
state_formula = list()
state_formula[[1]] = list(mu = ~ - x1, sigma = ~ 0.3) # for sigma positive
state_formula[[2]] = list(mu = ~ sin(x1)*5 -10, sigma = ~ abs(x2))

delta = rep(0.5, 2)
gamma = matrix(c(0.95, 0.05, 0.05, 0.95), ncol = 2)

d = data.frame(x1 = rnorm(200, sd =2 ))
d$x2 = d$x1*2
d$x2 = runif(200, max= 1.5, min = -1)

MS = MS.drawy(
  covariate = d,
  state = NULL,
  delta = delta,
  gamma = gamma,
  state_formula,
  yfamily = 'NO'
)

MS$param
ggplot(aes(x = x1, y = y, color=MS$state), data = MS$data) + geom_point()

# (yfamily) --------------------------------------------------------------------
state_formula = list()
state_formula[[1]] = list(mu = ~ x1, sigma = ~ 0.1) # for sigma positive
state_formula[[2]] = list(mu = ~ -x1*2 + 2, sigma = ~ 0.1)

d = data.frame(x1 = runif(200, min = 0 , max = 1))
MS= MS.drawy(
  covariate = d,
  state = NULL,
  delta = delta,
  gamma = gamma,
  state_formula,
  yfamily = 'GA'
)


ggplot(aes(x = x1, y = y, color=MS$state), data = MS$data) + geom_point()


# OVERLAP

# (CONFUSION)-------------------------------------------------------------------
fitevals$EM_MC_evals$n_200_tau_0.5$repetiton_1$confusion

# still issues with inf/nan with crossentropy
fitevals$EM_MC_evals$n_200_tau_0.5$repetiton_1$crossentropy_avg

# (fit on state_formula with different covariates) -----------------------------
state_formula = list()
state_formula[[1]] = list(mu = ~ exp(x1)- x1, sigma = ~ 0.3) # for sigma positive
state_formula[[2]] = list(mu = ~ sin(x2)*5 -10, sigma = ~ x2^2)

delta = rep(0.5, 2)
gamma = matrix(c(0.95, 0.05, 0.05, 0.95), ncol = 2)

d = data.frame(x1 = rnorm(200, sd =2 ))
d$x2 = runif(200, max= 1.5, min = -1)

fitevals = fiteval(
  data= d,
  state_formula,
  tau = c(0.5),
  n_obs = c(200),
  repetitions = 100,
  delta = delta,
  gamma = gamma,
  fixstatechain = FALSE
)

# ISSUE general family: link functions are not covered
# ISSUE is that plot methods all reference x1!
n_obs = 200
tau = 0.5
repetition = 1

MSobs = paste('n_', n_obs, sep = '')
fitconfig = paste('n_', n_obs, '_tau_', tau, sep = '')
reps = paste('repetiton_', repetition, sep = '')
MS.plot1realization(
  MS= fitevals$MSs[[MSobs]][[repetition]],
  fit=fitevals$fits[[fitconfig]][[repetition]],
  fiteval=fitevals$QR_evals[[fitconfig]][[repetition]],
  stateind = 1)









