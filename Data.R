library(ald) # for random draws / quantiles of ALD
library(gamlss) # for simulate data rfamily
library(MCMCpack) # for Dirichlet
library(Matrix) # for rankMatrix
library(expm) # for powers of matrices (stationariy)
library(dplyr)


#' @title Markov Switching (non-)linear Data generator
#' @description Dependent variable of Markov Switching ('y') is generated
#'   according to a state_formula list for each state, given a Dataframe Data,
#'   that contains the referenced covariates. The dependent is drawn from an
#'   Asymmetric Laplace Distribution or any gamlss family function according to
#'   the current state of the generated markov chain and to the parameters
#'   resulting from the stateregime applied to the covariates
#'   (=state-covariate-conditional distribution).
#'
#' @param covariate data.frame. Contains all covariates, which are referenced by
#'   state_formula. The provided length is also decisive for the number of
#'   generated observations.
#' @param state_formula list of lists. The Sublists contain named formula.
#'   Parameters are taken from gamlss families: mu, sigma, nu, tau. If not
#'   provided, they default to the rfamily values.
#' @param state vector. optional vector of a statechain. if not provided, a
#'   statechain is generated with delta & gamma arguments.
#' @param delta vector. Initial probabilities of the statechain. Semi-optional:
#'   required, if state is not supplied.
#' @param gamma Transition Probability matrix. To generate Statechain process.
#'   Semi-optional: required, if state is not supplied.
#' @param yfamily chr. Gamlss family name e.g. 'NO' from which y will be
#'   conditionally drawn. Also allows 'ALD'.
#'
#' @example
#' state_formula = list()
#' state_formula[[1]] = list(mu = ~x1, sigma = ~ 0.7, p = ~0.5)
#' state_formula[[2]] = list(mu = ~-1.5 * x1 - 10, sigma = ~ 0.3, p= ~0.5)
#'
#' MS = ms.y(
#'   covariate = data.frame(x1 = runif(2000, min = -2, max = 2)),
#'   state = NULL,
#'   delta = rep(0.5, 2),
#'   gamma = matrix(c(0.75, 0.25,
#'                    0.10, 0.90), byrow = T, ncol = 2),
#'   state_formula,
#'   yfamily = 'NO'
#' )
#'
#' # Or generate the markov chain with delta, gamma & state on the fly
#' set.seed(42)
#' r = rstateparam(s = 2)
#' MS = ms.y(
#'   covariate = data.frame(x1 = runif(2000, min = -2, max = 2)),
#'   state = state = stateprocess(n = 2000, r$delta, r$gamma),
#'   delta = NULL,
#'   gamma = NULL,
#'   state_formula,
#'   yfamily = 'NO'
#' )
#'
#' @return
#' data: data.frame: covriates augmented with y, conditional on the
#' state-dependent-parameters.
#' state: Markov chain; state process.
#' mcall: User call.
#' state_formula: user specified formula, that defines param, given the state.
#' gamma: transition matrix
#' delta: initial probability distribution
#' yfamily: distribution, of which y is drawn.
#' param: true parameters, from which each obs. y is drawn from (yfamily distrib.)
#' overlap: experimental measure. pairwise convolutions between
#'
#' @remark carefull state_formula may not be evaluated for a variable '~x' for
#'   some reason related to dplyr::transmute(data,
#'   eval(state_formula[[1]][[]][[2]]))) therefor: covariate = data.frame(x1 =
#'   rnorm(2000))
#' @export
MS.drawy = function(covariate = data.frame(x1 = rnorm(2000)),
                    state = NULL,
                    delta = rep(0.5, 2),
                    gamma = matrix(c(0.3, 0.7,
                                     0.6, 0.4), byrow = T, ncol = 2),
                    state_formula,
                    yfamily = 'NO') {

  mcall = match.call()

  # observations to draw
  n = nrow(covariate)

  # number of states
  if(is.null(state)){
    N = length(delta)
  } else {
    N = length(unique(state))
  }

  # (Error checking) --------------------------------------
  state_formula_error(state_formula, covariate)

  if (!is.null(delta) & !is.null(gamma)){
    stopifnot(nrow(gamma) == ncol(gamma))
    stopifnot(length(delta) == nrow(gamma))
  }

  # (Stateprocess & state-dependent parameters) -----------
  if (is.null(state)) {
    state = stateprocess(n, delta, gamma)
  }

  frames = extractStateframe(state_formula, covariate, state)
  param = frames$param # the state selected covariates/parameters


  # (Draw from conditional distribution) ------------------
  if (yfamily == 'ALD'){
    y = family_fun(family=yfamily,
                   func='r',
                   params=param,
                   n=1) # due to workaround
  } else {
    y = family_fun(family=yfamily,
                   func='r',
                   params=param,
                   n=n)
  }

  # (Experimental overlap measure on states) ------------------------
  # Density of y under both state
  conv = overlap(state_frame = frames$state_frame, yfamily=yfamily, n = length(y), N=nrow(gamma), y = y)


  return(
    list(
      data = cbind(y, covariate),
      state = state,
      mcall = mcall,
      state_formula = state_formula,
      gamma = gamma,
      delta = delta,
      yfamily = yfamily,
      param = param,
      overlap = conv
    )
  )
}


#' @title Family functions
#' @description This function is a convinience interface to handle the density,
#'   distribution, quantile and random generator function of a Gamlss family &
#'   ALD, conditional on a parameter frame.
#'
#' @param family family name e.g. NO, ALD, ...
#' @param func determines, which function is to be evaluated - d density, p
#'   distribution function, q quantile, r random generator.
#' @param params parameter frame containing the families' parameters
#' @param x see e.g. ?dnorm
#' @param q
#' @param p
#' @param n if NULL, for each observation, one is drawn - resulting in a vector
#'   of length==nrow(params). for n!=1, each observation (row) of params n obs.
#'   are drawn repeatedly.
#'
#' @remark func='n' in ALD case must be one, due to workaround caused by ALD
#'   funcs otherwise for each observations parameter vector n obs are evaluated.
#' @example
#' state_formula = list()
#' state_formula[[1]] = list(mu = ~0.5*x1 -0.5, sigma = ~ 0.1, p = ~0.5)
#' state_formula[[2]] = list(mu = ~ x1 + 0.5,   sigma = ~ 0.1, p= ~0.5)
#'
#'   MS = MS.drawy(
#'   covariate = data.frame(x1 = runif(2000, min = -2, max = 2)),
#'   state = NULL,
#'   delta = rep(0.5, 2),
#'   gamma = matrix(c(0.3, 0.7, 0.6, 0.4),byrow = T, ncol = 2),
#'   state_formula, yfamily = 'ALD' )
#'
#'   family_fun(family=MS$yfamily, func='d',params = MS$param, x=MS$data$y)
#'   family_fun(family=MS$yfamily, func='r',params = MS$param, n=1)
#' @return vector. Depending on the func argument.
#' @export
family_fun = function(family, func = c('d', 'p', 'q', 'r'), params,
                      p=NULL, q=NULL, x=NULL, n=NULL){

  func <- match.arg(func)
  fam_fun <- paste(func, family, sep = '')

  # Ensure correct argument pairs:
  if (any(!names(params) %in% names(formals(fam_fun)))) {
    stop("One of x,q,p,n,... arguments doesn't match with the distributional function
         (e.g. dNO, pNO, qNO, rNO). See the family's documentation for admissable arguments.")
  }

  if (family == 'ALD'){
    # workaround1, for ALD funcs, as vector valued parameter cause issues
    # workaround2, parameter p,q,x,n. Since 'p' is ALD skewness-parameter, and
    # the probab. argument 'p' is called 'prob' in qALD
    frame = list(p,q,x,n)
    frame = frame[!sapply(frame, is.null)]
    params$vec = frame[[1]]

    quantity = apply(
      params,
      MARGIN = 1,
      FUN = function(x)
        do.call(fam_fun,
                list(
                  x[4], # frame workaround
                  mu = x[1],
                  sigma  = x[2],
                  p = x[3]
                ))
    )

  } else { # gamlss case
    param <- list(
      mu = params$mu,
      sigma = params$sigma,
      nu = params$nu,
      tau = params$tau,
      x = x,
      q = q,
      p = p,
      n = n
    )

    # Kick out NULL parameters:
    param <- param[!sapply(param, is.null)]
    quantity = do.call(fam_fun, param)
  }

  return(quantity)
}




#(helper functions) -----------------------------------------------------------


#' @description generate parameterframe for each state from evaluated parameter
#'   formula with respect to the actual state of that observation. That is
#'   generate the parameters of the state-covariate-conditional distribution.
#' @param state_formula list with sublists. sublists containing the parameter
#'   formulas for one state each.
#' @param data data.frame supplying the covariates on which the parameter
#'   formulas are evaluated
#' @return list containing the set of state-dependent parameters derived from
#'   the state_formulas.
extractStateframe = function(state_formula, data, state) {

  state_frame = list()
  # iterate over all states to create parameterframes irrespective of statechain
  for (i in 1:length(state_formula)){
    state_frame[[i]] = as.data.frame(
      lapply(
        state_formula[[i]],
        FUN = function(x) dplyr::transmute(data, eval(x[[2]]))
    ))
    names(state_frame[[i]]) = names(state_formula[[1]])
  }

  # Select State-correct parameter frame
  param = data.frame(
    matrix(NA,  nrow = length(state),ncol = length(state_formula[[1]])))
  colnames(param) = names(state_frame[[1]])

  for (i in 1:length(state)) {
    param[i, ] = state_frame[[state[i]]][i, ]
  }

  return(list(param=param, state_frame=state_frame))
}





#' @description User error in state_param specification
state_formula_error = function(state_formula, covariate) {

  # Userinput for Data is sufficient to supply state_formula
  varnames = list()
  for (i in 1:length(state_formula)) {
    varnames[[i]] =  unique(unlist(sapply(state_formula[[i]], FUN = all.vars)))
  }
  varnames = unique(unlist(varnames))

  if (!any(varnames %in% names(covariate))) {
    stop('names(covariate) does not sufficie the used parameters in state_formula')
  }

  for (i in 1:length(state_formula)) {
    if (any(unlist(lapply(state_formula[[i]], FUN = class) != 'formula'))) {
      stop(
        'Check state_formula. Any of the equation systems is not of formula class.
        did you forget "~"? specify formulas e.g.: mu = ~exp(x1)'
      )
    }
  }
}





#' @description helper function for MS.drawy (experimental feature), estimates
#'   the state-pairwise convolution of the true densities.
#' @param state_frame list of the true parameter frames irrespectively of state
overlap <- function(state_frame, yfamily, n, N, y) {
  densities = data.frame(matrix(FALSE, nrow=n, ncol=N))
  for (s in 1:N){
    densities[,s] = family_fun(
      family=yfamily,
      func='d',
      params=state_frame[[s]],
      x=y)
  }

  # how much do the distributions overlap locally? each row of 'densities' df
  # describes the respective density value of an observation y evaluated under
  # the potential distribution regime of that state. (eventhough only one is
  # true) The product of the densities indicates the level of overlap (in
  # unimodal distributions)

  conv = data.frame(matrix(nrow = n))[-1]
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      conv[paste(i, 'with', j, sep = '')] = densities[, i] * densities[, j]
    }
  }
  return(conv)
}





# (markov chain simulation functions) ------------------------------------

#' @title alternative Markov chain generator
#'
#' @param k number of states
#' @param n length of the Markov chain
#' @param p power of TPM
#'
#' @example MC(4,10, 100)
MC = function(k, n, p=1){
  gamma = MCMCpack::rdirichlet(k, rep(1,k))
  delta = MCMCpack::rdirichlet(1, rep(1,k))

  TPM = gamma %^% p

  print('gamma')
  print(gamma)
  print('TPM')
  print(TPM)

  y = vector(length=n)
  y[1] = sample(seq(1,k),size=1, prob=delta)

  for (i in seq(2,n)){
    y[i] = sample(seq(1,k), size=1, prob=gamma[y[i-1],])
  }

  return(list(
    y=y,
    gamma=gamma,
    delta=delta
    )
  )
}




#' @title Markov state chain parameter generator
#' @description ensure that delta and each row are a probability distribution
#'   (as they are drawn from Dirichlet), but the properties of the Markov chain
#'   need not be fullfilled. CHECK ERDOGEDICITY!
#' @param s Integer. Number of states
#' @export
rstateparam <- function(s) {
  delta = c(MCMCpack::rdirichlet(1, rep(1, s)))
  gamma = matrix(NA, ncol = s, nrow = s)
  for(i in 1:nrow(gamma)){
    gamma[i,] = MCMCpack::rdirichlet(1, rep(1, s))
  }

  # check for stationary distribution
  # if (stationarity_flag(gamma)){
  #   gamma = rstateparam(s)$gamma # potential source of error??
  # }

  return(list(delta = delta, gamma = gamma))
}



#' @title Markov State Chain
#' @param delta initial probability vector
#' @param gamma transitional probability matrix
#' @param n length of chain
#' @export
stateprocess = function(n, delta, gamma){
  state = rep(NA, n)
  state[1] = sample(1:length(delta), 1, prob = delta)

  # if (stationarity_flag(gamma)){
  #   print('gamma has no stationary distribution')
  # }

  for (s in 2:length(state)) {
    # next state prob conditional on current state
    state[s] = sample(1:length(delta), 1, prob = gamma[state[s - 1], ])
  }

  # descriptive info
  for (i in sort(unique(state))) {
    cat('State', i , 'has', length(state[state == i]), 'observations\n')
  }
  return(state)
}


# (depreciated) ----------------------------------------------------------------

#' @title Markov-Switching-Process' Quantiles
#' @remark This function is depreciated, look at family_fun
#' @description given a Markov Switching process, this function allows to
#'   evaluate the correct quantiles of the state-dependent distribution. That
#'   is, given the markov regime, the distribution from which the current
#'   observation is drawn, can be evaluated at its quantiles. It thereby provides
#'   the expected value of a quantile regression. This may be of use in case one
#'   wants to evaluate the true realized bias.
#'
#' @param MS object resulted from ms.drawy call.
#' @param quantil quantiles at which the function is to be evaluated at
#'   MS.quantiles(MS, c(0.25, 0.5, 0.75))
MS.quantiles = function(MS, quantil){
  state = MS$state
  param = MS$param
  yfamily= MS$yfamily

  qfam <- paste('q', yfamily, sep = '')

  if (yfamily == 'ALD'){
    quantiles = t(apply(
      param,
      MARGIN = 1,
      FUN = function(x)
        do.call(ald::qALD,
                list(
                  quantil,
                  mu = x[1],
                  sigma  = x[2],
                  p = x[3]
                ))
    ))

    colnames(quantiles) = paste('q', quantil[1:length(quantil)], sep = '')

  } else { # gamlss case
    # estimate the corresponding quantiles of conditional family distribution
    # note the workaround, as do.call with vector valued p & mu fails. e.g.
    # qNO(p = c(0.25, 0.75), mu = 1:3, sigma = 1)
    quantiles = data.frame(matrix(nrow = length(state)))
    quantiles = quantiles[, FALSE]

    param = param[names(param) != 'n']
    for (i in quantil) {
      param$p = i
      quantiles[paste('q', i, sep = '')] = do.call(qfam, param)
    }
  }
  return(quantiles)
}


#' @description Check if Markov Chain has stationary distribution
#' @param gamma transition matrix.
#'
#' @remark needs remastering and decision upon shipping at all what are
#'   prequisites to check, whether or not gamma has a stationary distribution?
#'   check ergodecity
#'  https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
stationarity_flag = function(gamma){
  if (nrow(gamma) != ncol(gamma)){
    stop('gamma is not square')
  }

  if (!all.equal(apply(gamma, MARGIN = 1, sum),  rep(1,nrow(gamma)), tolerance=1e-3)){
    stop('stop gamma\'s rows are not probability distributions')
  }

  k = nrow(gamma)

  ones = rep(1, k)
  I = diag(rep(1, k))
  U = matrix(1, nrow=k, ncol=k)

  #stationary distribution
  PI = ones %*% solve(I - gamma + U)


  # -----------------------------
  # TODO stationarity_flag produces imagenary number in this approach
  # # Eigendecomposition approach:
  # # Get the eigenvectors of gamma, note: R returns right eigenvectors
  # r = eigen(gamma) # can produce imagenary numbers
  # rvec = r$vectors
  #
  # # left eigenvectors are the inverse of the right eigenvectors
  # lvec = ginv(r$vectors)
  #
  # # # The eigenvalues
  # # eigval<-r$values
  # # # second way of checking the spectral decomposition:
  # # rvec%*%diag(eigval)%*%ginv(rvec)
  #
  # pi_eig = Re(lvec[1,]/sum(lvec[1,]))


  # # second approach (failing)
  # s = nrow(gamma)
  # A = diag(s)-gamma
  # A = rbind(A, rep(1,s)) # linear constraint sum(pi)= 1 attached
  # b = c(rep(0, s), 1)
  # rank.A = as.numeric(Matrix::rankMatrix(A))
  # rank.Ab = as.numeric(Matrix::rankMatrix(cbind(A,b)))
  #
  # if(rank.A == rank.Ab){
  #   print('proposal is consistent')
  #   # approx. stationary distribution pi: solve normal equations
  #   pis = drop(solve(t(A)%*%A, t(A)%*% b))
  #   cat('pi' , pis , '\n')

  # third approach: Brute Force
  # # manual version of stationary distribution
  # g = gamma%^%40
  # cat('Check stationarity: Gamma to the power of 40:\n' )
  # print(g)
  # print(all.equal(g, matrix(rep(pis,s),byrow=T, nrow=s), tolerance=10e-4))
}

