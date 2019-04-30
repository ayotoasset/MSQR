library('ggplot2')
library('gridExtra')
library('fastDummies')
library('caret')
library('reshape2')
library('data.table')
library('mgcv')

#' @description This function is a wrapper around simulation, fit &
#'   fitevaluation, that evaluates MS.fitquantreg under different data and fit
#'   configurations.
#'
#' @param state_formula
#' @param tau vector of quantiles at which the fit is to be evaluated.
#' @param n_obs vector of numbers of observations for which the model is to be
#'   evaluated
#' @param repetitions number of repeated draws of the data (with n_obs
#'   observations configuration), to estimate bias and MSE over all repetitions.
#' @param delta vector. initial state distribution. Used to generate all statechains.
#' @param gamma transition probability matrix. Used to generate all statechains.
#' @param fixstatechain boolean, indicating whether or not to use the exact same
#'   state chain for all repetitions (of same n_obs length). If TRUE, one state
#'   chain per n_obs is generated and fix for all data simulation of same
#'   length. If FALSE, the state_chain is generated from delta & gamma for each
#'   repetition.
#'   Depending on delta and gamma to be NULL or not, the state chain
#'   is generated from the provided arguments or are drawn from rdirichlet. Note
#'   that in the latter case, the parameters are drawn only once and used to
#'   generate all statechains irrespectively of fixstatechain!
#'
#' @return output format of MSs, fits, EM_MC_evals, QR_evals: for each unique
#'   configuration one list, containig all realizations of this configuration.
#'   further, each realization carries its own parameters and return values
fiteval = function(state_formula,
                   data = data.frame(x1 = runif(n, min = -2, max = 2)),
                   tau =c(0.3,0.5),
                   n_obs = c(50, 100),
                   repetitions = 10,
                   delta = NULL,
                   gamma = NULL,
                   fixstatechain = TRUE) {


  # (generate statechain) ------------------
  # number of states
  N = length(state_formula)

  # generate MC parameters if not available
  if (is.null(delta) && is.null(gamma)){
    r = rstateparam(N)
    delta = r$delta
    gamma = r$gamma
  }

  # generate state and fix it for all repetitions
  if(fixstatechain){
    # bypassing MS.drawy's stateprocess generator (if fixstatechain is true!)
    state = list()
    for (n in n_obs){
      print(sprintf('chain with %s observations', n))
      state[[paste('n',n, sep='')]] = stateprocess(n, delta, gamma)
    }
  } else {
    # make use of MS.drawy stateprocess generator at each MS.drawy call
    state = NULL
  }

  # (generate data based on statechain)  -------------------------
  # nested structure: $config$repetitions
  MSs = list()
  for (n in n_obs) {
    MSs[[paste('n_', n, sep ='')]] = list()
    for (repet in 1:repetitions) {
      varname = paste('repetiton_', repet, sep = '')

      MSs[[paste('n_', n, sep='')]][[varname]] = MS.drawy(
        covariate = data,
        state = state[[paste('n', n, sep = '')]],
        delta = delta,
        gamma = gamma,
        state_formula,
        yfamily = 'NO'
      )
    }
  }

  # (fit on each of MSs and eval this fit) ----------------------------
  fits = list()
  EM_MC_evals = list()
  QR_evals = list()

  # deparse from state_chain, on what is to be fitted
  # TODO place this deparsing in MSfitquantreg!
  formulas = lapply(1:N,
    FUN=function(i) formula(
      paste('y',
            paste(all.vars(state_formula[[i]][[1]]),sep='+'),
            sep='~')
      )
  )
  # handling formulas as strings: deparse(formula)

  for (n in n_obs) {
    n_config = paste('n_', n, sep ='')
    for (t in tau) {
      # TODO refactor: ntau_config at the corresponding positions.
      ntau_config = paste('n_', n, '_tau_', t,sep='')


      fits[[paste('n_', n, '_tau_', t,sep='')]] =
        lapply(MSs[[n_config]],
               FUN = function(m) MS.fitquantreg(
                 data = m$data,
                 formula = formulas,
                 # FIXME what if they are changed or: general family: link functions?
                 N = length(state_formula),
                 stat = FALSE,
                 tau = t
               ))

      EM_MC_evals[[paste('n_', n, '_tau_', t,sep='')]] =
        mapply(FUN=function(m,f) EM_MC_eval(m,f),
               MSs[[n_config]],
               fits[[paste('n_', n, '_tau_', t,sep='')]],
               SIMPLIFY=FALSE)

      QR_evals[[paste('n_', n, '_tau_', t,sep='')]] =
        mapply(FUN=function(m,f) QR_eval(m,f),
               MSs[[n_config]],
               fits[[paste('n_', n, '_tau_', t,sep='')]],
               SIMPLIFY=FALSE)

      # #data access for some n config (developer tool)
      # lapply(MSs[[n_config]], names)
      # lapply(MSs[[n_config]], function(m) m$state)
      # lapply(fits[[paste('n_', n, '_tau_', t,sep='')]], names)
      # lapply(EM_MC_evals[[paste('n_', n, '_tau_', t,sep='')]], names)
      # lapply(QR_evals[[paste('n_', n, '_tau_', t,sep='')]], names)



      # (for plotting biases) ----------------------------------

      # biasvecs of all repetitions stacked (for each config)
      stackedbiasvec =
        unlist(
          lapply(
            QR_evals[[paste('n_', n, '_tau_', t,sep='')]],
            FUN= function(repet) repet[['biasvec']]),
          use.names=FALSE
        )

      stackedtruequantile =
        unlist(
          lapply(
            QR_evals[[paste('n_', n, '_tau_', t,sep='')]],
            FUN= function(repet) repet[['truequantile']]),
          use.names=FALSE
        )

      stackedframe = as.data.frame(data.table::rbindlist(
        lapply(MSs[[n_config]], FUN = function(repet) repet$data)))

      stackedstate = factor(
        unlist(
          lapply(
            MSs[[paste('n_', n,sep='')]],
            FUN= function(repet) repet[['state']]),
          use.names=FALSE),
        labels = sprintf('state_%s', 1:N)
      )

      stackedbiasframe = data.frame(
        stackedframe,
        state = stackedstate,
        biasvec = stackedbiasvec,
        truequantile = stackedtruequantile,
        realization = factor(rep(1:repetitions, each=n)))

      QR_evals[[paste('n_', n, '_tau_', t,sep='')]][['stackedbiasframe']] =
        stackedbiasframe


      # (for plotting EM Quantities) -----------------------------
      stackedgammaforb = unlist(
        lapply(
          EM_MC_evals[[paste('n_', n, '_tau_', t,sep='')]],
          FUN= function(repet) repet[['gammaforb']]),
        use.names=FALSE
      )

      stackeddeltaeucl = unlist(
        lapply(
          EM_MC_evals[[paste('n_', n, '_tau_', t,sep='')]],
          FUN= function(repet) repet[['deltaeucl']]),
        use.names=FALSE
      )

      EM_MC_evals[[paste('n_', n, '_tau_', t,sep='')]][['stackedgammaforb']] = stackedgammaforb
      EM_MC_evals[[paste('n_', n, '_tau_', t,sep='')]][['stackeddeltaeucl']] = stackeddeltaeucl

      EM_MC_evals[[paste('n_', n, '_tau_', t,sep='')]][['MCparameters']]=
        data.frame(
          gammaforb=stackedgammaforb,
          deltaeucl=stackeddeltaeucl)
    }
  }

  return(list(
    MSs=MSs,
    fits=fits,
    EM_MC_evals=EM_MC_evals,
    QR_evals=QR_evals
  ))
}









# (single fiteval helper) ------------------------------------------------------
#' @description Evaluation measures of Expectation Maximization part on Markov chain.
#' @return
#' # Markov Chain Measures -------
#'   correctstate: share of correctly predicted states.
#'   crossentropy_avg: average cross-entropy on EM classifiction on state chain.
#'   Note that sequence is ignored: only the individual classifications are
#'   considered.
#'   deltaeucl: euclidean distance of true delta & estimated
#'   gammadist: elementwise distance of estimated and true gamma.
#'   gammaforb: formbenius norm of true gamma and estimated gamma
#'   confusion: caret::confusionMatrix call
EM_MC_eval = function(MS, fit){
  n = nrow(MS$data)

  # In absolute decision: predicted state is most likely state i.e. stateprob max value
  pred_state = apply(fit$state_probs, MARGIN = 1, which.max)
  correctstate = sum(pred_state == MS$state) / n
  confusion = caret::confusionMatrix(data = factor(pred_state),
                                     reference = factor(MS$state))

  deltadist = MS$delta - fit$delta
  gammadist = MS$gamma - fit$gamma

  deltaeucl = sum(deltadist ^ 2)
  gammaforb = sqrt(sum(diag(gammadist %*% t(gammadist))))

  # (EM measure: Crossentropy on Markov-chain-state classification) ------------
  # IDEA: the EM algorithm faces a classification problem on the markov chain;
  # which observation originiated from which state?
  # (1) entropy is a measure of how uncertain the event is in the true distrib:
  # the average amount of information, that one gets from one sample of that
  # probability distribution, telling how unpredictable that probability
  # distribution is (on avg)
  # (2) cross entropy = entropy + Kullback Leibler div = - sum_i(q_i*log(p_i))
  # p true (one hot vector: e.g. N = 3, state = 2 ---> 0,1,0),
  # q predicted distribution (e.g. (0.1,0.7, 0.2))

  q = fastDummies::dummy_cols(data.frame(state = MS$state), select_columns = "state")[, -1]
  p = log(fit$state_probs) # carefull -Inf occurs due to stateprob = 0 exacly

  crossentropy = -1 * apply(q * p, MARGIN = 1, sum) # NaN occur 0 * -Inf
  crossentropy_avg = mean(crossentropy)

  return(list(
    correctstate = correctstate,
    crossentropy_avg = crossentropy_avg,
    deltaeucl = deltaeucl,
    gammadist = gammadist,
    gammaforb = gammaforb,
    confusion = confusion
    )
  )
}



#' @title quantile regression evaluation
#' @description
#'  CORE function to calculate BIAS: INPUT is a frame containing all qr fit
#'  predictions for the respective state (irrespectively of whether that state is
#'  true or not for the given observation) Based on these values, and the
#'  probability of each state (EM prediction) for a given observation, there are
#'  various ways how get a predition for one observation: pred_config
#'  distinguishes them.
#'
#' @param pred_config differentiates three calculations of
#'  MS.quantregs prediction:
#'  (1) select the predictions by true state chain
#'  (2) select by most probable state
#'  (3) weight the N predictions based on their states' probability
#' @remark ready for non-linear qr (RQSS), as long as the model has fitted()
#'   method
QR_eval = function(MS, fit , pred_config = 2){
  #TODO pred_config to fiteval argument
  N = nrow(MS$gamma)
  n = nrow(MS$data)

  # get fitted value irrespectively of statechain
  model_pred = sapply(1:N,
    FUN = function(s) fitted(fit$mod[[s]]))

  colnames(model_pred)=
    sprintf('state_%s',1:N)

  # (three prediction models) --------------------------------------------------
  # (1) select model_pred according to true statechain
  model_pred_true_state = sapply(1:n,
    FUN = function(j) model_pred[j,MS$state[j]])

  # (2) select model_pred according to max state_prob
  pred_state = apply(fit$state_probs, MARGIN = 1, which.max)
  model_pred_most_prob_state = sapply(1:n,
    FUN = function(j) model_pred[j,pred_state[j]])

  # # selection the matrix way
  # rowIndices = 1:n
  # colIndices = pred_state
  # mat = as.matrix(model_pred)
  # mat[rowIndices + nrow(mat) * (colIndices - 1)]

  # (3) weighted version of prediction according to state_prob
  weightedmodel_pred = apply(fit$state_probs * model_pred, MARGIN = 1, sum)

  model_pred_selected = list(
    model_pred_true_state,
    model_pred_most_prob_state,
    weightedmodel_pred)[[pred_config]]


  # ----------------------------------------------------------------------------

  # (true-state estimates)
  # [non-]linear quantile for true state is expectation in qr
  truequantile = family_fun(
    family = MS$yfamily,
    func = 'q',
    params = MS$param,
    p = fit$tau
  )

  # residual
  state_selected_true_residual = MS$data$y - truequantile
  state_selected_pred_residual = MS$data$y - model_pred_selected

  # difference between true prediction and true quantile (expectation) :
  biasvec = model_pred_selected - truequantile

  state_split_bias = lapply(1:N,
                            function(s)
                              biasvec[MS$state == s])


  return(list(
    truequantile= truequantile,
    model_pred = model_pred,
    residual_est = state_selected_pred_residual,
    residual_true = state_selected_true_residual,
    biasvec = biasvec
  ))
}






# (plotting) -------------------------------------------------------------------

#' @description plotting state predictions from EM part to their true value
#' only for the form y ~ x1,
#' @remark does not matter if rq- or rqss-fit. However this function should
#' be enabled to consider partials, if it is not y~x1
#' Also X1 as covariate is hard coded still!
#'
#' @param stateind. integer, indicating the coloring of state_probs on points:
#'  how certain is EM, that this point is in state == i
MS.plot1realization = function(MS, fit, fiteval, stateind = 1) {

  N = nrow(MS$gamma)

  d = data.frame(
    MS$data,
    model_pred=fiteval$model_pred,
    state_probs=fit$state_probs,
    true_state = factor(MS$state, labels = sprintf('state.%s',1:N)),
    true_quantile = fiteval$truequantile,
    bias = fiteval$biasvec)

  # convert to long format
  datlpred = reshape2::melt(
    data = d,
    id.vars = c("x1", 'true_state' , sprintf('state_probs.%s',1:N)),
    measure.vars = sprintf('model_pred.state_%s',1:N)
  )


  p1 = ggplot(aes(x = x1, alpha = 0.3), data = datlpred) +
    geom_point(data = d ,
               aes(y = y,
                   color = fit$state_probs[, stateind],
                   shape = true_state)) +
    geom_line(color = 'red',
              aes(y = value,
                  group = variable)) +
    geom_line(color = 'steelblue',
              data = d,
              aes(y = true_quantile,
                  group = true_state))+
    scale_color_continuous(name = paste("state_prob",stateind))

  # convert to long format
  datlstate_prob  = reshape2::melt(
    data = d,
    id.vars = c("x1", 'true_state', 'bias'),
    measure.vars = sprintf('state_probs.%s',1:N)
  )

  # plot the EM probability of each obs to be in state
  p2 = ggplot(aes(x = x1, color = true_state), data = datlstate_prob) +
    facet_wrap(~variable)+
    geom_point(aes(y = value))+
    ylab('state_probability')

  biases = ggplot(aes(x = x1, y = bias, color = true_state), data = datlstate_prob) +
    facet_wrap(true_state ~ .) +
    geom_point() +
    geom_hline(yintercept = 0)+
    ylab('bias realization')

  gridExtra::grid.arrange(p1, p2, biases, nrow = 3)
}



#' @description Plot true & predicted residuals seperated by true state. Further
#'   indicate the EM's certainty on that observation. As a result, overlaps are
#'   visable. Also biases are plotted colored according to the discrete
#'   convolution (hadamad product) of density values under both states.
fit.plot1resid <- function(MS, fit, fiteval, stateind=1) {
  N = nrow(MS$gamma)
  d = data.frame(
    residual_true = fiteval$residual_true,
    residual_est  = fiteval$residual_est,
    x1 = MS$data$x1,
    state = factor(MS$state, labels = sprintf('state_%s',1:N)),
    state_probs = fit$state_probs,
    bias = fiteval$biasvec
  )

  # reshape to long format
  datl = reshape2::melt(
    data = d,
    id.vars = c("x1", 'state' , sprintf('state_probs.%s', 1:N)),
    measure.vars = c('residual_true', 'residual_est')
  )

  # residual plot with user opportunity to choose the state_prob coloring,
  # tracking the EM's certainty about a residual to be in a certain state.
  # TODO if QR_evals pred_config is which.max, the actual decision could be
  # indicated by shape parameter!
  resids = ggplot(
    data = datl,
    aes(x = x1,
        y = value,
        alpha = 0.001,
        # workaround to select stateprob
        color = lapply(
          sprintf('state_probs.%s', 1:N),
          FUN=function(s) datl[[s]])[[stateind]])) +
    facet_grid(variable ~ state) +
    geom_point() +
    geom_hline(yintercept = 0) +
    ggtitle('Residuals based on true states') +
    ylab('')+
    scale_color_continuous(name = paste("state_prob",stateind))

  # convolution / overlap plot
  d2 = data.frame(MS$overlap, x1=MS$data$x1)
  datl2 = reshape2::melt(
    data = d2,
    id.vars = c("x1"),
    measure.vars =sapply(names(MS$overlap), FUN = function(n) paste('X', n , sep='' ))
  )
  overlap = ggplot(aes(x = x1, y = value), data = datl2)+ geom_point() +facet_grid(~variable)

  plot(overlap)
  plot(resids)

}


#' @description Plot all bias realization of one configuration, discriminated by
#'   true state and the gam /local kernel estimate of bias & MSE (on bias²)
#' @remark if QR_eval(pred_config=2), (==rowwise which.max state_prob is model's
#'   state prediction), the lines spike. The reason is that the state is
#'   misclassified - as consequence, the bias is calculated with the model's
#'   prediction (with faulty state) against true quantile (true state). However,
#'   the lines connect by true state. More or less, the same applies to
#'   (pred_config=3), where the prediction is the weighted average (according to
#'   states' prob) of all three model predictions
plotbiasrealizations = function(stackedbiasframe){

  # # biases with stat_smooth estimate (LOESS)
  # biases = ggplot(aes(x=x1, y=biasvec), data = stackedbiasframe)+
  #   facet_grid(~state)+
  #   geom_point(aes(color=realization, alpha = 0.01))+
  #   geom_line(aes(color=realization))+
  #   stat_smooth(formula = y~x, method = loess)+
  #   geom_hline(yintercept = 0)


  # (local bias/variance measure) -----------------------

  # (kernel-estimate) ------------------------------
  # FIXME: this measure is reaaaly slow for large datasets; as a huge amount of
  # distances must be estimated!! ---> GAM fit with B-splines

  # based on the x1 distance between two observations, build a weighted sum
  d = split(stackedbiasframe, stackedbiasframe$state)
  frame = lapply(d,FUN= function(s){
    # weighing sceme with normal kernel
    distances = dnorm(as.matrix(dist(s$x1, diag=TRUE, upper=TRUE)))
    weights = distances / apply(distances, 1, sum)

    localbias = weights %*% s$biasvec
    s$localbias = localbias

    localmse = weights %*% (s$biasvec)^2
    s$localmse = localmse

    return(s)
  })
  frame = as.data.frame(data.table::rbindlist(frame))


  localbias = ggplot(aes(x=x1),data= frame)+
    facet_grid(~state)+
    geom_point(aes(y=biasvec, color=realization))+
    geom_line( aes(y=biasvec, color=realization))+
    geom_line( aes(y=localbias), size = 1.1, color = 'red')+
    geom_hline(yintercept = 0)

  localMSE = ggplot(aes(x =x1, y = localmse), data = frame)+
    facet_grid(~state)+
    geom_line()+
    ylab('local MSE')

  cat('Number of bias observations over all realizations in states: \n',
      table(stackedbiasframe$state))


  # (gam) -----------------------------
  gbias = gam(formula = biasvec ~ s(x1, bs="ps", k=20 , by =state),
              data =stackedbiasframe)

  gmse = gam(formula = biasvec^2 ~ s(x1, bs="ps", k=20 , by =state),
             data =stackedbiasframe)

  stackedbiasframe$msefit = fitted(gmse)
  stackedbiasframe$biasfit = fitted(gbias)

  gamfit = ggplot(aes(x=x1),data= stackedbiasframe)+
    facet_grid(~state)+
    geom_point(aes(y=biasvec, color=realization))+
    geom_line( aes(y=biasvec, color=realization))+
    geom_line(aes(y=biasfit),linetype = 'dotdash', size = 1.1)+
    geom_line(aes(y=msefit),linetype = 'dotted')+
    geom_hline(yintercept = 0)+
    geom_line( aes(x= x1, y=localbias), size = 1, color = 'red',linetype = 'dotted', data = frame)

  # seperately plot MSE -----
  # ggplot(aes(x =x1 ,y=msefit),
  #        data = stackedbiasframe) +
  #   facet_grid(~state)+
  #   geom_line()

  # gammse = ggplot(aes(x=x1),data= stackedbiasframe)+
  #   facet_grid(~state)+
  #   geom_point(aes(y=biasvec, color=realization))+
  #   geom_line( aes(y=biasvec, color=realization))+
  #   geom_line(aes(y=msefit))+
  #   geom_hline(yintercept = 0)

  gridExtra::grid.arrange(localbias, localMSE, nrow = 2)
  plot(gamfit)
}



#' @description Plot evaluation of statechain parameter quality
#' plotting the forbenious norm of the estimated TPM to the true TPM
#' plotting the euclidean norm of the estimated delta to the true delta
#' @example plotgammas(stackedgammaforb=EM_MC_evals$n_50_tau_0.3[['stackedgammaforb']])
plotMCrealizations = function(stackedgammaforb, stackeddeltaeucl){
  gammas = qplot(y=stackedgammaforb)+
    scale_x_discrete(
      name ="Realization",
      limits=1:length(stackedgammaforb))+
    ylab('Forbeniusnorm')+
    ggtitle('Gamma')

  # TODO OPTIONAL informative feature expected share of observations in each
  # state based on gamma

  deltas = qplot(y=stackeddeltaeucl)+
    scale_x_discrete(
      name ="Realization",
      limits=1:length(stackeddeltaeucl))+
    ylab('Eucledian Norm')+
    ggtitle('Delta')

  gridExtra::grid.arrange(gammas, deltas, nrow = 2)
}

