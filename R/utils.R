.decomposition = function(
    data1, # a data frame
    data2, # a data frame
    analysis_formula, # a formula for OLS
    treatment_variable, # the name of the focal variable
    covariate_formula = NULL, # covariates to be balanced
    mediation_formula = NULL, # covariates + mediators to be balanced
    clusters_1 = NULL, # variables defining clusters.  must contain all variables in the covariate/mediation formulae
    clusters_2 = NULL
) {
  if (is.null(clusters_1)) {
    group1 = data1
    # group1$id = 1:nrow(group1)
    # data1$id = 1:nrow(data1)
  } else {
    group1 = unique(data1[,clusters_1])
  }
  if (is.null(clusters_2)) {
    group2 = data2
    # group2$idx = 1:nrow(group2)
    # data2$idx = 1:nrow(data2)
  } else {
    group2 = unique(data2[,clusters_2])
  }
  X = model.matrix(formula(covariate_formula), data=bind_rows(group1, group2))[,-1]
  S = c(rep(1, nrow(group1)), rep(0, nrow(group2)))
  group2$w = ebalance(Treatment=S, X=X, print.level=-1)$w # covariate weights

  if (!is.null(mediation_formula)){
    M = model.matrix(formula(mediation_formula), data=bind_rows(group1, group2))[,-1]
    group2$o = ebalance(Treatment=S, X=M, print.level=-1)$w # mediation weights
  }

  if (is.null(clusters_2)){
    data2weighted = group2
  }else{
    data2weighted = left_join(data2, group2, by=clusters_2)
  }

  theta1 = coef(lm(formula(analysis_formula), data=data1))[treatment_variable]
  theta2 = coef(lm(formula(analysis_formula), data=data2))[treatment_variable]
  theta2X= coef(lm(formula(analysis_formula), data=data2weighted, weights=w))[treatment_variable]

  if (!is.null(mediation_formula)){
    theta2M= coef(lm(formula(analysis_formula), data=data2weighted, weights=o))[treatment_variable]
    decomp = data.frame(
      "Original" = theta1,
      "Observed" = theta1 - theta2,
      "Covariates" = theta2X - theta2,
      "Mediators" = theta2M - theta2X,
      "Residual" = theta1 - theta2M
    )
  }else{ # mediation formula is null
    decomp = data.frame(
      "Original" = theta1,
      "Observed" = theta1 - theta2,
      "Covariates" = theta2X - theta2,
      "Residual" = theta1 - theta2X
    )
  }


  return(list("decomp" = decomp, "w" = group2$w, "o" = group2$o))
}





.selective_t = function(x, delta, sigma, threshold) {
  # Tests the null hypothesis mu[1] = delta under X ~ N(mu, Sigma)
  # Returns a P-value that remains valid conditional on |x2| > threshold
  tau = as.numeric(threshold)
  d = as.numeric(delta)
  stat = as.numeric(x[1])
  sd = as.numeric(sqrt(sigma[1,1]))
  beta = as.numeric(sigma[1,2]/sigma[1,1])
  R = as.numeric(x[2] - beta*x[1])
  low = min((tau-R)/beta, -(tau+R)/beta)
  high = max((tau-R)/beta, -(tau+R)/beta)
  num = pnorm((min(stat, low) - d)/sd) + pmax(0, pnorm((stat - d)/sd) - pnorm((high - d)/sd))
  denom = 1 + pnorm((low - d)/sd) - pnorm((high - d)/sd)
  if (denom == 0){
    return(Inf)
  }else{
    return(num/denom)
  }

}

.selective_ci = function(x, sigma, threshold, alpha = 0.1) {
  # Selective confidence interval for mu[1] conditional on |x[2]| > threshold
  lower = as.numeric(x[1] - 20*sqrt(sigma[1,1])) - 10
  upper = as.numeric(x[1] + 20*sqrt(sigma[1,1])) + 10

  delta.try = seq(lower, upper, by=0.01)
  all.val = sapply(delta.try, .selective_t, x=x, sigma=sigma, threshold=threshold)

  if (min(abs(all.val-alpha/2)) <= 0.1 * alpha){
    low.try = delta.try[which.min(abs(all.val-alpha/2))]
  }else{
    low.try = -Inf
  }

  if (min(abs(all.val-0.5)) <= 0.1 * alpha){
    med.try = delta.try[which.min(abs(all.val-0.5))]
  }else{
    med.try = -Inf
  }

  if (min(abs(all.val-1+alpha/2)) <= 0.1 * alpha){
    high.try = delta.try[which.min(abs(all.val-1+alpha/2))]
  }else{
    high.try = Inf
  }


  return(data.frame(
    "low" = min(low.try, high.try), "high" = max(low.try, high.try), "estimate" = med.try,
    "lowest" = lower, "highest" = upper
  ))
}

