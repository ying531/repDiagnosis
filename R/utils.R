.decomposition_helper = function(
    data1, # original data
    data2, # replication data
    analysis_formula, # a formula for OLS
    focal_variable, # the name of the focal variable
    covariate_formula = NULL, # covariates to be balanced
    mediation_formula = NULL, # covariates + mediators to be balanced
    cluster_id, # variable to cluster on. well-defined before calling this function
    selection_variable = NULL # variable on which the selection happens
) {


  group1 = data1 %>% group_by_at(cluster_id) %>% slice(1)
  group2 = data2 %>% group_by_at(cluster_id) %>% slice(1)

  no_covariate_formula = is.null(covariate_formula) || covariate_formula==""
  no_mediation_formula = is.null(mediation_formula) || mediation_formula==""


  # Create design matrices.  Care to drop unused factor levels, for the jackknife
  if (!no_covariate_formula) {
    X1 = model.matrix(formula(covariate_formula), data=group1)[,-1]
    X1 = X1[,apply(X1, 2, function(x) { !all(x==0) })]
    X2 = model.matrix(formula(covariate_formula), data=group2)[,colnames(X1)]
    X = rbind(X1, X2)
  }
  if (!no_mediation_formula) {
    M1 = model.matrix(formula(mediation_formula), data=group1)[,-1]
    M1 = M1[,apply(M1, 2, function(m) { !all(m==0) })]
    M2 = model.matrix(formula(mediation_formula), data=group2)[,colnames(M1)]
    M = rbind(M1, M2)
  }

  # Fitting covariate weights
  S = c(rep(1, nrow(group1)), rep(0, nrow(group2)))
  if (no_covariate_formula) {
    group2$w = rep(1, nrow(group2))
  } else {
    group2$w = ebalance(Treatment=S, X=X, print.level=-1)$w
  }

  # Fitting mediation weights
  if (no_mediation_formula) {
    group2$o = group2$w
  } else {
    group2$o = ebalance(Treatment=S, X=M, print.level=-1)$w
  }
  data2 = left_join(data2, group2[,c(cluster_id, "w", "o")], by = c(cluster_id))

  # Re-performing analysis
  selection = coef(lm(formula(analysis_formula), data=data1))[selection_variable]
  theta1 = coef(lm(formula(analysis_formula), data=data1))[focal_variable]
  theta2 = coef(lm(formula(analysis_formula), data=data2))[focal_variable]
  theta2X= coef(lm(formula(analysis_formula), data=data2, weights=data2$w))[focal_variable]
  theta2M= coef(lm(formula(analysis_formula), data=data2, weights=data2$o))[focal_variable]

  # Computing decomposition
  decomp = data.frame(
    "Selected" = selection,
    "Observed" = theta1 - theta2,
    "Covariates" = theta2X - theta2,
    "Mediators" = theta2M - theta2X,
    "Residual" = theta1 - theta2M
  )
  if (no_covariate_formula) {
    decomp = select(decomp, -Covariates)
  }
  if (no_mediation_formula) {
    decomp = select(decomp, -Mediators)
  }

  return(list("decomp" = decomp, "w" = group2$w, "o" = group2$o))

}



.selective_decomposition = function(decomp, Sigma, threshold, alpha=0.1) {
  # Selective confidence intervals and conditional MLEs
  x = decomp#$decomp
  S = Sigma#$Sigma
  observed = x[["Observed"]]

  # Compute selection-adjusted confidence intervals for each component
  components = setdiff(names(x), "Selected")
  bounds = lapply(
    X = components,
    FUN = function(v) { vars = c(v, "Selected"); .selective_ci(x[vars], S[vars, vars], threshold, alpha)}
  ) %>% bind_rows()
  rownames(bounds) = components

  # Compute selection-adjusted point estimates for the non-residual components
  vars = setdiff(components, "Residual")
  SX = S[vars, vars]
  SXY = S[vars, "Selected"]
  Omega = solve(SX)
  beta = as.numeric(Omega %*% SXY)
  sigma = as.numeric(sqrt(t(beta) %*% SX %*% beta))
  R = as.numeric(x["Selected"] - sum(beta*x[vars]))
  estimate = optim(
    par = as.numeric(unlist((bounds["high"] - bounds["low"])/2))[1:length(vars)], # Midpoint initialization
    fn = function(par) {
      0.5*t(as.numeric(x[vars]-par)) %*% Omega %*% as.numeric(x[vars]-par) +
        .truncprob(mean = sum(par*beta), sd = sigma, R = R, threshold = threshold, log.prob = T)
    }
  )$par # Maximum conditional likelihood
  names(estimate) = vars
  bounds$estimate = sapply(
    X = components,
    FUN = function(comp) {
      if (comp == "Observed") {
        return(observed)
      } else if (comp == "Residual") {
        return(as.numeric(2*estimate["Observed"] - sum(estimate)))
      } else {
        return(estimate[[comp]])
      }
    }
  )

  # Compute the estimate of sampling variability
  bounds = rbind(bounds, observed - c(bounds["Observed", "high"], bounds["Observed", "low"], estimate["Observed"]))

  # Return results
  bounds$component0 = c(components, "Sampling")
  clean_names = data.frame(
    component0 = c("Observed", "Covariates", "Mediators", "Residual", "Sampling"),
    component = c("Observed discrepancy", "+ Covariate shift", "+ Mediation shift", "+ Residual", "= Sampling variability")
  )
  bounds = left_join(data.frame(bounds), clean_names, by = c("component0")) %>% select(-component0)
  return(bounds)
}


.discrepancy_plot = function(bounds, alpha=0.1, base_size=28, digital=2) {
  # Bounds should be a data frame of the type returned by "selective_decomposition"
  df_plot = bounds %>% mutate(
    level = factor(component, levels=c("+ Residual", "+ Mediation shift", "+ Covariate shift", "= Sampling variability", "Observed discrepancy"))
  )
  # Location of the numbers and line
  xnums = min(df_plot$low) - 0.2*(max(df_plot$high) - min(df_plot$low))
  xline = min(df_plot$low) - 0.6*(max(df_plot$high) - min(df_plot$low)) - 0.008*(digital-2)
  balance_sheet_plot = df_plot %>% ggplot() +
    geom_crossbar(aes(y = level, x = estimate, col = level, fill = level,
                      xmin = low, xmax = high),
                  width = 0.5, alpha = 0.5) +
    geom_text(aes(x = xnums, y = level, label=paste(ifelse(estimate>=0, "+", ""),
                                                    as.character(round(estimate, digital)), sep="")),
              hjust = 1, size = base_size*0.27) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_vline(xintercept = df_plot$estimate[1], lty=2) +
    ylab("") +
    xlab("") +
    theme_bw(base_size = base_size) +
    theme(
      legend.position = "None",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(face = c(rep("italic", 4), "bold")),
      axis.text.x = element_text(size=base_size*0.5),
      plot.margin = unit(c(0.7, 0.7, -0.5, -0.5), "cm")
    ) +
    annotate(geom = "segment", x = xline, xend = xline, y = -Inf, yend = Inf) +
    scale_x_continuous(breaks = c(0, round(df_plot$estimate[1], 2)))
  return(balance_sheet_plot)
}


.truncprob = function(mean, sd, R, threshold, log.prob=TRUE) {
  # Probability that |N(mean, sd) + R| > threshold > 0
  upper = 1-pnorm((threshold-R-mean)/sd) # P{N(mu, sigma^2) > threshold - R}
  lower = pnorm(-(threshold+R+mean)/sd) # P{N(mu, sigma^2) < -threshold - R}
  return(ifelse(log.prob, log(upper+lower), upper+lower))
}

.selective_pval = function(mu, x, Sigma, threshold) {
  x = as.numeric(x)
  threshold = as.numeric(threshold)
  beta = as.numeric(Sigma[1,2]/Sigma[1,1])
  s1 = as.numeric(sqrt(Sigma[1,1]))
  R = as.numeric(x[2] - beta*x[1])
  if (beta >= 0) {
    num =
      pnorm((pmin(-(threshold+R)/beta, x[1]) - mu)/s1) +
      pmax(0, pnorm((x[1]-mu)/s1) - pnorm(((threshold-R)/beta - mu)/s1))
    denom =
      pnorm((-(threshold+R)/beta - mu)/s1) + pnorm(((threshold-R)/beta - mu)/s1, lower.tail=F)
    if (denom==0){return(1)}else{return(num/denom)}
    # return(num/denom)
  } else {
    num =
      pnorm((pmin((threshold-R)/beta, x[1]) - mu)/s1) +
      pmax(0, pnorm((x[1]-mu)/s1) - pnorm((-(threshold+R)/beta - mu)/s1))
    denom =
      pnorm(((threshold-R)/beta - mu)/s1) + pnorm((-(threshold+R)/beta - mu)/s1, lower.tail=F)
    if (denom==0){return(1)}else{return(num/denom)}
  }
}

.selective_ci = function(x, Sigma, threshold, alpha = 0.1) {
  # Confidence interval for mu[1], valid conditional on |x[2]| > threshold
  low = uniroot(
    f = function(mu) { .selective_pval(mu, x, Sigma, threshold) - alpha/2 },
    lower = as.numeric(x[1] - 10*sqrt(Sigma[1,1])),
    upper = as.numeric(x[1] + 10*sqrt(Sigma[1,1]))
  )$root
  high = uniroot(
    f = function(mu) { .selective_pval(mu, x, Sigma, threshold) - (1-alpha/2) },
    lower = as.numeric(x[1] - 10*sqrt(Sigma[1,1])),
    upper = as.numeric(x[1] + 10*sqrt(Sigma[1,1]))
  )$root
  return(data.frame("low" = min(low, high), "high" = max(low, high)))
}

















###
### old functions ###

#
# .decomposition = function(
#     data1, # a data frame
#     data2, # a data frame
#     analysis_formula, # a formula for OLS
#     treatment_variable, # the name of the focal variable
#     covariate_formula = NULL, # covariates to be balanced
#     mediation_formula = NULL, # covariates + mediators to be balanced
#     clusters_1 = NULL, # variables defining clusters.  must contain all variables in the covariate/mediation formulae
#     clusters_2 = NULL
# ) {
#   if (is.null(clusters_1)) {
#     group1 = data1
#     # group1$id = 1:nrow(group1)
#     # data1$id = 1:nrow(data1)
#   } else {
#     group1 = unique(data1[,clusters_1])
#   }
#   if (is.null(clusters_2)) {
#     group2 = data2
#     # group2$idx = 1:nrow(group2)
#     # data2$idx = 1:nrow(data2)
#   } else {
#     group2 = unique(data2[,clusters_2])
#   }
#   X = model.matrix(formula(covariate_formula), data=bind_rows(group1, group2))[,-1]
#   S = c(rep(1, nrow(group1)), rep(0, nrow(group2)))
#   group2$w = ebalance(Treatment=S, X=X, print.level=-1)$w # covariate weights
#
#   if (!is.null(mediation_formula)){
#     M = model.matrix(formula(mediation_formula), data=bind_rows(group1, group2))[,-1]
#     group2$o = ebalance(Treatment=S, X=M, print.level=-1)$w # mediation weights
#   }
#
#   if (is.null(clusters_2)){
#     data2weighted = group2
#   }else{
#     data2weighted = left_join(data2, group2, by=clusters_2)
#   }
#
#   theta1 = coef(lm(formula(analysis_formula), data=data1))[treatment_variable]
#   theta2 = coef(lm(formula(analysis_formula), data=data2))[treatment_variable]
#   theta2X= coef(lm(formula(analysis_formula), data=data2weighted, weights=w))[treatment_variable]
#
#   if (!is.null(mediation_formula)){
#     theta2M= coef(lm(formula(analysis_formula), data=data2weighted, weights=o))[treatment_variable]
#     decomp = data.frame(
#       "Original" = theta1,
#       "Observed" = theta1 - theta2,
#       "Covariates" = theta2X - theta2,
#       "Mediators" = theta2M - theta2X,
#       "Residual" = theta1 - theta2M
#     )
#   }else{ # mediation formula is null
#     decomp = data.frame(
#       "Original" = theta1,
#       "Observed" = theta1 - theta2,
#       "Covariates" = theta2X - theta2,
#       "Residual" = theta1 - theta2X
#     )
#   }
#
#
#   return(list("decomp" = decomp, "w" = group2$w, "o" = group2$o))
# }
#
#
#
#
#
# .selective_t = function(x, delta, sigma, threshold) {
#   # Tests the null hypothesis mu[1] = delta under X ~ N(mu, Sigma)
#   # Returns a P-value that remains valid conditional on |x2| > threshold
#   tau = as.numeric(threshold)
#   d = as.numeric(delta)
#   stat = as.numeric(x[1])
#   sd = as.numeric(sqrt(sigma[1,1]))
#   beta = as.numeric(sigma[1,2]/sigma[1,1])
#   R = as.numeric(x[2] - beta*x[1])
#   low = min((tau-R)/beta, -(tau+R)/beta)
#   high = max((tau-R)/beta, -(tau+R)/beta)
#   num = pnorm((min(stat, low) - d)/sd) + pmax(0, pnorm((stat - d)/sd) - pnorm((high - d)/sd))
#   denom = 1 + pnorm((low - d)/sd) - pnorm((high - d)/sd)
#   if (denom == 0){
#     return(Inf)
#   }else{
#     return(num/denom)
#   }
#
# }
#
# .selective_ci = function(x, sigma, threshold, alpha = 0.1) {
#   # Selective confidence interval for mu[1] conditional on |x[2]| > threshold
#   lower = as.numeric(x[1] - 20*sqrt(sigma[1,1])) - 10
#   upper = as.numeric(x[1] + 20*sqrt(sigma[1,1])) + 10
#
#   delta.try = seq(lower, upper, by=0.01)
#   all.val = sapply(delta.try, .selective_t, x=x, sigma=sigma, threshold=threshold)
#
#   if (min(abs(all.val-alpha/2)) <= 0.1 * alpha){
#     low.try = delta.try[which.min(abs(all.val-alpha/2))]
#   }else{
#     low.try = -Inf
#   }
#
#   if (min(abs(all.val-0.5)) <= 0.1 * alpha){
#     med.try = delta.try[which.min(abs(all.val-0.5))]
#   }else{
#     med.try = -Inf
#   }
#
#   if (min(abs(all.val-1+alpha/2)) <= 0.1 * alpha){
#     high.try = delta.try[which.min(abs(all.val-1+alpha/2))]
#   }else{
#     high.try = Inf
#   }
#
#
#   return(data.frame(
#     "low" = min(low.try, high.try), "high" = max(low.try, high.try), "estimate" = med.try,
#     "lowest" = lower, "highest" = upper
#   ))
# }

