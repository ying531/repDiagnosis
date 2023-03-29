#' @title Run diagnosis function
#' @description This function runs the diagnosis analysis.
#' @param df1 ...
#' @param df2 ...
#' @param treatment_name ...
#' @param id_name_1 ...
#' @param id_name_2 ...
#' @param analysis_formula ...
#' @param covariate_formula ...
#' @param mediation_formula ...
#' @param clusters_1 ...
#' @param clusters_2 ...
#' @return A list containing the result table, decomp.plot, and selective.plt.
#' @export


run_diagnosis <- function(
    df1,
    df2,
    treatment_name,
    id_name_1,
    id_name_2,
    analysis_formula, # = "pv ~ treatment + factor(delay)*factor(fv)"
    covariate_formula, # = "~ practice_religion + religion + race"
    mediation_formula, # = "~ (practice_religion + religion + race + panas + happiness_induced + mood_induced)*treatment"
    clusters_1, #= c("id", "treatment", "practice_religion", "happiness_induced", "mood_induced", "happiness", "religion", "panas", "race")
    clusters_2
){


  df1 = data.frame(df1)
  df2 = data.frame(df2)
  results = decomposition(
    data1 = df1,
    data2 = df2,
    analysis_formula = analysis_formula,
    treatment_variable = treatment_name,
    covariate_formula = covariate_formula,
    mediation_formula = mediation_formula,
    clusters_1 = clusters_1,
    clusters_2 = clusters_2
  )

  # choose the cluster id if it is not null
  if (!is.null(id_name_1) && !(id_name_1 == '')){
    id1 = unique(df1[,id_name_1])
  }else{
    id1 = 1:dim(df1)[1]
  }

  if (!is.null(id_name_2) && !(id_name_2 == '')){
    id2 = unique(df2[,id_name_2])
  }else{
    id2 = 1:dim(df2)[1]
  }

  N1 = length(id1)
  N2 = length(id2)
  phi1 = phi2 = data.frame()
  pb = progress_bar$new("Jackknifing [:bar] :percent", total=N1 + N2, width=80)

  # Add progress bar
  progress <- shiny::Progress$new()
  progress$set(message = "Running analysis...", value = 0)
  on.exit(progress$close())

  total_steps <- N1 + N2
  for (i in 1:N1) {
    pb$tick()
    # Update progress bar
    progress$set(value = i / total_steps)

    if (!is.null(id_name_1) && !(id_name_1 == '')){
      sub.df1 = filter(df1, !!sym(id_name_1) != id1[i])
    }else{
      sub.df1 = df1[-i,]
    }
    try({phi1 = bind_rows(phi1, decomposition(
      data1 = sub.df1,
      data2 = df2,
      analysis_formula = analysis_formula,
      treatment_variable = treatment_name,
      covariate_formula = covariate_formula,
      mediation_formula = mediation_formula,
      clusters_1 = clusters_1,
      clusters_2 = clusters_2
    )$decomp)}, silent = TRUE)
  }
  for (i in 1:N2) {
    pb$tick()
    # Update progress bar
    progress$set(value = (N1 + i) / total_steps)

    if (!is.null(id_name_2) && !(id_name_2 == '')){
      sub.df2 = df2[df2[,id_name_2]!=id2[i],] #filter(df2, !!sym(id_name_2) != id2[i])
    }else{
      sub.df2 = df2[-i,]
    }
    try({phi2 = bind_rows(phi2, decomposition(
      data1 = df1,
      data2 = sub.df2, #filter(df2, !!sym(id_name_2) != id2[i]),
      analysis_formula = analysis_formula,
      treatment_variable = treatment_name,
      covariate_formula = covariate_formula,
      mediation_formula = mediation_formula,
      clusters_1 = clusters_1,
      clusters_2 = clusters_2
    )$decomp)}, silent = TRUE)
  }
  SEs = sqrt(apply(phi1, 2, var)*(nrow(phi1)-1)^2/N1 + apply(phi2, 2, var)*(nrow(phi2)-1)^2/N2)
  Sigma = var(phi1)*(nrow(phi1)-1)^2/N1 + var(phi2)*(nrow(phi2)-1)^2/N2

  # ==================================================
  # Plotting results
  bounds = cbind(as.numeric(results$decomp),
                 as.numeric(results$decomp) - qnorm(0.95)*SEs,
                 as.numeric(results$decomp) + qnorm(0.95)*SEs) %>% data.frame()
  pvals = 2 - 2 * pnorm(abs(as.numeric(results$decomp))/SEs)
  colnames(bounds) = c("estimate", "low", "high")
  bounds$component = rownames(bounds)

  if (!is.null(mediation_formula)){
    decomp.plot = bounds %>% filter(component != "Original") %>%
      mutate(component = factor(component, levels = c("Observed", "Covariates", "Mediators", "Residual"))) %>%
      ggplot(aes(x = component, fill = component, y = estimate, ymin = low, ymax = high)) +
      geom_crossbar(aes(alpha = 0.8, col = component), width=0.5) +
      theme_bw() + xlab("") + ylab("") +
      # theme(legend.position = "None") +
      geom_hline(yintercept = 0, lty=2)+
      theme(text= element_text(family="Times", size=15),
            axis.text= element_text(family="Times", size=12),
            strip.text.x= element_text(family="Times", size=12),
            strip.text.y= element_text(family="Times", size=12),
            legend.title = element_text(family="Times", size=15),
            legend.position="None",
            plot.title = element_text(family="Times", size=15, hjust = 0.5))


    ret_table = cbind(as.numeric(results$decomp), SEs, as.numeric(results$decomp)/SEs, pvals)
    colnames(ret_table) = c("Estimate", "Std. Error", "t-stat", "Pr(>|z|)")
    rownames(ret_table) = c("Original", "Observed", "Covariates", "Mediators", "Residual")
    ret_table_df <- as.data.frame(ret_table)
    ret_table_df$Component <- rownames(ret_table)
    ret_table_df <- ret_table_df[, c("Component", "Estimate", "Std. Error", "t-stat", "Pr(>|z|)")]


    # ==============================================================
    # Accounting for publication bias at the p = 0.05 level

    threshold = qnorm(0.95)*summary(lm(analysis_formula, data = df1))$coefficients[treatment_name, "Std. Error"]
    components = c("Observed", "Covariates", "Mediators", "Residual")
    bounds = lapply(components, FUN = function(v) {
      selective_ci(
        x = results$decomp[c(v, "Original")],
        sigma = Sigma[c(v, "Original"), c(v, "Original")],
        threshold = threshold,
        alpha = 0.1
      )}) %>% bind_rows()

    bounds$component = components

    ret_table_sel = data.frame(bounds[,c("low", "high", "estimate")])
    rownames(ret_table_sel) = components


    # generate error message
    sel.message = "NOTE:"
    bounds.sel.tp = data.frame()
    for (v in 1:length(components)){

      if ((bounds$low[v] == -Inf) || (bounds$high[v] == Inf)){
        sel.message = paste(sel.message,
                            "No meaningful selective CI within [",
                            bounds$lowest[v], bounds$highest[v],"] for ",
                            components[v], "component! \n")
      }else{
        to.add = ret_table_sel[v,]
        if ((bounds$estimate[v] == -Inf) || (bounds$estimate[v] < bounds$low[v]) || (bounds$estimate[v] > bounds$high[v])){
          to.add$estimate = (to.add$low + to.add$high)/2
        }
        bounds.sel.tp = rbind(bounds.sel.tp, to.add)
      }
      #
      # if (bounds$estimate[v] == -Inf){
      #   ret_table_sel$estimate[v] = "NULL"
      # }
    }

    if (nrow(bounds.sel.tp)>0){
      bounds.sel.tp$components = rownames(bounds.sel.tp)
      selective.plt = bounds.sel.tp %>%
        mutate(component = factor(components, levels = c("Original", "Observed", "Covariates", "Mediators", "Residual"))) %>%
        ggplot(aes(x = component, fill = component, y = estimate, ymin = low, ymax = high)) +
        geom_crossbar(aes(col = component), alpha = 0.5, width=0.6) +
        theme_bw() + xlab("") + ylab("") +
        theme(legend.position = "None") +
        geom_hline(yintercept = 0, lty=2) +
        theme(text= element_text(family="Times", size=15),
              axis.text= element_text(family="Times", size=12),
              strip.text.x= element_text(family="Times", size=12),
              strip.text.y= element_text(family="Times", size=12),
              legend.title = element_text(family="Times", size=15),
              legend.position="None",
              plot.title = element_text(family="Times", size=15, hjust = 0.5))
    }else{
      selective.plt = NULL
    }


  }else{

    ###### mediation formula is null

    decomp.plot = bounds %>% filter(component != "Original") %>%
      mutate(component = factor(component, levels = c("Observed", "Covariates", "Residual"))) %>%
      ggplot(aes(x = component, fill = component, y = estimate, ymin = low, ymax = high)) +
      geom_crossbar(aes(alpha = 0.8, col = component), width=0.5) +
      theme_bw() + xlab("") + ylab("") +
      # theme(legend.position = "None") +
      geom_hline(yintercept = 0, lty=2)+
      theme(text= element_text(family="Times", size=15),
            axis.text= element_text(family="Times", size=12),
            strip.text.x= element_text(family="Times", size=12),
            strip.text.y= element_text(family="Times", size=12),
            legend.title = element_text(family="Times", size=15),
            legend.position="None",
            plot.title = element_text(family="Times", size=15, hjust = 0.5))


    ret_table = cbind(as.numeric(results$decomp), SEs, as.numeric(results$decomp)/SEs, pvals)
    colnames(ret_table) = c("Estimate", "Std. Error", "t-stat", "Pr(>|z|)")
    rownames(ret_table) = c("Original", "Observed", "Covariates", "Residual")
    ret_table_df <- as.data.frame(ret_table)
    ret_table_df$Component <- rownames(ret_table)
    ret_table_df <- ret_table_df[, c("Component", "Estimate", "Std. Error", "t-stat", "Pr(>|z|)")]


    # ==============================================================
    # Accounting for publication bias at the p = 0.05 level

    threshold = qnorm(0.95)*summary(lm(analysis_formula, data = df1))$coefficients[treatment_name, "Std. Error"]
    components = c("Observed", "Covariates", "Residual")
    bounds = lapply(components, FUN = function(v) {
      selective_ci(
        x = results$decomp[c(v, "Original")],
        sigma = Sigma[c(v, "Original"), c(v, "Original")],
        threshold = threshold,
        alpha = 0.1
      )}) %>% bind_rows()

    bounds$component = components

    ret_table_sel = data.frame(bounds[,c("low", "high", "estimate")])
    rownames(ret_table_sel) = components

    # generate error message
    sel.message = "NOTE:"
    bounds.sel.tp = data.frame()
    for (v in 1:length(components)){

      if ((bounds$low[v] == -Inf) || (bounds$high[v] == Inf)){
        sel.message = paste(sel.message,
                            "No meaningful selective CI within [",
                            bounds$lowest[v], bounds$highest[v],"] for ",
                            components[v], "component! <br>")
      }else{
        to.add = ret_table_sel[v,]
        if  ((bounds$estimate[v] == -Inf) || (bounds$estimate[v] < bounds$low[v]) || (bounds$estimate[v] > bounds$high[v])){
          to.add$estimate = (to.add$low + to.add$high)/2
        }
        bounds.sel.tp = rbind(bounds.sel.tp, to.add)
      }

    }



    if (nrow(bounds.sel.tp)>0){
      bounds.sel.tp$components = rownames(bounds.sel.tp)
      selective.plt = bounds.sel.tp %>%
        mutate(component = factor(components, levels = c("Original", "Observed", "Covariates", "Mediators", "Residual"))) %>%
        ggplot(aes(x = component, fill = component, y = estimate, ymin = low, ymax = high)) +
        geom_crossbar(aes(col = component), alpha = 0.5, width=0.6) +
        theme_bw() + xlab("") + ylab("") +
        theme(legend.position = "None") +
        geom_hline(yintercept = 0, lty=2) +
        theme(text= element_text(family="Times", size=15),
              axis.text= element_text(family="Times", size=12),
              strip.text.x= element_text(family="Times", size=12),
              strip.text.y= element_text(family="Times", size=12),
              legend.title = element_text(family="Times", size=15),
              legend.position="None",
              plot.title = element_text(family="Times", size=15, hjust = 0.5))
    }else{
      selective.plt = NULL
    }


  }


  # decomp.plot
  # print(ret_table_df)
  # selective.plt

  ret_table_sel$Component <- rownames(ret_table_sel)
  ret_table_sel <- ret_table_sel[, c("Component", "low", "high", "estimate")]
  colnames(ret_table_sel) = c("Component", "CI.Low", "CI.High", "Estimate")

  ret_table_sel$Estimate[ret_table_sel$Estimate==-Inf] = "NULL"

  return(list("plot" = decomp.plot, "table" = ret_table_df,
              "sel.plot" = selective.plt, "sel.table" = ret_table_sel,
              "sel.message" = sel.message))

}






#' @title Run reweight function
#' @description This function runs the reweight analysis.
#' @param new_means ...
#' @param df ...
#' @param regression_formula ...
#' @param covariates ...
#' @param group_id ...
#' @return A list containing the result table, coefficients, and weighted values.
#' @export
reweight_regression <- function(new_means, df, regression_formula, covariates, group_id) {
  group_flag = (!is.null(group_id) && group_id != "")
  if (group_flag){
    # compute means after clustering
    sliced_data <- df %>% group_by(!!sym(group_id)) %>% slice(1)
    means <- sapply(covariates, function(cov) mean(sliced_data[[cov]], na.rm = TRUE))
  }else{
    means <- sapply(covariates, function(cov) mean(df[[cov]], na.rm = TRUE))
  }


  if (group_flag) {
    cluster_formula = paste(regression_formula , "+", group_id, "+", paste(covariates,collapse="+"))
    # cluster_vars = c(covariates, outcome_name, regressors, group_id)
    group1 = unique(model.frame(as.formula(cluster_formula), data = df))
    cluster_vars = colnames(group1)
  } else {
    group1 = df
  }
  covariate_formula = paste("~ (", paste(covariates, collapse = "+"), ")")
  group2 = group1
  for (v in 1:length(covariates)){
    group2[[covariates[v]]] = group2[[covariates[v]]] - as.numeric(means[v]) + as.numeric(new_means[v])
  }
  X = model.matrix(formula(covariate_formula), data=bind_rows(group1, group2))[,-1]
  S = c(rep(0, nrow(group1)), rep(1, nrow(group2))) # reweight group 1 to group 2 means
  group1$w = ebalance(Treatment=S, X=X, print.level=-1)$w # covariate weights

  # point estimate
  if (!group_flag) {
    df$w = group1$w
    model_w <- lm(as.formula(regression_formula), data = df, weights=df$w)
    # model_summary = summary(model_w)
  } else {
    df = left_join(df, group1, by = cluster_vars)
    model_w <- lm(as.formula(regression_formula), data = df, weights=df$w)
  }
  coef_weighted = coefficients(model_w)

  df = df[,c("w", setdiff(colnames(df),c("w")))]
  return(list("coef" = coef_weighted, "weighted" = df))
}


run_reweight <- function(new_means, df, regression_formula, covariates, group_id){
  group_flag = (!is.null(group_id) && group_id != "")
  df = data.frame(df)
  results = reweight_regression(
    new_means = new_means,
    df = df,
    regression_formula = regression_formula,
    covariates = covariates,
    group_id = group_id
  )

  # choose the cluster id if it is not null
  if (group_flag){
    id1 = unique(df[,group_id])
  }else{
    id1 = 1:dim(df)[1]
  }

  N = length(id1)
  phi1 = data.frame()
  pb = progress_bar$new("Jackknifing [:bar] :percent", total=N, width=80)

  # Add progress bar
  progress <- shiny::Progress$new()
  progress$set(message = "Running analysis...", value = 0)
  on.exit(progress$close())

  total_steps <- N
  for (i in 1:N) {
    pb$tick()
    # Update progress bar
    progress$set(value = i / total_steps)

    if (group_flag){
      sub.df1 = filter(df, !!sym(group_id) != id1[i])
    }else{
      sub.df1 = df[-i,]
    }

    try({phi1 = bind_rows(phi1, reweight_regression(
      new_means = new_means,
      df = sub.df1,
      regression_formula = regression_formula,
      covariates = covariates,
      group_id = group_id
    )$coef)}, silent = TRUE)
  }

  # print(phi1)
  SEs = sqrt(apply(phi1, 2, var)*(nrow(phi1)-1)^2/N )

  # ==================================================
  # Plotting results
  bounds = cbind(as.numeric(results$coef),
                 as.numeric(results$coef) - qnorm(0.95)*SEs,
                 as.numeric(results$coef) + qnorm(0.95)*SEs) %>% data.frame()
  pvals = 2 - 2 * pnorm(abs(as.numeric(results$coef))/SEs)
  colnames(bounds) = c("estimate", "low", "high")
  bounds$component = rownames(bounds)

  print(bounds)

  ret_table = cbind(as.numeric(results$coef), SEs, as.numeric(results$coef)/SEs, pvals)
  colnames(ret_table) = c("Estimate", "Std. Error", "t-stat", "Pr(>|z|)")
  ret_table_df <- as.data.frame(ret_table)
  return(list("table" = ret_table_df, "coef" = results$coef, "weighted" = results$weighted))
}
