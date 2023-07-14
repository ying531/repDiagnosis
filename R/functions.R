#' @import ggplot2
#' @import ebal
#' @import tidyverse
#' @import progress
#' @import sandwich
#' @import lmtest
#' @title Run diagnosis function
#' @description This function runs the diagnosis analysis for replication studies. It returns a decomposition of the discrepancy between a pair of original and replication studies into pieces attributable to covariate shift, mediation shift, and residual factors. It also provides post-selective confidence intervals for these pieces that adjust for potential publication bias.
#' @param org_df The data frame for the original study
#' @param rep_df The data frame for the replication study
#' @param analysis_formula The analysis formula for the regression, such as "outcome ~ treatment"
#' @param treatment_variable String for the column name of the treatment indicator (1=treatment)
#' @param focal_variable Optional, string for the column name of the variable on which the decomposition is conducted; if not specified, it will be set as the treatment variable
#' @param covariates Optional, vector for column names or column number of covariates to balance
#' @param covariate_formula Optional, formula for covariates to balance, such as "~ age + gender" which means the covariates to balance are "age" and "gender", or "~(age+gender) * treatment" means they are balanced in each arm; the two choices should not differ significantly. If both covariates and covariate_formula are specified, covariate_formula will be taken.
#' @param mediators Optional, vector for column names or column number of mediators to balance
#' @param mediation_formula Optional, formula for mediators to balance, which should also include covariates, such as "~ (age + gender + mediator) * treatment", which means "mediator" is viewed as a mediator to balance. If both are specified, mediation_formula is taken
#' @param cluster_id Optional, string for the column name of the cluster id in the studies (specify if the study is clustered)
#' @param selection_variable Optional, string for the column name on which the publication bias happens; if not specified, it will be set as the focal variable
#' @param alpha Optional, default at 0.1, coverage level for the confidence intervals is 1-alpha
#' @param verbose Optional, default at TRUE, whether to show progress bars in analysis
#' @param if_selective Optional, default at TRUE, whether or not to run selective inference that adjusts for publication bias
#' @param pub_pvalue_threshold Optional, default at 0.05, the threshold for publication bias, i.e., the selective inference CIs are valid conditional on that the original study p-value is smaller than pub_pvalue_threshold

#' @return A list containing the summary table, decomposition plot, summary table for selective inference, decomposition plot for selective inference, and error message for selective inference in the case of infinite confidence intervals
#' @export
#'
#'

run_diagnosis <- function(
    org_df,
    rep_df,
    analysis_formula,
    treatment_variable, # column name for the treatment indicator (1 = treatment)
    focal_variable = NULL, # focal variable
    covariates = NULL,
    covariate_formula = NULL,
    mediators = NULL,
    mediation_formula = NULL,
    cluster_id = NULL,
    selection_variable = NULL,
    alpha = 0.1,
    verbose = TRUE,
    if_selective = TRUE,
    pub_pvalue_threshold = 0.05
){

  df1 = data.frame(org_df)
  df2 = data.frame(rep_df)



  #####################################################
  ######### sanity check for covariate inputs #########
  #####################################################

  # if covariate formula is not provided, check whether the covariates argument is meaningful
  # otherwise just use the covariate formula
  if (is.null(covariate_formula)){
    if (length(covariates) > 0){
      for (i.par in 1:length(covariates)){
        if (!is.na(suppressWarnings(as.integer(covariates[i.par])))){
          covariates[i.par] = colnames(org_df)[as.integer(covariates[i.par])]
        }
      }
      # if (is.null(treatment_variable)){
      #   stop("Please provide the treatment variable argument!")
      # }
      covariate_formula = paste("~(", paste(covariates, collapse = "+"), ") *", treatment_variable)
    }
  }

  # extract covariates and check whether it behaves normally
  if (!is.null(covariate_formula)){
    X1 = tryCatch({model.matrix(formula(covariate_formula), data=org_df)},
                  error = function(e) { return(NA) })
    X2 = tryCatch({model.matrix(formula(covariate_formula), data=rep_df)},
                  error = function(e) { return(NA) })
    if (is.na(X1) || is.na(X2)){ # any of the data fails to be extracted
      stop("Covariate input does not work! \n")
    }else{# extract covariate names into a vector
      covariates = setdiff(intersect(colnames(org_df), colnames(X1)[2:ncol(X1)]), c(treatment_variable))
      if (is.null(dim(X1[,2:ncol(X1)]))){
        if (sd(X1[,2:ncol(X1)])==0){stop("Constant covariate error! \n")}
      }
    }
  }

  #####################################################
  ######### sanity check for mediation inputs #########
  #####################################################

  # if mediation formula is not provided, check whether the mediators argument is meaningful
  if (is.null(mediation_formula)){
    if (length(mediators) > 0){
      # if (is.null(treatment_variable)){
      #   stop("Please provide the treatment variable argument!")
      # }

      for (i.par in 1:length(mediators)){
        if (!is.na(suppressWarnings(as.integer(mediators[i.par])))){
          mediators[i.par] = colnames(org_df)[as.integer(mediators[i.par])]
        }
      }

      cov_and_med = union(covariates, mediators)
      # generate mediation formula
      mediation_formula = paste("~(", paste(cov_and_med, collapse="+"), ")*", treatment_variable)

    }
  }

  # extract mediators and check whether it behaves normally
  if (!is.null(mediation_formula)){
    X1 = tryCatch({model.matrix(formula(mediation_formula), data=org_df)},
                  error = function(e) { return(NA) })
    X2 = tryCatch({model.matrix(formula(mediation_formula), data=rep_df)},
                  error = function(e) { return(NA) })
    if (is.na(X1) || is.na(X2)){ # any of the data fails to be extracted
      stop("Mediation input does not work! \n")
    }else{# extract covariate names into a vector
      mediators = setdiff(intersect(colnames(org_df), colnames(X1)[2:ncol(X1)]),
                          c(treatment_variable, covariates))
      if (is.null(dim(X1[,2:ncol(X1)]))){
        if (sd(X1[,2:ncol(X1)])==0){stop("Constant mediator error! \n")}
      }
    }
  }

  #####################################################
  ######### sanity check for other inputs #########
  #####################################################

  if (is.null(cluster_id)) {
    df1$idformyuse = 1:nrow(df1)
    df2$idformyuse = 1:nrow(df2)
    cluster_id = "idformyuse"
  }

  no_covariate_formula = is.null(covariate_formula) || covariate_formula==""
  no_mediation_formula = is.null(mediation_formula) || mediation_formula==""

  if (no_covariate_formula & no_mediation_formula) {
    warning("covariate_formula and mediation_formula cannot both be NULL.")
    return()
  }

  if (is.null(focal_variable)){
    focal_variable = treatment_variable
  }

  # selection variable equals the focal if not specified
  if (is.null(selection_variable)){
    selection_variable = focal_variable
  }

  fit = .decomposition_helper(
    data1 = df1,
    data2 = df2,
    analysis_formula = analysis_formula,
    focal_variable = focal_variable,
    covariate_formula = covariate_formula,
    mediation_formula = mediation_formula,
    cluster_id = cluster_id,
    selection_variable = selection_variable
  )

  results = fit$decomp
  clusters1 = unique(df1[[cluster_id]])
  clusters2 = unique(df2[[cluster_id]])


  N1 = length(clusters1)
  N2 = length(clusters2)
  phi1 = phi2 = data.frame()
  pb = progress_bar$new("Jackknifing [:bar] :percent", width=80, total=N1 + N2+1)

  # Add progress bar
  total_steps = N1+N2 + 1


  for (i in 1:N1) {
    if (verbose) { pb$tick() }

    phi1 = bind_rows(phi1, .decomposition_helper(
      data1 = df1[df1[[cluster_id]] != clusters1[i],],
      data2 = df2,
      analysis_formula = analysis_formula,
      focal_variable = focal_variable,
      covariate_formula = covariate_formula,
      mediation_formula = mediation_formula,
      cluster_id = cluster_id,
      selection_variable = selection_variable
    )$decomp)
  }
  for (i in 1:N2) {
    if (verbose) { pb$tick() }

    phi2 = bind_rows(phi2, .decomposition_helper(
      data1 = df1,
      data2 = df2[df2[[cluster_id]] != clusters2[i],],
      analysis_formula = analysis_formula,
      focal_variable = focal_variable,
      covariate_formula = covariate_formula,
      mediation_formula = mediation_formula,
      cluster_id = cluster_id,
      selection_variable = selection_variable
    )$decomp)
  }

  SEs = sqrt(apply(phi1, 2, var)*(nrow(phi1)-1)^2/N1 + apply(phi2, 2, var)*(nrow(phi2)-1)^2/N2)
  Sigma = var(phi1)*(nrow(phi1)-1)^2/N1 + var(phi2)*(nrow(phi2)-1)^2/N2

  # ==================================================
  # Plotting results

  # non-selective results
  pvals = 2 - 2 * pnorm(abs(as.numeric(fit$decomp))/SEs)
  ret_table = cbind(as.numeric(fit$decomp), SEs, as.numeric(fit$decomp)/SEs, pvals)
  colnames(ret_table) = c("Estimate", "Std. Error", "t-stat", "Pr(>|z|)")
  ret_table = ret_table[rownames(ret_table) != "Selected",]
  rownames(ret_table)[1] = "Observed discrepancy"
  ret_table_df <- as.data.frame(ret_table)
  ret_table_df$Component <- rownames(ret_table)
  ret_table_df <- ret_table_df[, c("Component", "Estimate", "Std. Error", "t-stat", "Pr(>|z|)")]

  clean_decomp = select(fit$decomp, -Selected)
  clean_SE = SEs[2:length(SEs)]
  clean_low = clean_decomp
  clean_low[2:length(clean_low)] = clean_low[2:length(clean_low)] - qnorm(1-alpha/2) * clean_SE[2:length(clean_SE)]
  clean_high = clean_decomp
  clean_high[2:length(clean_high)] = clean_high[2:length(clean_high)] + qnorm(1-alpha/2) * clean_SE[2:length(clean_SE)]

  component = c("Observed discrepancy")
  if (!no_covariate_formula){
    component = c(component, "+ Covariate shift")
  }
  if (!no_mediation_formula){
    component = c(component, "+ Mediation shift")
  }
  component = c(component, "+ Residual", "= Sampling variability")
  bounds_for_plot = data.frame("low" = c(as.numeric(clean_low),
                                         -qnorm(1-alpha/2) * (as.numeric(clean_SE[1]))),
                               "high" = c(as.numeric(clean_high),
                                          qnorm(1-alpha/2) * (as.numeric(clean_SE[1]))),
                               "estimate" = c(as.numeric(clean_decomp),0),
                               "component" = component)
  suppressWarnings({
    decomp.plot = .discrepancy_plot(bounds_for_plot, alpha=0.1, base_size=28, digital=2)
  })

  # selective part
  decomp.plot.sel = NULL
  sel.message = NULL
  ret_table_sel_df = NULL
  if (if_selective){
    sel.message = "NOTE:"

    if ((as.numeric(pvals['Selected']) > pub_pvalue_threshold)){
      sel.message = paste(sel.message,
                          "Original p-value > publication bias threshold, stop running selective inference! \n")
      pb$tick()
      return(list("plot" = decomp.plot, "table" = ret_table_df,
                  "sel.plot" = NULL, "sel.table" = NULL,
                  "sel.message" = sel.message))
    }else{
      sel.message = paste(sel.message,
                          paste("Running selective inference conditional on p <= ", pub_pvalue_threshold,"! \n",sep=''))
      bounds_sel = .selective_decomposition(fit$decomp, Sigma,
                                           qnorm(1-pub_pvalue_threshold/2) * sqrt(Sigma["Selected", "Selected"]),
                                           alpha=alpha)
    }


    # generate table
    ret_table_sel = bounds_sel[1:(nrow(bounds_sel)-1),]
    rownames(ret_table_sel) = ret_table_sel$component
    rownames(ret_table_sel)[2:nrow(ret_table_sel)] = sapply(rownames(ret_table_sel)[2:nrow(ret_table_sel)],
                                                            function(x) substr(x, 3, nchar(x)))
    ret_table_sel_df = data.frame(ret_table_sel)

    # selective table
    ret_table_sel_df$Component <- rownames(ret_table_sel_df)
    ret_table_sel_df <- ret_table_sel_df[, c("Component", "low", "high", "estimate")]
    colnames(ret_table_sel_df) = c("Component", "CI.Low", "CI.High", "Estimate")
    ret_table_sel_df$Estimate = as.numeric(ret_table_sel_df$Estimate)

    ret_table_sel_df$Estimate[ret_table_sel_df$Estimate==-Inf] = NULL
    ret_table_sel_df$Estimate[ret_table_sel_df$Estimate==Inf] = NULL

    # generate error message
    components = rownames(ret_table_sel_df)

    for (v in 1:length(components)){

      if ((ret_table_sel_df$CI.Low[v] == -Inf) || (ret_table_sel_df$CI.High[v] == Inf)){
        sel.message = paste(sel.message,
                            "No meaningful selective CI within est +- 10*se for ",
                            components[v], "component! \n")
      }

    }


    # generate selective decomposition plot
    bounds_for_plot_sel = data.frame(bounds_sel) %>% drop_na()
    if (nrow(bounds_for_plot_sel) >0 ){
      if (("Observed discrepancy" %in% bounds_for_plot_sel$component)){
        bounds_for_plot_sel$low[1] = bounds_for_plot_sel$estimate[1]
        bounds_for_plot_sel$high[1] = bounds_for_plot_sel$estimate[1]
      }
      suppressWarnings({
      decomp.plot.sel = .discrepancy_plot(bounds_for_plot_sel, alpha=alpha, base_size=28, digital=2)
      })

    }else{
      decomp.plot.sel = NULL
    }
  }

  if (verbose){
    cat("\n")
    cat("\n")
    cat("Message for selective inferece:\n")
    cat(sel.message)
  }


  pb$tick()

  ret_table_df = select(ret_table_df, -Component)
  ret_table_sel_df = select(ret_table_sel_df, -Component)

  ret_table_df = ret_table_df %>% mutate_if(is.numeric, function(x) round(x,4))
  ret_table_sel_df = ret_table_sel_df %>% mutate_if(is.numeric, function(x) round(x,4))

  return(list("plot" = decomp.plot, "table" = ret_table_df,
              "sel.plot" = decomp.plot.sel, "sel.table" = ret_table_sel_df,
              "sel.message" = sel.message))

}






