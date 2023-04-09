# repDiagnosis
Repository for R package `repDiagnosis` and R shiny app for live replication diagnosis tools.
 

**Reference.** Please use the following citation if you use our methods, datasets, or softwares for analyzing replication studies.

```

```

**Rshiny app.** To quickly get started with our method, we recommend you visit our Rshiny app and play with our preloaded  clean datasets. Our Rshiny app is also helpful if you want to quickly diagnose your replication study pairs with easy-to-use features, such as live displaying and visualizing data, specifying your study, automatically checking regression in two studies, etc., or if you want to probe the generalizability of your single study (this function is not in our package and paper).

**Dataset collection.** In a sibling Github repository [awesome-replicability-data](https://github.com/ying531/awesome-replicability-data), we collect individual-level data for 11 original-replication study pairs from publicly available sources, and 5 one-sided replication studies. Visit our data collection to explore these datasets and our pre-processed clean versions.

**Documentation website.** For more on the package and the Rshiny app, see the associated [documentation]() website.

## Install R Package `repDiagnosis`
  
1. Install the [devtools](https://github.com/hadley/devtools) package using `install.packages("devtools")`.
2. Install the latest development version using `devtools::install_github("ying531/repDiagnosis")`.
 

## Example 

We show how to conduct the diagnosis for the replication study on EMDR and misinformation, which decomposes the discrepancy between original and replication study estimates into several interpretable pieces (for details and cleaned datasets, see Data pair 3 in our dataset collection). 

**Load datasets.** First, the paired datasets can be directly loaded from the package. Use `data()` to load them as `org.dat` for original study and `rep.dat` for replication study. You can also use your own datasets, but make sure their column names are consistent. 

```R
data("emdr_misinfo")
```

**Specify regression.** In this study, the regression of interest is `totalcorrect` on the treatment indicator named `condition`. Specify the regression formula as you do for linear models, and specify the treatment variable name.

```R
analysis_formula = "totalcorrect ~ condition"
treatment_name = "condition"
```

**Specify covariates and mediators.** Suppose you are interested in the heterogeneity between two studies captured by covariates `age` and `gender`. Also, you wonder whether the mediators `postvividness` and `postemotionality` capture any more heterogeneity. Specify these quantities as vectors. 

```R
covariates = c("age", "gender")
mediators = c("postvividness", "postemotionality")
```

**Specify group IDs for clustered experiments.** This is not needed for this dataset. If you like, you can also specify group IDs as `X`, which is a redundant variable that indicates the participant ID. 
  

**Run the decomposition.** Finally, call our analysis function to run the diagnosis. It returns a summary table and a decomposition plot similar to our paper. By default, it also computes post-selective CIs that are valid at 1-`alpha` level after adjusting for the publication bias (i.e., only publishing p-values smaller than `pub_pvalue_threshold` whose default value is `0.05`). 

```R
results = run_diagnosis(
  org_df = org.dat, # data frame for original study
  rep_df = rep.dat, # data frame for replication study 
  analysis_formula = analysis_formula, # regression formula 
  treatment_name = treatment_name, # column name of treatment indicator
  id_name_1 = "X", # cluster ID for original study
  id_name_2 = "X", # cluster ID for replication study 
  covariates = covariates, # covariates to balance
  mediators = mediators, # mediators to balance
  alpha = 0.1, # confidence level for all CIs
  verbose = TRUE, # printing the tables and plots
  if_selective = TRUE, # conduct the selective inference CIs
  pub_pvalue_threshold = 0.05 # p-value threshold for publication bias
  )
```

If `verbose=TRUE`, it displays the following summary table:

```R
Summary of decomposition for replication diagnosis:
  
              Estimate Std. Error     t-stat     Pr(>|z|)
Original   -1.04878049  0.2634718 -3.9806177 6.873641e-05
Observed   -1.13211382  0.3361254 -3.3681297 7.567998e-04
Covariates  0.08133811  0.1487237  0.5469076 5.844422e-01
Mediators   0.07448551  0.1233981  0.6036197 5.460965e-01
Residual   -1.28793744  0.3984015 -3.2327622 1.225996e-03
```

and the following decomposition plot:
   
<img src="https://raw.githubusercontent.com/ying531/archiv/main/default_plot.png"  width="600" height="300">

If `if_selective=TRUE`, it runs selective inference procedures. It will additionally display the following summary table:

```R
Summary of post-selective decomposition for replication diagnosis:

                  low       high    estimate
Observed   -1.6846210 -0.5446210 -1.13462105
Covariates -0.1631354  0.3268646        -Inf
Mediators  -0.1334761  0.2765239  0.07652393
Residual   -1.9459682 -0.6159682 -1.28596819
```

and the following decomposition plot that adjusts for publication bias:
  

<img src="https://raw.githubusercontent.com/ying531/archiv/main/default_plot_selective.png"  width="600" height="300">

This plot does not differ much from the first because the original p-value is tiny, hence adjusting for publication bias does not change things too much.