{smcl}
{* *! version 1.0.0 12oct2017}{...}
{cmd:help car_sat}
{hline}

{title:Title}

    {hi:car_sat} {c -} Fully saturated linear regression with covariate-adaptive randomization standard error adjustment

{title:Syntax}
{p 8 17 2}
{cmd:car_sat} {it:dep_var} {it:treat_var} [{opt if}]{cmd:,}
 {opt strata:(strat_var)} [{it:options}]


{title:Description}

{pstd}
{cmd:car_sat} is a estimation command that estimates a "fully saturated" regression; that is, a regression of 
{it:dep_var} on all interaction terms between {it:treat_var} and {it:stratvar}. It then calculates
the estimated weighted coefficient of {it:treat_var} as well as its standard error adjusted to take into account 
potential bias from covariate-adaptive randomization, as per Bugni, Canay, and Shaikh (2017). 
{p_end}
{p 4 8}
Notes: {p_end}
{phang2}{cmd:.} The treatment variable must be such that the reference treatment (or control) takes on the value of zero.{p_end}
{phang2}{cmd:.} The strata variable must be numeric. {p_end}
{phang2}{cmd:.} Strata with missing cells are dropped in the estimation. {p_end}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt table}}display full regression table for the "fully saturated" regression, with robust (HC1) standard errors. {p_end}
{synopt:{opt covars(varlist)}}includes covariates into the regression.{p_end}

{title:Saved Results}

{phang}
{cmd:car_sat} saves the following in {cmd:e()}:{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}

{synopt:{cmd:e(df_r)}}degrees of freedom{p_end}
{synopt:{cmd:e(N)}}sample size{p_end}
{p2colreset}{...}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:car_sat}{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{p2colreset}{...}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(V_H)}}estimator of the heterogeneity term in the asymptotic variance (see Bugni, Canay, and Shaikh (2017), eqns. 11 and 15){p_end}
{synopt:{cmd:e(V_hc)}}HC-robust estimator of the asymptotic variance{p_end}

{p2colreset}{...}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{title:Author}

{phang}
Joe Long{p_end}
{phang}
jlong@u.northwestern.edu
{p_end}

{title:References}

{phang}
Bugni, Federico A., Ivan A. Canay, and Azeem M. Shaikh. (2017) "Inference under Covariate-Adaptive Randomization with Multiple Treatments." 
{p_end}
