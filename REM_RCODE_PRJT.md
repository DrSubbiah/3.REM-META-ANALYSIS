<style type="text/css">
body {
text-align: justify
}
body, td {
   font-size: 20px;
}
code.r{
  font-size: 14px;
}
pre {
  font-size: 12px
}
</style>

# <span style="color:maroon"> Meta Analytic Approach

Combining information from multiple studies is one of the wide spread
statistical practices in medical, social, statistical genetics, clinical
trials, and epidemiological research. Most of the prominent such
approaches are Multi-centric studies, randomized clinical trial, and
meta-analysis.

Objective and purpose of each method may be different, yet underlying
statistical model is mostly confined to fixed or random effects model
(FEM or REM); this notion is presented in the following figure

![*Studies from same population or derived from different
populations*](C:/Users/Subbi.DESKTOP-L4TAD5C/Pictures/REM_FIG01.png)

This notes explains the idea behind Random Effects Model (REM) and the
classical (frequentist) procedures for estimation of parameters in the
model

# <span style="color:maroon"> Random Effects Model - REM

Let there be *k* independent studies with an effect parameter
*θ*<sub>*i*</sub>   *i* = 1, 2, ⋯*k* which is subjected to have a
sampling error *ϵ*<sub>*i*</sub>   *i* = 1, 2, ⋯*k*

We can **estimate** *θ*<sub>*i*</sub>   *i* = 1, 2, ⋯*k* from the sample
data; let us denote these estimated values as
*y*<sub>*i*</sub>   *i* = 1, 2, ⋯*k*; that is
$\\hat{\\theta_i}=y_i \~\~\~i=1,2,\\cdots k$

Modelling assumptions are,
*y*<sub>*i*</sub> ∼ *N*(*θ*<sub>*i*</sub>,*σ*<sub>*i*</sub><sup>2</sup>)   *i* = 1, 2, ⋯*k*
  
*θ*<sub>*i*</sub> ∼ *N*(*μ*,*τ*<sup>2</sup>)

*σ*<sub>*i*</sub><sup>2</sup> amounts sampling ***variability*** in the
effect size *θ*<sub>*i*</sub> estimate

*τ*<sup>2</sup> amounts ***variability*** between the effect size (among
the grouping or stratifying variables)

In almost all cases of *y*<sub>*i*</sub>, *σ*<sub>*i*</sub><sup>2</sup>
is assumed to be known. This is achieved by estimating
*σ*<sub>*i*</sub><sup>2</sup> using asymptotic approaches for
corresponding statistic when exact method may be computationally
intensive or impossible to obtain closed form expressions.

Classical methods for point estimation of *τ*<sup>2</sup> can be broadly
grouped as

-   method of moments estimator (MOM): Cochran’s ANOVA (CA), DerSimonian
    and Laird (DL), and Paule and Mandel (PM)

-   variance component type estimator (VCT): Hedges (HE) or Hedges and
    Olkinand; Hunter and Schmidt (HS)

-   model error variance estimators (MEV): Sidik and Jonkman (SJ)

-   likelihood estimators (LIK): Maximum likelihood (ML) and Restricted
    Maximum Likelihood (REML)

These methods are further grouped on their nature of handling **zero
cases** and computational aspects which may be of lesser interest in
this era, few studies have mentioned the advantages of non-iterative
methods

This classification is presented in Figure 1.

![Fig 1: Classical Approaches for estimators for between variance
(*τ*<sup>2</sup>)](C:/Users/Subbi.DESKTOP-L4TAD5C/Pictures/REM_FIG1.png)

Except MEV, some extensions or modifications can be found in the
literature that are listed as follows; figure 2 refers to this list of
methods

1.  2SCA: Two-step Cochran’s ANOVA
2.  2SDL: Two-step DerSimonian-Laird
3.  HM: Hartung and Makambi
4.  PDL: Positive DerSimonian and Laid
5.  BDL: Bootstrap DerSimonian and Laird
6.  APM: Alternate Paule and Mandel
7.  ASJ: Alternate Sidik and Jonkman
8.  AREML: Approximate restricted maximum

![Fig 2: Alternate / Derived forms of Classical Approaches for
estimators for between variance
(*τ*<sup>2</sup>)](C:/Users/Subbi.DESKTOP-L4TAD5C/Pictures/REM_FIG2.png)

These methods can further be grouped based on the truncation limit or
the **aberration caused by negative estimates** of *τ*<sup>2</sup>, as
well as the estimation approaches, whether they are iterative or not.
This is depicted in Figure 3.

![Fig 3: Aberration/Truncation/Iterative/Non-Iterative in Classical
Approaches for estimators for between
variance(*τ*<sup>2</sup>)](C:/Users/Subbi.DESKTOP-L4TAD5C/Pictures/REM_FIG3.png)

It can be noted that BDL is different from other iterative procedures.
BDL relies on mean value of finite number of replication whereas other
iterative procedures are estimating the parameter in each step and stop
at the optimum stage.

# <span style="color:maroon">Estimating Confidence Interval for *τ*<sup>2</sup>

The following methods are available for obtaining confidence interval
(CI) estimates for the between-study variance *τ*<sup>2</sup>; this is
depicted in Figure 4

1.  LIKE: likelihood estimators
2.  QSTAT: Methods based on Q statistics
3.  BOST: Bootstrap methods
4.  Prof_Like: Profile Likelihood
5.  L_Ratio: Likelihood ratio and Chi square based
6.  ML: Maximum Likelihood
7.  REML: Restricted Maximum Likelihood
8.  BT: Biggerstaff and Tweedie
9.  SJ: Sidik and Jonkman

![Fig 4: Classical CI for estimators for between variance
(*τ*<sup>2</sup>)](C:/Users/Subbi.DESKTOP-L4TAD5C/Pictures/REM_FIG4.png)

Based on the above approaches three set of parameters can be estimated

1.  Individual study-effects *θ*<sub>*i*</sub>  *i* = 1, 2, 3⋯⋯*k*

2.  Overall study-effect *μ*

3.  Between variance (*τ*<sup>2</sup>)
