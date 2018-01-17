Meta Analysis\_Two Approaches
================

AIM
---

To explain Meta analytic approach in random effects model (REM) framework using 1. summary statistics (SSM) and 2. Fully Bayesian (FBM) method.

**Ref: Meta-analysis: formulating, evaluating, combining, and reporting Normand,ST Statist. Med. 1999**

REM here indicates combining multiple (k) studies each with a study effect *y*<sub>*i*</sub> where *i* = 1, 2, 3.....*k*

AIM 1 - SSM:
------------

Underlying model is <br><br>*y*<sub>*i*</sub> ∼ *N*(*μ*<sub>*i*</sub>, *v*<sub>*i*</sub><sup>2</sup>) <br><br> *μ*<sub>*i*</sub> ∼ *N*(*d*, *σ*<sub>*μ*</sub><sup>2</sup>). <br><br>Assumption is *v*<sub>*i*</sub><sup>2</sup> are known (estimated from data, most of the cases asymptotic variance) <br><br>quantities of interest are *μ* and *σ*<sub>*μ*</sub><sup>2</sup>.

Step 1:Summary Statistics
-------------------------

We use metafor pacakge to estimate the summary measure *y*<sub>*i*</sub> and asymptotic variance *v*<sub>*i*</sub><sup>2</sup>

In this note *y*<sub>*i*</sub> is log odds ratio of study i ∀1 ≤ *i* ≤ *k* and hence *v*<sub>*i*</sub><sup>2</sup> = *a*<sup>−1</sup> + *b*<sup>−1</sup> + *c*<sup>−1</sup> + *d*<sup>−1</sup> where *a*, *b*, *c*, *d* are cell counts of given 2 × 2 table

``` r
#Summary Statistics Approach
#Step 1
ami=c(0,    0,  0,  0,  1,  0,  0,  2,  0,  0,  0,  0,  0,  1,  0,  0,  1,  1,  0,  0,  0)
bmi=c(100,  100,    28, 14, 19, 80, 54, 148,    60, 15, 10, 15, 88, 141,    136,    204,    97, 23, 25, 41, 21)
cmi=c(0,    0,  1,  0,  1,  0,  0,  0,  2,  0,  0,  1,  2,  2,  0,  2,  2,  1,  0,  0,  0)
dmi=c(100,  101,    36, 20, 18, 80, 49, 150,    58, 16, 10, 19, 86, 137,    131,    182,    97, 25, 25, 67, 16)
require(metafor)
y=escalc(ai=ami,bi=bmi,ci=cmi,di=dmi,measure = "OR",add = 10^(-8), to="only0" ,  drop00=FALSE)
#End of step 1----------------------------------------------------------
```

Step 2: Stan code for REM Normal-Normal Model
=============================================

``` stan
data {
  int<lower=0> N; 
  real y[N]; 
  real<lower=0> v1[N]; 
} 
parameters {
  real d; 
  real<lower=0> sigma_mu; 
  vector[N] mu;
} 
transformed parameters {
  real<lower=0> sigmasq_mu; 
  sigmasq_mu=sigma_mu*sigma_mu;
} 
model {
  y ~ normal(mu, v1);
  mu ~ normal(d, sigma_mu);
  d ~ normal(0, 1E3); 
  sigmasq_mu ~ inv_gamma(1E-3, 1E-3); 
}
```

``` r
require(rstan)
```

    ## Loading required package: rstan

    ## Loading required package: ggplot2

    ## Loading required package: StanHeaders

    ## rstan (Version 2.17.2, GitRev: 2e1f913d3ca3)

    ## For execution on a local, multicore CPU with excess RAM we recommend calling
    ## options(mc.cores = parallel::detectCores()).
    ## To avoid recompilation of unchanged Stan programs, we recommend calling
    ## rstan_options(auto_write = TRUE)

``` r
yi=as.vector(y$yi);v1=sqrt(y$vi) #sq.root is for SD parameter in Normal in Stan
N=length(ami)
MA_SS_data <- list(y=yi,v=v1,N=N)
MA_SS_init <- function(){list(
  d = 0,
  mu = rnorm(N,0,1),
  sigma_mu = 1
)}
mod1 <- sampling(MA_SS_tbyt,data=MA_SS_data, init = MA_SS_init,
                 control = list(adapt_delta=0.99,max_treedepth=15),iter=2000,chains=4) 
```

    ## 
    ## SAMPLING FOR MODEL '9619ccf3b9cf03e6a4a776ff2a33354f' NOW (CHAIN 1).
    ## 
    ## Gradient evaluation took 0 seconds
    ## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Adjust your expectations accordingly!
    ## 
    ## 
    ## Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Iteration: 2000 / 2000 [100%]  (Sampling)
    ## 
    ##  Elapsed Time: 0.781 seconds (Warm-up)
    ##                1.882 seconds (Sampling)
    ##                2.663 seconds (Total)
    ## 
    ## 
    ## SAMPLING FOR MODEL '9619ccf3b9cf03e6a4a776ff2a33354f' NOW (CHAIN 2).
    ## 
    ## Gradient evaluation took 0 seconds
    ## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Adjust your expectations accordingly!
    ## 
    ## 
    ## Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Iteration: 2000 / 2000 [100%]  (Sampling)
    ## 
    ##  Elapsed Time: 0.619 seconds (Warm-up)
    ##                0.953 seconds (Sampling)
    ##                1.572 seconds (Total)
    ## 
    ## 
    ## SAMPLING FOR MODEL '9619ccf3b9cf03e6a4a776ff2a33354f' NOW (CHAIN 3).
    ## 
    ## Gradient evaluation took 0 seconds
    ## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Adjust your expectations accordingly!
    ## 
    ## 
    ## Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Iteration: 2000 / 2000 [100%]  (Sampling)
    ## 
    ##  Elapsed Time: 0.753 seconds (Warm-up)
    ##                0.539 seconds (Sampling)
    ##                1.292 seconds (Total)
    ## 
    ## 
    ## SAMPLING FOR MODEL '9619ccf3b9cf03e6a4a776ff2a33354f' NOW (CHAIN 4).
    ## 
    ## Gradient evaluation took 0 seconds
    ## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Adjust your expectations accordingly!
    ## 
    ## 
    ## Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Iteration: 2000 / 2000 [100%]  (Sampling)
    ## 
    ##  Elapsed Time: 0.517 seconds (Warm-up)
    ##                1.241 seconds (Sampling)
    ##                1.758 seconds (Total)

    ## Warning: There were 7 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See
    ## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

    ## Warning: There were 4 chains where the estimated Bayesian Fraction of Missing Information was low. See
    ## http://mc-stan.org/misc/warnings.html#bfmi-low

    ## Warning: Examine the pairs() plot to diagnose sampling problems

``` r
fitNor_Nor_summary<-as.data.frame(summary(mod1, pars = c("d","sigmasq_mu"), probs = c(0.025, 0.975))$summary)
ss_re=round(fitNor_Nor_summary,4)
```

AIM 2 - FBM:
------------

**Additional Reference:Smith et al 1995 Bayesian Approaches to Random-Effects Meta analysis: A Comparative Study**

Underlying model is <br>*r*<sub>*i*</sub><sup>*C*</sup> ∼ *B**i**n**o**m**i**a**l*(*p*<sub>*i*</sub><sup>*C*</sup>, *n*<sub>*i*</sub><sup>*C*</sup>) <br><br>*r*<sub>*i*</sub><sup>*T*</sup> ∼ *B**i**n**o**m**i**a**l*(*p*<sub>*i*</sub><sup>*T*</sup>, *n*<sub>*i*</sub><sup>*T*</sup>) <br><br>$logit(p^{C}\_{i}) = \\mu\_{i} - \\frac{\\delta\_{i}}{2}$<br><br> $logit(p^{T}\_{i}) = \\mu\_{i} + \\frac{\\delta\_{i}}{2}$<br><br>*δ*<sub>*i*</sub> ∼ *N**o**r**m**a**l*(*d*, *σ*<sup>2</sup>)

Here, *r*<sub>*i*</sub><sup>*C*</sup> indicates the number of successes in Control group arising from *n*<sub>*i*</sub><sup>*C*</sup> cases assumed to have probability of *p*<sub>*i*</sub><sup>*C*</sup> <br><br>Similarly *r*<sub>*i*</sub><sup>*T*</sup> can be defined for treatment group <br><br>Quantities of interest are *d* and *σ*<sup>2</sup>.

Fully Bayesian Method
=====================

**Stan Model**

``` stan
data {
  int<lower=0> N; 
  int<lower=0> nt[N]; 
  int<lower=0> rt[N]; 
  int<lower=0> nc[N]; 
  int<lower=0> rc[N]; 
} 
parameters {
  real d; 
  real<lower=0> sigma_delta; 
  vector[N] mu;
  vector[N] delta;
  } 
transformed parameters {
  real<lower=0> sigmasq_delta; 
  sigmasq_delta = sigma_delta*sigma_delta; 
 
} 
model {
  rt ~ binomial_logit(nt, mu+delta/2);
  rc ~ binomial_logit(nc, mu-delta/2);
  delta  ~ normal(d, sigma_delta); 
  mu ~ normal(0, 2);
  d ~ normal(0, sqrt(10)); 
  sigmasq_delta ~ inv_gamma(3, 1); 
  }
```

``` r
require(rstan)
nt=ami+bmi;nc=cmi+dmi
N=length(ami)

MA_FB_data <- list(nc=nc,rc=cmi,nt=nt,rt=ami,N=N)
MA_FB_init <- function(){list(
  d = 0,
  mu = rnorm(N,0,1),
  delta= rnorm(N,0,1),
  sigma_delta = 1
)}
fitMA_FB <- stan(file="MA_FB_tbyt.stan",
                   data=MA_FB_data, init = MA_FB_init,
                   control = list(adapt_delta=0.99,max_treedepth=15),iter=2000,chains=4)
```

    ## In file included from C:/Users/Lapilluz2/Documents/R/win-library/3.4/BH/include/boost/config.hpp:39:0,
    ##                  from C:/Users/Lapilluz2/Documents/R/win-library/3.4/BH/include/boost/math/tools/config.hpp:13,
    ##                  from C:/Users/Lapilluz2/Documents/R/win-library/3.4/StanHeaders/include/stan/math/rev/core/var.hpp:7,
    ##                  from C:/Users/Lapilluz2/Documents/R/win-library/3.4/StanHeaders/include/stan/math/rev/core/gevv_vvv_vari.hpp:5,
    ##                  from C:/Users/Lapilluz2/Documents/R/win-library/3.4/StanHeaders/include/stan/math/rev/core.hpp:12,
    ##                  from C:/Users/Lapilluz2/Documents/R/win-library/3.4/StanHeaders/include/stan/math/rev/mat.hpp:4,
    ##                  from C:/Users/Lapilluz2/Documents/R/win-library/3.4/StanHeaders/include/stan/math.hpp:4,
    ##                  from C:/Users/Lapilluz2/Documents/R/win-library/3.4/StanHeaders/include/src/stan/model/model_header.hpp:4,
    ##                  from file1218755f5b6a.cpp:8:
    ## C:/Users/Lapilluz2/Documents/R/win-library/3.4/BH/include/boost/config/compiler/gcc.hpp:186:0: warning: "BOOST_NO_CXX11_RVALUE_REFERENCES" redefined
    ##  #  define BOOST_NO_CXX11_RVALUE_REFERENCES
    ##  ^
    ## <command-line>:0:0: note: this is the location of the previous definition
    ## cc1plus.exe: warning: unrecognized command line option "-Wno-ignored-attributes"
    ## 
    ## SAMPLING FOR MODEL 'MA_FB_tbyt' NOW (CHAIN 1).
    ## 
    ## Gradient evaluation took 0 seconds
    ## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Adjust your expectations accordingly!
    ## 
    ## 
    ## Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Iteration: 2000 / 2000 [100%]  (Sampling)
    ## 
    ##  Elapsed Time: 12.824 seconds (Warm-up)
    ##                13.042 seconds (Sampling)
    ##                25.866 seconds (Total)
    ## 
    ## 
    ## SAMPLING FOR MODEL 'MA_FB_tbyt' NOW (CHAIN 2).
    ## 
    ## Gradient evaluation took 0 seconds
    ## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Adjust your expectations accordingly!
    ## 
    ## 
    ## Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Iteration: 2000 / 2000 [100%]  (Sampling)
    ## 
    ##  Elapsed Time: 13.07 seconds (Warm-up)
    ##                7.393 seconds (Sampling)
    ##                20.463 seconds (Total)
    ## 
    ## 
    ## SAMPLING FOR MODEL 'MA_FB_tbyt' NOW (CHAIN 3).
    ## 
    ## Gradient evaluation took 0 seconds
    ## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Adjust your expectations accordingly!
    ## 
    ## 
    ## Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Iteration: 2000 / 2000 [100%]  (Sampling)
    ## 
    ##  Elapsed Time: 14.447 seconds (Warm-up)
    ##                2.845 seconds (Sampling)
    ##                17.292 seconds (Total)
    ## 
    ## 
    ## SAMPLING FOR MODEL 'MA_FB_tbyt' NOW (CHAIN 4).
    ## 
    ## Gradient evaluation took 0 seconds
    ## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Adjust your expectations accordingly!
    ## 
    ## 
    ## Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Iteration: 2000 / 2000 [100%]  (Sampling)
    ## 
    ##  Elapsed Time: 18.834 seconds (Warm-up)
    ##                13.777 seconds (Sampling)
    ##                32.611 seconds (Total)

    ## Warning: There were 479 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See
    ## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

    ## Warning: There were 4 chains where the estimated Bayesian Fraction of Missing Information was low. See
    ## http://mc-stan.org/misc/warnings.html#bfmi-low

    ## Warning: Examine the pairs() plot to diagnose sampling problems

``` r
fitMA_FB_summary<-as.data.frame(summary(fitMA_FB, pars = c("d","sigmasq_delta"), probs = c(0.025, 0.975))$summary)
fb_re=round(fitMA_FB_summary,4)
```

``` r
library(knitr)
kable(ss_re,caption = "Summary Statistics Method")
```

|             |     mean|  se\_mean|      sd|     2.5%|   97.5%|    n\_eff|    Rhat|
|-------------|--------:|---------:|-------:|--------:|-------:|---------:|-------:|
| d           |  -0.4161|    0.0762|  0.6247|  -1.4972|  0.7776|   67.2033|  1.0312|
| sigmasq\_mu |   0.0401|    0.0247|  0.2514|   0.0003|  0.2953|  103.6459|  1.0272|

``` r
kable(fb_re,caption = "Fully Bayesian Method")
```

|                |     mean|  se\_mean|      sd|     2.5%|   97.5%|    n\_eff|    Rhat|
|----------------|--------:|---------:|-------:|--------:|-------:|---------:|-------:|
| d              |  -0.5447|    0.0601|  0.4162|  -1.3973|  0.1389|   47.9262|  1.0543|
| sigmasq\_delta |   0.0195|    0.0054|  0.0622|   0.0003|  0.1748|  132.8666|  1.0282|

This demonstrates the two MA approaches. Summary statistics requires continuity corrections for zero cells, which requires a careful investigation; where as this is completely alleviated when FBM is used.

Still a careful specification of priors (especially for between variance) is always a concern in Bayesian analysis
