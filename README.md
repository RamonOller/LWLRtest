A nonparametric test for the association between longitudinal covariates and censored survival data
================

This is a document explaining and providing the code used in the paper (Oller and Gómez 2019). The methodology is illustrated with a dataset from the study Epidemiology of Diabetes Interventions and Complications (Sparling et al. 2006). The dataset (in SAS format) can be downloaded from the web of the [GW Biostatistics Center](https://biostatcenter.gwu.edu/people/research-faculty/john-m-lachin){target="_blank"}.

You can import the data set into **R** like this

``` r
library(haven)
icetdcfit <- read_sas("icetdcfit.sas7bdat")
```

We focus on the time-dependent covariate *Edic\_Hba*, i.e. the updated current mean of glycemia during the EDIC study. The individual trajectories of this covariate are ploted below,

``` r
library(ggplot2)
ggplot(icetdcfit,aes(x=tau,y=Edic_Hba))+geom_line(aes(group=id),col="steelblue")+xlab("\nDays")
```

![unnamed-chunk-2-1](https://user-images.githubusercontent.com/45238159/52239449-cd0bfe00-28ce-11e9-9dbf-e631de93c0a2.png)

<br><br>The survival outcomes of interest are the times of progression of retinopathy. Of the N=1316 subjects, 1085 event times are right-censored and 231 are interval-censored. The left- and right- endpoints are obtained as follows,

``` r
icetdc.id <-icetdcfit[!duplicated(icetdcfit$id), ]
left<-icetdc.id$t1
right<-icetdc.id$t2
right[is.na(icetdc.id$t2)]<-Inf
```

We now plot the NPMLE of the survival function. The **R** estimation procedure is commented in order to reduce the computation times, instead we provide the resulting object `sFit.icetdc`,

``` r
library(interval)
# sFit.icetdc<-icfit(left,right,control=icfitControl(maxit =10^6))
load("sFit.icetdc.rda")
plot(sFit.icetdc,XLAB="Days")
```

![unnamed-chunk-4-1](https://user-images.githubusercontent.com/45238159/52239549-07759b00-28cf-11e9-9a98-07bdae4f8a3a.png)

<br><br>Our methodology is an extension of the log-rank test statistic to provide evidence of a plausible association between a time-to-event outcome *T* and a time-dependent covariate *z*(*t*). The computation of the test statistic *LWLR* requires the exact values of *z*<sub>*i*</sub>(*t*) (*i* = 1, …, *n*) at the jumping points of the NPMLE of the survival function. Since these values are usually not known, we propose to use predicted values based on a fitted linear mixed model for *z*<sub>*i*</sub>(*t*). In the case of the *Edic\_Hba*, we assume a linear evolution in time and propose the following linear mixed model with random intercept and slope,

``` r
library(JM)
lmeFit.icetdc <- lme(Edic_Hba ~ tau, random = ~ tau | id, data = icetdcfit)
```

The test statistic is implemented with the function `LWLRtest`.

``` r
load("LWLRtest.rda")
```

The arguments are the following:

-   *L*: Numeric vector of the left endpoints of the censoring intervals.
-   *R*: Numeric vector of the right endpoints of the censoring intervals.
-   *lmeObject*: A *lme* object to compute imputed values for *z*<sub>*i*</sub>(*t*).
-   *timeVar*: The name of the time variable in *lmeObject*.
-   *Lin*: Logical vector: should *L* be included in the interval?
-   *Rin*: Logical vector: should *R* be included in the interval?
-   *rho*: A scalar parameter that controls the type of test.
-   *lambda*: A scalar parameter that controls the type of test.
-   *sFit*: A precalculated *icfit* object for increased computation speed.

In this illustration, we use a longitudinal log-rank test (*rho=0* and *lambda=0*),

``` r
LWLRtest(left,right,lmeFit.icetdc,timeVar="tau",sFit=sFit.icetdc)
```

    ## $W.standardized
    ## [1] 9.757449
    ## 
    ## $p.value
    ## [1] 0

and a longitudinal Wilcoxon-type test (*rho=1* and *lambda=0*),

``` r
LWLRtest(left,right,lmeFit.icetdc,timeVar="tau",sFit=sFit.icetdc,rho=1)
```

    ## $W.standardized
    ## [1] 9.635641
    ## 
    ## $p.value
    ## [1] 0

Both results show a high and positive association between the risk of having diabetic retinopathy and the longitudinal values of glycemia over time.

References
==========

Oller, Ramon, and Guadalupe Gómez. 2019. “A Nonparametric Test for the Association Between Longitudinal Covariates and Censored Survival Data.” *To Appear in Biostatistics*.

Sparling, Yvonne H, Naji Younes, John M Lachin, and Oliver M Bautista. 2006. “Parametric Survival Models for Interval-Censored Data with Time-Dependent Covariates.” *Biostatistics* 7 (4): 599–614.

