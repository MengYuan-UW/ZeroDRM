# 
R packge for the inference under the density ratio model with semicontinous data
# R package ZeroDRM
This package is developed based on [Yuan et al. (2021a)](https://doi.org/10.1007/s10463-021-00804-4) and [Yuan et al. (2021b)](https://arxiv.org/abs/2106.02741). 
This package provides the estimators and confidence intervals for the functionals and Gini indices of two semicontibous populations under the density ratio models (DRMs).


### Table of Contents
**[Installation](#installation)**<br>
**[Functions](#functions)**<br>
**[Usage](#usage)**<br>
**[References](#references)**<br>
## Installation
Install this package from Github with 
```r
# If the package "devtool" has not been installed, install it fisrt with 
#install.packages("devtools")
library(devtools)
devtools::install_github("MengYuan-UW/ZeroDRM")
```
## Functions
This package contains the following functions:
- The Density Ratio Model (`DRM`): a function to fit the DRM.
- The inference on the functional (`ZeroFunc'): a function to estimate the functionals and construct the confidence interval under the DRM.
- The inference on the mean ratio (`MeanRatio'): a function to estimate the mean ratio of two populations and construct the confidence interval under the DRM.
- The inference on the Gini indices (`ZeroGini'): a function to estimate the Gini indices and construct the confidence intervals under the DRM.
- The inference on the difference of two Gini indices (`DiffGini'): a function to estimate the difference of Gini indices of two populations and construct the confidence intervals under the DRM.

## Usage
We provide two examples.
- Example 1: inference on population means and mean ratio
```r
library("YoudenDRM")
#Data generation function
DataGen = function(nu,para,size){
  sample = rbinom(size,1,nu)
  #since nu is the probability that x = 0
  n0 = sum(sample==1)
  n1 = sum(sample==0)
  data = c(rep(0,n0),
           rlnorm(n1,meanlog = para[1],sdlog = para[2]))
  data
}
# zero proportion
nu = c(0.3,0.5)
# mean for lognormal distribution in log scale
a = c(1/3,2/3)

set.seed(123)
x = DataGen(nu[1],c(a[1],1),50)
y = DataGen(nu[2],c(a[2],1),50)

# Now move to the inference part

# define the arguments in the 'ZeroFunc'
# (1) define the basis function q(x) = log(x)
basis = function(x){log(x)}
# (2) u function used to define the psi = c(mu_0, mu_1), where mu_0 and mu_1 are population means for x and y, respectively
Ufun = function(x,nu,theta){
  u1 = (1-nu[1])*x
  u2 = (1-nu[2])*x*exp(theta%*%c(1,log(x)))
  c(u1,u2)
}
# (3) the first derative of u functin with respect to parametres nu
Pnu = function(x,nu,theta){
  pu1 = c(-x,0)
  pu2 = c(0,-x*exp(theta%*%c(1,log(x))))
  rbind(pu1,pu2)
}
# (4)the first derative of u functin with respect to parametres theta
Ptheta = function(x,nu,theta){
  pu1 = c(0,0)
  pu2 = (1-nu[2])*x*c(exp(theta%*%c(1,log(x))))*c(1,log(x))
  rbind(pu1,pu2)
}

# get the estimate of psi and asymptotic variance-covaraince matrix
ZeroFunc(x,y,basis,method="optimal",Ufun,Pnu,Ptheta)

# otain the estimate of the mean ratio mu_1/mu_0 and asymptotic standard deviation
# construct the 95% confidence intervals using difference methods under the density ratio model
MeanRatio(x,y,basis,CItype = "None")
MeanRatio(x,y,basis,CItype = "NA-DRM")
MeanRatio(x,y,basis,CItype = "log-DRM")
```
- Example 2: inference on the Gini indices of two populations and their difference 
```r
# get the data 
data(Ilocos,package = "ineq")
x = Ilocos$income[Ilocos$urban == "urban"&Ilocos$pro == "Pangasinan"]
y = Ilocos$income[Ilocos$urban == "rural"&Ilocos$pro == "Pangasinan"]

# define the basis function q(x)
basis = function(x){log(x)}

# get the estimate of Gini indices and 95% Wald-type confidence intervals
ZeroGini(x,y,basis,CItype = "NA-DRM")

# get the estimate of the difference of Gini indices and 95% Wald-type confidence intervals
DiffGini(x,y,basis,CItype = "NA-DRM")
```

## References

[Yuan, M., Wang, C., Lin, B., and Li, P. (2021a). "Semiparametric inference on general functionals of two semicontinuous populations." Annals of the Institute of Statistical Mathematics. In press.](https://doi.org/10.1007/s10463-021-00804-4)

[Yuan, M., Li, P., and Wu, C. (2021b). "Semiparametric inference on Gini indices of two semicontinuous populations under density ratio models." arXiv:2106.02741.](https://arxiv.org/abs/2106.02741)
