#' @title Inference on Functionals
#' @description Estimate the functionals as well as construct the confidence interval for functionals based on the observed samples.
#'
#' @param x observed sample from \eqn{F_0}.
#' @param y observed sample from \eqn{F_1}.
#' @param basis pre-sepecified basis function \eqn{\boldsymbol{q}(x)} in the exponential term.
#' It should return a \eqn{d}-dimensional vector.
#' For example, \code{basis = function(x){c(x,log(x))}}.
#' @param method the method used to fit the DRM. The options include \code{"glm"}, \code{"multinom"} and \code{"optimal"}.
#' The default setting is \code{"optimal"}. See \code{help(DRM)} for `Details'.
#' @param Ufun the funtion \eqn{\boldsymbol{u}(x;\boldsymbol{\nu},\boldsymbol{\theta})} to specify the \eqn{p}-demensional functional \eqn{\boldsymbol{\psi}}
#' with zero proportion parameter \eqn{\boldsymbol{\nu} = (\nu_0,\nu_1)^\top} and DRM parameter \eqn{\boldsymbol{\theta} = (\alpha,\boldsymbol{\beta})^\top}.
#' It should return a \eqn{p}-dimensional vector. See `Details`.
#' @param Pnu a function which specifies the first derivative of \eqn{\boldsymbol{u}(x;\boldsymbol{\nu},\boldsymbol{\theta})} with respect to \eqn{\boldsymbol{\nu}}.
#' It should return a \eqn{p\times 2} matrix. See `Details`.
#' @param Ptheta a function which specifies the first derivative of function \eqn{\boldsymbol{u}(x;\boldsymbol{\nu},\boldsymbol{\theta})} with respect to \eqn{\boldsymbol{\theta}}.
#'It should return a \eqn{p\times (d+1)} matrix, where \eqn{d} is the dimensional of \eqn{\boldsymbol{\theta}}.
#'See `Details`.
#'
#' @details \insertCite{yuan2020zeros;textual}{ZeroDRM} considered a class of general functionals \eqn{\boldsymbol{\psi}} of dimension \eqn{p}, defined as
#' \deqn{\boldsymbol{\psi}= \int_0^\infty \boldsymbol{u}(x;\boldsymbol{\nu},\boldsymbol{\theta}) dG_0(x),}
#' where \eqn{\nu=(\nu_0,\nu_1)^\top}, \eqn{\boldsymbol{\theta} = (\alpha,\boldsymbol{\beta})^\top},
#' and  \eqn{\boldsymbol{u}(x;\boldsymbol{\nu},\boldsymbol{\theta})=\left(u_1(x;\boldsymbol{\nu},\boldsymbol{\theta}),\ldots,u_p(x;\boldsymbol{\nu},\boldsymbol{\theta})\right)^\top} is a given \eqn{(p\times 1)}-dimensional function.
#'
#' @details \code{Pnu} refers to the the first derivative of \eqn{\boldsymbol{u}(x;\boldsymbol{\nu},\boldsymbol{\theta})} with respect to \eqn{\boldsymbol{\nu}}:
#' \deqn{\frac{\partial \boldsymbol{u}(x;\boldsymbol{\nu},\boldsymbol{\theta})}{\partial \boldsymbol{\nu}}
#' = \left(\begin{array}{cc}
#' \frac{\partial u_1(x;\boldsymbol{\nu},\boldsymbol{\theta})}{\partial \nu_0}&\frac{\partial u_1(x;\boldsymbol{\nu},\boldsymbol{\theta})}{\partial \nu_1}\\
#' \vdots&\vdots\\
#' \frac{\partial u_p(x;\boldsymbol{\nu},\boldsymbol{\theta})}{\partial \nu_0}&\frac{\partial u_p(x;\boldsymbol{\nu},\boldsymbol{\theta})}{\partial \nu_1}
#' \end{array}\right).}
#'
#' @details \code{Ptheta} refers to the the first derivative of \eqn{\boldsymbol{u}(x;\boldsymbol{\nu},\boldsymbol{\theta})} with respect to \eqn{\boldsymbol{\theta}}:
#' \deqn{\frac{\partial \boldsymbol{u}(x;\boldsymbol{\nu},\boldsymbol{\theta})}{\partial \boldsymbol{\theta}}
#' = \left(\begin{array}{cc}
#' \frac{\partial u_1(x;\boldsymbol{\nu},\boldsymbol{\theta})}{\partial \alpha}&\frac{\partial u_1(x;\boldsymbol{\nu},\boldsymbol{\theta})}{\partial \boldsymbol{\beta}}\\
#' \vdots&\vdots\\
#' \frac{\partial u_p(x;\boldsymbol{\nu},\boldsymbol{\theta})}{\partial \alpha}&\frac{\partial u_p(x;\boldsymbol{\nu},\boldsymbol{\theta})}{\partial \boldsymbol{\beta}}
#' \end{array}\right).}
#' @return The function returns a list containing the following components:
#' \itemize{
#' \item \code{estimate}: the estimates of functional \eqn{\boldsymbol{\psi}}.
#' \item \code{AVar}: the estimate of asymptotic variance-covariance matrix \eqn{\Gamma}.
#' The explicit form of \eqn{\Gamma} van be found in \insertCite{yuan2020zeros;textual}{ZeroDRM}.
#' }
#'
#' @import nnet
#' @import stats
#' @importFrom Rdpack reprompt
#'
#'
#' @references
#'
#' \insertRef{yuan2020zeros}{ZeroDRM}
#'
#' @export
ZeroFunc = function(x,y,basis,method="optimal",Ufun,Pnu,Ptheta){
  n0 = length(x)
  n1 = length(y)
  n = n0+n1
  nu0 = sum(x == 0)/n0
  nu1 = sum(y == 0)/n1
  w = n0/n
  delta = w*(1-nu0) + (1-w)*(1-nu1)
  rho = (1-w)*(1-nu1)/delta

  S01 = x[x!=0]
  S11 = y[y!=0]
  coef = drmCoef(S01,S11,basis,method)

  t = c(S01,S11)
  Qt = cbind(1,matrix(basis(t),nrow = length(t),byrow=F))

  omega = c(exp(Qt%*%coef))
  P0 = 1/(length(S01)+length(S11)*omega)

  h = 1 + rho*(omega-1)
  h1 = rho*omega/h
  h0 = (1-rho)/h

  e = as.matrix(c(1,rep(0,ncol(Qt)-1)))

  U = sapply(1:length(t),function(i){Ufun(t[i],c(nu0,nu1),coef)})

  var1 = U%*%(t(U)*P0/h)/delta

  psi = U%*%P0
  var2 = psi%*%t(psi)/delta

  M1 = matrix(sapply(1:length(t), function(i){Pnu(t[i],c(nu0,nu1),coef)})%*%P0, nrow = 2)
  Utheta = matrix(sapply(1:length(t), function(i){Ptheta(t[i],c(nu0,nu1),coef)})%*%P0, nrow = length(coef))
  M2 = Utheta%*%e - rho*psi
  M3 = Utheta - U%*%(Qt*h1*P0)

  AnuInverse = matrix(c((nu0*(1-nu0))/w,0,0,(nu1*(1-nu1))/(1-w)),ncol = 2)
  Atheta = delta*(1-rho)*((t(Qt))%*%(Qt*P0*h1))

  var3 = M1%*%AnuInverse%*%t(M1) - M2%*%t(M2)/(delta*rho*(1-rho)) + M3%*%solve(Atheta)%*%t(M3)

  V = var1 - var2 + var3
  list(estimate = c(psi), AVar = V)
}

#' @title Inference on mean ratio
#' @description Estimate the mean ratio and construct the confidence interval based on the observed samples.
#'
#' @param x observed sample from \eqn{F_0}.
#' @param y observed sample from \eqn{F_1}.
#' @param basis pre-sepecified basis function \eqn{\boldsymbol{q}(x)} in the exponential term.
#' It should return a \eqn{d}-dimensional vector.
#' For example, \code{basis = function(x){c(x,log(x))}}.
#' @param method the method used to fit the DRM. The options include \code{"glm"}, \code{"multinom"} and \code{"optimal"}.
#' The default setting is \code{"optimal"}. See \code{help(DRM)} for `Details'.
#' @param CItype the method to be used for confidence interval construction. See `Details'.
#' @param tau the significance level. The defualt value is \eqn{0.05}.
#'
#' @details Let \eqn{\mu_i = \int_0^{\infty} x dF_i(x)} be the mean of populatio \eqn{i}.
#' Then we define the functional \eqn{\boldsymbol{\psi}} in \insertCite{yuan2020zeros;textual}{ZeroDRM} as
#'  \eqn{\boldsymbol{\psi} = (\mu_0,\mu_1)^\top}.
#'
#'Further let the function \eqn{g(\cdot)} in  \insertCite{yuan2020zeros;textual}{ZeroDRM} as
#'\eqn{g(x_1,x_2) = x_2/x_1}. We then have \eqn{g(\boldsymbol{\psi}) = \mu_1/\mu_0} is the mean ratio of two populations.
#'
#'@details The argument \code{CItype} refers to difference confidence intervals for the mean ratio \eqn{\mu_1/\mu_0}
#'under the DRMs in \insertCite{yuan2020zeros;textual}{ZeroDRM}.
#'\itemize{
#'\item  \code{"None"}: no confidence intervals for the difference is constructed;
#'\item  \code{"NA-DRM"}: the Wald-type confidence intervals based on the normal approximation of \eqn{\mu_1/\mu_0};
#'\item  \code{"log-DRM"}: the Wald-type confidence intervals based on the normal approximation of \eqn{\log(\mu_1/\mu_0)}.
#'}
#'
#'
#' @return The function returns a vector containing the following components:
#' \itemize{
#' \item \code{estimate}: the estimate of the difference \eqn{\mathcal{G}_0 - \mathcal{G}_1}.
#' \item \code{ASD}: the asymptotic standard deviation (ASD) of the estimator of the difference.
#' \item \code{lower bound}: the lower bound of the confidence interval for the difference.
#' \item \code{upper bound}: the upper bound of the confidence interval for the difference.
#' }
#' When \code{CItype = "None"}, it only returns the estimate and ASD of the difference.
#'
#'
#' @import nnet
#' @import stats
#' @importFrom Rdpack reprompt
#'
#'
#' @references
#'
#' \insertRef{yuan2020zeros}{ZeroDRM}
#'
#' @export
MeanRatio = function(x,y,basis,method="optimal",CItype, tau = 0.05){

  Ufun = function(x,nu,theta){
    u1 = (1-nu[1])*x
    u2 = (1-nu[2])*x*exp(theta%*%c(1,basis(x)))
    c(u1,u2)
  }

  Pnu = function(x,nu,theta){
    pu1 = c(-x,0)
    pu2 = c(0,-x*exp(theta%*%c(1,basis(x))))
    rbind(pu1,pu2)
  }

  Ptheta = function(x,nu,theta){
    pu1 = c(0,0)
    pu2 = (1-nu[2])*x*c(exp(theta%*%c(1,basis(x))))*c(1,basis(x))
    rbind(pu1,pu2)
  }

  model = ZeroFunc(x,y,basis,method,Ufun,Pnu,Ptheta)
  psi = model$estimate
  V = model$AVar
  muA = psi[1]
  muB = psi[2]
  ratio = muB/muA
  n = length(x)+length(y)

  if(CItype == "None"){
    #fisrt derivative of g = muB/muA
    grad = c(-muB/muA^2, 1/muA)
    sd = sqrt((t(grad)%*%V%*%grad)/n)
    val = c(ratio, sd)
    names(val) =c("estimate","ASD")
  }
  if(CItype == "NA-DRM"){
    #fisrt derivative of g = muB/muA
    grad = c(-muB/muA^2, 1/muA)
    sd = sqrt((t(grad)%*%V%*%grad)/n)
    ci = ratio - qnorm(c(1-tau/2,tau/2),0,1)*c(sd)
    val = c(ratio, sd, ci)
    names(val) =c("estimate","ASD", "lower bound","upper bound")
  }
  if(CItype == "log-DRM"){
    #fisrt derivative of g = log(muB)-log(muA)
    grad = c(-1/muA, 1/muB)
    sd = sqrt((t(grad)%*%V%*%grad)/n)
    logCI = log(ratio) - qnorm(c(1-tau/2,tau/2),0,1) *c(sd)
    val = c(ratio, sd, exp(logCI))
    names(val) =c("estimate","ASD", "lower bound","upper bound")
  }
  val
}


