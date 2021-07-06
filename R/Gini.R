estGini = function(x,y,basis,method = "optimal",prob = FALSE){
  n0 = length(x)
  n1 = length(y)
  n = n0+n1
  nu0 = sum(x == 0)/n0
  nu1 = sum(y == 0)/n1
  n01 = sum(x != 0)
  n11 = sum(y != 0)
  rho = n11/(n01+n11)
  delta = n0/n*(1-nu0) + n1/n*(1-nu1)

  S01 = x[x!=0]
  S11 = y[y!=0]
  coef = drmCoef(S01,S11,basis,method)

  t = c(S01,S11)
  Qt = cbind(1,matrix(basis(t),nrow = length(t),byrow=F))

  omega = c(exp(Qt%*%coef))
  P0 = 1/(n01+n11*omega)
  P1 = P0*omega

  h = 1 + rho*(omega-1)
  h1 = rho*omega/h

  mu0 = sum(P0*t)
  mu1 = sum(P1*t)
  G0 = function(a){sum(P0*ifelse(t<=a,1,0))}
  G1 = function(a){sum(P1*ifelse(t<=a,1,0))}
  psi0 = sum(P0*t*2*sapply(t,G0))
  psi1 = sum(P1*t*2*sapply(t,G1))
  est0 = (2*nu0-1)+(1-nu0)*psi0/mu0
  est1 = (2*nu1-1)+(1-nu1)*psi1/mu1

  #find the variance
  u0 = function(a){
    H0 =  2*( a*G0(a) + sum(P0*t*ifelse(t>=a,1,0)) ) - psi0
    (2*nu0-1)*a+(1-nu0)*H0
  }
  u1 = function(a){
    H1 =  2*( a*G1(a) + sum(P1*t*ifelse(t>=a,1,0)) ) - psi1
    (2*nu1-1)*a+(1-nu1)*H1
  }
  U = rbind(t,sapply(t,u0),omega*t,omega*sapply(t,u1))
  tildU = rbind(-rho*t,-rho*sapply(t,u0),(1-rho)*t,(1-rho)*sapply(t,u1))

  UUh = (U)%*%(t(U)/h*P0)
  B0 = (tildU)%*%(Qt*h1*P0)
  Atheta = delta*(1-rho)*((t(Qt))%*%(Qt*P0*h1))
  B = B0%*%solve(Atheta)%*%t(B0)
  gDer = matrix(c(-est0/mu0, 1/mu0,rep(0,4),-est1/mu1, 1/mu1),
                nrow = 2,byrow = T)
  Sigma0 =  gDer%*%(UUh/delta + B/(delta*rho^2))%*%t(gDer)
  Sigma1 =matrix(c(nu0*(1-est0)^2/(delta*(1-rho)), rep(0,2),
                   nu1*(1-est1)^2/(delta*rho)),
                 nrow = 2,byrow = T)
  Sigma = (Sigma0 + Sigma1)/n
  if(prob == FALSE){
    list(est =c(est0,est1), var = Sigma)
  }else{
    list(est =c(est0,est1), var = Sigma,
         prob = cbind(c(0,t),c(nu0,(1-nu0)*P0),c(nu1,(1-nu1)*P1)))
  }
}


#' @title Inference on Gini Indices
#' @description Estimate the Gini indices as well as construct the confidence intervals for the Gini indices
#' under the density ratio model (DRM) based on the observed samples.
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
#' @param B the number of bootstap samples. The defualt value is \eqn{1000}.
#'
#' @details The Gini index of the semicontinuous population \eqn{i} can be expessed by
#' \deqn{\mathcal{G}_i=(2\nu_i-1)+(1-\nu_i)\frac{\int_0^{\infty} \{2xG_i(x)\}dG_i(x)}{ \int_0^{\infty} x dG_i(x)}.}
#' \insertCite{yuan2021Gini;textual}{ZeroDRM} estimated the Gini indices \eqn{\mathcal{G}_0} and \eqn{\mathcal{G}_1} based on the above expression.
#' When there is no excess of zero values in the data, the proposed method is also applicable by setting \eqn{\nu_i = 0}.
#'
#'@details The argument \code{CItype} refers to difference confidence intervals of Gini indices
#'under the DRMs in \insertCite{yuan2021Gini;textual}{ZeroDRM}.
#'\itemize{
#'\item  \code{"None"}: no confidence intervals for Gini indices is constructed;
#'\item  \code{"NA-DRM"}: the Wald-type confidence intervals based on the normal approximation;
#'\item  \code{"NL-DRM"}: the Wald-type confidence intervals using logit transformation;
#'\item  \code{"BT-DRM"}: the bootstrap-t confidence intervals;
#'\item  \code{"BL-DRM"}: the bootstrap-t confidence intervals using logit transformation.
#'}
#'
#' @return The function returns a list containing the following components:
#' \itemize{
#' \item \code{Gini0}: the estimate and the confidence interval (lower bound and uppper bound) of \eqn{\mathcal{G}_0}.
#' When \code{CItype = "None"}, it only returns the estimate of \eqn{\mathcal{G}_0}.
#' \item \code{Gini1}: the estimate and the confidence interval (lower bound and uppper bound) of \eqn{\mathcal{G}_1}.
#' When \code{CItype = "None"}, it only returns the estimate of \eqn{\mathcal{G}_1}.
#' \item \code{var}: the estimate of asymptotic variance-covariance matrix \eqn{\Sigma}.
#'  The explicit form of \eqn{\Sigma} van be found in \insertCite{yuan2021Gini;textual}{ZeroDRM}.
#' }
#'
#' @import nnet
#' @import stats
#' @import ineq
#' @importFrom Rdpack reprompt
#'
#'@examples
#'data(Ilocos,package = "ineq")
#'x = Ilocos$income[Ilocos$urban == "urban"&Ilocos$pro == "Pangasinan"]
#'y = Ilocos$income[Ilocos$urban == "rural"&Ilocos$pro == "Pangasinan"]
#'basis = function(x){log(x)}
#'ZeroGini(x,y,basis,CItype = "NA-DRM")
#'
#' @references
#'
#'\insertRef{yuan2021Gini}{ZeroDRM}
#'
#' @export
ZeroGini = function(x,y,basis,method = "optimal",CItype,tau = 0.05, B = 1000){
  model = estGini(x,y,basis,method,prob = TRUE)
  est = model$est
  sd = sqrt(diag(model$var))

  if(missing(CItype)){CItype == "None"}
  if(CItype == "None"){
    val0 = est[1]
    val1 = est[2]
  }

  if(CItype == "NA-DRM"){
    ci0 = est[1] - qnorm(c(1-tau/2,tau/2),0,1)*sd[1]
    ci1 = est[2] - qnorm(c(1-tau/2,tau/2),0,1)*sd[2]
    val0 = c(est[1],ci0)
    names(val0) =c("estimate","lower bound","upper bound")
    val1 = c(est[2],ci1)
    names(val1) =c("estimate","lower bound","upper bound")
  }

  if(CItype == "NL-DRM"){
    estlogit = log(est/(1-est))
    sdlogit = sd/(est*(1-est))
    ci0 = estlogit[1]- qnorm(c(1-tau/2,tau/2),0,1)*sdlogit[1]
    ci1 = estlogit[2]- qnorm(c(1-tau/2,tau/2),0,1)*sdlogit[2]
    val0 = c(est[1],exp(ci0)/(1+exp(ci0)))
    names(val0) =c("estimate","lower bound","upper bound")
    val1 = c(est[2],exp(ci1)/(1+exp(ci1)))
    names(val1) =c("estimate","lower bound","upper bound")
  }

  if(CItype == "BT-DRM"){
    prob = model$prob
    valt = c()
    for(i in 1:B){
      Bootx = sample(prob[,1],length(x),replace = T,prob = prob[,2])
      Booty = sample(prob[,1],length(y),replace = T,prob = prob[,3])
      BootModel = estGini(Bootx,Booty,basis,method)
      Bootest = BootModel$est
      Bootsd = sqrt(diag(BootModel$var))
      valt = rbind(valt,(Bootest- est)/Bootsd)
    }
    ci0 = est[1] - quantile(valt[,1],prob = c(1-tau/2,tau/2))*sd[1]
    ci1 = est[2] - quantile(valt[,2],prob = c(1-tau/2,tau/2))*sd[2]
    val0 = c(est[1],ci0)
    names(val0) =c("estimate","lower bound","upper bound")
    val1 = c(est[2],ci1)
    names(val1) =c("estimate","lower bound","upper bound")
  }

  if(CItype == "BL-DRM"){
    estlogit = log(est/(1-est))
    sdlogit = sd/(est*(1-est))
    prob = model$prob
    valt = c()
    for(i in 1:B){
      Bootx = sample(prob[,1],length(x),replace = T,prob = prob[,2])
      Booty = sample(prob[,1],length(y),replace = T,prob = prob[,3])
      BootModel = estGini(Bootx,Booty,basis,method)
      Bootest = BootModel$est
      Bootsd = sqrt(diag(BootModel$var))
      valt = rbind(valt,(log(Bootest/(1-Bootest))- estlogit)/(Bootsd/(Bootest*(1-Bootest))))
    }
    ci0 = estlogit[1]- quantile(valt[,1],prob = c(1-tau/2,tau/2))*sdlogit[1]
    ci1 = estlogit[2]- quantile(valt[,2],prob = c(1-tau/2,tau/2))*sdlogit[2]
    val0 = c(est[1],exp(ci0)/(1+exp(ci0)))
    names(val0) =c("estimate","lower bound","upper bound")
    val1 = c(est[2],exp(ci1)/(1+exp(ci1)))
    names(val1) =c("estimate","lower bound","upper bound")
  }
  list(Gini0 = val0,Gini1 = val1,AVar = model$var)
}


#' @title Inference on the difference of Gini Indices
#' @description Estimate the difference of Gini indices as well as construct the confidence interval for the different
#' under the density ratio model (DRM) based on the observed samples.
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
#' @param B the number of bootstap samples. The defualt value is \eqn{1000}.
#'
#'@details
#'Let \eqn{\mathcal{G}_i} be the Gini index of population \eqn{i}.
#'The argument \code{CItype} refers to difference confidence intervals for the difference \eqn{\mathcal{G}_0 - \mathcal{G}_1}
#'under the DRMs in \insertCite{yuan2021Gini;textual}{ZeroDRM}.
#'\itemize{
#'\item  \code{"None"}: no confidence intervals for the difference is constructed;
#'\item  \code{"NA-DRM"}: the Wald-type confidence intervals for the difference based on the normal approximation;
#'\item  \code{"BT-DRM"}: the bootstrap-t confidence intervals for the difference ;
#'}
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
#' @import nnet
#' @import stats
#' @import ineq
#' @importFrom Rdpack reprompt
#'
#'@examples
#'data(Ilocos,package = "ineq")
#'x = Ilocos$income[Ilocos$urban == "urban"&Ilocos$pro == "Pangasinan"]
#'y = Ilocos$income[Ilocos$urban == "rural"&Ilocos$pro == "Pangasinan"]
#'basis = function(x){log(x)}
#'DiffGini(x,y,basis,CItype = "NA-DRM")
#'
#' @references
#'
#'\insertRef{yuan2021Gini}{ZeroDRM}
#'
#' @export
DiffGini = function(x,y,basis,method = "optimal",CItype,tau = 0.05, B = 1000){
  model = estGini(x,y,basis,method,prob = TRUE)
  est = model$est[1] - model$est[2]
  sd = c(sqrt(c(1,-1)%*%(model$var)%*%c(1,-1)))

  if(missing(CItype)){CItype == "None"}
  if(CItype == "None"){
    val = c(est,sd)
    names(val) = c("estimate","ASD")
  }

  if(CItype == "NA-DRM"){
    ci = est - qnorm(c(1-tau/2,tau/2),0,1)*sd
    val = c(est,sd,ci)
    names(val) =c("estimate","ASD", "lower bound","upper bound")
  }

  if(CItype == "BT-DRM"){
    prob = model$prob
    valt = c()
    for(i in 1:B){
      Bootx = sample(prob[,1],length(x),replace = T,prob = prob[,2])
      Booty = sample(prob[,1],length(y),replace = T,prob = prob[,3])
      BootModel = estGini(Bootx,Booty,basis,method)
      Bootest = BootModel$est[1] - BootModel$est[2]
      Bootsd = sqrt(c(1,-1)%*%(model$var)%*%c(1,-1))
      valt = c(valt,(Bootest- est)/Bootsd)
    }
    ci = est - quantile(valt,prob = c(1-tau/2,tau/2))*sd
    val = c(est,sd,ci)
    names(val) =c("estimate","ASD", "lower bound","upper bound")
  }
  val
}

