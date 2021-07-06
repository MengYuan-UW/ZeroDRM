#' @title Density Ratio Model (DRM)
#' @description Fit a semiparametric two-sample DRM based on positive data.
#'
#' @param x observed sample from \eqn{F_0}.
#' @param y observed sample from \eqn{F_1}.
#' @param basis pre-sepecified basis function \eqn{\boldsymbol{q}(x)} in the exponential term.
#' It should return a \eqn{d}-dimensional vector.
#' For example, \code{basis = function(x){c(x,log(x))}}.
#' @param method the method used to fit the DRM. The options include \code{"glm"}, \code{"multinom"} and \code{"optimal"}.
#' The default setting is \code{"optimal"}. See `Details'.
#'
#' @details Suppose that two independent samples are from following mixture models:
#' \deqn{X_{i1},\ldots, X_{in_i} \sim F_i(x)= \nu_iI(x \geq 0)+(1-\nu_i)I(x>0)G_i(x), \quad \mbox{for~}~i=0,1,}
#' where \eqn{\nu_i\in (0,1)}, \eqn{n_i} is the sample size for the \eqn{i}th sample,
#' \eqn{I(\cdot)} is an indicator function, and the \eqn{G_i}'s are the cumulative distribution functions (CDFs) of
#' the positive observations in the \eqn{i}th sample.
#'
#' \insertCite{yuan2020zeros;textual}{ZeroDRM} and \insertCite{yuan2021Gini;textual}{ZeroDRM} linked \eqn{G_0} and \eqn{G_1} via a DRM:
#' \deqn{dG_1(x) = \exp\{\alpha + \boldsymbol{\beta}^\top \boldsymbol{q}(x)\}dG_0(x) = \exp\{\boldsymbol{\theta}^\top \boldsymbol{Q}(x)\}dG_0(x),}
#' where \eqn{dG_i} denote the density of \eqn{G_i}; \eqn{\boldsymbol{\theta} = (\alpha,\boldsymbol{\beta}^\top)^\top} are unknown parameters for the DRM;
#' \eqn{\boldsymbol{Q}(x) = (1,\boldsymbol{q}(x)^\top)^\top} with \eqn{\boldsymbol{q}(x)} being a pre-specified, non-trivial function of dimension \eqn{d};
#'the baseline distribution \eqn{F_0} is unspecified.
#'
#' Let \eqn{n_{i0}} and \eqn{n_{i1}} be the (random) numbers of zero observations and positive observations, respectively, in each sample \eqn{i=0,1}.
#'  Without loss of generality, we assume that the first \eqn{n_{i1}} observations in group \eqn{i}, \eqn{X_{i1},\ldots, X_{in_{i1}}}, are positive,
#'  and the remaining \eqn{n_{i0}} observations are 0.
#' According to \insertCite{yuan2020zeros;textual}{ZeroDRM},
#' the estimators of \eqn{\boldsymbol{\theta}} maximize the following empirical log-likelihood function
#' \deqn{\ell_{n}(\boldsymbol{\theta} ) = \sum_{j = 1}^{n_{11}} \{ \boldsymbol{\theta}^\top \boldsymbol{Q}(X_{1j})\} - \sum_{i = 0}^1\sum_{j = 1}^{n_{i1}}\log\left[1+\frac{n_{11}}{n_{01}+n_{11}} \exp\{\boldsymbol{\theta}^\top \boldsymbol{Q}(X_{ij})\}\right].}
#'
#' @details Since the DRM is equivalent to the logistic regression model with some justification of \eqn{\alpha},
#' we can use function \code{glm} or \code{multinom} to solve for parameter \eqn{\boldsymbol{\theta}} instead of maximizing the above likelihood function.
#' The method \code{optimal} estimates the parameter \eqn{\boldsymbol{\theta}} by \code{glm} or \code{multinom}, whichever gives the larger likelihood \eqn{\ell_{n}(\boldsymbol{\theta})}.
#'
#' @return The function returns the estimates of parameter \eqn{\boldsymbol{\theta}}.
#'
#' @import nnet
#' @import stats
#' @importFrom Rdpack reprompt
#' @references
#'
#' \insertRef{yuan2020zeros}{ZeroDRM}
#'
#' \insertRef{yuan2021Gini}{ZeroDRM}
#' @export
#'
drmCoef = function(x,y,basis,method = "optimal"){
  # the size of sample from f_0
  n0=length(x)
  # the size of sample from f_1
  n1=length(y)

  t = c(x,y)
  group=c(rep(0,length(x)),rep(1,length(y)))
  Qt = cbind(1,matrix(basis(t),nrow = length(t),byrow=F))
  glm_frame =data.frame(cbind(group,Qt[,-1]))
  formula = as.formula(paste("group~",
                             paste(names(glm_frame)[-1],collapse = "+"),collapse = ""))
  Qy =cbind(1,matrix(basis(y),nrow = length(y),byrow=F))

  dual = function(coef){
    pi = 1/(1+n1/n0*exp(Qt%*%coef))
    loglik = sum(log(pi)) + sum(Qy%*%coef)
    loglik
  }

  if(method == "optimal"){
    out1 = tryCatch(multinom(formula,data = glm_frame, trace = F),
                    error = function(e) NULL)
    out2 = tryCatch(glm(formula,data = glm_frame,family = binomial(link = "logit")),
                    error = function(e) NULL)
    if (is.null(out1) & is.null(out2) ) {stop("parameters of the density ratio model can not been solved")}
    if (is.null(out1)){
      dual1 = -1e5
      outcoef1 = NULL
    }else{
      outcoef1=coef(out1)
      outcoef1[1]=outcoef1[1]-log(n1/n0)
      dual1 = dual(outcoef1)
    }
    if (is.null(out2)){
      dual2 = -1e5
      outcoef2 = NULL
    }else{
      outcoef2 = coef(out2)
      outcoef2[1]=outcoef2[1]-log(n1/n0)
      dual2 = dual(outcoef2)
    }
    # the estimated coefficients
    if (dual1 > dual2){
      outcoef = outcoef1
    }else{
      outcoef = outcoef2
    }
  }
  if(method == "glm"){
    out = tryCatch(glm(formula,data = glm_frame,family = binomial(link = "logit")),
                   error = function(e) NULL)
    if(is.null(out)){
      stop("parameters of the density ratio model can not been solved")
    }else{
      outcoef=coef(out)
      outcoef[1]=outcoef[1]-log(n1/n0)
    }
  }
  if(method == "multinom"){
    out = tryCatch(multinom(formula,data = glm_frame, trace = F),
                   error = function(e) NULL)
    if(is.null(out)){
      stop("parameters of the density ratio model can not been solved")
    }else{
      outcoef=coef(out)
      outcoef[1]=outcoef[1]-log(n1/n0)
    }
  }
  c(outcoef)
}
