% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DRM.R
\name{drmCoef}
\alias{drmCoef}
\title{Density Ratio Model (DRM)}
\usage{
drmCoef(x, y, basis, method = "optimal")
}
\arguments{
\item{x}{observed sample from \eqn{F_0}.}

\item{y}{observed sample from \eqn{F_1}.}

\item{basis}{pre-sepecified basis function \eqn{\boldsymbol{q}(x)} in the exponential term.
It should return a \eqn{d}-dimensional vector.
For example, \code{basis = function(x){c(x,log(x))}}.}

\item{method}{the method used to fit the DRM. The options include \code{"glm"}, \code{"multinom"} and \code{"optimal"}.
The default setting is \code{"optimal"}. See `Details'.}
}
\value{
The function returns the estimates of parameter \eqn{\boldsymbol{\theta}}.
}
\description{
Fit a semiparametric two-sample DRM based on positive data.
}
\details{
Suppose that two independent samples are from following mixture models:
\deqn{X_{i1},\ldots, X_{in_i} \sim F_i(x)= \nu_iI(x \geq 0)+(1-\nu_i)I(x>0)G_i(x), \quad \mbox{for~}~i=0,1,}
where \eqn{\nu_i\in (0,1)}, \eqn{n_i} is the sample size for the \eqn{i}th sample,
\eqn{I(\cdot)} is an indicator function, and the \eqn{G_i}'s are the cumulative distribution functions (CDFs) of
the positive observations in the \eqn{i}th sample.

\insertCite{yuan2020zeros;textual}{ZeroDRM} and \insertCite{yuan2021Gini;textual}{ZeroDRM} linked \eqn{G_0} and \eqn{G_1} via a DRM:
\deqn{dG_1(x) = \exp\{\alpha + \boldsymbol{\beta}^\top \boldsymbol{q}(x)\}dG_0(x) = \exp\{\boldsymbol{\theta}^\top \boldsymbol{Q}(x)\}dG_0(x),}
where \eqn{dG_i} denote the density of \eqn{G_i}; \eqn{\boldsymbol{\theta} = (\alpha,\boldsymbol{\beta}^\top)^\top} are unknown parameters for the DRM;
\eqn{\boldsymbol{Q}(x) = (1,\boldsymbol{q}(x)^\top)^\top} with \eqn{\boldsymbol{q}(x)} being a pre-specified, non-trivial function of dimension \eqn{d};
the baseline distribution \eqn{F_0} is unspecified.

Let \eqn{n_{i0}} and \eqn{n_{i1}} be the (random) numbers of zero observations and positive observations, respectively, in each sample \eqn{i=0,1}.
 Without loss of generality, we assume that the first \eqn{n_{i1}} observations in group \eqn{i}, \eqn{X_{i1},\ldots, X_{in_{i1}}}, are positive,
 and the remaining \eqn{n_{i0}} observations are 0.
According to \insertCite{yuan2020zeros;textual}{ZeroDRM},
the estimators of \eqn{\boldsymbol{\theta}} maximize the following empirical log-likelihood function
\deqn{\ell_{n}(\boldsymbol{\theta} ) = \sum_{j = 1}^{n_{11}} \{ \boldsymbol{\theta}^\top \boldsymbol{Q}(X_{1j})\} - \sum_{i = 0}^1\sum_{j = 1}^{n_{i1}}\log\left[1+\frac{n_{11}}{n_{01}+n_{11}} \exp\{\boldsymbol{\theta}^\top \boldsymbol{Q}(X_{ij})\}\right].}

Since the DRM is equivalent to the logistic regression model with some justification of \eqn{\alpha},
we can use function \code{glm} or \code{multinom} to solve for parameter \eqn{\boldsymbol{\theta}} instead of maximizing the above likelihood function.
The method \code{optimal} estimates the parameter \eqn{\boldsymbol{\theta}} by \code{glm} or \code{multinom}, whichever gives the larger likelihood \eqn{\ell_{n}(\boldsymbol{\theta})}.
}
\references{
\insertRef{yuan2020zeros}{ZeroDRM}

\insertRef{yuan2021Gini}{ZeroDRM}
}
