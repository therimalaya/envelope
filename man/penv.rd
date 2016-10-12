\name{penv}
\alias{penv}
\title{Fit the partial envelope model}
\description{
 Fit the partial envelope model in multivariate linear regression with dimension u. 
}
\usage{
penv(X1, X2, Y, u, asy = TRUE) 
}
\arguments{
  \item{X1}{Predictors of main interest. An n by p1 matrix, n is the number of observations, and p1 is the number of main predictors. The predictors can be univariate or multivariate, discrete or continuous.}
  \item{X2}{Covariates, or predictors not of main interest.  An n by p2 matrix, p2 is the number of covariates.}
  \item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is number of observations. The responses must be continuous variables.}
  \item{u}{Dimension of the partial envelope. An integer between 0 and r.}
  \item{asy}{Flag for computing the asymptotic variance of the partial envelope estimator.  The default is \code{TRUE}. When p and r are large, computing the asymptotic variance can take much time and memory.  If only the partial envelope estimators are needed, the flag can be set to \code{asy = FALSE}.}
}
\details{
This function fits the partial envelope model to the responses Y and predictors X1 and X2, \deqn{
 Y = \alpha + \Gamma\eta X_{1} + \beta_{2}X_{2} +\varepsilon, \Sigma=\Gamma\Omega\Gamma'+\Gamma_{0}\Omega_{0}\Gamma'_{0}
}
using the maximum likelihood estimation.  When the dimension of the envelope is between 1 and r - 1, we implemented the algorithm in Su and Cook (2011), but the partial envelope subspace is estimated using the blockwise coordinate descent algorithm in Cook et al. (2015).  When the dimension is r, then the partial envelope model degenerates to the standard multivariate linear regression with Y as the responses and both X1 and X2 as predictors.  When the dimension is 0, X1 and Y are uncorrelated, and the fitting is the standard multivariate linear regression with Y as the responses and X2 as the predictors.
}
\value{The output is a list that contains the following components:
\item{beta1}{The partial envelope estimator of beta1, which is the regression coefficients for X1.} 
\item{beta2}{The partial envelope estimator of beta2, which is the regression coefficients for X2.} 
\item{Sigma}{The partial envelope estimator of the error covariance matrix.}
\item{Gamma}{An orthogonal basis of the partial envelope subspace.}
\item{Gamma0}{An orthogonal basis of the complement of the partial envelope subspace.}
\item{eta}{The coordinates of beta1 with respect to Gamma.}
\item{Omega}{The coordinates of Sigma with respect to Gamma.}
\item{Omega0}{The coordinates of Sigma with respect to Gamma0.}
\item{alpha}{The estimated intercept in the partial envelope model.}
\item{loglik}{The maximized log likelihood function.}
\item{covMatrix}{The asymptotic covariance of vec(beta), while beta = (beta1, beta2). The covariance matrix returned are asymptotic.  For the actual standard errors, multiply by 1 / n.}
\item{asySE}{The asymptotic standard error for elements in beta1 and beta2 under the partial envelope model.  The standard errors returned are asymptotic, for actual standard errors, multiply by 1 / sqrt(n).}
\item{ratio}{The asymptotic standard error ratio of the stanard multivariate linear regression estimator over the partial envelope estimator, for each element in beta1.}
\item{n}{The number of observations in the data.}
}
\references{
Su, Z. and Cook, R.D. (2011). Partial envelopes for efficient estimation in multivariate linear regression. Biometrika 98, 133 - 146.
}

\examples{
data(fiberpaper)
X1 <- fiberpaper[, 7]
X2 <- fiberpaper[, 5:6]
Y <- fiberpaper[, 1:4]
u <- u.penv(X1, X2, Y)
u

m <- penv(X1, X2, Y, 1)
m
m$beta1
}

