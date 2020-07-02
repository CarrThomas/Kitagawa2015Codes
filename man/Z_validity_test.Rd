\name{Z_validity_test}
\alias{Z_validity_test}
\title{A Test for Instrument Validity}
\usage{
Z_validity_test(Y, D, Z, Z_order, xis, B = 500, alpha = c(0.1, 0.05, 0.01))
}
\description{
\code{Z_validity_test} performs the test for instrument validity described in Kitagawa (2015) in the case of a binary treatment, an instrument that takes at least two values, and no conditioning covariates.
}

\arguments{
\item{Y}{A vector containing observations of the outcome variable.}
\item{D}{A vector containing observations of the treatment variable. The treatment must be a binary variable. The test procedure assumes that a value of \eqn{1} denotes receipt of the treatment.}
\item{Z}{A vector containing observations of the instrument variable.}
\item{Z_order}{A vector containing the unique values of the instrument variable ordered so that the propensity score is increasing with the value of the instrument. For example, if the instrument takes values \eqn{1, 2, 3} and we belive that \eqn{P(D=1|Z=1) < P(D=1|Z=2) < P(D=1|Z=3)}, then \code{Z_order} would be \code{c(1, 2, 3)}.}
\item{xis}{A vector containing values of the trimming constant. Individual entries should be \eqn{>0}. For example, \code{c(0.05, 0.1, 0.3, 0.5)}.}
\item{B}{The number of bootstrap iterations when constructing standard errors. Defaults to \code{500}.}
\item{alpha}{A vector of significance levels for which critical values will be calculated. Each entry should be between \eqn{0} and \eqn{1}. By default, critical values are returned for the 10\%, 5\% and 1\% levels of significance.}
}

\details{
The test procedure is described in detail in Kitagawa (2015). Calculating the test statistic and p-values requires maximising a function by searching over intervals of the outcome variable. This is achieved by first restricting the set of intervals under consideration using the properties of the function and then evaluating the function at all remaining intervals. In the worst case, the number of intervals for which the function is evaluated is determined by the number of unique values of the outcome variable for which assignment violations are observed (in the binary instrument - binary treatment case where \eqn{D=1} denotes receipt of the treatment and it is believed that \eqn{P(D=1|Z=1) \ge P(D=1|Z0)}, these would be points where \eqn{(Z,D)=(0,1)} and \eqn{(Z,D)=(1,0)}).

It is generally quicker to run the function once with a \code{xis} vector containing many values than many times with a singleton \code{xis}. Holding the bootstrap samples fixed, p-values will be the same for all values of the trimming constant \eqn{\ge 1/2} so, to avoid unnecessary calculations, it is recommend that all entries in \code{xis} be between \eqn{0} and \eqn{1/2}.
}

\value{
Returns a list with entries
\item{teststat}{A \code{length(xis)} vector with  entries containing the value of the test statistic corresponding to each value of the trimming constant.}
\item{pvals}{A \code{length(xis)} vector of containing the p-value corresponding to each value of the trimming constant.}
\item{cvs}{A matrix with \code{length(alpha)} rows and \code{length(xis)} columns. Entry \eqn{[i,j]} is the critical value for significance level \code{alpha[i]} when the trimming constant is \code{xis[j]}.}

}

\examples{
## Example where the null is valid

set.seed(1234);
Z <- c(rep(1, 1000), rep(0, 1000));
D <- rbinom(2000, 1, 1 / 2);
Y <- D + rnorm(2000);

xis <- c(sqrt(0.005 * (1 - 0.005)), sqrt(0.05 * (1 - 0.05)), sqrt(0.1 * (1 - 0.1)), 1);
Z_validity_test(Y, D, Z, c(0, 1), xis)

## Example where the null is not valid

set.seed(5678);
Z <- c(rep(1, 1000), rep(0, 1000));
D <- rbinom(2000, 1, 0.55 * Z + 0.45 * (1 - Z));
Y <- (1 - D * (1 - Z)) * rnorm(2000) + D * (1 - Z) * rnorm(2000, -0.7, 1);

xis <- c(sqrt(0.005 * (1 - 0.005)), sqrt(0.05 * (1 - 0.05)), sqrt(0.1 * (1 - 0.1)), 1);
Z_validity_test(Y, D, Z, c(0, 1), xis)
}
