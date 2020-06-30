\name{Z_validity_test}
\alias{Z_validity_test}
\title{A Test for Instrument Validity}
\usage{
Z_validity_test(Y, D, Z, Z_order, xis, B = 500, alpha = c(0.1, 0.05, 0.01))
}
\description{
\code{Z_validity_test} performs the test for instrument validity described in Kitagawa (2015) in the case of a binary treatment, an instrument that can take many values, and no conditioning covariates.
}

\arguments{
\item{Y}{A vector containing observations of the outcome variable.}
\item{D}{A vector containing observations of the treatment variable.}
\item{Z}{A vector containing observations of the instrument variable.}
\item{Z_order}{A vector containing the unique values of the instrument variable ordered so that the propensity score is increasing with the value of the instrument. For example, if the instrument takes values \eqn{1, 2, 3} and we belive that \eqn{P(D=1|Z=1) < P(D=1|Z=2) < P(D=1|Z=3)}, then \code{Z_order} would be \code{c(1, 2, 3)}.}
\item{xis}{A vector containing values of the trimming constant. Individual entries should be \eqn{>0}. For example, \code{c(0.05, 0.1, 0.3, 0.5)}.}
\item{B}{The number of bootstrap iterations when constructing standard errors. Defaults to \code{500}.}
\item{alpha}{A vector of significance levels for which critical values will be calculated. Each entry should be between \eqn{0} and \eqn{1}. By default, critical values are returned for the 10\%, 5\% and 1\% levels of significance.}
}

\details{
a
}

\value{
Returns a list with entries
\item{teststat}{A vector with \code{length(xis)} entries containing the value of the test statistic corresponding to each value of the trimming constant.}
\item{pvals}{A vector of with \code{length(xis)} containing the p-value corresponding to each value of the trimming constant.}
\item{cvs}{A matrix with \code{length(alpha)} rows and \code{length(xis)} columns. Entry \eqn{(i,j)} is the critical value for significance level \code{alpha[i]} when the trimming constant is \code{xis[j]}.}

}

\examples{
## Example where the null is valid

set.seed(1234);
Z <- c(rep(1, 1000), rep(0, 1000));
D <- rbinom(2000, 1, 1 / 2);
Y <- D + rnorm(2000);

xis <- c(sqrt(0.005 * (1 - 0.005)), sqrt(0.05 * (1 - 0.05)), sqrt(0.1 * (1 - 0.1)), 1);
Z_validity_test(Y, D, Z, c(0, 1), xis)

# Returns
$teststat
[1] 3.1702131 2.1545545 1.5652476 0.4695743

$pvals
[1] 0.710 0.912 0.978 0.984

$cvs
        [,1]     [,2]     [,3]     [,4]
90\% 3.855944 3.449081 3.218573 1.274559
95\% 4.036757 3.693522 3.532302 1.408723
99\% 4.245746 4.039350 3.953901 1.677051

## Example where the null is not valid

set.seed(5678);
Z <- c(rep(1, 1000), rep(0, 1000));
D <- rbinom(2000, 1, 0.55 * Z + 0.45 * (1 - Z));
Y <- (1 - D * (1 - Z)) * rnorm(2000) + D * (1 - Z) * rnorm(2000, -0.7, 1);

xis <- c(sqrt(0.005 * (1 - 0.005)), sqrt(0.05 * (1 - 0.05)), sqrt(0.1 * (1 - 0.1)), 1);
Z_validity_test(Y, D, Z, c(0, 1), xis)

# Returns
$teststat
[1] 7.261042 7.261042 7.261042 2.593839

$pvals
[1] 0 0 0 0

$cvs
        [,1]     [,2]     [,3]     [,4]
90\% 3.868514 3.466813 3.282482 1.299155
95\% 4.059778 3.718391 3.624888 1.431084
99\% 4.331277 4.254834 4.018666 1.744133


}
