This script contains functions implementing the Kitagawa (2015) test without covariates. The instrument may be binary or multivalued. 

After running the script, use the function "mergedZtest" to perform the test

mergedZtest(Y, D, Z, Z_order, xis, B = 500, alpha = c(0.1, 0.05, 0.01))

Inputs are
* Y: A vector containing observations of the outcome variable.
* D: A vector containing observations of the treatment variable. The treatment must be a binary variable. The test procedure assumes that a value of 1 denotes receipt of the treatment.
* Z: A vector containing observations of the instrumental variable.
* Z_order: A vector containing the unique values of the instrument variable ordered so that the propensity score is increasing with the value of the instrument. For example, if the instrument takes values {1, 2, 3} and we believe that P(D=1|Z=1) < P(D=1|Z=2) < P(D=1|Z=3), then Z_order would be c(1, 2, 3)}.
* xis: A vector containing values of the trimming constant. Individual entries should be >0. For example, c(0.05, 0.1, 0.3, 0.5) are the values used in the simulations in Kitagawa (2015).
* B: The number of bootstrap iterations when constructing standard errors. Defaults to 500.
* alpha: A vector of significance levels for which critical values will be calculated. Each entry should be between 0 and 1. By default, critical values are returned for the 10\%, 5\% and 1\% levels of significance.

The function returns a list with entries
* teststat: A length(xis) vector with  entries containing the value of the test statistic corresponding to each value of the trimming constant.
* pvals: A length(xis) vector containing the p-value corresponding to each value of the trimming constant.
* cvs: A matrix with length(alpha) rows and length(xis) columns. Entry [i,j] is the critical value for significance level alpha[i] when the trimming constant is xis[j].

Notes

The test procedure is described in detail in Kitagawa (2015). Calculating the test statistic and p-values requires maximising a function by searching over intervals of the outcome variable. This is achieved by first restricting the set of intervals under consideration using the properties of the function and then evaluating the function at all remaining intervals. In the worst case, the number of intervals for which the function is evaluated is determined by the number of unique values of the outcome variable for which assignment violations are observed. In the binary instrument - binary treatment case where D=1 denotes receipt of the treatment and it is believed that P(D=1|Z=1) >= P(D=1|Z0), these would be points where (Z,D)=(0,1) and (Z,D)=(1,0).

It is generally quicker to run the function once with a xis vector containing many values than many times with a singleton xis. Holding the bootstrap samples fixed, p-values will be the same for all values of the trimming constant >= 1/2 so, to avoid unnecessary calculations, it is recommend that all entries in xis are between 0 and 1/2.


Examples:
 
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