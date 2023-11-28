################################# FRONT MATTER #################################

library(dplyr)

# A script for testing functions with K > 1 covariates
# Data is assumed to be a data frame with columns 1:K containing covariates and
# column K + 1 containing the outcome variable.

rm(list=ls())
gc()

############################### KADANE FUNCTIONS ###############################

# Simple Kadane algorithm: only tracks and returns the minimum subsequence sum
kadane <- function(kappa_vec){

  # Takes a sorted vector and returns the minimum subsequence sum

  val <- 0
  min_val <- 0

  for (i in 1:length(kappa_vec)){

    val <- kappa_vec[i] + min(val, 0)

    if (val < min_val){
      min_val <- val
    }
  }

  return(min_val)
}

# Kadane's algorithm that tracks and returns both the minimum subsequence sum
# and the interval which contains it.
kadane_2 <- function(kappa_vec){

  # Takes a sorted vector and returns the minimum subsequence sum plus the
  # start and end indices of this sum

  val <- 0
  min_val <- 0

  start_index <- 1
  end_index <- 1
  min_start <- 1
  min_end <- 1

  for (i in 1:length(kappa_vec)){

    if (val > 0){

      start_index <- i
      val <- kappa_vec[i]

    }else{

      val <- kappa_vec[i] + val

    }

    end_index <- i

    if (val < min_val){

      min_val <- val
      min_start <- start_index
      min_end <- end_index
    }

  }

  return(list(val = min_val, start = min_start, end = min_end))

}

############################ BOTTOM LAYER FUNCTIONS ############################

# Takes data sorted on a single dimension and aggregates it into intervals which
# contain points of the same sign and columns containing row indices of the
# first and last points in the interval
get_intervals <- function(data, k, K){

  # k is the current layer (on which data is ordered), K is the total number
  # of covariates

  # get an indicator for a change in sign in included points
  temp <- data[as.logical(data[, K + 4]), ]
  ind <- c(1, abs(temp[2:nrow(temp), K + 3] - temp[1:(nrow(temp) - 1), K + 3]) == 2)
  bot <- subset(temp[, k], ind == 1)
  top <- subset(temp[, k], c(ind[2:length(ind)], 1) == 1)
  ind <- cumsum(ind)
  intervals <- aggregate(temp[, K + 2], list(ind), sum)[, 2]
  intervals <- cbind(intervals, bot, top)

  #Remove the first and last intervals if the kappa are positive
  if (intervals[1, 1] >= 0){
    intervals <- intervals[2:nrow(intervals), , drop = F]
  }
  if (intervals[nrow(intervals), 1] >= 0){
    intervals <- intervals[1:(nrow(intervals) - 1), , drop = F]
  }

  return(intervals)

}


# final layer 1: start from the top row (i.e. the lower bound is 1) and work
# down going row by row. At each stage, call Kadane for the interval in the
# current row with the minimum value (if it is below the incumbent) or move on
# to the next row
final_layer_1 <- function(data, incumbent, k, K){

  data[data[, k], ] <- data
  intervals <- get_intervals(data, k, K)

  if (nrow(intervals) == 1){
    incumbent <- min(incumbent, intervals[1, 1])
  }else{

    # Get indices for the negative intervals, and running sums of positive and
    # negative intervals
    neg_sum <- cumsum(intervals[seq.int(1, nrow(intervals), 2), 1])
    pos_sum <- cumsum(intervals[seq.int(nrow(intervals) - 1, 2, -2), 1])
    pos_sum <- c(pos_sum[length(pos_sum):1], 0)

    data <- cbind(data, rep(0, nrow(data)))
    n <- length(neg_sum)
    bounds <- neg_sum

    flag <- 0
    kadane_count <- 0
    j_max <- n
    j <- 1 # current row

    # loop over rows:
    while (flag == 0){

      i <- which.min(bounds)

      if (bounds[i] < incumbent){
        # need to check interval on current row

        kappa_vec <- data
        # exclude all observations outside the current interval
        if (intervals[(2 * j) - 1, 2] > 1){
          kappa_vec[1:intervals[(2 * j) - 1, 2], K + 4] <- 0
        }
        if (intervals[(2 * i) - 1, 3] < nrow(kappa_vec)){
          kappa_vec[intervals[(2 * i) - 1, 3]:nrow(kappa_vec), K + 4] <- 0
        }
        kappa_vec[kappa_vec[, K + 1], ] <- kappa_vec
        kappa_vec <- subset(kappa_vec, kappa_vec[, K + 4] == 1)
        candidate <- kadane(kappa_vec[, K + 2])

        # update bounds
        if (i > 1){
          new_bounds <- candidate + c(- pos_sum[1:(i - 1)] + pos_sum[i],
                                      neg_sum[i:n] - neg_sum[i])
        }else{
          new_bounds <- candidate + neg_sum[i:n] - neg_sum[i]
        }
        bounds <- pmax(bounds, new_bounds)

        # update incumbent
        incumbent <- min(incumbent, candidate)

      }else{

        # move to a new row
        new_row <- neg_sum - neg_sum[j]
        bounds <- pmax(new_row, bounds - pos_sum[j] + pos_sum[j + 1])
        j <- j + 1

      }

      if (j > j_max){

        flag <- 1

      }

    }

  }

  return(incumbent)

}

# final layer 2: consider all possible intervals at once. At each iteration,
# find the minimum bound. If this is below the incumbent, call Kadane. Otherwise,
# terminate and return results.
final_layer_2 <- function(data, incumbent, k, K){

  data[data[, k], ] <- data
  intervals <- get_intervals(data, k, K)

  if (nrow(intervals) == 1){

    incumbent <- min(incumbent, intervals[1, 1])

  }else{

    # Get indices for the negative intervals, and running sums of positive and
    # negative intervals
    neg_sum <- cumsum(intervals[seq.int(1, nrow(intervals), 2), 1])
    pos_sum <- cumsum(intervals[seq.int(nrow(intervals) - 1, 2, -2), 1])
    pos_sum <- c(pos_sum[length(pos_sum):1], 0)

    alt_neg_sum <- cumsum(intervals[seq.int(nrow(intervals), 1, -2), 1])
    alt_neg_sum <- rev(alt_neg_sum)

    data <- cbind(data, rep(0, nrow(data)))

    n <- length(neg_sum)
    bounds <- matrix(neg_sum, n, n, byrow = T) -
      matrix(c(0, neg_sum[1:(n - 1)]), n, n)
    bounds[lower.tri(bounds)] <- 0

    flag <- 0

    while(flag == 0){

      # get the current minimum
      i <- which.min(bounds)
      col <- floor(i / n) + 1
      row <- i - (col - 1) * n

      if(bounds[row, col] < incumbent){
        # Pass onto Kadane's algorithm
        kadane_count <- kadane_count + 1

        kappa_vec <- data
        # exclude all observations outside the current interval
        if (intervals[2 * row - 1, 2] > 1){
          kappa_vec[1:intervals[2 * row - 1, 2], K + 4] <- 0
        }
        if (intervals[(2 * col) - 1, 3] < nrow(kappa_vec)){
          kappa_vec[intervals[(2 * col) - 1, 3]:nrow(kappa_vec), K + 4] <- 0
        }
        kappa_vec[kappa_vec[, K + 1], ] <- kappa_vec
        kappa_vec <- subset(kappa_vec, kappa_vec[, K + 4] == 1)
        candidate <- kadane(kappa_vec[, K + 2])

        # update incumbent
        incumbent <- min(incumbent, candidate)

        # update bounds
        new_bounds <- matrix(0, n, n)

        # update row containing this interval
        if (col > 1){
          new_bounds[row, ] <- candidate + c(- pos_sum[1:(col - 1)] + pos_sum[col],
                                             neg_sum[col:n] - neg_sum[col])
        }else{
          new_bounds[row, ] <- candidate + neg_sum[col:n] - neg_sum[col]
        }
        # update other rows
        if (row > 1){
          new_bounds[1:(row - 1), ] <-
            matrix(new_bounds[row, ], row - 1, n, byrow = T) +
            matrix(alt_neg_sum[1:(row - 1)], row - 1, n) - alt_neg_sum[row]
        }
        if (row < n){
          new_bounds[(row + 1):n, ] <-
            matrix(new_bounds[row, ], n - row, n, byrow = T) -
            pos_sum[row] + matrix(pos_sum[(row + 1):n], (n - row), n)
        }

        bounds <- pmax(bounds, new_bounds)

      }else{
        flag <- 1
      }

    }

  }

  return(incumbent)

}

############################ UPPER LAYER FUNCTIONS #############################

# upper layer 1: go row-by-row, and at each stage evaluate only the minimum for
# the current row
upper_layer_1 <- function(data, incumbent, k, K, final_layer_option){

  data[data[, k], ] <- data
  intervals <- get_intervals(data, k, K)

  if (nrow(intervals) == 1){
    incumbent <- min(incumbent, intervals[1, 1])
  }else{

    # Get indices for the negative intervals, and running sums of positive and
    # negative intervals
    neg_sum <- cumsum(intervals[seq.int(1, nrow(intervals), 2), 1])
    pos_sum <- cumsum(intervals[seq.int(nrow(intervals) - 1, 2, -2), 1])
    pos_sum <- c(pos_sum[length(pos_sum):1], 0)

    n <- length(neg_sum)
    bounds <- neg_sum

    flag <- 0
    j_max <- n
    j <- 1 # current row

    # loop over rows:
    while (flag == 0){

      i <- which.min(bounds)

      if (bounds[i] < incumbent){
        # need to check interval on current row

        kappa_vec <- data
        # exclude all observations outside the current interval
        if (intervals[(2 * j) - 1, 2] > 1){
          kappa_vec[1:intervals[(2 * j) - 1, 2], K + 4] <- 0
        }
        if (intervals[(2 * i) - 1, 3] < nrow(kappa_vec)){
          kappa_vec[intervals[(2 * i) - 1, 3]:nrow(kappa_vec), K + 4] <- 0
        }

        # pass to the next layer
        if (k + 1 < K){
          # intermediate layer
          candidate <- upper_layer_1(kappa_vec, 0, k + 1, K, final_layer_option)
        }else{
          # final layer
          if(final_layer_option == 1){
            candidate <- final_layer_1(kappa_vec, incumbent, k + 1, K)
          }else{
            candidate <- final_layer_2(kappa_vec, incumbent, k + 1, K)
          }

        }

        # update bounds
        if (i > 1){
          new_bounds <- candidate + c(- pos_sum[1:(i - 1)] + pos_sum[i],
                                      neg_sum[i:n] - neg_sum[i])
        }else{
          new_bounds <- candidate + neg_sum[i:n] - neg_sum[i]
        }
        bounds <- pmax(bounds, new_bounds)

        # update incumbent
        incumbent <- min(incumbent, candidate)

      }else{

        # move to a new row
        new_row <- neg_sum - neg_sum[j]
        bounds <- pmax(new_row, bounds - pos_sum[j] + pos_sum[j + 1])
        j <- j + 1

      }

      if (j > j_max){

        flag <- 1

      }

    }
    }

  return(incumbent)

}

# upper layer 2: consider all possible intervals simultaneously. At each stage,
# look at the value for the minimum bound
upper_layer_2 <- function(data, incumbent, k, K, final_layer_option){

  data[data[, k], ] <- data
  intervals <- get_intervals(data, k, K)

  if (nrow(intervals) == 1){
    incumbent <- min(incumbent, intervals[1, 1])
  }else{

    # Get indices for the negative intervals, and running sums of positive and
    # negative intervals
    neg_sum <- cumsum(intervals[seq.int(1, nrow(intervals), 2), 1])
    pos_sum <- cumsum(intervals[seq.int(nrow(intervals) - 1, 2, -2), 1])
    pos_sum <- c(pos_sum[length(pos_sum):1], 0)

    alt_neg_sum <- cumsum(intervals[seq.int(nrow(intervals), 1, -2), 1])
    alt_neg_sum <- rev(alt_neg_sum)

    data <- cbind(data, rep(0, nrow(data)))

    n <- length(neg_sum)
    bounds <- matrix(neg_sum, n, n, byrow = T) -
      matrix(c(0, neg_sum[1:(n - 1)]), n, n)
    bounds[lower.tri(bounds)] <- 0

    flag <- 0

    # loop over rows:
    while (flag == 0){

      i <- which.min(bounds)
      col <- floor(i / n) + 1
      row <- i - (col - 1) * n

      if (bounds[row, col] < incumbent){
        # need to check interval on current row

        kappa_vec <- data

        # exclude all observations outside the current interval
        if (intervals[2 * row - 1, 2] > 1){
          kappa_vec[1:intervals[2 * row - 1, 2], K + 4] <- 0
        }
        if (intervals[(2 * col) - 1, 3] < nrow(kappa_vec)){
          kappa_vec[intervals[(2 * col) - 1, 3]:nrow(kappa_vec), K + 4] <- 0
        }

        # pass to the next layer
        if (k + 1 < K){
          # intermediate layer
          candidate <- upper_layer_2(kappa_vec, 0, k + 1, K, final_layer_option)
        }else{
          # final layer
          if(final_layer_option == 1){
            candidate <- final_layer_1(kappa_vec, incumbent, k + 1, K)
          }else{
            candidate <- final_layer_2(kappa_vec, incumbent, k + 1, K)
          }

        }

        # update bounds
        new_bounds <- matrix(0, n, n)

        # update row containing this interval
        if (col > 1){
          new_bounds[row, ] <- candidate + c(- pos_sum[1:(col - 1)] + pos_sum[col],
                                             neg_sum[col:n] - neg_sum[col])
        }else{
          new_bounds[row, ] <- candidate + neg_sum[col:n] - neg_sum[col]
        }
        # update other rows
        if (row > 1){
          new_bounds[1:(row - 1), ] <-
            matrix(new_bounds[row, ], row - 1, n, byrow = T) +
            matrix(alt_neg_sum[1:(row - 1)], row - 1, n) - alt_neg_sum[row]
        }
        if (row < n){
          new_bounds[(row + 1):n, ] <-
            matrix(new_bounds[row, ], n - row, n, byrow = T) -
            pos_sum[row] + matrix(pos_sum[(row + 1):n], (n - row), n)
        }

        bounds <- pmax(bounds, new_bounds)

        # update incumbent
        incumbent <- min(incumbent, candidate)

      }else{
        # terminate

        flag <- 1
      }

    }

  }

  return(incumbent)

}

################################ GENERATE DATA #################################

# number of observations
N <- 500

# number of conditioning covariates
K <- 3

# draw parameters
B <- runif(K, -1, 1)
gamma <- runif(K, -1, 1)

# draw variables
X <- matrix(rnorm(N * K), N, K)

PZ1 <- pnorm(X %*% gamma)
Z <- rbinom(N, 1, PZ1)

D <- rbinom(N, 1, 0.5)
Y <- X %*% B + D + rnorm(N)

# get the kappa functions
K1 <- D * (Z - PZ1) / (PZ1 * (1 - PZ1))
K0 <- (1 - D) * (PZ1 - Z) / (PZ1 * (1 - PZ1))

# form data
data <- cbind(X, Y, K0, K1)

################################## TIDY UP DATA ################################

# Construct subsets for D = 0 and D = 1 such that
# 1) data is trimmed so that extreme points (i.e. max and min for each X and Y)
# have negative kappa
# 2) X and Y values are replaced with ranks
# 3) data is sorted by the first X

# D = 0 subsample
data0 <- subset(data, data[, K + 2] != 0)
data0 <- data0[, 1:(K + 2)]
data0 <- cbind(data0, sign(data0[, K + 2]), rep(1, nrow(data0)))
data0 <- data.frame(data0)
names(data0) <- c(paste("X", 1:K, sep = ""), "Y", "K0", "sign", "incl_ind")

# trim extreme observations where kappa is positive
temp <- dplyr::filter(data0, sign == -1)
max_Y <- max(temp$Y)
min_Y <- min(temp$Y)
data0 <- dplyr::filter(data0, Y <= max_Y & Y >= min_Y)
for (k in 1:K){
  max_X <- max(temp[, k])
  min_X <- min(temp[, k])
  data0 <- dplyr::filter(data0, data0[, k] <= max_X & data0[, k] >= min_X)
}

# replace values with rankings
for (k in 1:(K + 1)){
  data0[, k] <- rank(data0[, k])
}
data0[data0[, 1], ] <- data0

# D = 1 subsample
data1 <- subset(data, data[, K + 3] != 0)
data1 <- data1[, c(1:(K + 1), K + 3)]
data1 <- cbind(data1, sign(data1[, K + 2]), rep(1, nrow(data1)))
data1 <- data.frame(data1)
names(data1) <- c(paste("X", 1:K, sep = ""), "Y", "K1", "sign", "incl_ind")

# trim extreme observations where kappa is positive
temp <- dplyr::filter(data1, sign == -1)
max_Y <- max(temp$Y)
min_Y <- min(temp$Y)
data1 <- dplyr::filter(data1, Y <= max_Y & Y >= min_Y)
for (k in 1:K){
  max_X <- max(temp[, k])
  min_X <- min(temp[, k])
  data1 <- dplyr::filter(data1, data1[, k] <= max_X & data1[, k] >= min_X)
}

# replace values with rankings
for (k in 1:(K + 1)){
  data1[, k] <- rank(data1[, k])
}
data1[data1[, 1], ] <- data1

############################## EXHAUSTIVE SOLUTION #############################

# # exhaustive solution
# exhaustive0 <- exhaustive_solution(data0)
# exhaustive1 <- exhaustive_solution(data1)
#
# true_solution <- min(exhaustive0$true_vals, exhaustive1$true_vals)

############################# STARTING VALUE 1 RUNS ############################

# get starting value by optimisation over a single dimension
incumbent <- rep(0, 2 * (K + 1))
for (k in 1:(K + 1)){
  data0[data0[, K + 2 - k], ] <- data0
  incumbent[2 * (k - 1) + 1] <- kadane(data0[, K + 2])
  data1[data1[, K + 2 - k], ] <- data1
  incumbent[2 * (k - 1) + 2] <- kadane(data1[, K + 2])
}
incumbent <- min(incumbent)

# different optimisers
start_time <- Sys.time()
upper_layer_1(data0, incumbent, 1, K, 1)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
upper_layer_2(data0, incumbent, 1, K, 1)
end_time <- Sys.time()
end_time - start_time



