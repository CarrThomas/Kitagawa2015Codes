################################# FRONT MATTER #################################

library(dplyr)

#A script for testing functions with a single covariate.

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

############################## SOLUTION FUNCTIONS ##############################

# Takes data sorted on a single dimension and aggregates it into intervals which
# contain points of the same sign and columns containing row indices of the
# first and last points in the interval
get_intervals <- function(data){

  ind <- c(1, data[2:nrow(data), 4] != data[1:(nrow(data) - 1), 4])
  bot <- subset(data[, 1], ind == 1)
  top <- subset(data[, 1], c(ind[2:length(ind)], 1) == 1)
  ind <- cumsum(ind)
  intervals <- aggregate(data[, 3], list(ind), sum)[, 2]
  intervals <- cbind(intervals, bot, top)

  # #Remove the first and last intervals if the kappa are positive
  # if (intervals[1, 1] >= 0){
  #   intervals <- intervals[2:nrow(intervals), ]
  # }
  # if (intervals[nrow(intervals), 1] >= 0){
  #   intervals <- intervals[1:(nrow(intervals) - 1), ]
  # }

  return(intervals)

}

# Exhaustive solution (apply Kadane's algorithm at each possible interval)
exhaustive_solution <- function(data){

  intervals <- get_intervals(data)

  data <- cbind(data, rep(0, nrow(data)))

  # exhaustive solution at each point
  n <- floor(nrow(intervals) / 2) + 1
  vals <- matrix(0, n, n)

  kadane_count <- 0
  start_time <- Sys.time()
  for (i in 1:n){
    for (j in i:n){

      kadane_count <- kadane_count + 1

      kappa_vec <- data
      kappa_vec[intervals[(2 * i) - 1, 2]:intervals[(2 * j) - 1, 3], 5] <- 1
      kappa_vec[kappa_vec[, 2], ] <- kappa_vec
      kappa_vec <- subset(kappa_vec, kappa_vec[, 5] == 1)
      vals[i, j] <- kadane(kappa_vec[, 3])

    }
  }
  end_time <- Sys.time()

  return(list(true_vals = vals,
              kadane_count = kadane_count,
              run_time = end_time - start_time))

}

################################ GENERATE DATA #################################

# number of observations
N <- 200

# draw parameters
B <- runif(1, -1, 1)
gamma <- runif(1, -1, 1)

# draw variables
X <- rnorm(N)

PZ1 <- pnorm(X * gamma)
Z <- rbinom(N, 1, PZ1)

D <- rbinom(N, 1, 0.5)
Y <- X * B + D + rnorm(N)

# get the kappa functions
K1 <- D * (Z - PZ1) / (PZ1 * (1 - PZ1))
K0 <- (1 - D) * (PZ1 - Z) / (PZ1 * (1 - PZ1))

# form data
data <- cbind(X, Y, K0, K1)

################################## TIDY UP DATA ################################

# Construct subsets for D = 0 and D = 1 such that
# 1) data is trimmed so that extreme points (i.e. max and min X and Y) have
# negative kappa
# 2) X and Y values are replaced with ranks
# 3) data is sorted by X

# D = 0 subsample
data0 <- subset(data, data[, 3] != 0)
data0 <- data0[, 1:3]
data0 <- cbind(data0, data0[, 3] < 0)
data0 <- data.frame(data0)
names(data0)[4] <- "ind"

# trim extreme observations where kappa is positive
temp <- dplyr::filter(data0, ind == 1)
max_Y <- max(temp$Y)
min_Y <- min(temp$Y)
max_X <- max(temp$X)
min_X <- min(temp$X)
data0 <- dplyr::filter(data0, Y <= max_Y & Y >= min_Y & X <= max_X & X >= min_X)

# replace values with rankings
for (k in 1:2){
  data0[, k] <- rank(data0[, k])
}
data0[data0[, 1], ] <- data0

# D = 1 subsample
data1 <- subset(data, data[, 4] != 0)
data1 <- data1[, c(1,2,4)]
data1 <- cbind(data1, data1[, 3] < 0)
data1 <- data.frame(data1)
names(data1)[4] <- "ind"

# trim extreme observations where kappa is positive
temp <- dplyr::filter(data1, ind == 1)
max_Y <- max(temp$Y)
min_Y <- min(temp$Y)
max_X <- max(temp$X)
min_X <- min(temp$X)
data1 <- dplyr::filter(data1, Y <= max_Y & Y >= min_Y & X <= max_X & X >= min_X)

# replace values with rankings
for (k in 1:2){
  data1[, k] <- rank(data1[, k])
}
data1[data1[, 1], ] <- data1

############################## EXHAUSTIVE SOLUTION #############################

# exhaustive solution
exhaustive0 <- exhaustive_solution(data0)
exhaustive1 <- exhaustive_solution(data1)

############################### CHECK PROPERTIES ###############################

temp_data <- data0
temp_data <- cbind(temp_data, rep(0, nrow(temp_data)))

# Find the minimum
vals <- exhaustive0$true_vals
n <- nrow(vals)
i <- which.min(vals)

col <- floor(i / n) + 1
row <- i - n * (col - 1)

min(vals)
vals[row, col]

min_val <- min(vals)

# Find the points included in the final square
intervals <- get_intervals(data0)
temp_data[intervals[(2 * row) - 1, 2]:intervals[(2 * col) - 1, 3], 5] <- 1
kappa_vec <- temp_data
kappa_vec[kappa_vec[, 2], ] <- kappa_vec
kappa_vec <- subset(kappa_vec, kappa_vec[, 5] == 1)
second_round <- kadane_2(kappa_vec[, 3])
kappa_vec <- kappa_vec[second_round$start:second_round$end, ]

# check properties
kappa_vec[, 1] <- rank(kappa_vec[, 1])
kappa_vec[, 2] <- rank(kappa_vec[, 2])

# 1) cumulative sum over dimension 2 is always negative
cumsum(kappa_vec[, 3])
sum(cumsum(kappa_vec[, 3]) > 0)
cumsum(rev(kappa_vec[, 3]))
sum(cumsum(rev(kappa_vec[, 3])) > 0)

# 2) cumulative sum over dimension 1 is always negative
kappa_vec[kappa_vec[, 1], ] <- kappa_vec
cumsum(kappa_vec[, 3])
sum(cumsum(kappa_vec[, 3]) > 0)
cumsum(rev(kappa_vec[, 3]))
sum(cumsum(rev(kappa_vec[, 3])) > 0)

######################## REPEATED APPLICATION OF KADANE ########################

temp_data <- data0
temp_data <- cbind(temp_data, rep(1, nrow(temp_data)))

# apply kadane until there is no more trimming
kadane_2(temp_data[])

