################################# FRONT MATTER #################################

library(dplyr)
library(tidyr)
library(ggplot2)

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
N <- 1000

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

#### SOME FUNCTIONS ####

update_bounds <- function(bounds, candidate, pos_sum, neg_sum, alt_neg_sum, row, col, n){

  # function to update bounds
  new_bounds <- matrix(0, n, n)

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

  return(bounds)

}

solve_given_row_col <- function(data, intervals, row, col){

  data[intervals[2 * row - 1, 2]:intervals[(2 * col) - 1, 3], 5] <- 1
  data[data[, 2], ] <- data
  results <- kadane_2(subset(data, data[, 5] == 1)[, 3])
  start <- data[results$start, 2]
  end <- data[results$end, 2]

  return(list(val = results$val, start = start, end = end))

}

get_next_row_col <- function(data, results){

  data[data[, 2], ] <- data
  data[results$start:results$end, 5] <- 1
  data[data[, 1], ] <- data
  data <- subset(data, data[, 5] == 1)
  temp <- kadane_2(data[, 3])

  row <- 1 + (data[temp$start, 6] - 1) / 2
  col <- 1 + (data[temp$end, 6] - 1) / 2

  return(list(row = row, col = col))

}

# iterates by one dimenionsal kadane through to convergence
oneD_iteration <- function(data, intervals, pos_sum, neg_sum, alt_neg_sum,
                           row, col, n, record, kadane_count, incumbent){

  flag <- 0

  while(flag == 0){

    # make a record of the current iteration
    if (record[row, col] != 0){
      flag <- 1
    }else{
      record[row, col] <- kadane_count
    }

    kadane_count <- kadane_count + 1

    # given this choice of row and column, solve for the y interval
    results <- solve_given_row_col(data, intervals, row, col)

    # update bounds
    bounds <-
      update_bounds(bounds, results$val, pos_sum, neg_sum, alt_neg_sum, row, col, n)

    if (results$val > incumbent){

      flag <- 1

    }else{

      incumbent <- results$val

      # get the next interval for which we want to solve
      new <- get_next_row_col(data, results)
      new_row <- new$row
      new_col <- new$col

      if(row == new_row & col == new_col){
        flag <- 1
      }else{
        row <- new_row
        col <- new_col
      }
    }

  }

  return(list(bounds = bounds,
              incumbent = incumbent,
              record = record,
              kadane_count = kadane_count))

}





#### MAIN BODY ####

data <- data0
data <- cbind(data, rep(0, nrow(data)))
names(data)[5] <- "inclusion"

# get intervals for x
data[data[, 1], ] <- data
ind <- c(1, data[2:nrow(data), 4] != data[1:(nrow(data) - 1), 4])
bot <- subset(data[, 1], ind == 1)
top <- subset(data[, 1], c(ind[2:length(ind)], 1) == 1)
ind <- cumsum(ind)
data <- cbind(data, ind)
names(data)[6] <- "x_interval"

# get intervals for y
data[data[, 2], ] <- data
ind <- c(1, data[2:nrow(data), 4] != data[1:(nrow(data) - 1), 4])
bot <- subset(data[, 1], ind == 1)
top <- subset(data[, 1], c(ind[2:length(ind)], 1) == 1)
ind <- cumsum(ind)
data <- cbind(data, ind)
names(data)[7] <- "y_interval"

data[data[, 1], ] <- data
data_store <- data

# Similar method to new_method_2, but points are chosen by repeated application
# of Kadane's algorithm.
# 1) no trimming over x, pass to kadane
# 2) given y interval returned by kadane, choose x interval
# 3) given this x interval, use kadane to get y interval
# 4) repeat until you get to a local optimum.
# Update bounds at each stage. After convergence, choose a new starting point
# by looking for the minimum bound

# setup:
intervals <- get_intervals(data)

neg_sum <- cumsum(intervals[seq.int(1, nrow(intervals), 2), 1])
pos_sum <- cumsum(intervals[seq.int(nrow(intervals) - 1, 2, -2), 1])
pos_sum <- c(pos_sum[length(pos_sum):1], 0)

alt_neg_sum <- cumsum(intervals[seq.int(nrow(intervals), 1, -2), 1])
alt_neg_sum <- rev(alt_neg_sum)

n <- length(neg_sum)

bounds <- matrix(neg_sum, n, n, byrow = T) -
  matrix(c(0, neg_sum[1:(n - 1)]), n, n)
bounds[lower.tri(bounds)] <- 0

# matrix to store which intervals we check at which point in time
record <- matrix(0, n, n)

# initial loop: start from top right (no trimming in x)

# initialize
kadane_count <- 1
incumbent <- 0

flag <- 0

start_time <- Sys.time()

while(flag == 0){

  i <- which.min(bounds)
  col <- floor(i / n) + 1
  row <- i - (col - 1) * n

  if (bounds[row, col] < incumbent){

    stage <- oneD_iteration(data, intervals, pos_sum, neg_sum, alt_neg_sum, row,
                       col, n, record, kadane_count, incumbent)

    bounds <- stage$bounds
    record <- stage$record
    kadane_count <- stage$kadane_count
    incumbent <- stage$incumbent

  }else{

    flag <- 1

  }

}

end_time <- Sys.time()

end_time - start_time

temp <- data.frame(x = rep(1:n, n), y = rep(1:n, each = n), record = c(record))

temp$ind <- temp$record != 0

temp %>% ggplot() +
  geom_point(aes(x = x, y = y, col = as.factor(ind)))


new_method_1(data0, 0)



# check that all running sums are negative
data[data[, 1], ] <- data
data[intervals[2 * row - 1, 2]:intervals[(2 * col) - 1, 3], 5] <- 1
data[data[, 2], ] <- data
data <- subset(data, data[, 5] == 1)
data[, 5] <- 0
data[results$start:results$end, 5] <- 1
data <- subset(data, data[, 5] == 1)

cumsum(data[, 3])
sum(cumsum(data[, 3]) > 0)
cumsum(rev(data[, 3]))
sum(cumsum(rev(data[, 3])) > 0)

data <- data[order(data[, 1]), ]
cumsum(data[, 3])
sum(cumsum(data[, 3]) > 0)
cumsum(rev(data[, 3]))
sum(cumsum(rev(data[, 3])) > 0)

# plot of the data
x_min <- min(data[, 1])
x_max <- max(data[, 1])
y_min <- min(data[, 2])
y_max <- max(data[, 2])

data <- data_store

data %>% ggplot() +
  geom_point(aes(x = X, y = Y, col = as.factor(ind), size = abs(K0))) +
  geom_hline(yintercept = y_max) +
  geom_hline(yintercept = y_min) +
  geom_vline(xintercept = x_min) +
  geom_vline(xintercept = x_max)


############################## FIGURING OUT BOUNDS #############################

rm(list=ls())
gc()

# number of observations
N <- 1000

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

#### Different way of calculating bounds

# Use all the information obtained so far to update bounds?

data <- data0
data <- cbind(data, rep(0, nrow(data)))
names(data)[5] <- "inclusion"

# get intervals for x
data[data[, 1], ] <- data
ind <- c(1, data[2:nrow(data), 4] != data[1:(nrow(data) - 1), 4])
bot <- subset(data[, 1], ind == 1)
top <- subset(data[, 1], c(ind[2:length(ind)], 1) == 1)
ind <- cumsum(ind)
data <- cbind(data, ind)
names(data)[6] <- "x_interval"

# get intervals for y
data[data[, 2], ] <- data
ind <- c(1, data[2:nrow(data), 4] != data[1:(nrow(data) - 1), 4])
bot <- subset(data[, 1], ind == 1)
top <- subset(data[, 1], c(ind[2:length(ind)], 1) == 1)
ind <- cumsum(ind)
data <- cbind(data, ind)
names(data)[7] <- "y_interval"

data[data[, 1], ] <- data
data_store <- data


intervals <- get_intervals(data)

neg_sum <- cumsum(intervals[seq.int(1, nrow(intervals), 2), 1])
pos_sum <- cumsum(intervals[seq.int(nrow(intervals) - 1, 2, -2), 1])
pos_sum <- c(pos_sum[length(pos_sum):1], 0)

alt_neg_sum <- cumsum(intervals[seq.int(nrow(intervals), 1, -2), 1])
alt_neg_sum <- rev(alt_neg_sum)

n <- length(neg_sum)

bounds <- matrix(neg_sum, n, n, byrow = T) -
  matrix(c(0, neg_sum[1:(n - 1)]), n, n)
bounds[lower.tri(bounds)] <- 0

# matrix to store which intervals we check at which point in time
record <- matrix(0, n, n)

# initial loop: start from top right (no trimming in x)

# initialize
kadane_count <- 1
incumbent <- 0

# initial pass over smaller intervals
row <- c(1, 6, 11)
col <- c(5, 10, 15)

results <- vector(mode = "list", length = 3)

for (i in 1:3){

  results[[i]] <- solve_given_row_col(data, intervals, row[i], col[i])

}





new_method_2 <- function(data, incumbent){

  intervals <- get_intervals(data)

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
  kadane_count <- 0

  start_time <- Sys.time()
  while(flag == 0){

    # get the current minimum
    i <- which.min(bounds)
    col <- floor(i / n) + 1
    row <- i - (col - 1) * n

    if(bounds[row, col] < incumbent){
      # Pass onto Kadane's algorithm
      kadane_count <- kadane_count + 1

      kappa_vec <- data
      kappa_vec[intervals[2 * row - 1, 2]:intervals[(2 * col) - 1, 3], 5] <- 1
      kappa_vec[kappa_vec[, 2], ] <- kappa_vec
      kappa_vec <- subset(kappa_vec, kappa_vec[, 5] == 1)
      candidate <- kadane(kappa_vec[, 3])

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
  end_time <- Sys.time()

  return(list(incumbent = incumbent,
              bounds = bounds,
              kadane_count = kadane_count,
              run_time = end_time - start_time))

}
