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

# August 2019 algorithm: starts from the bottom (bounds N, N) and works up
# (i.e. next start value considered is N - 1, then N - 2, and so on)
august2019algo <- function(data, incumbent){

  intervals <- get_intervals(data)

  # Get indices for the negative intervals, and running sums of positive and
  # negative intervals
  neg_sum <- cumsum(intervals[seq.int(1, nrow(intervals), 2), 1])
  pos_sum <- cumsum(intervals[seq.int(nrow(intervals) - 1, 2, -2), 1])
  pos_sum <- c(pos_sum[length(pos_sum):1], 0)
  data <- cbind(data, rep(0, nrow(data)))

  bounds <- neg_sum - neg_sum[length(neg_sum)]
  kadane_count <- 0

  start_time <- Sys.time()
  for (i in length(neg_sum):1){
    bounds <- bounds + intervals[(2 * i) - 1, 1];
    for (j in i:length(neg_sum)){
      if (bounds[j] < incumbent){
        kadane_count <- kadane_count + 1
        kappa_vec <- data
        kappa_vec[intervals[(2 * i) - 1, 2]:intervals[(2 * j) - 1, 3], 5] <- 1
        kappa_vec[kappa_vec[, 2], ] <- kappa_vec
        kappa_vec <- subset(kappa_vec, kappa_vec[, 5] == 1)
        candidate <- kadane(kappa_vec[, 3]); #minimum for this box

        if (j > i + 1){#case where points to the left and right can be updated
          new_bounds <- candidate +
            c(- pos_sum[(i + 1):(j - 1)] + pos_sum[j],
              neg_sum[j:length(neg_sum)] - neg_sum[j])
          bounds[(i + 1):length(bounds)] <-
            pmax(bounds[(i + 1):length(bounds)], new_bounds)
        }
        else{#case where only points to the right can be updated
          new_bounds <- candidate + neg_sum[j:length(neg_sum)] - neg_sum[j]
          bounds[j:length(bounds)] <- pmax(bounds[j:length(bounds)], new_bounds)
        }

        incumbent <- min(incumbent, candidate)
      }
    }
  }
  end_time <- Sys.time()

  return(list(incumbent = incumbent,
              kadane_count = kadane_count,
              run_time = end_time - start_time))

}

# New method 1: start from the top row (i.e. the lower bound is 1) and work down
# going row by row. At each stage, call Kadane for the interval in the current
# row with the minimum value (if it is below the incumbent) or move on to the
# next row
new_method_1 <- function(data, incumbent){

  intervals <- get_intervals(data)

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
  start_time <- Sys.time()
  while (flag == 0){

    i <- which.min(bounds)

    if (bounds[i] < incumbent){
      # need to check interval on current row

      kadane_count <- kadane_count + 1

      kappa_vec <- data
      kappa_vec[intervals[(2 * j) - 1, 2]:intervals[(2 * i) - 1, 3], 5] <- 1
      kappa_vec[kappa_vec[, 2], ] <- kappa_vec
      kappa_vec <- subset(kappa_vec, kappa_vec[, 5] == 1)
      candidate <- kadane(kappa_vec[, 3])

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
  end_time <- Sys.time()

  return(list(incumbent = incumbent,
              kadane_count = kadane_count,
              run_time = end_time - start_time))

}

# New method 2: consider all possible intervals at once. At each iteration, find
# the minimum bound. If this is below the incumbent, call Kadane. Otherwise,
# terminate and return results.
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

# New method 3: consider all possible intervals at once. At each iteration, find
# the minimum bound. If this is above the incumbent, terminate. Otherwise, call
# Kadane. If the fitst-state result is above the current minimum, apply Kadane
# iteratively until convergence, then move onto another point
new_method_3 <- function(data, incumbent){

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

  # initialize
  kadane_count <- 1
  flag <- 0

  start_time <- Sys.time()

  while(flag == 0){

    i <- which.min(bounds)
    col <- floor(i / n) + 1
    row <- i - (col - 1) * n

    if (bounds[row, col] < incumbent){

      stage <- oneD_iteration(data, intervals, pos_sum, neg_sum, alt_neg_sum, row,
                              col, n, kadane_count, incumbent, bounds)

      bounds <- stage$bounds
      kadane_count <- stage$kadane_count
      incumbent <- stage$incumbent

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
                           row, col, n, kadane_count, incumbent, bounds){

  flag <- 0

  while(flag == 0){

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
              kadane_count = kadane_count))

}

data_prep <- function(data){

  # indicator for negative values
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

  return(data)

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

true_solution <- min(exhaustive0$true_vals, exhaustive1$true_vals)

############################# STARTING VALUE 1 RUNS ############################

# get starting value by optimisation over a single dimension
# incumbent <- rep(0, 4)
# for (k in 1:2){
#   data0[data0[, 3 - k], ] <- data0
#   incumbent[k] <- kadane(data0[, 3])
# }
# for (k in 1:2){
#   data1[data1[, 3 - k], ] <- data1
#   incumbent[k + 2] <- kadane(data1[, 3])
# }
# incumbent <- min(incumbent)
incumbent <- 0

# different optimisers
august2019_0 <- august2019algo(data0, incumbent)
august2019_1 <- august2019algo(data1, august2019_0$incumbent)

august2019_solution <- august2019_1$incumbent

new_method_1_D0 <- new_method_1(data0, incumbent)
new_method_1_D1 <- new_method_1(data1, new_method_1_D0$incumbent)

new_method_1_solution <- new_method_1_D1$incumbent

new_method_2_D0 <- new_method_2(data0, incumbent)
new_method_2_D1 <- new_method_2(data1, new_method_2_D0$incumbent)

new_method_2_solution <- new_method_2_D1$incumbent

data0 <- data_prep(data0)
data1 <- data_prep(data1)

new_method_3_D0 <- new_method_3(data0, incumbent)
new_method_3_D1 <- new_method_3(data1, new_method_3_D0$incumbent)

new_method_3_solution <- new_method_3_D1$incumbent

################################## COMPARISONS #################################

# comparison of solutions
solutions <- c(true_solution, august2019_solution, new_method_1_solution,
               new_method_2_solution, new_method_3_solution)
solutions

# comparison of run times
run_times_D0 <- c(exhaustive0$run_time, august2019_0$run_time,
                  new_method_1_D0$run_time, new_method_2_D0$run_time,
                  new_method_3_D0$run_time)

run_times_D1 <- c(exhaustive1$run_time, august2019_1$run_time,
                  new_method_1_D1$run_time, new_method_2_D1$run_time,
                  new_method_3_D1$run_time)


# comparison of Kadane calls
kadane_calls_D0 <- c(exhaustive0$kadane_count, august2019_0$kadane_count,
                     new_method_1_D0$kadane_count, new_method_2_D0$kadane_count,
                     new_method_3_D0$kadane_count)

kadane_calls_D1 <- c(exhaustive1$kadane_count, august2019_1$kadane_count,
                     new_method_1_D1$kadane_count, new_method_2_D1$kadane_count,
                     new_method_3_D1$kadane_count)

# comparison of bounds for new method 2 and the true values in each interval
bounds_check_D0 <-
  sum(round(exhaustive0$true_vals[upper.tri(exhaustive0$true_vals)], 8) <
        round(new_method_2_D0$bounds[upper.tri(new_method_2_D0$bounds)], 8))
bounds_check_D1 <-
  sum(round(exhaustive1$true_vals[upper.tri(exhaustive1$true_vals)], 8) <
      round(new_method_2_D1$bounds[upper.tri(new_method_2_D1$bounds)], 8))

bounds_check_3_D0 <-
  sum(round(exhaustive0$true_vals[upper.tri(exhaustive0$true_vals)], 8) <
        round(new_method_3_D0$bounds[upper.tri(new_method_3_D0$bounds)], 8))
bounds_check_3_D1 <-
  sum(round(exhaustive1$true_vals[upper.tri(exhaustive1$true_vals)], 8) <
        round(new_method_3_D1$bounds[upper.tri(new_method_3_D1$bounds)], 8))


comparisons_starting_1 <- cbind(solutions, run_times_D0, run_times_D1,
                                kadane_calls_D0, kadane_calls_D1,
                                c(rep(NA, 3), bounds_check_D0, bounds_check_3_D0),
                                c(rep(NA, 3), bounds_check_D1, bounds_check_3_D1))


row.names(comparisons_starting_1) <- c("Exhaustive", "2019", "New 1", "New 2", "New 3")
colnames(comparisons_starting_1) <- c("Minimum",
                                      "Run Time D = 0",
                                      "Run Time D = 1",
                                      "Kadane Calls D = 0",
                                      "Kadane Calls D = 1",
                                      "Bounds Check D = 0",
                                      "Bounds Check D = 1")


comparisons_starting_1

############################# STARTING VALUE 2 RUNS ############################

# get a starting guess for the alternative algorithms by "greedy" optimisation
# over both dimensions
# "greedy solution" over both dimensions to find the incumbent
incumbent <- rep(0, 4)

data0[data0[, 2], ] <- data0
temp <- kadane_2(data0[, 3])
temp_data <- data0[temp$start:temp$end, ]
temp_data <- temp_data[order(temp_data[, 1]), ]
incumbent[1] <- kadane(temp_data[, 3])

data0[data0[, 1], ] <- data0
temp <- kadane_2(data0[, 3])
temp_data <- data0[temp$start:temp$end, ]
temp_data <- temp_data[order(temp_data[, 2]), ]
incumbent[2] <- kadane(temp_data[, 3])

data1[data1[, 2], ] <- data1
temp <- kadane_2(data1[, 3])
temp_data <- data1[temp$start:temp$end, ]
temp_data <- temp_data[order(temp_data[, 1]), ]
incumbent[3] <- kadane(temp_data[, 3])

data1[data1[, 1], ] <- data1
temp <- kadane_2(data1[, 3])
temp_data <- data1[temp$start:temp$end, ]
temp_data <- temp_data[order(temp_data[, 2]), ]
incumbent[4] <- kadane(temp_data[, 3])

incumbent <- min(incumbent)

# different optimisers
august2019_0 <- august2019algo(data0, incumbent)
august2019_1 <- august2019algo(data1, august2019_0$incumbent)

august2019_solution <- august2019_1$incumbent

new_method_1_D0 <- new_method_1(data0, incumbent)
new_method_1_D1 <- new_method_1(data1, new_method_1_D0$incumbent)

new_method_1_solution <- new_method_1_D1$incumbent

new_method_2_D0 <- new_method_2(data0, incumbent)
new_method_2_D1 <- new_method_2(data1, new_method_2_D0$incumbent)

new_method_2_solution <- new_method_2_D1$incumbent

################################## COMPARISONS #################################

# comparison of solutions
solutions <- c(true_solution, august2019_solution, new_method_1_solution,
               new_method_2_solution)

# comparison of run times
run_times_D0 <- c(exhaustive0$run_time, august2019_0$run_time,
                  new_method_1_D0$run_time, new_method_2_D0$run_time)

run_times_D1 <- c(exhaustive1$run_time, august2019_1$run_time,
                  new_method_1_D1$run_time, new_method_2_D1$run_time)

# comparison of Kadane calls
kadane_calls_D0 <- c(exhaustive0$kadane_count, august2019_0$kadane_count,
                     new_method_1_D0$kadane_count, new_method_2_D0$kadane_count)

kadane_calls_D1 <- c(exhaustive1$kadane_count, august2019_1$kadane_count,
                     new_method_1_D1$kadane_count, new_method_2_D1$kadane_count)

# comparison of bounds for new method 2 and the true values in each interval
bounds_check_D0 <-
  sum(round(exhaustive0$true_vals[upper.tri(exhaustive0$true_vals)], 8) <
        round(new_method_2_D0$bounds[upper.tri(new_method_2_D0$bounds)], 8))
bounds_check_D1 <-
  sum(round(exhaustive1$true_vals[upper.tri(exhaustive1$true_vals)], 8) <
        round(new_method_2_D1$bounds[upper.tri(new_method_2_D1$bounds)], 8))


comparisons_starting_2 <- cbind(solutions, run_times_D0, run_times_D1,
                                kadane_calls_D0, kadane_calls_D1,
                                c(rep(NA, 3), bounds_check_D0),
                                c(rep(NA, 3), bounds_check_D1))


row.names(comparisons_starting_2) <- c("Exhaustive", "2019", "New 1", "New 2")
colnames(comparisons_starting_2) <- c("Minimum",
                                      "Run Time D = 0",
                                      "Run Time D = 1",
                                      "Kadane Calls D = 0",
                                      "Kadane Calls D = 1",
                                      "Bounds Check D = 0",
                                      "Bounds Check D = 1")

comparisons_starting_1
comparisons_starting_2

############################### OVERALL FUNCTION ###############################

#Construct a function that solves the 2D case.
#Inputs will be data and the number of bootstrap iterations. The output will be a P-value.

solve <- function(data, B){
  #Data assumed to be a N * 4 matrix where column 1 is x, column 2 is y, column 3 is kappa0 and column 4 is kappa1.

  #Re-order the data
  data0 <- subset(data, data[, 3] != 0);
  for (k in 1:2){
    data0[, k] <- rank(data0[, k]);
  }
  data1 <- subset(data, data[, 4] != 0);
  for (k in 1:2){
    data1[, k] <- rank(data1[, k]);
  }
  data0[data0[, 1], ] <- data0;
  data1[data1[, 1], ] <- data1;
  data <- rbind(data0, data1);
  N0 <- nrow(data0); #comes in handy later
  #Note that data is now always initially ordered by x

  data0 <- data0[, 1:3];
  data0 <- cbind(data0, data0[, 3] < 0);
  data1 <- data1[, c(1,2,4)];
  data1 <- cbind(data1, data1[, 3] < 0);

  #Solve the unidimensional problem to find the incumbent
  incumbent <- rep(0, 4);
  for (k in 1:2){
    data0[data0[, 3 - k], ] <- data0;
    incumbent[k] <- kadane(data0[, 3]);
  }
  for (k in 1:2){
    data1[data1[, 3 - k], ] <- data1;
    incumbent[k + 2] <- kadane(data1[, 3]);
  }
  incumbent <- min(incumbent)

  #Solve the 2D problem
  incumbent <- enumerate(data0, incumbent);
  incumbent <- enumerate(data1, incumbent);

  #Bootstrap P-values
  #Data contains data0 followed by data1, both sub-matrices are sorted by x
  flag <- rep(0, B); #stores whether the bootstrap statistic exceeds the original or not
  for (b in 1:B){
    boot_sample <- sample.int(nrow(data), nrow(data), replace = TRUE);
    boot_sample <- tabulate(boot_sample, nrow(data));
    boot_data <- data[, 3:4] * (boot_sample - 1);
    boot_data <- cbind(data[, 1:2], boot_data);

    data0 <- boot_data[1:N0, ]; #points that could potentially have kappa0 != 0
    data0 <- subset(data0, data0[, 3] != 0); #remove the zeroes then re-rank
    data0[, 1] <- seq.int(1, nrow(data0));
    data0[, 2] <- rank(data0[, 2]);

    data1 <- boot_data[(N0 + 1):nrow(data), ]; #points that could potentially have kappa0 != 0
    data1 <- subset(data1, data1[, 4] != 0); #remove the zeroes then re-rank
    data1[, 1] <- seq.int(1, nrow(data1));
    data1[, 2] <- rank(data1[, 2]);

    data0 <- data0[, 1:3];
    data0 <- cbind(data0, data0[, 3] < 0);
    data1 <- data1[, c(1,2,4)];
    data1 <- cbind(data1, data1[, 3] < 0);

    flag[b] <- boot_enumerate(data0, incumbent);
    if (flag[b] == 0){
      flag[b] <- boot_enumerate(data1, incumbent);
    }
  }
  return(mean(flag));
}

############################## KADANE'S ALGORITHM ##############################

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

############################################### ENUMERATE INTERVALS ###############################################

enumerate <- function(data, incumbent){
  #Data assumed to be a matrix with values ordered by the x variable

  #Collapse to a set of intervals composed of the sums of neighbouring points of the same sign
  ind <- c(1, data[2:nrow(data), 4] != data[1:(nrow(data) - 1), 4]);
  bot <- subset(data[, 1], ind == 1);
  top <- subset(data[, 1], c(ind[2:length(ind)], 1) == 1);
  ind <- cumsum(ind);
  intervals <- aggregate(data[, 3], list(ind), sum)[, 2];
  intervals <- cbind(intervals, bot, top);

  #Remove the first and last intervals if the kappa are positive
  if (intervals[1, 1] >= 0){
    intervals <- intervals[2:nrow(intervals), ];
  }
  if (intervals[nrow(intervals), 1] >= 0){
    intervals <- intervals[1:(nrow(intervals) - 1), ];
  }

  #Get indices for the negative intervals, and kappa sums for positive and negative intervals
  neg_sum <- cumsum(intervals[seq.int(1, nrow(intervals), 2), 1]);
  pos_sum <- cumsum(intervals[seq.int(nrow(intervals) - 1, 2, -2), 1]);
  pos_sum <- c(pos_sum[length(pos_sum):1], 0);
  bounds <- neg_sum - neg_sum[length(neg_sum)];
  data <- cbind(data, rep(0, nrow(data)));

  for (i in length(neg_sum):1){
    bounds <- bounds + intervals[(2 * i) - 1, 1];
    for (j in i:length(neg_sum)){
      if (bounds[j] < incumbent){
        kappa_vec <- data;
        kappa_vec[intervals[(2 * i) - 1, 2]:intervals[(2 * j) - 1, 3], 5] <- 1;
        kappa_vec[kappa_vec[, 2], ] <- kappa_vec;
        kappa_vec <- subset(kappa_vec, kappa_vec[, 5] == 1);
        candidate <- kadane(kappa_vec[, 3]); #minimum for this box

        if (j > i + 1){#case where points to the left and right can be updated
          new_bounds <- candidate + c(- pos_sum[(i + 1):(j - 1)] + pos_sum[j], neg_sum[j:length(neg_sum)] - neg_sum[j])
          bounds[(i + 1):length(bounds)] <- pmax(bounds[(i + 1):length(bounds)], new_bounds);
        }
        else{#case where only points to the right can be updated
          new_bounds <- candidate + neg_sum[j:length(neg_sum)] - neg_sum[j]
          bounds[j:length(bounds)] <- pmax(bounds[j:length(bounds)], new_bounds);
        }

        incumbent <- min(incumbent, candidate);
      }
    }
  }
  return(incumbent);
}

############################################### BOOTSTRAP ENUMERATE ###############################################

boot_enumerate <- function(data, incumbent){
  #Data assumed to be a matrix with values ordered by the x variable
  flag <- 0;

  #Collapse to a set of intervals composed of the sums of neighbouring points of the same sign
  ind <- c(1, data[2:nrow(data), 4] != data[1:(nrow(data) - 1), 4]);
  bot <- subset(data[, 1], ind == 1);
  top <- subset(data[, 1], c(ind[2:length(ind)], 1) == 1);
  ind <- cumsum(ind);
  intervals <- aggregate(data[, 3], list(ind), sum)[, 2];
  intervals <- cbind(intervals, bot, top);

  #Remove the first and last intervals if the kappa are positive
  if (intervals[1, 1] >= 0){
    intervals <- intervals[2:nrow(intervals), ];
  }
  if (intervals[nrow(intervals), 1] >= 0){
    intervals <- intervals[1:(nrow(intervals) - 1), ];
  }

  #Get indices for the negative intervals, and kappa sums for positive and negative intervals
  neg_sum <- cumsum(intervals[seq.int(1, nrow(intervals), 2), 1]);
  pos_sum <- cumsum(intervals[seq.int(nrow(intervals) - 1, 2, -2), 1]);
  pos_sum <- c(pos_sum[length(pos_sum):1], 0);
  bounds <- neg_sum - neg_sum[length(neg_sum)];
  data <- cbind(data, rep(0, nrow(data)));

  for (i in length(neg_sum):1){
    if (flag == 1){
      break
    }
    bounds <- bounds + intervals[(2 * i) - 1, 1];
    for (j in i:length(neg_sum)){
      if (bounds[j] < incumbent){
        kappa_vec <- data;
        kappa_vec[intervals[(2 * i) - 1, 2]:intervals[(2 * j) - 1, 3], 5] <- 1;
        kappa_vec[kappa_vec[, 2], ] <- kappa_vec;
        kappa_vec <- subset(kappa_vec[, 3], kappa_vec[, 5] == 1);
        candidate <- kadane(kappa_vec); #minimum for this box
        if (candidate < incumbent){
          flag <- 1;
          break;
        }
        if (j > i + 1){#case where points to the left and right can be updated
          new_bounds <- candidate + c(- pos_sum[(i + 1):(j - 1)] + pos_sum[j], neg_sum[j:length(neg_sum)] - neg_sum[j]);
          bounds[(i + 1):length(bounds)] <- pmax(bounds[(i + 1):length(bounds)], new_bounds);
        }
        else{#case where only points to the right can be updated
          new_bounds <- candidate + neg_sum[j:length(neg_sum)] - neg_sum[j]
          bounds[j:length(bounds)] <- pmax(bounds[j:length(bounds)], new_bounds);
        }
      }
    }
  }
  return(flag);
}

############################################### MONTE CARLO SETUP ###############################################

set.seed(1);
MC_testsize <- list();

J <- 20; #Number of Monte Carlo iterations
B <- 300; #Number of bootstap iterations
a <- c(0.1, 0.05, 0.01); #Nominal test size
N_a <- length(a);

############################################### SMALL SAMPLE SIZE ###############################################

N <- 200;
data.array <- array(0, dim = c(N, 4, J));

#Generate J Datasets
for (j in 1:J){
  parameter <- runif(2, 0, 1); #B and p parameters
  data.array[, 1, j] <- rnorm(N); #X
  D <- rbinom(N, 1, 1 / 2); #D
  data.array[, 2, j] <- D +  parameter[1] * data.array[, 1, j] + rnorm(N); #Y
  PZ1 <- pnorm(data.array[, 1, j] * parameter[2]);
  Z <- rbinom(N, 1, PZ1); #Z
  kap <-  (PZ1 - Z) / (PZ1 *(1 - PZ1));
  data.array[, 3, j] <- (1 - D) * kap; #kappa0
  data.array[, 4, j] <- - D * kap; #kappa1
}

pvals <- apply(data.array, 3, solve, B = B);
reject <- matrix(0, J, N_a);
for (i in 1:N_a){
  reject[, i] <- pvals <= a[i];
}
reject_rate <- colMeans(reject);

MC_testsize$Small$rejection <- reject_rate;
#save(MC_testsize, file="MCsize.Rdata");

############################################ SMALL MEDIUM SAMPLE SIZE ############################################

N <- 600;
data.array <- array(0, dim = c(N, 4, J));
rho <- qnorm(1/6);

#Generate J Datasets
for (j in 1:J){
  parameter <- runif(2, 0, 1); #B and p parameters
  data.array[, 1, j] <- rnorm(N); #X
  D <- rbinom(N, 1, 1 / 2); #D
  data.array[, 2, j] <- D +  parameter[1] * data.array[, 1, j] + rnorm(N); #Y
  PZ1 <- pnorm(rho + data.array[, 1, j] * parameter[2]);
  Z <- rbinom(N, 1, PZ1); #Z
  kap <-  (PZ1 - Z) / (PZ1 *(1 - PZ1));
  data.array[, 3, j] <- (1 - D) * kap; #kappa0
  data.array[, 4, j] <- - D * kap; #kappa1
}

pvals <- apply(data.array, 3, solve, B = B);
reject <- matrix(0, J, N_a);
for (i in 1:N_a){
  reject[, i] <- pvals <= a[i];
}
reject_rate <- colMeans(reject);

MC_testsize$SmallMedium$rejection <- reject_rate;
#save(MC_testsize, file="MCsize.Rdata");

############################################ MEDIUM SAMPLE SIZE ############################################

N <- 1000;
data.array <- array(0, dim = c(N, 4, J));

#Generate J Datasets
for (j in 1:J){
  parameter <- runif(2, 0, 1); #B and p parameters
  data.array[, 1, j] <- rnorm(N); #X
  D <- rbinom(N, 1, 1 / 2); #D
  data.array[, 2, j] <- D +  parameter[1] * data.array[, 1, j] + rnorm(N); #Y
  PZ1 <- pnorm(data.array[, 1, j] * parameter[2]);
  Z <- rbinom(N, 1, PZ1); #Z
  kap <-  (PZ1 - Z) / (PZ1 *(1 - PZ1));
  data.array[, 3, j] <- (1 - D) * kap; #kappa0
  data.array[, 4, j] <- - D * kap; #kappa1
}

pvals <- apply(data.array, 3, solve, B = B);
reject <- matrix(0, J, N_a);
for (i in 1:N_a){
  reject[, i] <- pvals <= a[i];
}
reject_rate <- colMeans(reject);

MC_testsize$Medium$rejection <- reject_rate;
#save(MC_testsize, file="MCsize.Rdata");

############################################ SMALL LARGE SAMPLE SIZE ############################################

N <- 1100;
data.array <- array(0, dim = c(N, 4, J));
rho <- qnorm(1/11);

#Generate J Datasets
for (j in 1:J){
  parameter <- runif(2, 0, 1); #B and p parameters
  data.array[, 1, j] <- rnorm(N); #X
  D <- rbinom(N, 1, 1 / 2); #D
  data.array[, 2, j] <- D +  parameter[1] * data.array[, 1, j] + rnorm(N); #Y
  PZ1 <- pnorm(rho + data.array[, 1, j] * parameter[2]);
  Z <- rbinom(N, 1, PZ1); #Z
  kap <-  (PZ1 - Z) / (PZ1 *(1 - PZ1));
  data.array[, 3, j] <- (1 - D) * kap; #kappa0
  data.array[, 4, j] <- - D * kap; #kappa1
}

pvals <- apply(data.array, 3, solve, B = B);
reject <- matrix(0, J, N_a);
for (i in 1:N_a){
  reject[, i] <- pvals <= a[i];
}
reject_rate <- colMeans(reject);

MC_testsize$SmallLarge$rejection <- reject_rate;
#save(MC_testsize, file="MCsize.Rdata");

############################################ LARGE SAMPLE SIZE ############################################

N <- 2000;
data.array <- array(0, dim = c(N, 4, J));

#Generate J Datasets
for (j in 1:J){
  parameter <- runif(2, 0, 1); #B and p parameters
  data.array[, 1, j] <- rnorm(N); #X
  D <- rbinom(N, 1, 1 / 2); #D
  data.array[, 2, j] <- D +  parameter[1] * data.array[, 1, j] + rnorm(N); #Y
  PZ1 <- pnorm(data.array[, 1, j] * parameter[2]);
  Z <- rbinom(N, 1, PZ1); #Z
  kap <-  (PZ1 - Z) / (PZ1 *(1 - PZ1));
  data.array[, 3, j] <- (1 - D) * kap; #kappa0
  data.array[, 4, j] <- - D * kap; #kappa1
}

pvals <- apply(data.array, 3, solve, B = B);
reject <- matrix(0, J, N_a);
for (i in 1:N_a){
  reject[, i] <- pvals <= a[i];
}
reject_rate <- colMeans(reject);

MC_testsize$Large$rejection <- reject_rate;
#save(MC_testsize, file="MCsize.Rdata");
