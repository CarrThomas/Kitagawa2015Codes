################################# FRONT MATTER #################################

# Fourth attempt: main idea is to get tighter bounds by exploiting properties

rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)

################################### FUNCTIONS ##################################

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

############################# GENERATE + TIDY DATA #############################

# number of observations
N <- 400

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

################################### MAIN BODY ##################################

# use the code where at each point the minimum is calculated for the interval
# with the minimum bound but
# 1) track which intervals we get the true minimum for
# 2) track which points whave their bound updated at each stage

data <- data0

true_vals <- exhaustive_solution(data)$true_vals

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
incumbent <- 0

flag <- 0
kadane_count <- 0

point_record <- matrix(0, 1000, 2)
update_record <- list()

while(flag == 0){

  # get the current minimum
  i <- which.min(bounds)
  col <- floor(i / n) + 1
  row <- i - (col - 1) * n

  if(bounds[row, col] < incumbent){
    # Pass onto Kadane's algorithm
    kadane_count <- kadane_count + 1
    point_record[kadane_count, ] <- c(row, col)

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
    update_record[[kadane_count]] <- new_bounds > bounds
    bounds <- pmax(bounds, new_bounds)


  }else{
      flag <- 1
  }
}

# Make plots. In each plot show
# 1) points where the bound was updated
# 2) most recent point where the exact minimum was obtained
# 3) all past points where the minimum was calculated

point_record <- point_record[1:kadane_count, ]
point_record <- cbind(point_record, 1:kadane_count)
point_record <- data.frame(point_record)
names(point_record) <- c("row", "col", "label")

for (i in 1:kadane_count){

  temp_data <- data.frame(ind = c(update_record[[i]]),
                        row = rep(1:n, n),
                        col = rep(1:n, each = n))
  print(temp_data %>% ggplot() +
    geom_tile(aes(x = col, y = n + 1 - row, fill = as.factor(ind))) +
    geom_text(aes(x =  col, y = n + 1 - row, label = label),
              data = point_record[1:i, ]))
}

############################# ALTERNATIVE APPROACH #############################

# initial split + updating of bounds

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
og_bounds <- bounds

incumbent <- 0

kadane_count <- 0

point_record <- matrix(0, 1000, 2)
update_record <- list()

# split into three intervals
a <- floor(n / 3)
b <- floor(2 * n / 3)

# solve for the initial areas
init_rows <- 1 + c(0, a, b)
init_cols <- c(a, b, n)

candidates <- rep(0, 3)

for (i in 1:3){

  col <- init_cols[i]
  row <- init_rows[i]

  kadane_count <- kadane_count + 1
  point_record[kadane_count, ] <- c(row, col)

  kappa_vec <- data
  kappa_vec[intervals[2 * row - 1, 2]:intervals[2 * col - 1, 3], 5] <- 1
  kappa_vec[kappa_vec[, 2], ] <- kappa_vec
  kappa_vec <- subset(kappa_vec, kappa_vec[, 5] == 1)
  candidate <- kadane(kappa_vec[, 3])
  candidates[i] <- candidate

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

  update_record[[kadane_count]] <- new_bounds > bounds
  bounds <- pmax(bounds, new_bounds)

}

#### additional updating ####

# 1) region 1 + region 2
row <- 1
col <- b

kadane_count <- kadane_count + 1
point_record[kadane_count, ] <- c(row, col)

candidate <- candidates[1] + candidates[2]

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

update_record[[kadane_count]] <- new_bounds > bounds
bounds <- pmax(bounds, new_bounds)

# 2) region 2 + region 3
row <- a + 1
col <- n

kadane_count <- kadane_count + 1
point_record[kadane_count, ] <- c(row, col)

candidate <- candidates[2] + candidates[3]

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

update_record[[kadane_count]] <- new_bounds > bounds
bounds <- pmax(bounds, new_bounds)

# 3) region 1 + region 2 + region 3
row <- 1
col <- n

kadane_count <- kadane_count + 1
point_record[kadane_count, ] <- c(row, col)

candidate <- candidates[1] + candidates[2] + candidates[3]

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

update_record[[kadane_count]] <- new_bounds > bounds
bounds <- pmax(bounds, new_bounds)

#### Plot ####

point_record_store <- point_record

point_record <- point_record[1:kadane_count, ]
point_record <- cbind(point_record, 1:kadane_count)
point_record <- data.frame(point_record)
names(point_record) <- c("row", "col", "label")

for (i in 1:kadane_count){

  temp_data <- data.frame(ind = c(update_record[[i]]),
                          row = rep(1:n, n),
                          col = rep(1:n, each = n))
  print(temp_data %>% ggplot() +
          geom_tile(aes(x = col, y = n + 1 - row, fill = as.factor(ind))) +
          geom_text(aes(x =  col, y = n + 1 - row, label = label),
                    data = point_record[1:i, ]))
}


point_record <- point_record_store
flag <- 0

while(flag == 0){

  # get the current minimum
  i <- which.min(bounds)
  col <- floor(i / n) + 1
  row <- i - (col - 1) * n

  if(bounds[row, col] < incumbent){
    # Pass onto Kadane's algorithm
    kadane_count <- kadane_count + 1
    point_record[kadane_count, ] <- c(row, col)

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
    update_record[[kadane_count]] <- new_bounds > bounds
    bounds <- pmax(bounds, new_bounds)


  }else{
    flag <- 1
  }
}


#### Plot

point_record <- point_record[1:kadane_count, ]
point_record <- cbind(point_record, 1:kadane_count)
point_record <- data.frame(point_record)
names(point_record) <- c("row", "col", "label")

for (i in 1:kadane_count){

  temp_data <- data.frame(ind = c(update_record[[i]]),
                          row = rep(1:n, n),
                          col = rep(1:n, each = n))
  print(temp_data %>% ggplot() +
          geom_tile(aes(x = col, y = n + 1 - row, fill = as.factor(ind))) +
          geom_text(aes(x =  col, y = n + 1 - row, label = label),
                    data = point_record[1:i, ]))
}
