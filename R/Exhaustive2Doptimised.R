################################# FRONT MATTER #################################

library(dplyr)

#A script for running Monte Carlo exercises with the 2D algorithm plus data.

rm(list=ls())
gc()


################################### FUNCTIONS ##################################


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

############################## CHECK OPTIMISATION ##############################

data0 <- subset(data, data[, 3] != 0)
for (k in 1:2){
  data0[, k] <- rank(data0[, k])
}
data1 <- subset(data, data[, 4] != 0)
for (k in 1:2){
  data1[, k] <- rank(data1[, k])
}
data0[data0[, 1], ] <- data0
data1[data1[, 1], ] <- data1
data <- rbind(data0, data1)
N0 <- nrow(data0); #comes in handy later
#Note that data is now always initially ordered by x

data0 <- data0[, 1:3]
data0 <- cbind(data0, data0[, 3] < 0)
data1 <- data1[, c(1,2,4)]
data1 <- cbind(data1, data1[, 3] < 0)

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
# incumbent <- min(incumbent)
incumbent

# "greedy solution" over both dimensions to find the incumbent
alt_incumbent <- rep(0, 4)

data0[data0[, 2], ] <- data0
temp <- kadane_2(data0[, 3])
temp_data <- data0[temp$start:temp$end, ]
temp_data <- temp_data[order(temp_data[, 1]), ]
alt_incumbent[1] <- kadane(temp_data[, 3])

data0[data0[, 1], ] <- data0
temp <- kadane_2(data0[, 3])
temp_data <- data0[temp$start:temp$end, ]
temp_data <- temp_data[order(temp_data[, 2]), ]
alt_incumbent[2] <- kadane(temp_data[, 3])

data1[data1[, 2], ] <- data1
temp <- kadane_2(data1[, 3])
temp_data <- data1[temp$start:temp$end, ]
temp_data <- temp_data[order(temp_data[, 1]), ]
alt_incumbent[3] <- kadane(temp_data[, 3])

data1[data1[, 1], ] <- data1
temp <- kadane_2(data1[, 3])
temp_data <- data1[temp$start:temp$end, ]
temp_data <- temp_data[order(temp_data[, 2]), ]
alt_incumbent[4] <- kadane(temp_data[, 3])

incumbent
alt_incumbent

incumbent_store <- list(og = incumbent, alt = alt_incumbent)

############################# COMPLETELY EXHAUSTIVE ############################

# Get the solution for each possible interval

#Collapse to a set of intervals composed of points of the same sign
ind <- c(1, data0[2:nrow(data0), 4] != data0[1:(nrow(data0) - 1), 4])
bot <- subset(data0[, 1], ind == 1)
top <- subset(data0[, 1], c(ind[2:length(ind)], 1) == 1)
ind <- cumsum(ind)
intervals <- aggregate(data0[, 3], list(ind), sum)[, 2]
intervals <- cbind(intervals, bot, top)

#Remove the first and last intervals if the kappa are positive
if (intervals[1, 1] >= 0){
  intervals <- intervals[2:nrow(intervals), ];
}
if (intervals[nrow(intervals), 1] >= 0){
  intervals <- intervals[1:(nrow(intervals) - 1), ];
}

# Get indices for the negative intervals, and sums for positive and negative intervals
neg_sum <- cumsum(intervals[seq.int(1, nrow(intervals), 2), 1])
pos_sum <- cumsum(intervals[seq.int(nrow(intervals) - 1, 2, -2), 1])
pos_sum <- c(pos_sum[length(pos_sum):1], 0)
data0 <- cbind(data0, rep(0, nrow(data0)))

# Simple bounds (sum of all negative points)
n <- length(neg_sum)
bounds <- matrix(neg_sum, n, n, byrow = T) -
  matrix(c(0, neg_sum[- 1 + 1:n]), n, n)
bounds <- bounds * (bounds < 0)
bounds_mat <- bounds

# exhaustive solution at each point
vals <- matrix(0, n, n)

start_time <- Sys.time()
for (i in 1:n){
  for (j in i:n){

    kappa_vec <- data0
    kappa_vec[intervals[(2 * i) - 1, 2]:intervals[(2 * j) - 1, 3], 5] <- 1
    kappa_vec[kappa_vec[, 2], ] <- kappa_vec
    kappa_vec <- subset(kappa_vec, kappa_vec[, 5] == 1)
    vals[i, j] <- kadane(kappa_vec[, 3])

  }
}
end_time <- Sys.time()
start_time - end_time

# current alternative
bounds <- neg_sum - neg_sum[length(neg_sum)]
incumbent <- min(incumbent)

kadane_calls <- 0

start_time <- Sys.time()
for (i in length(neg_sum):1){
  bounds <- bounds + intervals[(2 * i) - 1, 1];
  for (j in i:length(neg_sum)){
    if (bounds[j] < incumbent){
      kadane_calls <- kadane_calls + 1
      kappa_vec <- data0
      kappa_vec[intervals[(2 * i) - 1, 2]:intervals[(2 * j) - 1, 3], 5] <- 1
      kappa_vec[kappa_vec[, 2], ] <- kappa_vec
      kappa_vec <- subset(kappa_vec, kappa_vec[, 5] == 1)
      candidate <- kadane(kappa_vec[, 3]); #minimum for this box

      if (j > i + 1){#case where points to the left and right can be updated
        new_bounds <- candidate + c(- pos_sum[(i + 1):(j - 1)] + pos_sum[j], neg_sum[j:length(neg_sum)] - neg_sum[j])
        bounds[(i + 1):length(bounds)] <- pmax(bounds[(i + 1):length(bounds)], new_bounds)
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
start_time - end_time

kadane_calls
(n + 1) * n / 2

# different incumbent
bounds <- neg_sum - neg_sum[length(neg_sum)]
incumbent <- min(alt_incumbent)

kadane_calls <- 0

start_time <- Sys.time()
for (i in length(neg_sum):1){
  bounds <- bounds + intervals[(2 * i) - 1, 1];
  for (j in i:length(neg_sum)){
    if (bounds[j] < incumbent){
      kadane_calls <- kadane_calls + 1
      kappa_vec <- data0
      kappa_vec[intervals[(2 * i) - 1, 2]:intervals[(2 * j) - 1, 3], 5] <- 1
      kappa_vec[kappa_vec[, 2], ] <- kappa_vec
      kappa_vec <- subset(kappa_vec, kappa_vec[, 5] == 1)
      candidate <- kadane(kappa_vec[, 3]); #minimum for this box

      if (j > i + 1){#case where points to the left and right can be updated
        new_bounds <- candidate + c(- pos_sum[(i + 1):(j - 1)] + pos_sum[j], neg_sum[j:length(neg_sum)] - neg_sum[j])
        bounds[(i + 1):length(bounds)] <- pmax(bounds[(i + 1):length(bounds)], new_bounds)
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
start_time - end_time

kadane_calls
(n + 1) * n / 2



############################### MODIFIED FUNCTION ##############################

# split into D = 0 and D = 1 subsamples + trim off any positive observations
# at the end

# D = 0
data0 <- subset(data, data[, 3] != 0)[, 1:3]
data0 <- cbind(data0, data0[, 3] < 0)

# trimming
data0 <- data.frame(data0)
temp <- dplyr::filter(data0, V4 == 1)

Ymax <- max(temp$Y)
Ymin <- min(temp$Y)

Xmax <- max(temp$X)
Xmin <- min(temp$X)

data0 <- dplyr::filter(data0, Y <= Ymax & Y >= Ymin & X <= Xmax & X >= Xmin)

# replace values with rank
for (k in 1:2){
  data0[, k] <- rank(data0[, k])
}

# order by X
data0[data0[, 1], ] <- data0

# D = 1
data1 <- subset(data, data[, 4] != 0)[, c(1:2, 4)]
data1 <- cbind(data1, data1[, 3] < 0)

# trimming
data1 <- data.frame(data1)
temp <- dplyr::filter(data1, V4 == 1)

Ymax <- max(temp$Y)
Ymin <- min(temp$Y)

Xmax <- max(temp$X)
Xmin <- min(temp$X)

data1 <- dplyr::filter(data1, Y <= Ymax & Y >= Ymin & X <= Xmax & X >= Xmin)

# replace values with rank
for (k in 1:2){
  data1[, k] <- rank(data1[, k])
}

# order by X
data1[data1[, 1], ] <- data1

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

incumbent_store <- incumbent
incumbent <- min(incumbent)

# Solution for D = 0 data

# Collapse to intervals
ind <- c(1, data0[2:nrow(data0), 4] != data0[1:(nrow(data0) - 1), 4])
bot <- subset(data0[, 1], ind == 1)
top <- subset(data0[, 1], c(ind[2:length(ind)], 1) == 1)
ind <- cumsum(ind)
intervals <- aggregate(data0[, 3], list(ind), sum)[, 2]
intervals <- cbind(intervals, bot, top)

#Remove the first and last intervals if the kappa are positive
if (intervals[1, 1] >= 0){
  intervals <- intervals[2:nrow(intervals), ];
}
if (intervals[nrow(intervals), 1] >= 0){
  intervals <- intervals[1:(nrow(intervals) - 1), ];
}

#Get indices for the negative intervals, and kappa sums for positive and negative intervals
neg_sum <- cumsum(intervals[seq.int(1, nrow(intervals), 2), 1])
pos_sum <- cumsum(intervals[seq.int(nrow(intervals) - 1, 2, -2), 1])
pos_sum <- c(pos_sum[length(pos_sum):1], 0)
bounds <- neg_sum - neg_sum[length(neg_sum)]
data <- cbind(data, rep(0, nrow(data)))

# Construct initial bounds
n <- length(neg_sum)
bounds <- matrix(neg_sum, n, n, byrow = T) -
  matrix(c(0, neg_sum[-1 + 1:n]), n, n)
bounds <- bounds * upper.tri(bounds, diag = T)

sum(bounds < incumbent)


# Main loop:

# loop over rows, starting with row 1
# - get the minimum possible value for the row
# - if minimum < incumbent, do Kadane and update
# - if minumum > incumbent, move on to the next row


# top row only
flag <- 0
kadane_calls <- 0
bounds <- neg_sum
n <- length(bounds)
# loop over row 1:
while (flag == 0){

  i <- which.min(bounds)

  if (bounds[i] < incumbent){

    kadane_calls <- kadane_calls + 1

    kappa_vec <- data
    kappa_vec[1:intervals[(2 * i) - 1, 3], 5] <- 1
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

    flag <- 1

  }

}
# all rows
flag <- 0
kadane_calls <- 0
bounds <- neg_sum
n <- length(bounds)
j_max <- n

j <- 1 # current row

# loop over row 1:
while (flag == 0){

  i <- which.min(bounds)

  if (bounds[i] < incumbent){
    # need to check interval on current row

    kadane_calls <- kadane_calls + 1

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

kadane_calls

################################ GLOBAL MINIMUM ################################

# At each stage update all bounds
# Find the minimum over all current bounds
# If this is below the incumbent, do kadane + update
# otherwise we are done

# all rows
flag <- 0
kadane_calls <- 0

alt_neg_sum <- cumsum(intervals[seq.int(nrow(intervals), 1, -2), 1])
alt_neg_sum <- rev(alt_neg_sum)

n <- length(neg_sum)
bounds <- matrix(neg_sum, n, n, byrow = T) -
  matrix(c(0, neg_sum[1:(n - 1)]), n, n)
bounds[lower.tri(bounds)] <- 0

while(flag == 0){

  # get the current minimum
  i <- which.min(bounds)
  col <- floor(i / n) + 1
  row <- i - (col - 1) * n

  if(bounds[row, col] < incumbent){
    # Pass onto Kadane's algorithm
    kadane_calls <- kadane_calls + 1

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

kadane_calls


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
