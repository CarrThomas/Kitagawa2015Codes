#A script for running Monte Carlo exercises with the 2D algorithm plus data.

rm(list=ls());


############################################### OVERALL FUNCTION ###############################################

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

############################################### KADANE'S ALGORITHM ###############################################

kadane <- function(kappa_vec){
  #Kadane's algorithm. Takes a sorted vector and returns the minimum subsequence sum.
  val <- 0;
  min_val <- 0;

  for (i in 1:length(kappa_vec)){

    val <- kappa_vec[i] + min(val, 0);

    if (val < min_val){
      min_val <- val;
    }
  }

  return(min_val);
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
