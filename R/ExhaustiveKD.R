#Exhaustive search method for an arbitrary number of dimensions

rm(list=ls());
set.seed(1);

##############################################  GENERATE DATA   ##############################################

N <- 1000; #number of observations
K <- 2; #number of x variables (total dimensions will be K + 1)

parameters <- runif(2 * K, 0, 1);
data <- matrix(0, N, K + 3);

data[, 1:K] <- rnorm((K*N)); #x variables
PZ1 <- pnorm(data[, 1:K] %*% parameters[1:K]);
Z <- rbinom(N, 1, PZ1);
D <- rbinom(N, 1, 1 / 2);
data[, K + 1] <- D + data[, 1:K] %*% parameters[(K + 1):(2 * K)] + rnorm(N);#Y
kap <- (PZ1 - Z) / (PZ1 * (1 - PZ1));
data[, K + 2] <- (1 - D) * kap; #k0
data[, K + 3] <- - D * kap; #k1


############################################# SUBSET AND ADD COLUMNs #############################################

data0 <- data[, 1:(K + 2)]; #remove the kappa1 column
data0 <- subset(data0, data0[, K + 2] != 0); #subset of relevant points
for (k in 1:(K + 1)){
  data0[, k] <- rank(data0[, k]); #rank by each x variable and by y
}
data0 <- cbind(data0, data0[, K + 2] < 0); #indicator for kappa0 < 0

data1 <- data[, c(1:(K + 1), K + 3)]; #remove the kappa0 column
data1 <- subset(data1, data1[, K + 2] != 0); #subset of relevant points
for (k in 1:(K + 1)){
  data1[, k] <- rank(data1[, k]); #rank by each x variable and by y
}
data1 <- cbind(data1, data1[, K + 2] < 0); #indicator for kappa1 < 0

################################################# BUILD A FUNCTION #################################################

#Sove the 1D problem using Kadane's algorithm. In Kadane's algorithm, we move along an ordered vector updating a sum
#as we go. At point i, the new sum is given by the minimum of i + sum(i - 1) and i. That is, if the sum up to this
#point is negative, it is included. If it is positive, we essentially start again.

kadane <- function(data){
  #Assumes inputs is a vector sorted along the relevant dimension

  #Kadane's algorithm
  val <- data[1];
  bot <- 1;

  min_val <- min(0, val);

  for (i in 1:nrow(data)){
    val <- data[i] + min(0, val);
    if (val < min_val){
      min_val <- val;
    }
  }
  return(min_val);
}

############################################# UNIDIMENSIONAL SOLUTION #############################################

#Solve each problem `as if' you were only allowed to choose bounds for one variable.
#Although we solve over one dimension only, we check that there are no points where kappa is positive that can be
#excluded by modifying the bounds in the other dimension, while not excluding any points where kappa is negative.

#Start with the subset where kappa0 is not 0.
incumbent <- rep(0, (2 * (K + 1)));
for (k in 1:(K + 1)){
  data0[data0[, k], ] <- data0;
  incumbent[k] <- kadane(data0[, K + 2]);
  data1[data1[, k], ] <- data1;
  incumbent[K + 1 + k] <- kadane(data1[, K + 2]);
}
incumbent <- min(incumbent);

############################################## ENUMERATE FUNCTION ##############################################


#Rewritten version of the enumerate function that is hopefully more efficient.
#Makes two major changes
#i) merge neighbouring positive points at all stages of the enumeration process. Can do this as including the points
#together does not introduce any additional positive points
#ii) more efficient revision of minimum bounds when a new solution is obtained (that is, all bounds are updated
#whenever a new solution is obtained, this keeps the minimum for an area more tightly bounded)

enumerate <- function(data, k, K, incumbent){

  data[data[, k], ] <- data; #order along dimension k

  #pass to another function to format the data
  tempdata <- intervals(data, K);


  #Generate an initial matrix of lower bounds equal to the sum of negative points in a region
  bounds <- cumsum(tempdata$neg);
  bounds <- matrix(bounds, length(bounds), length(bounds), T) - matrix(c(0, bounds[1:(length(bounds) - 1)]), length(bounds), length(bounds));
  #keep track of indices
  top_ind <- matrix(max_ind, length(max_ind), length(max_ind), T);
  bot_ind <- matrix(min_ind, length(min_ind), length(min_ind));

  #loop over rows of the bounds matrix, starting at the bottom
  for (i in 1:(nrow(bounds) - 1)){#loop over rows
    for (j in (nrow(bounds) - i + 1):ncol(bounds)){ #loop over columns within a row
      #check whether we can beat the incumbent
      if (bounds[nrow(bounds) - i, j] < incumbent[1]){ #interval which can beat the incumbent
        #"take out" all points outside of the current interval
        bot <- bot_ind[nrow(bounds) - i, j];
        top <- top_ind[nrow(bounds) - i, j];
        inputdata <- data;
        if (bot > 1){
          inputdata[1:(bot - 1), c(K + 2, K + 3)] <- 0;
        }
        if (top < nrow(data)){
          inputdata[(top + 1):nrow(data), c(K + 2, K + 3)] <- 0;
        }
        #either pass to the minimisation algorithm, or to another round of enumeration
        if (k < K){
          incumbent <- enumerate(inputdata, k + 1, K, incumbent);


          #tempsum[(j - i + 1):length(tempsum)] <- tempsum[(j - i + 1):length(tempsum)] - tempsum[j - i + 1] + incumbent[2]; #update bounds
        }
        else{
          incumbent[3] <- searchfinal(inputdata, k + 1, K);
          #tempsum[(j - i + 1):length(tempsum)] <- tempsum[(j - i + 1):length(tempsum)] - tempsum[j - i + 1] + incumbent[3]; #update bounds
          #incumbent[2] <- min(incumbent[2], tempsum[length(tempsum)]);
          #incumbent[1] <- min(incumbent[1], incumbent[3]);
        }
      }
    }
  }
  return(incumbent)
}

############################################## INTERVALS FUNCTION ##############################################

#A child function of enumerate. Formats the data to aid searching.

intervals <- function(data, K){

  #Aggregate to intervals of the same sign
  tempdata <- data[, c(k, K + 2, K + 3)];
  ind <- tempdata[, 3] != c(0, tempdata[1:(nrow(tempdata) - 1), 3]);
  ind <- cumsum(ind);
  max_ind <- aggregate(tempdata[, 1], list(ind), max)[, 2];
  min_ind <- aggregate(tempdata[, 1], list(ind), min)[, 2];
  tempdata <- aggregate(tempdata[, 2:3], list(ind), sum)[, 2:3];

  #Remove the first and last intervals if they contain positive kappa
  if (tempdata[1, 2] == 0){
    tempdata <- tempdata[2:nrow(tempdata), ];
    max_ind <- max_ind[2:length(max_ind)];
    min_ind <- min_ind[2:length(min_ind)];
  }
  if (tempdata[nrow(tempdata), 2] == 0){
    tempdata <- tempdata[1:(nrow(tempdata) - 1), ];
    max_ind <- max_ind[1:(length(max_ind) - 1)];
    min_ind <- min_ind[2:(length(min_ind) - 1)];
  }

  #Seperate out the positive and negative intervals
  temp_neg <- tempdata[seq(1,nrow(tempdata), 2), 1];
  max_ind <- max_ind[c(seq(1,length(max_ind), 2))];
  min_ind <- min_ind[c(seq(1,length(min_ind), 2))];
  temp_pos <- tempdata[seq(2, nrow(tempdata) - 1, 2), 1];
  #entry k in temp pos corresponds to the positive points between intervals k and k + 1 in the negative vector

  return(list('pos' = temp_pos, 'neg' = temp_neg, 'max_ind' = max_ind, 'min_ind' = min_ind));
}

############################################## UPDATE BOUNDS FUNCTION ##############################################

#A function used to update bounds after we have solved for the minimum at a particular point.

bound_update <- function(bounds, candidate, i, j,){

  #Assume we have solved at point (i,j).

  #Every point (l,m) with l<=i and m>=j can be easily updated
  bounds[1:i,j:ncol(bounds)] <- bounds[1:i,j:ncol(bounds)] - bounds[i,j] + candidate;

  #Points to with l<=i and m<j can sometimes be updated
  #For the bound for an interval to be updated, it must that the interval overlaps with the solution in some dimension
  #Need to consider the elements of row i starting from i+1 and ending at j-1.
  left_update <- bounds[i, (i + 1):(j - 1)];
  for (k in 1:length(left_updata)){
    #loop over each element and see if updating makes sense

  }

}


############################################## K DIMENSIONAL SOLUTION ##############################################

#For now, we solve over dimensions in an arbitrary order (just loop over them in column order)
#Need to make a two element vector, first element should be the incumbent value, second can be the candidate value
incumbent <- c(incumbent, 0, 0);
incumbent <- enumerate(data0, 1, K, incumbent);
incumbent <- enumerate(data1, 1, K, incumbent);


###################################################### NOTES ######################################################

#GUROBI solutions are -89.0871 for kappa0 and -105.353 for kappa1 when N = 1000 K = 2.
#Takes about 20 minutes to solve this way.
#GUROBI solutions are -17.571 for kappa0 and -16.7292 for kappa 1 when N = 100 K = 2

#N = 500 K = 2 solution take about 8 minutes and 10 seconds arrives at -64.3
#Gurobi takes about a minute, but arrives at the same solution
