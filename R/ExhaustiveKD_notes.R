#Need to revise the values for the bounds as much as possible. 

#Assume we have the value at point (i,j)

#For any point (l,m) with l <= i and m >= j, we can get a new bound by subtracting the current value
bounds[1:i,j:ncol(bounds)] <- bounds[1:i,j:ncol(bounds)] - bounds[i, j] + solution

#For points to the left
#Extract points to the left whose values can still be revised
left_replace <- bounds[i, (i + 1):(j - 1)];
#Subtract this from all rows above
bounds[1:i, (i + 1):(j - 1)] <- bounds[1:i, (i + 1):(j - 1)] - left_replace;

#initial conditions for the replacement value
left_edge <- j; #the left-side boundary of the values which will be replaced
replace_value <- solution; #the replacement value to be recursively updated

#Get the replacement values
for (l in 1:length(left_replace)){
  replace_value <- replace_value + positive; #add the additional positive point
  if (left_replace[] < replace_value){
    left_replace <- replace_value
    left_edge <- left_edge - 1;
  }
  else{
    break;
  }
}
#Put in the replacement values
bounds[1:i, (i + 1):(j - 1)] <- bounds[1:i, (i + 1):(j - 1)] + left_replace;

#For points below

