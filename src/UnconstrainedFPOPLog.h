int UnconstrainedFPOPLog
(int *data_vec, double *weight_vec, int data_count,
 double penalty,
 // the following vectors are for output.
 // cost_vec and intervals_vec store the optimal cost and number of intervals
 // at each time point, for the up and down cost models.
 // end_vec and mean_vec store the best model up to and including the
 // last data point. 
 double *cost_vec, //data_count
 int *end_vec, //data_count
 double *mean_vec,//data_count
 int *intervals_vec);
