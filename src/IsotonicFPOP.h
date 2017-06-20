void IsotonicFPOP
(double *data_vec,
  int data_count, double penalty,
  // the following matrices are for output.
  double *cost_mat,
  int *end_vec,
  double *mean_vec,
  int *intervals_mat, 
  bool *constraint);