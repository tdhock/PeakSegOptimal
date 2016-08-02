void PeakSegFPOPLog
(int *data_vec, double *weight_vec,
 int data_count, double penalty,
 // the following matrices are for output.
 double *cost_mat,
 int *end_vec,
 double *mean_vec,
 int *intervals_mat);

void PeakSegFPOPall
(int *data_vec, double *weight_vec,
 int data_count, double penalty,
 // the following matrices are for output.
 double *cost_mat,
 int *end_vec,
 double *mean_vec,
 int *intervals_mat,
  int *label_vec);
