#define ERROR_NEGATIVE_DATA 2

int PoissonFPOPunconstrainedLog
(int *data_vec, double *weight_vec,
 int data_count, double penalty,
 double *cost_vec,
 int *end_vec,
 double *mean_vec,
 int *intervals_vec);
