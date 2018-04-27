/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "ARFPOP.h"

extern "C" {
  void ARFPOP_interface
  (double *data_ptr,
   int *data_count, double *penalty,
   double *gam,
   double *cost_mat, int *end_vec,
   double *mean_vec, int *intervals_mat, bool *constraint, int *success, double EPS){
    ARFPOP(data_ptr,
                 *data_count, *penalty, *gam,
                 cost_mat, end_vec, mean_vec, intervals_mat, constraint, success, EPS);
  }
}
    
