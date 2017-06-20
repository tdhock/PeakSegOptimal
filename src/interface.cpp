/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "IsotonicFPOP.h"

extern "C" {
 void IsotonicFPOP_interface
  (double *data_ptr,
   int *data_count, double *penalty,
   double *cost_mat, int *end_vec,
   double *mean_vec, int *intervals_mat, bool *constraint){
    IsotonicFPOP(data_ptr,
                  *data_count, *penalty,
                  cost_mat, end_vec, mean_vec, intervals_mat, constraint);
 }
}
    
