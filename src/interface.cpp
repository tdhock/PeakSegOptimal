/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "ARFPOP.h"
#include "FitSegmentModel.h"

extern "C" {
  void ARFPOP_interface
  (double *data_ptr,
   int *data_count, double *penalty,
   double *gam,
   double *cost_mat, int *end_vec,
   double *mean_vec, int *intervals_mat, bool *constraint, int *success, bool *compute_fitted_values, double *EPS){
    ARFPOP(data_ptr,
                 *data_count, *penalty, *gam,
                 cost_mat, end_vec, mean_vec, intervals_mat, constraint, success, compute_fitted_values, *EPS);
  }
  
  void FitSegmentModel_interface
  (double *data_vec, int *data_count,
   double *gam, // decay parameter
   int *end_vec, //input changepts
   double *mean_vec,//data_count
   double *EPS){
    FitSegmentModel(data_vec, *data_count, *gam, end_vec, mean_vec, *EPS);
}
}
    
