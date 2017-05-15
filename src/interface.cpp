/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "PeakSegPDPALog.h"
#include "PeakSegFPOPLog.h"

#include "IsotonicFPOP.h"

extern "C" {

  void PeakSegPDPALog_interface
  (int *data_ptr, double *weight_ptr,
   int *data_count, int *maxSegments,
   double *cost_mat, int *end_mat,
   double *mean_mat, int *intervals_mat
   ){
    PeakSegPDPALog(data_ptr, weight_ptr, *data_count, *maxSegments,
		   cost_mat, end_mat, mean_mat, intervals_mat);
  }
  
  void PeakSegPDPAInf_interface
  (int *data_ptr, double *weight_ptr,
   int *data_count, int *maxSegments,
   double *cost_mat, int *end_mat,
   double *mean_mat, int *intervals_mat
   ){
    PeakSegPDPAInf(data_ptr, weight_ptr, *data_count, *maxSegments,
		   cost_mat, end_mat, mean_mat, intervals_mat);
  }
  
  void PeakSegFPOPLog_interface
  (int *data_ptr, double *weight_ptr,
   int *data_count, double *penalty,
   double *cost_mat, int *end_vec,
   double *mean_vec, int *intervals_mat){
    PeakSegFPOPLog(data_ptr, weight_ptr,
		   *data_count, *penalty,
		   cost_mat, end_vec, mean_vec, intervals_mat);
  }
 
 void IsotonicFPOP_interface
  (int *data_ptr, double *weight_ptr,
   int *data_count, double *penalty,
   double *cost_mat, int *end_vec,
   double *mean_vec, int *intervals_mat){
    IsotonicFPOP(data_ptr, weight_ptr,
                  *data_count, *penalty,
                  cost_mat, end_vec, mean_vec, intervals_mat);
 }
  
}
    
