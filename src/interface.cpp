/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "PeakSegPDPA.h"

extern "C" {

  void PeakSegPDPALog_interface
  (int *data_ptr, int *weight_ptr,
   int *data_count, int *maxSegments,
   double *cost_mat, int *end_mat,
   double *mean_mat, int *intervals_mat
   ){
    PeakSegPDPALog(data_ptr, weight_ptr, *data_count, *maxSegments,
		   cost_mat, end_mat, mean_mat, intervals_mat);
  }
  
  void PeakSegPDPA_interface
  (int *data_ptr, int *weight_ptr,
   int *data_count, int *maxSegments,
   double *cost_mat, int *end_mat,
   double *mean_mat, int *intervals_mat
   ){
    PeakSegPDPA(data_ptr, weight_ptr, *data_count, *maxSegments,
		   cost_mat, end_mat, mean_mat, intervals_mat);
  }
  
}
    
