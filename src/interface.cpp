/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "PeakSegPDPA.h"

extern "C" {

  void PeakSegPDPA_interface
  (double *data_ptr, double *weight_ptr,
   int *data_count, int *maxSegments,
   double *cost_mat, int *end_mat, double *mean_mat
   ){
    PeakSegPDPA(data_ptr, weight_ptr, *data_count, *maxSegments,
		cost_mat, end_mat, mean_mat);
  }
  
}
    
