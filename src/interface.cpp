/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "PeakSegPDPALog.h"
#include "PeakSegFPOPLog.h"
#include <R.h>
#include <R_ext/Rdynload.h>

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

R_CMethodDef cMethods[] = {
  {"PeakSegPDPALog_interface",
   (DL_FUNC) &PeakSegPDPALog_interface, 8
   //,{REALSXP, REALSXP, INTSXP, INTSXP, REALSXP}
  },
  {"PeakSegPDPAInf_interface",
   (DL_FUNC) &PeakSegPDPAInf_interface, 8
   //,{REALSXP, REALSXP, INTSXP, INTSXP, REALSXP}
  },
  {"PeakSegFPOPLog_interface",
   (DL_FUNC) &PeakSegFPOPLog_interface, 8
   //,{INTSXP, REALSXP, REALSXP, INTSXP}
  },
  {NULL, NULL, 0}
};

extern "C" {
  void R_init_coseg(DllInfo *info) {
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    //R_useDynamicSymbols call says the DLL is not to be searched for
    //entry points specified by character strings so .C etc calls will
    //only find registered symbols.
    R_useDynamicSymbols(info, FALSE);
  }

}
    
