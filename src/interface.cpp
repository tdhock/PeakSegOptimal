/* -*- compile-command: "R CMD INSTALL .." -*- */
// R DLL adaptions from https://github.com/tdhock/PeakSegOptimal/blob/master/src/interface.cpp
#include "ARFPOP.h"          // to estimate spikes
#include "FitSegmentModel.h" // to estimate calcium concentration
#include <R.h>
#include <R_ext/Rdynload.h>

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


R_CMethodDef cMethods[] = {
        {"ARFPOP_interface",
                (DL_FUNC) &ARFPOP_interface, 12
        },
        {"FitSegmentModel_interface",
                (DL_FUNC) &FitSegmentModel_interface, 6
        },
        {NULL, NULL}
};

extern "C" {
void R_init_FastLZeroSpikeInference(DllInfo *info) {
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

}