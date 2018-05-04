#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cmath>


double regression_coef(double *data_vec, int segment_start, int segment_end, int start_i, double gam, double EPS, double *ss) { 
  int length_y = start_i - segment_start + 1; 
  double prefactor = (gam * gam - 1) / (pow(gam, 2 * (1 - length_y)) * (pow(gam, 2 * length_y) - 1));
  *ss = *ss / gam + data_vec[segment_start + length_y - 1];
  double coef = *ss * prefactor;
  if (coef > (EPS / gam) && start_i < (segment_end - 1)) {  // not at the last segment and must inforce continunity constraint
    coef = EPS;
  } else if (coef < EPS && start_i == (segment_end - 1)) { // basically unconstrained regression on the whole segment 
    coef = EPS;
  }
  
  return(coef);
}

double rss(double *data_vec, int segment_start, int segment_end, double *fitted_values) { 
  int length_y = segment_end - segment_start; 
  double rss = 0;
  for(int data_i=0; data_i< length_y; data_i++) {
    rss += 0.5 * pow((data_vec[data_i + segment_start] - fitted_values[data_i]) , 2);
  }
  return(rss);
}

void update_fitted_values(double *mean_vec, int segment_start, int segment_end, double *fitted_values) {
  int length_y = segment_end - segment_start; 
  for(int data_i=0; data_i< length_y; data_i++) {
    mean_vec[data_i + segment_start] = fitted_values[data_i];
  }
}

void fit_from_regression(double end_value, double *fitted_values, int length_y, int length_sub, double gam, double EPS) {
  fitted_values[length_sub - 1] = end_value; 
  double current_value;
  
  for (int data_i= (length_sub - 2); data_i > -1; data_i--) {
    fitted_values[data_i] = fitted_values[data_i + 1] / gam;
  }
  for (int data_i = length_sub; data_i < length_y; data_i++) {
    fitted_values[data_i] = EPS;
  }
  
}

void printfitted(double *mean_vec, int dat_count) {
  for (int data_i = 0; data_i < dat_count; data_i++) {
  }
}


void FitSegmentModel
  (double *data_vec, int data_count,
   double gam, // decay parameter
   int *end_vec, //input changepts
   double *mean_vec,//data_count
   double EPS){
  
  int segment_start = end_vec[0];
  int segment_end = data_count;
  double rss_best;
  int length_segment; 
  
  for(int data_i=0; data_i< data_count; data_i++){
    segment_start = end_vec[data_i] + 1;
    if (segment_start < segment_end) { 
      // do the regression 
      rss_best = std::numeric_limits<double>::infinity();
      length_segment = segment_end - segment_start;
      double ss = 0; 
      double fitted_values[length_segment];
      for (int start_i = segment_start; start_i < segment_end; start_i++){
        double end_value = regression_coef(data_vec, segment_start, segment_end, start_i, gam, EPS, &ss);
        int length_sub = start_i - segment_start + 1;
        fit_from_regression(end_value, fitted_values, length_segment, length_sub, gam, EPS);
        double rss_now = rss(data_vec, segment_start, segment_end, fitted_values);
        if (rss_now < rss_best) {
          update_fitted_values(mean_vec, segment_start, segment_end, fitted_values);
          rss_best = rss_now; 
        }
      }
    }
     // update the pointers
     segment_end = segment_start; 
    } // increment segment_end until start
  }
  
