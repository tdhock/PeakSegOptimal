#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cmath>


double regression_coef(double *data_vec, int segment_start, int start_i, double gam, double EPS) { 
  int length_y = start_i - segment_start + 1; 
  // printf("len y %d\n", length_y);
  double prefactor = (1 - gam * gam) / (1 - pow(gam, (2 * length_y)));
  // printf("prefactor %f\n", prefactor);
  double ss = 0;
  for(int data_i=0; data_i< length_y; data_i++) {
    // printf("data pt %f\n", data_vec[segment_start + data_i]);
    ss += data_vec[segment_start + data_i] * pow(gam, data_i);
  }
  double coef = ss * prefactor;
  if (coef < EPS) { coef = EPS;}
  return(coef);
}

double rss(double *data_vec, int segment_start, int segment_end, double *fitted_values) { 
  int length_y = segment_end - segment_start; 
  // printf("length inside rss %d\n", length_y);
  double rss = 0;
  for(int data_i=0; data_i< length_y; data_i++) {
    // printf("data value %f \t fitted value %f \n", data_vec[data_i + segment_start], fitted_values[data_i]);
    rss += 0.5 * pow((data_vec[data_i + segment_start] - fitted_values[data_i]) , 2);
  }
  return(rss);
}

void update_fitted_values(double *mean_vec, int segment_start, int segment_end, double *fitted_values) {
  int length_y = segment_end - segment_start; 
  for(int data_i=0; data_i< length_y; data_i++) {
    // printf("assigning val \t %f \n", fitted_values[data_i]);
    mean_vec[data_i + segment_start] = fitted_values[data_i];
  }
}

void fit_from_regression(double initial_value, double *fitted_values, int length_y, int length_sub, double gam, double EPS) {
  fitted_values[0] = initial_value; 
  double current_value;
  
  for (int data_i=1; data_i < length_y; data_i++) {
    if (data_i < length_sub) {
      current_value = gam * fitted_values[data_i - 1];  
    } else { 
      current_value = EPS; 
    }
    
    if (current_value < EPS) {
      current_value = EPS;
    }
    fitted_values[data_i] = current_value; 
  }
}

bool check_feasible(double *fitted_values, int length_y, double gam, double EPS) {
  double nxt_val;
  if (length_y == 1) {
    if (fitted_values[0] > EPS) 
      return true; 
  } else { 
    for (int data_i = 1; data_i < length_y; data_i++) {
      nxt_val = fitted_values[data_i - 1] * gam;
      if (nxt_val < EPS && fitted_values[data_i] == EPS) {
        // nothing
      } else if (fitted_values[data_i] == nxt_val) {
        // nothing
      } else { 
        return false; 
      }
    }
    return true; 
  }
}

void printfitted(double *mean_vec, int dat_count) {
  for (int data_i = 0; data_i < dat_count; data_i++) {
    // printf("data at %d \t %f\n", data_i, mean_vec[data_i]);
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
    // printf("data_i, segment_start, segment_end %d, %d, %d \n", data_i, segment_start, segment_end);
    segment_start = end_vec[data_i] + 1;
    // printf("current data_i %d \t current segment_start %d\n", data_i, segment_start);
    if (segment_start < segment_end) { 
      // printf("--- hit a segment (start, end) \t (%d, %d)\n", segment_start, segment_end);
      // do the regression 
      rss_best = std::numeric_limits<double>::infinity();;
      length_segment = segment_end - segment_start;
      // printf("full seg length ---- %d \n", length_segment);
      double fitted_values[length_segment];
      for (int start_i = segment_start; start_i < segment_end; start_i++){
        double initial_value = regression_coef(data_vec, segment_start, start_i, gam, EPS);
         // printf("initial value %f\n", initial_value);
        int length_sub = start_i - segment_start + 1; 
        // printf("len of seg. %d\n", length_sub);
        fit_from_regression(initial_value, fitted_values, length_segment, length_sub, gam, EPS);
        // printfitted(fitted_values, length_segment);
        double rss_now = rss(data_vec, segment_start, segment_end, fitted_values);
        // printf("current rss %f \t best rss %f \n", rss_now, rss_best);
        // printf("--**--\n");
        if (rss_now < rss_best && check_feasible(fitted_values, length_segment, gam, EPS)) {
          // printf("update fits now \n");
          update_fitted_values(mean_vec, segment_start, segment_end, fitted_values);
          rss_best = rss_now; 
        }
      }
    }
     // update the pointers
     segment_end = segment_start; 
    } // increment segment_end until start
  
  // printf("----- fitted values \n");
  // printfitted(mean_vec, data_count);
  }
  
