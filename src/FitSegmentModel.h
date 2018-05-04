#include <vector>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cmath>

void FitSegmentModel
(double *data_vec, int data_count,
  double gam, // decay parameter
  int *end_vec, //input changepts
  double *mean_vec,//data_count
  double EPS);

double regression_coef(double *data_vec, int segment_start, int segment_end, int start_i, double gam, double EPS, double *ss);

double rss(double *data_vec, int segment_start, int segment_end, double *fitted_values);

void update_fitted_values(double *mean_vec, int segment_start, int segment_end, double *fitted_values);

void fit_from_regression(double end_value, double *fitted_values, int length_y, int length_sub, double gam, double EPS);