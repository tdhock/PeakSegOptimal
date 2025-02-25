/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <vector>
#include <stdio.h>
#include "funPieceListLog.h"
#include <math.h>
#include <R.h> // for Rprintf

int UnconstrainedFPOPLog
(int *data_vec, double *weight_vec, int data_count,
 double penalty,
 // the following vectors are for output.
 // cost_vec and intervals_vec store the optimal cost and number of intervals
 // at each time point, for the up and down cost models.
 // end_vec and mean_vec store the best model up to and including the
 // last data point. 
 double *cost_vec, //data_count
 int *end_vec, //data_count
 double *mean_vec,//data_count
 int *intervals_vec){//data_count
  double min_log_mean=INFINITY, max_log_mean=-INFINITY;
  for(int data_i=0; data_i<data_count; data_i++){
    double log_data = log( (double)(data_vec[data_i]) );
    if(log_data < min_log_mean){
      min_log_mean = log_data;
    }
    if(max_log_mean < log_data){
      max_log_mean = log_data;
    }
  }
  if(min_log_mean == max_log_mean){
    return ERROR_MIN_MAX_SAME;
  }
  std::vector<PiecewisePoissonLossLog> cost_model_vec(data_count);
  PiecewisePoissonLossLog *this_cost, *prev_cost;
  PiecewisePoissonLossLog min_prev_cost;
  int verbose=0;
  double cum_weight_i = 0.0, cum_weight_prev_i;
  for(int data_i=0; data_i<data_count; data_i++){
    cum_weight_i += weight_vec[data_i];
    this_cost = &cost_model_vec[data_i];
    if(data_i==0){
      // initialization C_1(m)=gamma_1(m)/w_1
      this_cost->piece_list.emplace_back
	(1.0, -data_vec[0], 0.0,
	 min_log_mean, max_log_mean, -1, false);
    }else{
      min_prev_cost.set_to_unconstrained_min_of(prev_cost, verbose);
      min_prev_cost.set_prev_seg_end(data_i-1);
      min_prev_cost.add(0.0, 0.0, penalty/cum_weight_prev_i);
      this_cost->set_to_min_env_of(&min_prev_cost, prev_cost, verbose);
      this_cost->multiply(cum_weight_prev_i);
      this_cost->add
	(weight_vec[data_i],
	 -data_vec[data_i]*weight_vec[data_i],
	 0.0);
      this_cost->multiply(1/cum_weight_i);
    }//if(data_i initialization else update
    cum_weight_prev_i = cum_weight_i;
    prev_cost = this_cost;
  }
  // Decoding the cost_model_vec, and writing to the output.
  double best_cost, best_log_mean, prev_log_mean;
  int prev_seg_end=data_count;
  for(int i=0; i<data_count; i++){
    mean_vec[i] = INFINITY;
    end_vec[i] = -2;
  }
  for(int i=0; i<data_count; i++){
    this_cost = &cost_model_vec[i];
    intervals_vec[i] = this_cost->piece_list.size();
    this_cost->Minimize
      (cost_vec+i, &best_log_mean,
       &prev_seg_end, &prev_log_mean);
  }
  int out_i=0;
  prev_seg_end = data_count-1;
  while(0 <= prev_seg_end){
    this_cost = &cost_model_vec[prev_seg_end];
    this_cost->Minimize
      (&best_cost, &best_log_mean,
       &prev_seg_end, &prev_log_mean);
    mean_vec[out_i] = exp(best_log_mean);
    end_vec[out_i] = prev_seg_end;
    out_i++;
  }//for(data_i
  return 0;
}
