/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <stdio.h>
#include "funPieceListLog.h"
#include <math.h>
#include <R.h> // for Rprintf

int UnconstrainedFPOPLog
(const int *data_vec, const double *weight_vec, const int data_count,
 const double penalty, const char* verbose_file,
 // the following vectors are for output.
 // cost_vec and intervals_vec store the optimal cost and number of intervals
 // at each time point, for the up and down cost models.
 // end_vec and mean_vec store the best model up to and including the
 // last data point. 
 double *cost_vec, //data_count
 int *end_vec, //data_count
 double *mean_vec,//data_count
 int *intervals_vec){//data_count
  bool verbose = strcmp(verbose_file, "") != 0;
  std::ofstream verbose_fstream;
  if(verbose){
    verbose_fstream.open(verbose_file);
    verbose_fstream << "data_index" << "\t" << "change" << "\t" << "cost" << "\n";
  }
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
  PiecewisePoissonLossLog this_cost, prev_cost, min_prev_cost;
  int all_end_vec[data_count];
  double all_log_mean_vec[data_count];
  double cum_weight_i = 0.0, cum_weight_prev_i;
  for(int data_i=0; data_i<data_count; data_i++){
    cum_weight_i += weight_vec[data_i];
    if(data_i==0){
      // initialization C_1(m)=gamma_1(m)/w_1
      this_cost.piece_list.emplace_back
	(1.0, -data_vec[0], 0.0,
	 min_log_mean, max_log_mean, -1, false);
    }else{
      min_prev_cost.set_to_unconstrained_min_of(&prev_cost, 0);
      min_prev_cost.set_prev_seg_end(data_i-1);
      min_prev_cost.add(0.0, 0.0, penalty/cum_weight_prev_i);
      this_cost.set_to_min_env_of(&min_prev_cost, &prev_cost, 0);
      this_cost.multiply(cum_weight_prev_i);
      this_cost.add
	(weight_vec[data_i],
	 -data_vec[data_i]*weight_vec[data_i],
	 0.0);
      this_cost.multiply(1/cum_weight_i);
    }//if(data_i initialization else update
    cum_weight_prev_i = cum_weight_i;
    prev_cost = this_cost;
    intervals_vec[data_i] = this_cost.piece_list.size();
    this_cost.Min_maybe_verbose
      (cost_vec+data_i,
       all_log_mean_vec+data_i,
       all_end_vec+data_i,
       data_i,
       verbose_fstream);
  }
  // Decoding, and writing to the output.
  for(int i=0; i<data_count; i++){
    mean_vec[i] = INFINITY;
    end_vec[i] = -2;
  }
  int last_i=data_count-1;
  int out_i=0;
  while(0 <= last_i){
    end_vec[out_i] = last_i;
    double best_log_mean = all_log_mean_vec[last_i];
    mean_vec[out_i] = exp(best_log_mean);
    out_i++;
    last_i = all_end_vec[last_i];
  }
  return 0;
}
