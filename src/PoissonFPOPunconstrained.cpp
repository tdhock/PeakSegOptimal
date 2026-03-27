/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <vector>
#include <math.h>
#include "funPieceListLog.h"
#include "PoissonFPOPunconstrained.h"

int PoissonFPOPunconstrainedLog
(int *data_vec, double *weight_vec, int data_count,
 double penalty,
 double *cost_vec,
 int *end_vec,
 double *mean_vec,
 int *intervals_vec){
  double min_log_mean=INFINITY, max_log_mean=-INFINITY;
  for(int data_i=0; data_i<data_count; data_i++){
    if(data_vec[data_i] < 0) return ERROR_NEGATIVE_DATA;
    double log_data = log((double)data_vec[data_i]);
    if(log_data < min_log_mean) min_log_mean = log_data;
    if(max_log_mean < log_data) max_log_mean = log_data;
  }
  if(min_log_mean == max_log_mean){
    return ERROR_MIN_MAX_SAME;
  }

  std::vector<PiecewisePoissonLossLog> cost_model_mat(data_count);
  PiecewisePoissonLossLog *cost, *cost_prev;
  PiecewisePoissonLossLog min_prev_cost;
  int verbose = 0;
  double cum_weight_i = 0.0, cum_weight_prev_i;

  for(int data_i=0; data_i<data_count; data_i++){
    cost = &cost_model_mat[data_i];
    cum_weight_i += weight_vec[data_i];
    if(data_i==0){
      // C_1(m) = gamma_1(m)/w_1
      cost->piece_list.emplace_back(
        1.0, -data_vec[0], 0.0,
        min_log_mean, max_log_mean, -1, INFINITY);
    }else{
      // C_t(m) = (gamma_t + w_{1:t-1} * M_t(m)) / w_{1:t}, where
      // M_t(m) = min{ C_{t-1}(m), min_m' C_{t-1}(m') + lambda/w_{1:t-1} }
      min_prev_cost.set_to_unconstrained_min_of(cost_prev);
      min_prev_cost.set_prev_seg_end(data_i-1);
      min_prev_cost.add(0.0, 0.0, penalty/cum_weight_prev_i);
      cost->set_to_min_env_of(&min_prev_cost, cost_prev, verbose);
      cost->multiply(cum_weight_prev_i);
      cost->add(
        weight_vec[data_i],
        -data_vec[data_i]*weight_vec[data_i],
        0.0);
      cost->multiply(1.0/cum_weight_i);
    }
    cum_weight_prev_i = cum_weight_i;
    cost_prev = cost;
  }

  double best_log_mean, prev_log_mean;
  int prev_seg_end;
  for(int i=0; i<data_count; i++){
    cost = &cost_model_mat[i];
    intervals_vec[i] = cost->piece_list.size();
    cost->Minimize(cost_vec+i, &best_log_mean, &prev_seg_end, &prev_log_mean);
  }

  for(int i=0; i<data_count; i++){
    mean_vec[i] = INFINITY;
    end_vec[i] = -2;
  }
  mean_vec[0] = exp(best_log_mean);
  end_vec[0] = prev_seg_end;
  int out_i=1;
  while(0 <= prev_seg_end && out_i < data_count){
    cost = &cost_model_mat[prev_seg_end];
    if(prev_log_mean != INFINITY){
      best_log_mean = prev_log_mean;
    }
    cost->findMean(best_log_mean, &prev_seg_end, &prev_log_mean);
    mean_vec[out_i] = exp(best_log_mean);
    end_vec[out_i] = prev_seg_end;
    out_i++;
  }
  return 0;
}
