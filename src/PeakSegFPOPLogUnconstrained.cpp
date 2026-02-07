/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <vector>
#include <stdio.h>
#include "funPieceListLog.h"
#include "PeakSegFPOPLog.h"
#include <math.h>
#include <R.h>

int PeakSegFPOPLogUnconstrained(
    int *data_vec, double *weight_vec, int data_count,
    double penalty,
    double *cost_mat,
    int *end_vec,
    double *mean_vec,
    int *intervals_mat)
{

    double min_log_mean = INFINITY, max_log_mean = -INFINITY;
    for (int data_i = 0; data_i < data_count; data_i++)
    {
        double log_data = log((double)(data_vec[data_i]));
        if (log_data < min_log_mean)
        {
            min_log_mean = log_data;
        }
        if (max_log_mean < log_data)
        {
            max_log_mean = log_data;
        }
    }

    if (min_log_mean == max_log_mean)
    {
        return ERROR_MIN_MAX_SAME;
    }

    std::vector<PiecewisePoissonLossLog> cost_model_mat(data_count);
    PiecewisePoissonLossLog *cost, *cost_prev;
    PiecewisePoissonLossLog min_prev_cost;
    int verbose = 0;
    double cum_weight_i = 0.0, cum_weight_prev_i;

    for(int data_i = 0; data_i < data_count; data_i++){
      cum_weight_i += weight_vec[data_i];
      cost = &cost_model_mat[data_i];
      
      if(data_i == 0){
        // initialization
        cost->piece_list.emplace_back(
            1.0, -data_vec[0], 0.0,
            min_log_mean, max_log_mean, -1, INFINITY
        );
      } else {
        // get unconstrained minimum of previous cost
        min_prev_cost.set_to_unconstrained_min_of(cost_prev, verbose);
        int status = min_prev_cost.check_min_of(cost_prev, cost_prev);
        if(status){
          Rprintf("BAD UNCONSTRAINED MIN CHECK data_i=%d status=%d\n", data_i, status);
          throw status;
        }
        
        // penalty is added here for creating new segment
        min_prev_cost.set_prev_seg_end(data_i - 1);
        min_prev_cost.add(0.0, 0.0, penalty / cum_weight_prev_i);
        
        // Compare: add changepoint vs continue segment
        // This works for ALL data_i >= 1, including data_i == 1
        cost->set_to_min_env_of(&min_prev_cost, cost_prev, verbose);
        status = cost->check_min_of(&min_prev_cost, cost_prev);
        if(status){
          Rprintf("BAD MIN ENV CHECK data_i=%d status=%d\n", data_i, status);
          throw status;
        }
        
        // Add current data point
        cost->multiply(cum_weight_prev_i);
        cost->add(weight_vec[data_i], -data_vec[data_i] * weight_vec[data_i], 0.0);
        cost->multiply(1.0 / cum_weight_i);
      }
      
      cum_weight_prev_i = cum_weight_i;
      cost_prev = cost;
    }

// Decoding
double best_cost, best_log_mean, prev_log_mean;
int prev_seg_end = data_count;

for (int i = 0; i < data_count; i++)
{
    mean_vec[i] = INFINITY;
    end_vec[i] = -2;
}

for (int i = 0; i < data_count; i++)
{
    cost = &cost_model_mat[i];
    intervals_mat[i] = cost->piece_list.size();
    cost->Minimize(&cost_mat[i], &best_log_mean, &prev_seg_end, &prev_log_mean);
}

cost = &cost_model_mat[data_count - 1];
cost->Minimize(&best_cost, &best_log_mean, &prev_seg_end, &prev_log_mean);

mean_vec[0] = exp(best_log_mean);
end_vec[0] = prev_seg_end;

int out_i = 1;
while (0 <= prev_seg_end)
{
    cost = &cost_model_mat[prev_seg_end];

    if (prev_log_mean != INFINITY)
    {
        best_log_mean = prev_log_mean;
    }

    cost->findMean(best_log_mean, &prev_seg_end, &prev_log_mean);
    mean_vec[out_i] = exp(best_log_mean);
    end_vec[out_i] = prev_seg_end;
    out_i++;
}

return 0;
}