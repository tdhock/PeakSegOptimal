/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <vector>
#include <stdio.h>
#include "funPieceListLog.h"
#include <math.h>

void IsotonicFPOP
  (double *data_vec, int data_count,
   double penalty,
   // the following matrices are for output.
   // cost_mat and intervals_mat store the optimal cost and number of intervals
   // at each time point, for the up and down cost models.
   // end_vec and mean_vec store the best model up to and including the
   // last data point. 
   double *cost_mat,
   int *end_vec, //data_count
   double *mean_vec,//data_count
   int *intervals_mat, //data_count
   bool *constraint){
  double min_mean=data_vec[0], max_mean=data_vec[0];
  double EPS = 1e-40;
  for(int data_i=1; data_i<data_count; data_i++){
    double data = data_vec[data_i];
    if(data < min_mean){
      min_mean = data;
    }
    if(max_mean < data){
      max_mean = data;
    }
  }
  
  std::vector<PiecewiseSquareLoss> cost_model_mat(data_count);
  PiecewiseSquareLoss *cost, *cost_prev;
  PiecewiseSquareLoss min_prev_cost;
  int verbose=0;
  for(int data_i=0; data_i<data_count; data_i++){
    cost = &cost_model_mat[data_i];
    // Alg 3, ln 4
    if(data_i==0){
      cost->piece_list.emplace_back
      (1.0, -2 * data_vec[0], data_vec[0] * data_vec[0],
       min_mean, max_mean, -1, false);
    }else{ // Alg 3 ln 6 - 8 
      
      if (*constraint) {
        min_prev_cost.set_to_min_less_of(cost_prev, verbose);  
      } else {
        min_prev_cost.set_to_unconstrained_min_of(cost_prev, EPS, verbose);  
      }

      min_prev_cost.set_prev_seg_end(data_i-1, EPS);
      min_prev_cost.add(0.0, 0.0, penalty);
      cost->set_to_min_env_of(&min_prev_cost, cost_prev, verbose);
      
      if (verbose) {
      printf("new cost - likelihood \n");
      cost -> print();
      }
      
      cost->add
        (1.0,
         -2 * data_vec[data_i], data_vec[data_i] * data_vec[data_i]);
    }
    cost_prev = cost;
  }

  if (verbose) {
  // Print mu functions 
  printf("Printing Cost*_T(mu) functions to get handle on tau*_T(mu) \n");
  cost -> print();
  }
  
  // Decoding the cost_model_vec, and writing to the output matrices.
  double best_cost, best_mean, prev_mean;
  int prev_seg_end=data_count;
  for(int i=0; i<data_count; i++){
    mean_vec[i] = INFINITY;
    end_vec[i] = -2;
  }
  for(int i=0; i< data_count; i++){
    cost = &cost_model_mat[i];
    intervals_mat[i] = cost->piece_list.size();
    cost->Minimize
      (cost_mat+i, &best_mean,
       &prev_seg_end, &prev_mean);
  }
  cost = &cost_model_mat[data_count - 1];
  cost->Minimize
    (&best_cost, &best_mean,
     &prev_seg_end, &prev_mean);
  mean_vec[0] = best_mean;
  end_vec[0] = prev_seg_end;
  int out_i=1;
  while(0 <= prev_seg_end){
    cost = &cost_model_mat[prev_seg_end];
    if(prev_mean < INFINITY){
      //equality constraint inactive
      best_mean = prev_mean;
    }
    cost->findMean
      (best_mean, &prev_seg_end, &prev_mean);
    mean_vec[out_i] = best_mean;
    end_vec[out_i] = prev_seg_end;
    out_i++;
  }
}

