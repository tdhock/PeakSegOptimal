/* -*- compile-command: "R CMD INSTALL .." -*- */
  
  #include <vector>
  #include <stdio.h>
  #include "funPieceListLog.h"
  #include <math.h>
  
  void IsotonicFPOP
(int *data_vec, double *weight_vec, int data_count,
  double penalty,
  // the following matrices are for output.
  // cost_mat and intervals_mat store the optimal cost and number of intervals
  // at each time point, for the up and down cost models.
  // end_vec and mean_vec store the best model up to and including the
  // last data point. 
  double *cost_mat, //data_count x 2.
  int *end_vec, //data_count
  double *mean_vec,//data_count
  int *intervals_mat){//data_count x 2
  double min_mean=data_vec[0], max_mean=data_vec[0];
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
  PiecewiseSquareLoss *up_cost, *up_cost_prev;
  PiecewiseSquareLoss min_prev_cost;
  int verbose=1;
  double cum_weight_i = 0.0, cum_weight_prev_i;
  for(int data_i=0; data_i<data_count; data_i++){
    cum_weight_i += weight_vec[data_i];
    up_cost = &cost_model_mat[data_i];
    // Alg 3, ln 4
    if(data_i==0){
      // initialization Cdown_1(m)=gamma_1(m)/w_1
      up_cost->piece_list.emplace_back
      (1.0, -data_vec[0], 0.0,
        min_mean, max_mean, -1, false);
    }else{ // Alg 3 ln 6 - 8 
      // if data_i is up, it could have come from down_cost_prev.
      min_prev_cost.set_to_min_less_of(up_cost_prev, verbose);
      int status = min_prev_cost.check_min_of(up_cost_prev, up_cost_prev);
      if(status){
        printf("BAD MIN LESS CHECK data_i=%d status=%d\n", data_i, status);
        min_prev_cost.set_to_min_less_of(up_cost_prev, true);
        printf("=prev down cost\n");
        up_cost_prev->print();
        printf("=min less(prev down cost)\n");
        min_prev_cost.print();
        throw status;
      }
      // C^up_t(m) = (gamma_t + w_{1:t-1} * M^up_t(m))/w_{1:t}, where
      // M^up_t(m) = min{
        //   C^up_{t-1}(m),
        //   C^{<=}_down_{t-1}(m) + lambda/w_{1:t-1}
        // in other words, we need to divide the penalty by the previous cumsum,
        // and add that to the min-less-ified function, before applying the min-env.
        min_prev_cost.set_prev_seg_end(data_i-1);
        min_prev_cost.add(0.0, 0.0, penalty/cum_weight_prev_i);
        if(data_i==1){
          *up_cost = min_prev_cost;
        }else{
          up_cost->set_to_min_env_of(&min_prev_cost, up_cost_prev, verbose);
          status = up_cost->check_min_of(&min_prev_cost, up_cost_prev);
          if(status){
            printf("BAD MIN ENV CHECK data_i=%d status=%d\n", data_i, status);
            up_cost->set_to_min_env_of(&min_prev_cost, up_cost_prev, true);
            printf("=min less(prev down cost) + %f\n", penalty);
            min_prev_cost.print();
            printf("=prev up cost\n");
            up_cost_prev->print();
            printf("=new up cost model\n");
            up_cost->print();
            throw status;
          }
        }
        up_cost->multiply(cum_weight_prev_i);
        up_cost->add
        (weight_vec[data_i],
          -data_vec[data_i]*weight_vec[data_i],
          0.0);
        up_cost->multiply(1/cum_weight_i);
    }
    cum_weight_prev_i = cum_weight_i;
    up_cost_prev = up_cost;
    // Decoding the cost_model_vec, and writing to the output matrices.
    double best_cost, best_mean, prev_mean;
    int prev_seg_end=data_count;
    for(int i=0; i<data_count; i++){
      mean_vec[i] = INFINITY;
      end_vec[i] = -2;
    }
    for(int i=0; i<2*data_count; i++){
      up_cost = &cost_model_mat[i];
      intervals_mat[i] = up_cost->piece_list.size();
      up_cost->Minimize
      (cost_mat+i, &best_mean,
        &prev_seg_end, &prev_mean);
    }
    // last segment is down (offset N) so the second to last segment is
    // up (offset 0).
    // int prev_seg_offset = 0;
    // down_cost = &cost_model_mat[data_count*2-1];
    // down_cost->Minimize
    //   (&best_cost, &best_log_mean,
    //    &prev_seg_end, &prev_log_mean);
    // mean_vec[0] = exp(best_log_mean);
    // end_vec[0] = prev_seg_end;
    // int out_i=1;
    // while(0 <= prev_seg_end){
    //   // up_cost is actually either an up or down cost.
    //   up_cost = &cost_model_mat[prev_seg_offset + prev_seg_end];
    //   //printf("decoding out_i=%d prev_seg_end=%d prev_seg_offset=%d\n", out_i, prev_seg_end, prev_seg_offset);
    //   //up_cost->print();
    //   if(prev_log_mean != INFINITY){
    //     //equality constraint inactive
    //     best_log_mean = prev_log_mean;
    //   }
    //   up_cost->findMean
    //     (best_log_mean, &prev_seg_end, &prev_log_mean);
    //   mean_vec[out_i] = exp(best_log_mean);
    //   end_vec[out_i] = prev_seg_end;
    //   // change prev_seg_offset and out_i for next iteration.
    //   if(prev_seg_offset==0){
    //     //up_cost is actually up
    //     prev_seg_offset = data_count;
    //   }else{
    //     //up_cost is actually down
    //     prev_seg_offset = 0;
    //   }
    //   out_i++;
    }
  }
  