/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <vector>
#include <stdio.h>
#include "funPieceListLog.h"
#include <math.h>

void PeakSegFPOPLog
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
  double min_log_mean=log(data_vec[0]), max_log_mean=log(data_vec[0]);
  for(int data_i=1; data_i<data_count; data_i++){
    double log_data = log(data_vec[data_i]);
    if(log_data < min_log_mean){
      min_log_mean = log_data;
    }
    if(max_log_mean < log_data){
      max_log_mean = log_data;
    }
  }
  std::vector<PiecewisePoissonLossLog> cost_model_mat(data_count * 2);
  PiecewisePoissonLossLog *up_cost, *down_cost, *up_cost_prev, *down_cost_prev;
  PiecewisePoissonLossLog min_prev_cost;
  int verbose=0;
  double cum_weight_i = 0.0, cum_weight_prev_i;
  for(int data_i=0; data_i<data_count; data_i++){
    cum_weight_i += weight_vec[data_i];
    up_cost = &cost_model_mat[data_i];
    down_cost = &cost_model_mat[data_i + data_count];
    if(data_i==0){
      // initialization Cdown_1(m)=gamma_1(m)/w_1
      down_cost->piece_list.emplace_back
	(1.0, -data_vec[0], 0.0,
	 min_log_mean, max_log_mean, -1, false);
    }else{
      // if data_i is up, it could have come from down_cost_prev.
      // if(data_i==8323){
      // 	printf("computing cost data_i=%d\n", data_i);
      // 	verbose=1;
      // }else{
      // 	verbose=0;
      // }
      min_prev_cost.set_to_min_less_of(down_cost_prev, verbose);
      int status = min_prev_cost.check_min_of(down_cost_prev, down_cost_prev);
      if(status){
	printf("BAD MIN LESS CHECK data_i=%d status=%d\n", data_i, status);
	printf("=prev down cost\n");
	down_cost_prev->print();
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
	  printf("=prev down cost\n");
	  down_cost_prev->print();
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
      // compute down_cost.
      if(data_i==1){
	//for second data point, the cost is only a function of the
	//previous down cost (there is no first up cost).
	*down_cost = *down_cost_prev;
      }else{
	// if data_i is down, it could have come from up_cost_prev.
	min_prev_cost.set_to_min_more_of(up_cost_prev, verbose);
	status = min_prev_cost.check_min_of(up_cost_prev, up_cost_prev);
	if(status){
	  printf("BAD MIN MORE CHECK data_i=%d status=%d\n", data_i, status);
	  printf("=prev up cost\n");
	  up_cost_prev->print();
	  printf("=min more(prev up cost)\n");
	  min_prev_cost.print();
	  throw status;
	}
	min_prev_cost.set_prev_seg_end(data_i-1);
	//NO PENALTY FOR DOWN CHANGE
	down_cost->set_to_min_env_of(&min_prev_cost, down_cost_prev, verbose);
	status = down_cost->check_min_of(&min_prev_cost, down_cost_prev);
	if(status){
	  printf("BAD MIN ENV CHECK data_i=%d status=%d\n", data_i, status);
	  printf("=prev up cost\n");
	  up_cost_prev->print();
	  printf("=min more(prev up cost)\n");
	  min_prev_cost.print();
	  printf("=prev down cost\n");
	  down_cost_prev->print();
	  printf("=new down cost model\n");
	  down_cost->print();
	  throw status;
	}
      }
      down_cost->multiply(cum_weight_prev_i);
      down_cost->add
	(weight_vec[data_i],
	 -data_vec[data_i]*weight_vec[data_i],
	 0.0);
      down_cost->multiply(1/cum_weight_i);
    }//if(data_i initialization else update
    cum_weight_prev_i = cum_weight_i;
    up_cost_prev = up_cost;
    down_cost_prev = down_cost;
  }
  // Decoding the cost_model_vec, and writing to the output matrices.
  double best_cost, best_log_mean, prev_log_mean;
  int prev_seg_end=data_count;
  for(int i=0; i<data_count; i++){
    mean_vec[i] = INFINITY;
    end_vec[i] = -2;
  }
  for(int i=0; i<2*data_count; i++){
    up_cost = &cost_model_mat[i];
    intervals_mat[i] = up_cost->piece_list.size();
    up_cost->Minimize
      (cost_mat+i, &best_log_mean,
       &prev_seg_end, &prev_log_mean);
  }
  // last segment is down (offset N) so the second to last segment is
  // up (offset 0).
  int prev_seg_offset = 0;
  down_cost = &cost_model_mat[data_count*2-1];
  down_cost->Minimize
    (&best_cost, &best_log_mean,
     &prev_seg_end, &prev_log_mean);
  mean_vec[0] = exp(best_log_mean);
  end_vec[0] = prev_seg_end;
  int out_i=1;
  while(0 <= prev_seg_end){
    // up_cost is actually either an up or down cost.
    up_cost = &cost_model_mat[prev_seg_offset + prev_seg_end];
    //printf("decoding out_i=%d prev_seg_end=%d prev_seg_offset=%d\n", out_i, prev_seg_end, prev_seg_offset);
    //up_cost->print();
    if(prev_log_mean != INFINITY){
      //equality constraint inactive
      best_log_mean = prev_log_mean;
    }
    up_cost->findMean
      (best_log_mean, &prev_seg_end, &prev_log_mean);
    mean_vec[out_i] = exp(best_log_mean);
    end_vec[out_i] = prev_seg_end;
    // change prev_seg_offset and out_i for next iteration.
    if(prev_seg_offset==0){
      //up_cost is actually up
      prev_seg_offset = data_count;
    }else{
      //up_cost is actually down
      prev_seg_offset = 0;
    }
    out_i++;
  }//for(data_i
}

// void PeakSegFPOPall
// (int *data_vec, double *weight_vec, int data_count,
//  double penalty,
//  // the following matrices are for output.
//  // cost_mat and intervals_mat store the optimal cost and number of intervals
//  // at each time point, for the up and down cost models.
//  // end_vec and mean_vec store the best model up to and including the
//  // last data point. 
//  double *cost_mat, //data_count x 2.
//  int *end_vec, //data_count
//  double *mean_vec,//data_count
//  int *intervals_mat,//data_count x 2
//  int *label_vec){//data_count
//   double min_log_mean=log(data_vec[0]), max_log_mean=log(data_vec[0]);
//   for(int data_i=1; data_i<data_count; data_i++){
//     double log_data = log(data_vec[data_i]);
//     if(log_data < min_log_mean){
//       min_log_mean = log_data;
//     }
//     if(max_log_mean < log_data){
//       max_log_mean = log_data;
//     }
//   }
//   std::vector<PiecewisePoissonLossLog> cost_model_mat(data_count * 2);
//   PiecewisePoissonLossLog *up_cost, *down_cost, *up_cost_prev, *down_cost_prev;
//   PiecewisePoissonLossLog min_prev_cost;
//   int verbose=0;
//   for(int data_i=0; data_i<data_count; data_i++){
//     //printf("computing cost data_i=%d\n", data_i);
//     up_cost = &cost_model_mat[data_i];
//     down_cost = &cost_model_mat[data_i + data_count]; 
//     if(data_i==0){
//       // initialization C_1(m)=gamma_1(m)
//       up_cost->piece_list.emplace_back
// 	(weight_vec[0], -data_vec[0]*weight_vec[0], 0.0,
// 	 min_log_mean, max_log_mean, -1, false);
//       down_cost->piece_list.emplace_back
// 	(weight_vec[0], -data_vec[0]*weight_vec[0], 0.0,
// 	 min_log_mean, max_log_mean, -1, false);
//     }else{
//       // if data_i is up, it could have come from down_cost_prev.
//       min_prev_cost.set_to_min_less_of(down_cost_prev, verbose);
//       min_prev_cost.set_prev_seg_end(data_i-1);
//       min_prev_cost.add(0.0, 0.0, penalty);
//       up_cost->set_to_min_env_of(&min_prev_cost, up_cost_prev, verbose);
//       up_cost->add
// 	(weight_vec[data_i],
// 	 -data_vec[data_i]*weight_vec[data_i],
// 	 0.0);
//       // if data_i is down, it could have come from up_cost_prev.
//       min_prev_cost.set_to_min_more_of(up_cost_prev, verbose);
//       min_prev_cost.set_prev_seg_end(data_i-1);
//       min_prev_cost.add(0.0, 0.0, penalty);
//       down_cost->set_to_min_env_of(&min_prev_cost, down_cost_prev, verbose);
//       down_cost->add
// 	(weight_vec[data_i],
// 	 -data_vec[data_i]*weight_vec[data_i],
// 	 0.0);
//     }
//     up_cost_prev = up_cost;
//     down_cost_prev = down_cost;
//   }
//   // Decoding the cost_model_vec, and writing to the output matrices.
//   double best_cost, best_log_mean, candidate_cost, candidate_log_mean;
//   bool equality_constraint_active, candidate_constraint;
//   int prev_seg_end=data_count, candidate_end;
//   for(int i=0; i<data_count; i++){
//     mean_vec[i] = INFINITY;
//     end_vec[i] = -2;
//     label_vec[i] = -1;
//   }
//   for(int i=0; i<2*data_count; i++){
//     up_cost = &cost_model_mat[i];
//     intervals_mat[i] = up_cost->piece_list.size();
//     up_cost->Minimize
//       (cost_mat+i, &best_log_mean,
//        &prev_seg_end, &equality_constraint_active);
//   }
//   int prev_seg_offset;
//   up_cost = &cost_model_mat[data_count-1];
//   up_cost->Minimize
//     (&candidate_cost, &candidate_log_mean,
//      &candidate_end, &candidate_constraint);
//   //printf("final up cost=%f at log_mean=%f\n", candidate_cost, candidate_log_mean);
//   //up_cost->print();
//   down_cost = &cost_model_mat[data_count*2-1];
//   down_cost->Minimize
//     (&best_cost, &best_log_mean,
//      &prev_seg_end, &equality_constraint_active);
//   //printf("final down cost=%f at log_mean=%f\n", best_cost, best_log_mean);
//   //down_cost->print();
//   if(candidate_cost < best_cost){
//     //more likely to end up, so previous segment was down.
//     best_cost = candidate_cost;
//     best_log_mean = candidate_log_mean;
//     prev_seg_end = candidate_end;
//     equality_constraint_active = candidate_constraint;
//     prev_seg_offset = data_count;
//     label_vec[0] = 1;
//   }else{
//     prev_seg_offset = 0;
//     label_vec[0] = 0;
//   }
//   mean_vec[0] = exp(best_log_mean);
//   end_vec[0] = prev_seg_end;
//   int out_i=1;
//   while(0 <= prev_seg_end){
//     //printf("decoding out_i=%d prev_seg_end=%d prev_seg_offset=%d\n", out_i, prev_seg_end, prev_seg_offset);
//     up_cost = &cost_model_mat[prev_seg_offset + prev_seg_end];
//     //up_cost->print();
//     if(equality_constraint_active){
//       up_cost->findMean
// 	(best_log_mean, &prev_seg_end, &equality_constraint_active);
//     }else{
//       up_cost->Minimize
// 	(&best_cost, &best_log_mean,
// 	 &prev_seg_end, &equality_constraint_active);
//     }
//     if(prev_seg_offset==0){
//       prev_seg_offset = data_count;
//       label_vec[out_i] = 1;
//     }else{
//       prev_seg_offset = 0;
//       label_vec[out_i] = 0;
//     }
//     mean_vec[out_i] = exp(best_log_mean);
//     end_vec[out_i] = prev_seg_end;
//     out_i++;
//   }//for(data_i
// }
