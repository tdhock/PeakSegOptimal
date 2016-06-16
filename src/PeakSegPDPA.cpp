/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <vector>
#include "funPieceList.h"

void PeakSegPDPA
(double *data_vec, double *weight_vec, int data_count,
 int maxSegments,
 // the following matrices are for output, data_count x maxSegments.
 double *cost_mat,
 int *end_mat,
 double *mean_mat){
  double min_mean=data_vec[0], max_mean=data_vec[0];
  for(int data_i=1; data_i<data_count; data_i++){
    if(data_vec[data_i] < min_mean){
      min_mean = data_vec[data_i];
    }
    if(max_mean < data_vec[data_i]){
      max_mean = data_vec[data_i];
    }
  }
  std::vector<PiecewisePoissonLoss> cost_model_vec(data_count * maxSegments);
  double Log_cumsum = 0;
  double Linear_cumsum = 0;
  PoissonLossPiece piece;
  PiecewisePoissonLoss cost_model;
  for(int data_i=0; data_i<data_count; data_i++){
    Linear_cumsum += weight_vec[data_i];
    Log_cumsum += -data_vec[data_i]*weight_vec[data_i];
    cost_model.piece_list.clear();
    cost_model.piece_list.emplace_back
      (Linear_cumsum, Log_cumsum, 0.0, min_mean, max_mean, -1, false);
    cost_model_vec[data_i] = cost_model;
  }

  // DP.
  PiecewisePoissonLoss *prev_cost_model;
  PiecewisePoissonLoss candidate;
  for(int total_changes=1; total_changes<maxSegments; total_changes++){
    int prev_i = total_changes-1;
    prev_cost_model = &cost_model_vec[prev_i + prev_i*data_count];
    if(total_changes % 2){
      prev_cost_model->min_less(&candidate);
    }else{
      prev_cost_model->min_more(&candidate);
    }
    candidate.add
      (weight_vec[total_changes],
       -data_vec[total_changes]*weight_vec[total_changes],
       0.0);
    candidate.set_prev_seg_end(prev_i);
    cost_model_vec[total_changes + total_changes*data_count] = candidate;
  }

  double best_cost, best_mean;
  double *best_mean_vec;
  int *prev_seg_vec;
  bool equality_constraint_active;
  int prev_seg_end;
  // for(int total_changes=0; total_changes<maxSegments;total_changes++){
  //   for(int data_i=0; data_i<data_count; data_i++){
  //     cost_mat[data_i + total_changes*data_count] = INFINITY;
  //   }
  // }
  
  // Decoding the cost_model_vec, and writing to the output matrices.
  maxSegments=1;//TODO change.
  for(int total_changes=0; total_changes<maxSegments;total_changes++){
    for(int data_i=0; data_i<data_count; data_i++){
      cost_model = cost_model_vec[data_i + total_changes*data_count];
      cost_model.Minimize(&best_cost, &best_mean,
			  &prev_seg_end, &equality_constraint_active);
      // for the models up to any data point, we store the best cost.
      cost_mat[data_i + total_changes*data_count] = best_cost;
      if(data_i == data_count-1){
	// for the models up to the last data point, we also store the
	// segment means and the change-points. we have already
	// computed the values for the last segment so store them.
	best_mean_vec = mean_mat + total_changes*maxSegments;
	prev_seg_vec = end_mat + total_changes*maxSegments;
	best_mean_vec[total_changes] = best_mean;
	prev_seg_vec[total_changes] = prev_seg_end;
	for(int seg_i=total_changes-1; 0 <= seg_i; seg_i--){
	  cost_model = cost_model_vec[prev_seg_end + seg_i*data_count];
	  if(equality_constraint_active){
	    cost_model.findMean
	      (best_mean, &prev_seg_end, &equality_constraint_active);
	  }else{
	    cost_model.Minimize(&best_cost, &best_mean,
				&prev_seg_end, &equality_constraint_active);
	  }
	  best_mean_vec[seg_i] = best_mean;
	  prev_seg_vec[seg_i] = prev_seg_end;
	}//for(seg_i
      }//if(data_i
    }//for(data_i
  }//for(total_changes
}
