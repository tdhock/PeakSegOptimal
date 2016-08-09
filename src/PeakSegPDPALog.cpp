/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <vector>
#include <stdio.h>
#include "funPieceListLog.h"
#include <math.h>

#define IFPRINT(arg) if(data_i==370 && total_changes==-308) (arg)

void PeakSegPDPALog
(int *data_vec, double *weight_vec, int data_count,
 int maxSegments,
 // the following matrices are for output:
 double *cost_mat, //data_count x maxSegments.
 int *end_mat,//maxSegments x maxSegments(model up to last data point).
 double *mean_mat,
 int *intervals_mat){//maxSegments x maxSegments(model up to last data point).
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
  std::vector<PiecewisePoissonLossLog> cost_model_vec(data_count * maxSegments);
  double data_weight_cumsum = 0.0;
  double weight_cumsum = 0.0;
  std::vector<double> weight_cumsum_vec(data_count);
  for(int data_i=0; data_i<data_count; data_i++){
    weight_cumsum += weight_vec[data_i];
    weight_cumsum_vec[data_i] = weight_cumsum;
    data_weight_cumsum += data_vec[data_i]*weight_vec[data_i];
    PiecewisePoissonLossLog *cost_model = &cost_model_vec[data_i];
    cost_model->piece_list.emplace_back
      (1.0, -data_weight_cumsum/weight_cumsum, 0.0, min_log_mean, max_log_mean, -1, false);
  }

  // DP: compute functional model of best cost in S segments up to
  // data point N.
  PiecewisePoissonLossLog *prev_cost_model, *new_cost_model;
  PiecewisePoissonLossLog min_prev_cost, cost_model;
  for(int total_changes=1; total_changes<maxSegments; total_changes++){
    for(int data_i=total_changes; data_i<data_count; data_i++){
      int prev_i = data_i-1;
      prev_cost_model = &cost_model_vec[prev_i + (total_changes-1)*data_count];
      IFPRINT(printf("DP changes=%d data_i=%d\n", total_changes, data_i));
      IFPRINT(printf("=prev cost model\n"));
      IFPRINT(prev_cost_model->print());
      int verbose = 0, status;
      IFPRINT(verbose=1);
      if(total_changes % 2){
	min_prev_cost.set_to_min_less_of(prev_cost_model, verbose);
      }else{
	min_prev_cost.set_to_min_more_of(prev_cost_model, verbose);
      }
      status = min_prev_cost.check_min_of(prev_cost_model, prev_cost_model);
      if(status){
	printf("BAD MIN LESS/MORE CHECK status=%d changes=%d data_i=%d\n",
	       status, total_changes, data_i);
	printf("prev cost\n");
	prev_cost_model->print();
	printf("min less/more(prev cost)\n");
	min_prev_cost.print();
	throw status;
      }
      min_prev_cost.set_prev_seg_end(prev_i);
      new_cost_model = &cost_model_vec[data_i + total_changes*data_count];
      if(data_i==total_changes){//first cost model, only one candidate.
	IFPRINT(printf("=new cost model = min prev cost\n"));
	IFPRINT(min_prev_cost.print());
	*new_cost_model = min_prev_cost;
      }else{
	IFPRINT(printf("=min prev cost\n"));
        IFPRINT(min_prev_cost.print());
	IFPRINT(printf("=cost model\n"));
	IFPRINT(cost_model.print());
	new_cost_model->set_to_min_env_of
	  (&min_prev_cost, &cost_model, verbose);
	status = new_cost_model->check_min_of(&min_prev_cost, &cost_model);
	if(status){
	  printf("DP changes=%d data_i=%d BAD CHECK status=%d\n", total_changes, data_i, status);
	  printf("=prev cost model\n");
	  prev_cost_model->print();
	  printf("=min prev cost\n");
	  min_prev_cost.print();
	  printf("=cost model\n");
	  cost_model.print();
	  printf("=new cost model\n");
	  new_cost_model->print();
	  throw status;
	}
      }
      IFPRINT(printf("=new cost model\n"));
      IFPRINT(new_cost_model->print());
      new_cost_model->multiply(weight_cumsum_vec[prev_i]);
      new_cost_model->add
	(weight_vec[data_i],
	 -data_vec[data_i]*weight_vec[data_i],
	 0.0);
      new_cost_model->multiply(1/weight_cumsum_vec[data_i]);
      IFPRINT(new_cost_model->print());
      cost_model = *new_cost_model;
    }
  }

  double best_cost, best_log_mean, prev_log_mean;
  double *best_mean_vec;
  int *prev_seg_vec;
  int prev_seg_end;
  
  // Decoding the cost_model_vec, and writing to the output matrices.
  for(int i=0; i<maxSegments*maxSegments; i++){
    mean_mat[i] = INFINITY;
    end_mat[i] = -1;
  }
  for(int i=0; i<maxSegments*data_count; i++){
    intervals_mat[i] = -1;
    cost_mat[i] = INFINITY;
  }
  for(int total_changes=0; total_changes<maxSegments;total_changes++){
    for(int data_i=total_changes; data_i<data_count; data_i++){
      IFPRINT(printf("decoding changes=%d data_i=%d\n", total_changes, data_i));
      PiecewisePoissonLossLog *cost_model =
	&cost_model_vec[data_i + total_changes*data_count];
      IFPRINT(cost_model->print());
      cost_model->Minimize
	(&best_cost, &best_log_mean,
	 &prev_seg_end, &prev_log_mean);
      IFPRINT(printf("cost=%f log_mean=%f prev_end=%d prev_log_mean=%f\n", best_cost, best_log_mean, prev_seg_end, prev_log_mean));
      // for the models up to any data point, we store the best cost
      // and the total number of intervals.
      cost_mat[data_i + total_changes*data_count] = best_cost;
      intervals_mat[data_i + total_changes*data_count] =
	cost_model->piece_list.size();
      if(data_i == data_count-1){
	// for the models up to the last data point, we also store the
	// segment means and the change-points. we have already
	// computed the values for the last segment so store them.
	best_mean_vec = mean_mat + total_changes*maxSegments;
	prev_seg_vec = end_mat + total_changes*maxSegments;
	best_mean_vec[total_changes] = exp(best_log_mean);
	prev_seg_vec[total_changes] = prev_seg_end;
	for(int seg_i=total_changes-1; 0 <= seg_i; seg_i--){
	  //printf("seg_i=%d prev_seg_end=%d\n", seg_i, prev_seg_end);
	  cost_model = &cost_model_vec[prev_seg_end + seg_i*data_count];
	  if(prev_log_mean != INFINITY){
	    //equality constraint inactive
	    best_log_mean = prev_log_mean;
	  }
	  cost_model->findMean
	    (best_log_mean, &prev_seg_end, &prev_log_mean);
	  best_mean_vec[seg_i] = exp(best_log_mean);
	  prev_seg_vec[seg_i] = prev_seg_end;
	}//for(seg_i
      }//if(data_i
    }//for(data_i
  }//for(total_changes
}
