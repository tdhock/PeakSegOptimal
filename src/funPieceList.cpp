/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "funPieceList.h"
#include <list>
#include <math.h>
#include <stdio.h>

#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_errno.h>

// s//printf("%a",-1/exp(1))
// "-0.367879441171442334024277442949824035167694091796875000"
#define POISSON_THRESH -0.367879441171442334024277442949824035167694091796875

PoissonLossPiece::PoissonLossPiece
(int li, int lo, double co, double m, double M, int i, bool a){
  Linear = li;
  Log = lo;
  Constant = co;
  min_mean = m;
  max_mean = M;
  data_i = i;
  equality_constraint_active = a;
}

double PoissonLossPiece::getDiscriminant(double add_constant){
  return (double)Linear/(double)Log*exp((add_constant-Constant)/(double)Log);
}

double PoissonLossPiece::discriminant2mean_principal(double discriminant){
  gsl_sf_result result;
  int status = gsl_sf_lambert_W0_e(discriminant, &result);
  if(status != GSL_SUCCESS){
    throw status;
  }
  return (double)Log/(double)Linear*result.val;
}

double PoissonLossPiece::discriminant2mean_secondary(double discriminant){
  gsl_sf_result result;
  int status = gsl_sf_lambert_Wm1_e(discriminant, &result);
  if(status != GSL_SUCCESS){
    throw status;
  }
  return (double)Log/(double)Linear*result.val;
}

double PoissonLossPiece::discriminant2mean_larger(double discriminant){
  if(0 < (double)Log/(double)Linear){
    return discriminant2mean_principal(discriminant);
  }else{
    return discriminant2mean_secondary(discriminant);
  }
}

double PoissonLossPiece::discriminant2mean_smaller(double discriminant){
  if(0 < (double)Log/(double)Linear){
    return discriminant2mean_secondary(discriminant);
  }else{
    return discriminant2mean_principal(discriminant);
  }
}

double PoissonLossPiece::getMinMean(){
  return - (double)Log / (double)Linear;
}

double PoissonLossPiece::PoissonLoss(double mean, int verbose){
  double loss_without_log_term = (double)Linear*mean + Constant;
  //verbose=0;
  //if(verbose)printf("loss before adding log term %a\n", loss_without_log_term);
  if(Log==0){
    return loss_without_log_term;
  }
  //printf(""); //does not work here.
  double log_mean_only = log(mean);
  //printf(""); //does work here.
  //if(log_mean_only < 0)return INFINITY;
  //printf("log_mean_only=%a\n", log_mean_only);
  double log_coef_only = (double)Log;
  //if(verbose)printf("log_coef_only=%a\n", log_coef_only);
  double product = log_mean_only * log_coef_only;
  //if(verbose)printf("product=%a\n", product);
  return loss_without_log_term + product;
}

double PoissonLossPiece::PoissonDeriv(double mean){
  return Linear + (double)Log/mean;
}

void PiecewisePoissonLoss::set_to_min_less_of
(PiecewisePoissonLoss *input, int verbose){
  int prev_data_i = -2;
  piece_list.clear();
  PoissonLossPieceList::iterator it = input->piece_list.begin();
  double prev_min_cost = INFINITY, prev_min_mean = it->min_mean;
  while(it != input->piece_list.end()){
    if(prev_min_cost == INFINITY){
      // Look for min achieved in this interval.
      if(it->Log==0){
	// degenerate linear function. since the Linear coef is never
	// negative, we know that this function must be increasing or
	// numerically constant on this interval.
	double cost_left = it->PoissonLoss(it->min_mean, verbose);
	double cost_right = it->PoissonLoss(it->max_mean, verbose);
	//printf("DEGENERATE LINEAR FUNCTION IN MIN LESS\n");
	if(cost_left==cost_right){
	  //printf("NUMERICALLY EQUAL\n");
	  // store this numerically constant interval.
	  piece_list.emplace_back
	    (it->Linear, it->Log, it->Constant, prev_min_mean, it->max_mean,
	     it->data_i, true); // equality constraint active on convex piece.
	  prev_min_mean = it->max_mean;
	}else{
	  //printf("NUMERICALLY INCREASING\n");
	  // don't store this interval, but store its min cost as a
	  // constant.
	  prev_min_cost = cost_left;
	  prev_data_i = it->data_i;
	}
      }else{
	double mu = it->getMinMean();
	//printf("getMinMean=%a\n", mu);
	if(mu <= it->min_mean){
	  /* The minimum is achieved on the left or before this
	     interval, so this function is always increasing in this
	     interval. We don't need to store it, but we do need to keep
	     track of the minimum cost, which occurs at the min mean
	     value in this interval. */
	  prev_min_cost = it->PoissonLoss(it->min_mean, verbose);
	  prev_data_i = it->data_i;
	}else if(mu < it->max_mean){
	  // Minimum in this interval, so add a convex piece up to the
	  // min, and keep track of the min cost to create a constant
	  // piece later.
	  piece_list.emplace_back
	    (it->Linear, it->Log, it->Constant, prev_min_mean, mu,
	     it->data_i, true); // equality constraint active on convex piece.
	  prev_min_mean = mu;
	  prev_min_cost = it->PoissonLoss(mu, verbose);
	  if(verbose)printf("prev_min_cost=%a\n", prev_min_cost);
	  prev_data_i = it->data_i;
	}else{
	  // Minimum after this interval, so this function is
	  // decreasing on this entire interval, and so we can just
	  // store it as is.
	  piece_list.emplace_back
	    (it->Linear, it->Log, it->Constant, prev_min_mean, it->max_mean,
	     it->data_i, true); // equality constraint active on convex piece.
	  prev_min_mean = it->max_mean;
	}//if(non-degenerate mu in interval
      }//if(degenerate linear cost.
    }else{//prev_min_cost is finite
      // Look for a function with prev_min_cost in its interval.
      if(it->Log==0){
	//degenerate Linear case
	if(it->Linear < 0){
	  //decreasing linear function.
	  throw(500);//this should never happen.
	}else{
	  //increasing linear function, so will not intersect the
	  //constant below.
	}
      }else{// Log is not zero.
	// Theoretically there can be zero, one, or two intersection
	// points between the constant function prev_min_mean and a
	// non-degenerate Poisson loss function piece. To check what
	// case we are in, we compare the discriminant to
	// POISSON_THRESH=-1/e (the point x at which LambertW(x)=-1
	// for both branches).
	double discriminant = it->getDiscriminant(prev_min_cost);
	if(POISSON_THRESH < discriminant){
	  // There are two mean values where the Poisson loss piece
	  // intersects the constant function prev_min_mean, but we
	  // are only concerned with the first mean value (the
	  // lesser of the two). Since the Log coef of all Poisson
	  // loss pieces is negative, the principal branch results
	  // in the smaller of the two mean values.
	  double mu = it->discriminant2mean_principal(discriminant);
	  if(it->min_mean < mu && mu < it->max_mean){
	    // The smaller intersection point occurs within the
	    // interval, so the constant interval ends here, and we
	    // can store it immediately.
	    piece_list.emplace_back
	      (0, 0, prev_min_cost, prev_min_mean, mu, prev_data_i,
	       false);// equality constraint inactive on constant piece.
	    prev_min_cost = INFINITY;
	    prev_min_mean = mu;
	    it--;
	  }//if(mu occurs in interval
	}//if(there are two roots
      }//if(Log is zero
    }//if(prev_min_cost is finite
    it++;
  }//while(it
  if(prev_data_i != -2){
    // ending on a constant piece -- we never have a convex piece at
    // the end, because the end is the maximum observed data.
    it--;
    piece_list.emplace_back
      (0, 0, prev_min_cost, prev_min_mean, it->max_mean, prev_data_i,
       false);//equality constraint inactive on constant piece.
  }
}

void PiecewisePoissonLoss::set_to_min_more_of(PiecewisePoissonLoss *input){
  int verbose=false;
  int prev_data_i = -2;
  piece_list.clear();
  PoissonLossPieceList::iterator it = input->piece_list.end();
  it--;
  double prev_min_cost = INFINITY, prev_max_mean = it->max_mean;
  it++;
  while(it != input->piece_list.begin()){
    it--;
    if(prev_min_cost == INFINITY){
      // Look for min achieved in this interval.
      if(it->Log==0){
	//degenerate Linear function. since the Linear coef is never
	//negative, we know that this function must be increasing or
	//numerically constant on this interval. In both cases we
	//should just store this interval.
	piece_list.emplace_front
	  (it->Linear, it->Log, it->Constant, it->min_mean, prev_max_mean,
	   it->data_i, true); // equality constraint active on convex piece.
	prev_max_mean = it->min_mean;
      }else{
	double mu = it->getMinMean();
	if(it->max_mean <= mu){
	  /* The minimum is achieved after this interval, so this
	     function is always decreasing in this interval. We don't
	     need to store it. */
	  prev_min_cost = it->PoissonLoss(it->max_mean, verbose);
	  prev_data_i = it->data_i;
	}else if(it->min_mean < mu){
	  // Minimum in this interval, so add a convex piece up to the
	  // min, and keep track of the min cost to create a constant
	  // piece later.
	  piece_list.emplace_front
	    (it->Linear, it->Log, it->Constant, mu, prev_max_mean, 
	     it->data_i, true); // equality constraint active on convex piece.
	  prev_max_mean = mu;
	  prev_min_cost = it->PoissonLoss(mu, verbose);
	  prev_data_i = it->data_i;
	}else{
	  // Minimum before this interval, so this function is
	  // increasing on this entire interval, and so we can just
	  // store it as is.
	  piece_list.emplace_front
	    (it->Linear, it->Log, it->Constant, it->min_mean, prev_max_mean,
	     it->data_i, true); // equality constraint active on convex piece.
	  prev_max_mean = it->min_mean;
	}
      }//if(degenerate linear)else
    }else{//prev_min_cost is finite
      // Look for a function with prev_min_cost in its interval.
      double mu = INFINITY;
      if(it->Log==0){
	//degenerate Linear case, there is one intersection point.
	mu = (prev_min_cost - it->Constant)/it->Linear;
      }else{// Log is not zero.
	// Theoretically there can be zero, one, or two intersection
	// points between the constant function prev_min_mean and a
	// non-degenerate Poisson loss function piece. To check what
	// case we are in, we compare the discriminant to
	// POISSON_THRESH=-1/e (the point x at which LambertW(x)=-1
	// for both branches).
	double discriminant = it->getDiscriminant(prev_min_cost);
	if(POISSON_THRESH < discriminant){
	  // There are two mean values where the Poisson loss piece
	  // intersects the constant function prev_min_mean, but we
	  // are only concerned with the second mean value (the
	  // greater of the two). Since the Log coef of all Poisson
	  // loss pieces is negative, the non-principal branch results in
	  // the larger of the two mean values.
	  mu = it->discriminant2mean_secondary(discriminant);
	}//if(there are two roots
      }//if(Log is zero
      if(it->min_mean < mu && mu < it->max_mean){
	// The smaller intersection point occurs within the
	// interval, so the constant interval ends here, and we
	// can store it immediately.
	piece_list.emplace_front
	  (0, 0, prev_min_cost,
	   mu, prev_max_mean,
	   prev_data_i,
	   false);// equality constraint inactive on constant piece.
	prev_min_cost = INFINITY;
	prev_max_mean = mu;
	it++;
      }//if(mu occurs in interval
    }//if(prev_min_cost is finite
  }//while(it
  if(prev_data_i != -2){
    // ending on a constant piece -- we never have a convex piece at
    // the start, because the end is the min observed data.
    piece_list.emplace_front
      (0, 0, prev_min_cost,
       it->min_mean, prev_max_mean,
       prev_data_i,
       false);//equality constraint inactive on constant piece.
  }
}

void PiecewisePoissonLoss::add(int Linear, int Log, double Constant){
  PoissonLossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->Linear += Linear;
    it->Log += Log;
    it->Constant += Constant;
  }
}

void PiecewisePoissonLoss::set_prev_seg_end(int prev_seg_end){
  PoissonLossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->data_i = prev_seg_end;
  }
}

void PiecewisePoissonLoss::findMean(double mean, int *seg_end, bool *equality_constraint_active){
  PoissonLossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    if(it->min_mean <= mean && mean <= it->max_mean){
      *seg_end = it->data_i;
      *equality_constraint_active = it->equality_constraint_active;
      return;
    }
  }
}

double PiecewisePoissonLoss::findCost(double mean){
  PoissonLossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    if(it->min_mean <= mean && mean <= it->max_mean){
      int verbose = 0;
      return it->PoissonLoss(mean, verbose);
    }
  }
}

void PiecewisePoissonLoss::print(){
  PoissonLossPieceList::iterator it;
  printf("%10s %10s %10s %10s %10s %10s\n",
	 "Linear", "Log", "Constant", "min_mean", "max_mean", "data_i");
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    printf("%10d %10d %20f %60.55f %60.55f %d\n",
	   it->Linear, it->Log, it->Constant,
	   it->min_mean, it->max_mean, it->data_i);
  }
}
  

void PiecewisePoissonLoss::Minimize(double *best_cost,
	      double *best_mean,
	      int *data_i,
	      bool *equality_constraint_active){
  double candidate_cost, candidate_mean;
  int verbose=false;
  PoissonLossPieceList::iterator it;
  *best_cost = INFINITY;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    candidate_mean = it->getMinMean();
    if(it->min_mean <= candidate_mean && candidate_mean <= it->max_mean){
      candidate_cost = it->PoissonLoss(candidate_mean, verbose);
      if(candidate_cost < *best_cost){
	*best_cost = candidate_cost;
	*best_mean = candidate_mean;
	*data_i = it->data_i;
	*equality_constraint_active = it->equality_constraint_active;
      }
    }
  }
}

// check that this function is the minimum on all pieces.
int PiecewisePoissonLoss::check_min_of
(PiecewisePoissonLoss *prev, PiecewisePoissonLoss *model){
  PoissonLossPieceList::iterator it;
  int verbose = 0;
  for(it = piece_list.begin(); it != piece_list.end(); it++){
    if(it != piece_list.begin()){
      PoissonLossPieceList::iterator pit = it;
      pit--;
      if(pit->max_mean != it->min_mean){
	printf("prev->max_mean != it->min_mean min\n");
	return 3;
      }
    }
    if(it->max_mean <= it->min_mean){
      printf("max_mean<=min_mean=%15.10f min\n", it->min_mean);
      return 2;
    }
    double mid_mean = (it->min_mean + it->max_mean)/2;
    double cost_min = it->PoissonLoss(mid_mean, verbose);
    double cost_prev = prev->findCost(mid_mean);
    if(cost_prev+1e-6 < cost_min){
      printf("prev(%.55f)=%f %a\n", mid_mean, cost_prev, cost_prev);
      prev->print();
      printf("min(%.55f)=%f %a\n", mid_mean, cost_min, cost_min);
      print();
      return 1;
    }
    double cost_model = model->findCost(mid_mean);
    if(cost_model+1e-6 < cost_min){
      printf("model(%.55f)=%f %a\n", mid_mean, cost_model, cost_model);
      model->print();
      printf("min(%.55f)=%f %a\n", mid_mean, cost_min, cost_min);
      print();
      return 1;
    }
  }
  for(it = prev->piece_list.begin(); it != prev->piece_list.end(); it++){
    if(it != prev->piece_list.begin()){
      PoissonLossPieceList::iterator pit = it;
      pit--;
      if(pit->max_mean != it->min_mean){
	printf("prev->max_mean != it->min_mean prev\n");
	return 3;
      }
    }
    if(it->max_mean <= it->min_mean){
      printf("max_mean<=min_mean=%15.10f prev\n", it->min_mean);
      return 2;
    }
    double mid_mean = (it->min_mean + it->max_mean)/2;
    double cost_prev = it->PoissonLoss(mid_mean, verbose);
    double cost_min = findCost(mid_mean);
    if(cost_prev+1e-6 < cost_min){
      printf("prev(%.55f)=%f %a\n", mid_mean, cost_prev, cost_prev);
      prev->print();
      printf("min(%.55f)=%f %a\n", mid_mean, cost_min, cost_min);
      print();
      return 1;
    }
  }
  for(it = model->piece_list.begin(); it != model->piece_list.end(); it++){
    if(it != model->piece_list.begin()){
      PoissonLossPieceList::iterator pit = it;
      pit--;
      if(pit->max_mean != it->min_mean){
	printf("prev->max_mean != it->min_mean model\n");
	return 3;
      }
    }
    if(it->max_mean <= it->min_mean){
      printf("max_mean<=min_mean=%15.10f model\n", it->min_mean);
      return 2;
    }
    double mid_mean = (it->min_mean + it->max_mean)/2;
    double cost_model = it->PoissonLoss(mid_mean, verbose);
    double cost_min = findCost(mid_mean);
    if(cost_model+1e-6 < cost_min){
      printf("model(%.55f)=%f %a\n", mid_mean, cost_model, cost_model);
      model->print();
      printf("min(%.55f)=%f %a\n", mid_mean, cost_min, cost_min);
      print();
      return 1;
    }
  }
  return 0;
}

void PiecewisePoissonLoss::set_to_min_env_of
(PiecewisePoissonLoss *fun1, PiecewisePoissonLoss *fun2, int verbose){
  PoissonLossPieceList::iterator
    it1 = fun1->piece_list.begin(),
    it2 = fun2->piece_list.begin();
  piece_list.clear();
  while(it1 != fun1->piece_list.end() &&
	it2 != fun2->piece_list.end()){
    push_min_pieces(fun1, fun2, it1, it2, verbose);
    if(verbose){
      print();
      printf("------\n");
    }
    double last_max_mean = piece_list.back().max_mean;
    if(it1->max_mean == last_max_mean){
      it1++;
    }
    if(it2->max_mean == last_max_mean){
      it2++;
    }
  }
}

bool sameFuns
(PoissonLossPieceList::iterator it1,
 PoissonLossPieceList::iterator it2){
  return it1->Linear == it2->Linear &&
    it1->Log == it2->Log &&
    it1->Constant == it2->Constant;
}

bool sameSigns(double x, double y){
  return
    (x  < 0 && y < 0) ||
    (0  < x && 0 < y) ||
    (x == 0 && y == 0);
}

void PiecewisePoissonLoss::push_min_pieces
(PiecewisePoissonLoss *fun1,
 PiecewisePoissonLoss *fun2,
 PoissonLossPieceList::iterator it1,
 PoissonLossPieceList::iterator it2,
 int verbose){
  bool same_at_left;
  double last_min_mean;
  PoissonLossPieceList::iterator prev2 = it2;
  prev2--;
  PoissonLossPieceList::iterator prev1 = it1;
  prev1--;
  if(it1->min_mean < it2->min_mean){
    //it1 function piece starts to the left of it2.
    same_at_left = sameFuns(prev2, it1);
    last_min_mean = it2->min_mean;
  }else{
    //it1 function piece DOES NOT start to the left of it2.
    last_min_mean = it1->min_mean;
    if(it2->min_mean < it1->min_mean){
      //it2 function piece starts to the left of it1.
      same_at_left = sameFuns(prev1, it2);
    }else{
      //it1 and it2 start at the same min_mean value.
      if(it1==fun1->piece_list.begin() &&
	 it2==fun2->piece_list.begin()){
	same_at_left = false;
      }else{
	same_at_left = sameFuns(prev1, prev2);
      }
    }
  }
  PoissonLossPieceList::iterator next2 = it2;
  next2++;
  PoissonLossPieceList::iterator next1 = it1;
  next1++;
  bool same_at_right;
  double first_max_mean;
  if(it1->max_mean < it2->max_mean){
    // it2 function piece continues to the right of it1.
    same_at_right = sameFuns(next1, it2);
    first_max_mean = it1->max_mean;
  }else{
    first_max_mean = it2->max_mean;
    if(it2->max_mean < it1->max_mean){
      // it2 function piece ends before it1.
      same_at_right = sameFuns(it1, next2);
    }else{
      if(next1==fun1->piece_list.end() &&
	 next2==fun2->piece_list.end()){
	same_at_right = false;
      }else{
	same_at_right = sameFuns(next1, next2);
      }
    }
  }
  if(last_min_mean == first_max_mean){
    // we should probably never get here, but if we do, no need to
    // store this interval.
    if(verbose){
      printf("prev\n");
      fun1->print();
      printf("model\n");
      fun2->print();
      printf("interval size 0!-----------------\n");
    }
    return;
  }
  if(sameFuns(it1, it2)){
    // The functions are exactly equal over the entire interval so we
    // can push either of them.
    push_piece(it1, last_min_mean, first_max_mean);
    if(verbose)printf("exactly equal over entire interval\n");
    return;
  }
  PoissonLossPiece diff_piece
    (it1->Linear - it2->Linear,
     it1->Log - it2->Log,
     it1->Constant - it2->Constant,
     last_min_mean, first_max_mean,
     -5, false);
  // printf("it1->Constant=%a\nit2->Constant=%a\n",
  // 	 it1->Constant, it2->Constant);
  double mid_mean = (first_max_mean + last_min_mean)/2;
  double cost_diff_mid = diff_piece.PoissonLoss(mid_mean, false);
  if(diff_piece.Log == 0){
    if(diff_piece.Linear == 0){
      // They are offset by a Constant.
      if(diff_piece.Constant < 0){
	push_piece(it1, last_min_mean, first_max_mean);
      }else{
	push_piece(it2, last_min_mean, first_max_mean);
      }
      if(verbose)printf("offset by a constant\n");
      return;
    }
    if(diff_piece.Constant == 0){
      // The only difference is the Linear coef.
      if(diff_piece.Linear < 0){
	push_piece(it1, last_min_mean, first_max_mean);
      }else{
	push_piece(it2, last_min_mean, first_max_mean);
      }
      if(verbose)printf("only diff is linear coef\n");
      return;
    }
    double mean_at_equal_cost = -diff_piece.Constant/(double)diff_piece.Linear;
    if(last_min_mean < mean_at_equal_cost &&
       mean_at_equal_cost < first_max_mean){
      // the root is in the interval, so we need to add two intervals.
      if(0 < diff_piece.Linear){
	push_piece(it1, last_min_mean, mean_at_equal_cost);
	push_piece(it2, mean_at_equal_cost, first_max_mean);
      }else{
	push_piece(it2, last_min_mean, mean_at_equal_cost);
	push_piece(it1, mean_at_equal_cost, first_max_mean);
      }
      if(verbose)printf("Log zero with one root in interval\n");
      return;
    }
    // the root is outside the interval, so one is completely above
    // the other over this entire interval.
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    if(verbose)printf("Log zero with no roots in interval\n");
    return;
  }//if(diff->Log == 0
  double cost_diff_left = diff_piece.PoissonLoss(last_min_mean, false);
  double cost_diff_right = diff_piece.PoissonLoss(first_max_mean, false);
  double discriminant = diff_piece.getDiscriminant(0.0);
  bool two_roots = false;
  double larger_mean, smaller_mean;
  if(POISSON_THRESH < discriminant){
    larger_mean = diff_piece.discriminant2mean_larger(discriminant);
    smaller_mean = diff_piece.discriminant2mean_smaller(discriminant);
    if(smaller_mean < larger_mean){
      two_roots = true;
    }
  }
  if(verbose)printf("discriminant=%a two_roots=%d Linear=%d Log=%d\nConstant=%a\n", discriminant, two_roots, diff_piece.Linear, diff_piece.Log, diff_piece.Constant);
  if(same_at_right){
    // they are equal on the right, but we don't know if there is
    // another crossing point somewhere to the left.
    if(two_roots){
      // there could be a crossing point to the left.
      double mean_at_equal_cost = smaller_mean;
      if(verbose)printf("smaller_mean=%a\n", mean_at_equal_cost);
      if(last_min_mean < mean_at_equal_cost &&
      	 mean_at_equal_cost < first_max_mean){
	//the cross point is in the interval.
	if(cost_diff_left < 0){
	  push_piece(it1, last_min_mean, mean_at_equal_cost);
	  push_piece(it2, mean_at_equal_cost, first_max_mean);
	}else{
	  push_piece(it2, last_min_mean, mean_at_equal_cost);
	  push_piece(it1, mean_at_equal_cost, first_max_mean);
	}
	if(verbose)printf("equal on the right with one crossing in interval\n");
	return;
      }
    }//if(two_roots
    // Test the cost at the midpoint, since the cost may be equal on
    // both the left and the right.
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    if(verbose)printf("equal on the right with no crossing in interval\n");
    return;
  }
  if(same_at_left){
    // equal on the left.
    if(two_roots){
      // There could be a crossing point to the right.
      double mean_at_equal_cost = larger_mean;
      if(verbose)printf("larger_mean=%a\n", mean_at_equal_cost);
      if(last_min_mean < mean_at_equal_cost &&
      	 mean_at_equal_cost < first_max_mean){
	// the crossing point is in this interval.
	if(cost_diff_right < 0){
	  push_piece(it2, last_min_mean, mean_at_equal_cost);
	  push_piece(it1, mean_at_equal_cost, first_max_mean);
	}else{
	  push_piece(it1, last_min_mean, mean_at_equal_cost);
	  push_piece(it2, mean_at_equal_cost, first_max_mean);
	}
	if(verbose)printf("equal on the left with crossing in interval\n");
	return;
      }
    }//if(there may be crossing
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    if(verbose)printf("equal on the left with no crossing in interval\n");
    return;
  }
  // The only remaining case is that the curves are equal neither on
  // the left nor on the right of the interval. However they may be
  // equal inside the interval, so let's check for that.
  double first_mean = INFINITY, second_mean = INFINITY;
  if(two_roots){
    bool larger_inside =
      last_min_mean < larger_mean && larger_mean < first_max_mean;
    if(verbose)printf("smaller_mean=%f %a\nlarger_mean=%f %a\n",
		      smaller_mean, smaller_mean,
		      larger_mean, larger_mean);
    bool smaller_inside =
      last_min_mean < smaller_mean && smaller_mean < first_max_mean;
    if(larger_inside){
      if(smaller_inside){
	// both are in the interval.
	first_mean = smaller_mean;
	second_mean = larger_mean;
	if(verbose){
	  printf("%f and %f in [%f,%f]\n",
		 smaller_mean, larger_mean,
		 last_min_mean, first_max_mean);
	}
      }else{
	// smaller mean is not in the interval, but the larger is.
	first_mean = larger_mean;
	if(verbose){
	  printf("%f in [%f,%f]\n",
		 first_mean,
		 last_min_mean, first_max_mean);
	}
      }
    }else{
      // larger mean is not in the interval
      if(smaller_inside){
	// smaller mean is in the interval, but not the larger.
	first_mean = smaller_mean;
	if(verbose){
	  printf("%f in [%f,%f]\n",
		 first_mean,
		 last_min_mean, first_max_mean);
	}
      }
    }
  }//if(two_roots
  if(second_mean != INFINITY){
    // two crossing points.
    double before_mean = (last_min_mean + first_mean)/2;
    double cost_diff_before = diff_piece.PoissonLoss(before_mean, false);
    if(cost_diff_before < 0){
      push_piece(it1, last_min_mean, first_mean);
      push_piece(it2, first_mean, second_mean);
      push_piece(it1, second_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_mean);
      push_piece(it1, first_mean, second_mean);
      push_piece(it2, second_mean, first_max_mean);
    }
    if(verbose)printf("not equal on the sides, 2 crossing points\n");
  }else if(first_mean != INFINITY){
    // "one" crossing point. actually sometimes we have last_min_mean
    // < first_mean < first_max_mean but cost_diff_before and
    // cost_diff_after have the same sign! In that case we need to
    // just push one piece.
    double before_mean = (last_min_mean + first_mean)/2;
    double cost_diff_before = diff_piece.PoissonLoss(before_mean, false);
    if(verbose)printf("cost_diff_before(%.55f)=%f\n", before_mean, cost_diff_before);
    double after_mean = (first_max_mean + first_mean)/2;
    double cost_diff_after = diff_piece.PoissonLoss(after_mean, false);
    if(verbose)printf("cost_diff_after(%.55f)=%f\n", after_mean, cost_diff_after);
    if(cost_diff_before < 0){
      if(cost_diff_after < 0){
	// f1-f2<0 meaning f1<f2 on the entire interval, so just push it1.
	push_piece(it1, last_min_mean, first_max_mean);
      }else{
	push_piece(it1, last_min_mean, first_mean);
	push_piece(it2, first_mean, first_max_mean);
      }
    }else{//f1(before)-f2(before)>=0 meaning f1(before)>=f2(before)
      if(cost_diff_after < 0){
	//f1(after)-f2(after)<0 meaning f1(after)<f2(after)
	push_piece(it2, last_min_mean, first_mean);
	push_piece(it1, first_mean, first_max_mean);
      }else{
	//f1(after)-f2(after)>=0 meaning f1(after)>=f2(after)
	push_piece(it2, last_min_mean, first_max_mean);
      }
    }
    if(verbose)printf("not equal on the sides, 1 crossing point\n");
  }else{
    // "zero" crossing points. actually there may be a crossing point
    // in the interval that is numerically so close as to be identical
    // with last_min_mean or first_max_mean.
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    if(verbose)printf("not equal on the sides, zero crossing points\n");
  }
}

void PiecewisePoissonLoss::push_piece
(PoissonLossPieceList::iterator it, double min_mean, double max_mean){
  PoissonLossPieceList::iterator last=piece_list.end();
  if(piece_list.size() && sameFuns(--last, it)){
    //it is the same function as last, so just make last extend
    //further to the right.
    last->max_mean = max_mean;
  }else{
    //it is a different function than last, so push it on to the end
    //of the list.
    piece_list.emplace_back
      (it->Linear, it->Log, it->Constant,
       min_mean, max_mean,
       it->data_i, it->equality_constraint_active);
  }
}
