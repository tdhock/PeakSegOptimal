/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "funPieceList.h"
#include <list>
#include <math.h>
#include <stdio.h>

#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_errno.h>

// s//printf("%.90f",-1/exp(1))
// "-0.367879441171442334024277442949824035167694091796875000"
#define POISSON_THRESH -0.367879441171442334024277442949824035167694091796875

PoissonLossPiece::PoissonLossPiece
(double li, double lo, double co, double m, double M, int i, bool a){
  Linear = li;
  Log = lo;
  Constant = co;
  min_mean = m;
  max_mean = M;
  data_i = i;
  equality_constraint_active = a;
}

double PoissonLossPiece::getDiscriminant(double add_constant){
  return Linear/Log*exp((add_constant-Constant)/Log);
}

double PoissonLossPiece::discriminant2mean_principal(double discriminant){
  gsl_sf_result result;
  int status = gsl_sf_lambert_W0_e(discriminant, &result);
  if(status != GSL_SUCCESS){
    throw status;
  }
  return Log/Linear*result.val;
}

double PoissonLossPiece::discriminant2mean_secondary(double discriminant){
  gsl_sf_result result;
  int status = gsl_sf_lambert_Wm1_e(discriminant, &result);
  if(status != GSL_SUCCESS){
    throw status;
  }
  return Log/Linear*result.val;
}

double PoissonLossPiece::discriminant2mean_larger(double discriminant){
  if(0 < Log/Linear){
    return discriminant2mean_principal(discriminant);
  }else{
    return discriminant2mean_secondary(discriminant);
  }
}

double PoissonLossPiece::discriminant2mean_smaller(double discriminant){
  if(0 < Log/Linear){
    return discriminant2mean_secondary(discriminant);
  }else{
    return discriminant2mean_principal(discriminant);
  }
}

double PoissonLossPiece::getMinMean(){
  return - Log / Linear;
}

double PoissonLossPiece::PoissonLoss(double mean){
  double loss = Linear*mean + Constant;
  if(Log!=0){
    loss += Log * log(mean);
  }
  return loss;
}

double PoissonLossPiece::PoissonDeriv(double mean){
  return Linear + Log/mean;
}

void PiecewisePoissonLoss::set_to_min_less_of(PiecewisePoissonLoss *input){
  int prev_data_i = -2;
  piece_list.clear();
  PoissonLossPieceList::iterator it = input->piece_list.begin();
  double prev_min_cost = INFINITY, prev_min_mean = it->min_mean;
  while(it != input->piece_list.end()){
    if(prev_min_cost == INFINITY){
      // Look for min achieved in this interval.
      double mu = it->getMinMean();
      if(mu <= it->min_mean){
	/* The minimum is achieved before this interval, so this
	   function is always increasing in this interval. We don't
	   need to store it, but we do need to keep track of the
	   minimum cost, which occurs at the min mean value in this
	   interval. */
	prev_min_cost = it->PoissonLoss(it->min_mean);
	prev_data_i = it->data_i;
      }else if(mu < it->max_mean){
	// Minimum in this interval, so add a convex piece up to the
	// min, and keep track of the min cost to create a constant
	// piece later.
	piece_list.emplace_back
	  (it->Linear, it->Log, it->Constant, prev_min_mean, mu,
	   it->data_i, true); // equality constraint active on convex piece.
	prev_min_mean = mu;
	prev_min_cost = it->PoissonLoss(mu);
	prev_data_i = it->data_i;
      }else{
	// Minimum after this interval, so this function is
	// decreasing on this entire interval, and so we can just
	// store it as is.
	piece_list.emplace_back
	  (it->Linear, it->Log, it->Constant, prev_min_mean, it->max_mean,
	   it->data_i, true); // equality constraint active on convex piece.
	prev_min_mean = it->max_mean;
      }
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
	      (0.0, 0.0, prev_min_cost, prev_min_mean, mu, prev_data_i,
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
      (0.0, 0.0, prev_min_cost, prev_min_mean, it->max_mean, prev_data_i,
       false);//equality constraint inactive on constant piece.
  }
}

void PiecewisePoissonLoss::set_to_min_more_of(PiecewisePoissonLoss *input){
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
      double mu = it->getMinMean();
      if(it->max_mean <= mu){
	/* The minimum is achieved after this interval, so this
	   function is always decreasing in this interval. We don't
	   need to store it. */
	prev_min_cost = it->PoissonLoss(it->max_mean);
	prev_data_i = it->data_i;
      }else if(it->min_mean < mu){
	// Minimum in this interval, so add a convex piece up to the
	// min, and keep track of the min cost to create a constant
	// piece later.
	piece_list.emplace_front
	  (it->Linear, it->Log, it->Constant, mu, prev_max_mean, 
	   it->data_i, true); // equality constraint active on convex piece.
	prev_max_mean = mu;
	prev_min_cost = it->PoissonLoss(mu);
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
	  (0.0, 0.0, prev_min_cost,
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
      (0.0, 0.0, prev_min_cost,
       it->min_mean, prev_max_mean,
       prev_data_i,
       false);//equality constraint inactive on constant piece.
  }
}

void PiecewisePoissonLoss::add(double Linear, double Log, double Constant){
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
    }
  }
}

void PiecewisePoissonLoss::print(){
  PoissonLossPieceList::iterator it;
  printf("%d pieces\n", piece_list.size());
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    printf("Piece: [%f,%f] data_i=%d Linear=%f Log=%f Constant=%f\n",
	   it->min_mean, it->max_mean, it->data_i,
	   it->Linear, it->Log, it->Constant);
  }
}
  

void PiecewisePoissonLoss::Minimize(double *best_cost,
	      double *best_mean,
	      int *data_i,
	      bool *equality_constraint_active){
  double candidate_cost, candidate_mean;
  PoissonLossPieceList::iterator it;
  *best_cost = INFINITY;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    candidate_mean = it->getMinMean();
    if(it->min_mean <= candidate_mean && candidate_mean <= it->max_mean){
      candidate_cost = it->PoissonLoss(candidate_mean);
      if(candidate_cost < *best_cost){
	*best_cost = candidate_cost;
	*best_mean = candidate_mean;
	*data_i = it->data_i;
	*equality_constraint_active = it->equality_constraint_active;
      }
    }
  }
}

void PiecewisePoissonLoss::set_to_min_env_of
(PiecewisePoissonLoss *fun1, PiecewisePoissonLoss *fun2){
  PoissonLossPieceList::iterator
    it1 = fun1->piece_list.begin(),
    it2 = fun2->piece_list.begin();
  piece_list.clear();
  while(it1 != fun1->piece_list.end() &&
	it2 != fun2->piece_list.end()){
    push_min_pieces(fun1, fun2, it1, it2);
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
  

void PiecewisePoissonLoss::push_min_pieces
(PiecewisePoissonLoss *fun1,
 PiecewisePoissonLoss *fun2,
 PoissonLossPieceList::iterator it1,
 PoissonLossPieceList::iterator it2){
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
  if(sameFuns(it1, it2)){
    // The functions are exactly equal over the entire interval so we
    // can push either of them.
    push_piece(it1, last_min_mean, first_max_mean);
    //printf("exactly equal over entire interval\n");
    return;
  }
  PoissonLossPiece diff_piece
    (it1->Linear - it2->Linear,
     it1->Log - it2->Log,
     it1->Constant - it2->Constant,
     last_min_mean, first_max_mean,
     -5, false);
  if(diff_piece.Log == 0){
    if(diff_piece.Linear == 0){
      // They are offset by a Constant.
      if(diff_piece.Constant < 0){
	push_piece(it1, last_min_mean, first_max_mean);
      }else{
	push_piece(it2, last_min_mean, first_max_mean);
      }
      //printf("offset by a constant\n");
      return;
    }
    if(diff_piece.Constant == 0){
      // The only difference is the Linear coef.
      if(diff_piece.Linear < 0){
	push_piece(it1, last_min_mean, first_max_mean);
      }else{
	push_piece(it2, last_min_mean, first_max_mean);
      }
      //printf("only diff is linear coef\n");
      return;
    }
    double mean_at_equal_cost = -diff_piece.Constant/diff_piece.Linear;
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
      //printf("Log zero with one root in interval\n");
      return;
    }
    // the root is outside the interval, so one is completely above
    // the other over this entire interval.
    if(mean_at_equal_cost < last_min_mean){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    //printf("Log zero with no roots in interval\n");
    return;
  }//if(diff->Log == 0
  double cost_diff_left = diff_piece.PoissonLoss(last_min_mean);
  double cost_diff_right = diff_piece.PoissonLoss(first_max_mean);
  double discriminant = diff_piece.getDiscriminant(0.0);
  bool two_roots = POISSON_THRESH < discriminant;
  bool fun1_min_on_right;
  if(!same_at_right){
    fun1_min_on_right = cost_diff_right < 0;
  }else{
    // They are equal on the right limit, so use the first and second
    // derivatives to see which is minimal just before the right
    // limit. Do we need to check if they intersect before the right
    // limit? Only if the sign of the slope of the more curved
    // function is positive. And in that case we just need to check
    // the smaller root.
    double deriv1_right = it1->PoissonDeriv(first_max_mean);
    double deriv2_right = it2->PoissonDeriv(first_max_mean);
    bool maybe_cross = (it2->Log < it1->Log && 0 < deriv2_right) ||
				    (it1->Log < it2->Log && 0 < deriv1_right);
    if(two_roots && maybe_cross){
      // there could be a crossing point to the left.
      double mean_at_equal_cost =
	diff_piece.discriminant2mean_smaller(discriminant);
      if(last_min_mean < mean_at_equal_cost &&
	 mean_at_equal_cost && first_max_mean){
	//the cross point is in the interval.
	if(cost_diff_left < 0){
	  push_piece(it1, last_min_mean, mean_at_equal_cost);
	  push_piece(it2, mean_at_equal_cost, first_max_mean);
	}else{
	  push_piece(it2, last_min_mean, mean_at_equal_cost);
	  push_piece(it1, mean_at_equal_cost, first_max_mean);
	}
	//printf("equal on the right with one crossing in interval\n");
	return;
      }
    }//if(two_roots && maybe cross
    bool fun1_min_before_right;
    if((deriv1_right < 0 && deriv2_right < 0) ||
       (0 < deriv1_right && 0 < deriv2_right)){
      // the derivatives have the same sign on the right, so are
      // certainly not equal on the left, and we can just compare
      // their cost values on the left to determine which is larger.
      fun1_min_before_right = cost_diff_left < 0;
    }else{
      // the derivatives have the same sign, so we can compare them to
      // determine which function is larger.
      fun1_min_before_right = deriv2_right < deriv1_right;
    }
    if(fun1_min_before_right){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    //printf("equal on the right with no crossing in interval\n");
    return;
  }
  bool fun1_min_on_left;
  if(!same_at_left){
    fun1_min_on_left = cost_diff_left < 0;
  }else{
    // equal on the left.
    double deriv1_left = it1->PoissonDeriv(last_min_mean);
    double deriv2_left = it2->PoissonDeriv(last_min_mean);
    bool maybe_cross = (it2->Log < it1->Log && deriv2_left < 0) ||
      (it1->Log < it2->Log && deriv1_left < 0);
    if(two_roots && maybe_cross){
      // There could be a crossing point to the right.
      double mean_at_equal_cost =
	diff_piece.discriminant2mean_larger(discriminant);
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
	//printf("equal on the left with crossing in interval\n");
	return;
      }
    }//if(there may be crossing
    bool fun1_min_after_left;
    if((deriv1_left < 0 && deriv2_left < 0) ||
       (0 < deriv1_left && 0 < deriv2_left)){
      fun1_min_after_left = cost_diff_right < 0;
    }else{
      fun1_min_after_left = deriv1_left < deriv2_left;
    }
    if(fun1_min_after_left){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    //printf("equal on the left with no crossing in interval\n");
    return;
  }
  // The only remaining case is that the curves are equal neither on
  // the left nor on the right of the interval. However they may be
  // equal inside the interval, so let's check for that.
  double first_mean = INFINITY, second_mean = INFINITY;
  if(two_roots){
    double larger_mean = diff_piece.discriminant2mean_larger(discriminant);
    bool larger_inside =
      last_min_mean < larger_mean && larger_mean < first_max_mean;
    double smaller_mean = diff_piece.discriminant2mean_smaller(discriminant);
    bool smaller_inside =
      last_min_mean < smaller_mean && smaller_mean < first_max_mean;
    if(larger_inside){
      if(smaller_inside){
	// both are in the interval.
	first_mean = smaller_mean;
	second_mean = larger_mean;
	//printf("%f and %f in [%f,%f]\n",
	//smaller_mean, larger_mean,
	//last_min_mean, first_max_mean);
      }else{
	// smaller mean is not in the interval, but the larger is.
	first_mean = larger_mean;
	//printf("%f in [%f,%f]\n",
	//first_mean,
	//last_min_mean, first_max_mean);
      }
    }else{
      // larger mean is not in the interval
      if(smaller_inside){
	// smaller mean is in the interval, but not the larger.
	first_mean = smaller_mean;
	//printf("%f in [%f,%f]\n",
	//first_mean,
	//last_min_mean, first_max_mean);
      }
    }
  }
  if(second_mean != INFINITY){
    // two crossing points.
    if(fun1_min_on_left){
      push_piece(it1, last_min_mean, first_mean);
      push_piece(it2, first_mean, second_mean);
      push_piece(it1, second_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_mean);
      push_piece(it1, first_mean, second_mean);
      push_piece(it2, second_mean, first_max_mean);
    }
    //printf("not equal on the sides, 2 crossing points\n");
  }else if(first_mean != INFINITY){
    // one crossing point.
    if(fun1_min_on_left){
      push_piece(it1, last_min_mean, first_mean);
      push_piece(it2, first_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_mean);
      push_piece(it1, first_mean, first_max_mean);
    }
    //printf("not equal on the sides, 1 crossing point\n");
  }else{
    // zero crossing points.
    if(fun1_min_on_left){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    //printf("not equal on the sides, zero crossing points\n");
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
