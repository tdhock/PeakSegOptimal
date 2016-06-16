/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "funPieceList.h"
#include <list>
#include <math.h>

#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_errno.h>

// sprintf("%.90f",-1/exp(1))
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


void PiecewisePoissonLoss::min_less(PiecewisePoissonLoss *output){
  double prev_min_cost = INFINITY, prev_min_mean;
  int prev_data_i = -2;
  output->piece_list.clear();
  PoissonLossPieceList::iterator it = piece_list.begin();
  while(it != piece_list.end()){
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
	output->piece_list.emplace_back
	  (it->Linear, it->Log, it->Constant, prev_min_mean, mu,
	   it->data_i, true); // equality constraint active on convex piece.
	prev_min_mean = mu;
	prev_min_cost = it->PoissonLoss(mu);
	prev_data_i = it->data_i;
      }else{
	// Minimum after this interval, so this function is
	// decreasing on this entire interval, and so we can just
	// store it as is.
	output->piece_list.emplace_back
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
	    output->piece_list.emplace_back
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
    output->piece_list.emplace_back
      (0.0, 0.0, prev_min_cost, prev_min_mean, it->max_mean, prev_data_i,
       false);//equality constraint inactive on constant piece.
  }
}

void PiecewisePoissonLoss::min_more(PiecewisePoissonLoss *output){
  double prev_min_cost = INFINITY, prev_max_mean;
  int prev_data_i = -2;
  output->piece_list.clear();
  PoissonLossPieceList::iterator it = piece_list.end();
  while(it != piece_list.begin()){
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
	output->piece_list.emplace_front
	  (it->Linear, it->Log, it->Constant, prev_max_mean, mu,
	   it->data_i, true); // equality constraint active on convex piece.
	prev_max_mean = mu;
	prev_min_cost = it->PoissonLoss(mu);
	prev_data_i = it->data_i;
      }else{
	// Minimum before this interval, so this function is
	// increasing on this entire interval, and so we can just
	// store it as is.
	output->piece_list.emplace_front
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
	output->piece_list.emplace_front
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
    output->piece_list.emplace_back
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
