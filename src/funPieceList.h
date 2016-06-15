/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <list>
#include <math.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_errno.h>


// sprintf("%.90f",-1/exp(1))
// "-0.367879441171442334024277442949824035167694091796875000"
#define POISSON_THRESH -0.367879441171442334024277442949824035167694091796875

class PoissonLossPiece {
 public:
  double Linear;
  double Log;
  double Constant;
  double min_mean;
  double max_mean;
  int data_i;
  bool equality_constraint_active;
  PoissonLossPiece(){
  }
  double getDiscriminant(double add_constant){
    return Linear/Log*exp((add_constant-Constant)/Log);
  }
  double discriminant2mean_principal(double discriminant){
    gsl_sf_result result;
    int status = gsl_sf_lambert_W0_e(discriminant, &result);
    if(status != GSL_SUCCESS){
      throw status;
    }
    return Log/Linear*result.val;
  }
  PoissonLossPiece
    (double li, double lo, double co, double m, double M, int i, bool a){
    Linear = li;
    Log = lo;
    Constant = co;
    min_mean = m;
    max_mean = M;
    data_i = i;
    equality_constraint_active = a;
  }
  double getMinMean(){
    return - Log / Linear;
  }
  double PoissonLoss(double mean){
    double loss = Linear*mean + Constant;
    if(Log!=0){
      loss += Log * log(mean);
    }
    return loss;
  }
};

typedef std::list<PoissonLossPiece> PoissonLossPieceList;

class PiecewisePoissonLoss {
 public:
  PoissonLossPieceList piece_list;
  void add(double Linear, double Log, double Constant){
    PoissonLossPieceList::iterator it;
    for(it=piece_list.begin(); it != piece_list.end(); it++){
      it->Linear += Linear;
      it->Log += Log;
      it->Constant += Constant;
    }
  }
  void set_prev_seg_end(int prev_seg_end){
    PoissonLossPieceList::iterator it;
    for(it=piece_list.begin(); it != piece_list.end(); it++){
      it->data_i = prev_seg_end;
    }
  }
  void findMean(double mean, int *seg_end, bool *equality_constraint_active){
    PoissonLossPieceList::iterator it;
    for(it=piece_list.begin(); it != piece_list.end(); it++){
      if(it->min_mean <= mean && mean <= it->max_mean){
	*seg_end = it->data_i;
	*equality_constraint_active = it->equality_constraint_active;
      }
    }
  }
  void min_less(PiecewisePoissonLoss *output){
    double prev_min_cost = INFINITY, prev_min_mean;
    int prev_data_i;
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
	  PoissonLossPiece new_piece
	    (it->Linear, it->Log, it->Constant, prev_min_mean, mu,
	     it->data_i, true); // equality constraint active on convex piece.
	  output->piece_list.push_back(new_piece);
	}else{
	  // Minimum after this interval, so this function is
	  // decreasing on this entire interval, and so we can just
	  // store it as is.
	  PoissonLossPiece new_piece
	    (it->Linear, it->Log, it->Constant, prev_min_mean, it->max_mean,
	     it->data_i, true); // equality constraint active on convex piece.
	  output->piece_list.push_back(new_piece);
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
	}else{
	  double discriminant = it->getDiscriminant(prev_min_mean);
	  if(POISSON_THRESH < discriminant){
	    // since the Log coef is negative, the principal branch
	    // results in the smaller of the two mean values.
	    double mu = it->discriminant2mean_principal(discriminant);
	  }
	}
      }
    }
  }
  void min_more(PiecewisePoissonLoss *output){
  }
  void Minimize(double *best_cost,
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
};
