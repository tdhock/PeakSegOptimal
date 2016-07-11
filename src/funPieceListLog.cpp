/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "funPieceListLog.h"
#include <list>
#include <math.h>
#include <stdio.h>

#define NEWTON_EPSILON 1e-8
#define NEWTON_STEPS 100

#define ABS(x) ((x)<0 ? -(x) : (x))

PoissonLossPieceLog::PoissonLossPieceLog
(int li, int lo, double co, double m, double M, int i, bool a){
  Linear = li;
  Log = lo;
  Constant = co;
  min_log_mean = m;
  max_log_mean = M;
  data_i = i;
  equality_constraint_active = a;
}

bool PoissonLossPieceLog::has_two_roots(double equals){
  // are there two solutions to the equation Linear*e^x + Log*x +
  // Constant = equals ?
  if(Log == 0){//degenerate linear function.
    // f(mean) = Linear*mean + Constant - equals
    throw "problem";
    return false;
  }
  double optimal_log_mean = argmin(); //min or max!
  double optimal_cost = getCost(optimal_log_mean);
  double optimal_mean = - (double)Log / (double)Linear; //min or max!
  double optimal_cost2 = PoissonLoss(optimal_mean);
  // does g(x) = Linear*e^x + Log*x + Constant = equals? compare the
  // cost at optimum to equals.
  // g'(x)= Linear*e^x + Log,
  // g''(x)= Linear*e^x.
  if(0 < Linear){//convex.
    return optimal_cost < equals && optimal_cost2 < equals;
  }
  //concave.
  return equals < optimal_cost && equals < optimal_cost2;
}

double PoissonLossPieceLog::PoissonLoss(double mean){
  double loss_without_log_term = (double)Linear*mean + Constant;
  if(Log==0){
    return loss_without_log_term;
  }
  double log_mean_only = log(mean);
  double log_coef_only = (double)Log;
  double product = log_mean_only * log_coef_only;
  return loss_without_log_term + product;
}

double PoissonLossPieceLog::PoissonDeriv(double mean){
  return Linear + (double)Log/mean;
}

// This function performs the root finding in the positive (mean)
// space, but needs to return a value in the log(mean) space.
double PoissonLossPieceLog::get_larger_root(double equals){
  double optimal_mean = - (double)Log / (double)Linear; //min or max!
  double optimal_cost = PoissonLoss(optimal_mean);
  // Approximate the solution by the line through
  // (optimal_mean,optimal_cost) with the asymptotic slope. As m tends
  // to Inf, f'(m)=Linear+Log/m tends to Linear.
  //double candidate_root = optimal_mean + (equals-optimal_cost)/(double)Linear;
  double candidate_root = optimal_mean + 1;
  // find the larger root of f(m) = Linear*m + Log*log(m) + Constant -
  // equals = 0.
  double candidate_cost, possibly_outside, deriv;
  double closest_positive_cost = INFINITY, closest_positive_mean;
  double closest_negative_cost = -INFINITY, closest_negative_mean;
  if(optimal_cost < 0){
    closest_negative_cost = optimal_cost;
    closest_negative_mean = optimal_mean;
  }else{
    closest_positive_cost = optimal_cost;
    closest_positive_mean = optimal_mean;
  }
  int step=0;
  do{
     candidate_cost = PoissonLoss(candidate_root) - equals;
     if(0 < candidate_cost && candidate_cost < closest_positive_cost){
       closest_positive_cost = candidate_cost;
       closest_positive_mean = candidate_root;
     }
     if(closest_negative_cost < candidate_cost && candidate_cost < 0){
       closest_negative_cost = candidate_cost;
       closest_negative_mean = candidate_root;
     }
     if(NEWTON_STEPS <= ++step){
       printf("larger root MAXSTEPS with equals=%e\n", equals);
       print();
       printf("step=%d mean=%e cost=%e\n", step, candidate_root, candidate_cost);
       return log((closest_positive_mean + closest_negative_mean)/2);
     }
     deriv = PoissonDeriv(candidate_root);
     possibly_outside = candidate_root - candidate_cost/deriv;
     if(possibly_outside < optimal_mean){
       //it overshot to the left of the optimum, so the root is
       //probably very close to the optimum, and we have probably
       //already explored very close to the zero.
       printf("larger root WRONG SIDE equals=%e\n", equals);
       print();
       printf("neg_cost=%e neg_mean=%e pos_cost=%e pos_mean=%e\n", closest_negative_cost, closest_negative_mean, closest_positive_cost, closest_positive_mean);
       if(closest_negative_cost==-INFINITY){
	 double optimal_log_mean = argmin(); //min or max!
	 double optimal_cost2 = getCost(optimal_log_mean);
	 printf("optimal_mean=%e=%e=exp(%e) optimal_cost=%e=%e=\n", optimal_mean, exp(optimal_log_mean), optimal_log_mean, optimal_cost, optimal_cost2);
	 throw 1;
       }
       return log((closest_positive_mean + closest_negative_mean)/2);
     }else{
       // it is to the left so that is fine.
       candidate_root = possibly_outside;
     }
  }while(NEWTON_EPSILON < ABS(candidate_cost));
  //printf("found root %e in %d steps!\n", candidate_root, step);
  return log(candidate_root);
}

double PoissonLossPieceLog::get_smaller_root(double equals){
  double optimal_log_mean = argmin(); //min or max!
  double optimal_cost = getCost(optimal_log_mean);
  // Approximate the solution by the line through
  // (optimal_mean,optimal_cost) with the asymptotic slope. As x tends
  // to -Inf, g'(x)=Linear*e^x+Log tends to Log.
  //double candidate_root = optimal_log_mean + (equals-optimal_cost)/(double)Log;
  double candidate_root = optimal_log_mean - 1;
  // find the smaller root of g(x) = Linear*e^x + Log*x + Constant -
  // equals = 0.
  double candidate_cost, possibly_outside, deriv;
  // as we search we will store bounds on the left and the right of
  // the zero point.
  double closest_positive_cost = INFINITY, closest_positive_log_mean;
  double closest_negative_cost = -INFINITY, closest_negative_log_mean;
  if(optimal_cost < 0){
    closest_negative_cost = optimal_cost;
    closest_negative_log_mean = optimal_log_mean;
  }else{
    closest_positive_cost = optimal_cost;
    closest_positive_log_mean = optimal_log_mean;
  }
  int step=0;
  do{
     candidate_cost = getCost(candidate_root) - equals;
     if(0 < candidate_cost && candidate_cost < closest_positive_cost){
       closest_positive_cost = candidate_cost;
       closest_positive_log_mean = candidate_root;
     }
     if(closest_negative_cost < candidate_cost && candidate_cost < 0){
       closest_negative_cost = candidate_cost;
       closest_negative_log_mean = candidate_root;
     }
     if(NEWTON_STEPS <= ++step){
       printf("smaller root MAXSTEPS equals=%e\n", equals);
       print();
       printf("step=%d log_mean=%e cost=%e\n", step, candidate_root, candidate_cost);
       return (closest_positive_log_mean + closest_negative_log_mean)/2;
     }
     deriv = getDeriv(candidate_root);
     possibly_outside = candidate_root - candidate_cost/deriv;
     if(possibly_outside < optimal_log_mean){
       // it's to the left of the optimum, no problem.
       candidate_root = possibly_outside;
     }else{
       // it's on the right of the optimum, so the root is probably
       //very close to the optimum, and we have probably already
       //explored very close to the zero.
       printf("smaller root WRONG SIDE equals=%e\n", equals);
       print();
       printf("neg_cost=%e neg_log_mean=%e pos_cost=%e pos_log_mean=%e\n", closest_negative_cost, closest_negative_log_mean, closest_positive_cost, closest_positive_log_mean);
       return (closest_positive_log_mean + closest_negative_log_mean)/2;
     }
  }while(NEWTON_EPSILON < ABS(candidate_cost));
  return candidate_root;
}

double PoissonLossPieceLog::argmin(){
  // g(x) = Linear*e^x + Log*x + Constant,
  // g'(x)= Linear*e^x + Log = 0 means
  // x = log(-Log/Linear).
  return log(- (double)Log / (double)Linear);
}

double PoissonLossPieceLog::getCost(double log_mean){
  // f(m) = Linear*m + Log*log(m) + Constant,
  // x = log(m),
  // g(x) = Linear*e^x + Log*x + Constant.
  double linear_term, log_term;
  if(log_mean == -INFINITY){
    linear_term = 0.0;
  }else{
    linear_term = (double)Linear*exp(log_mean);
  }
  if(Log==0){
    log_term = 0.0;
  }else{
    log_term = (double)Log*log_mean;
  }
  return linear_term + log_term + Constant;
}

double PoissonLossPieceLog::getDeriv(double log_mean){
  // g(x) = Linear*e^x + Log*x + Constant,
  // g'(x)= Linear*e^x + Log.
  double linear_term;
  if(log_mean == -INFINITY){
    linear_term = 0.0;
  }else{
    linear_term = (double)Linear*exp(log_mean);
  }
  return linear_term + (double)Log;
}

void PiecewisePoissonLossLog::set_to_min_less_of
(PiecewisePoissonLossLog *input, int verbose){
  int prev_data_i = -2;
  piece_list.clear();
  PoissonLossPieceListLog::iterator it = input->piece_list.begin();
  double prev_min_cost = INFINITY, prev_min_log_mean = it->min_log_mean;
  while(it != input->piece_list.end()){
    if(prev_min_cost == INFINITY){
      // Look for min achieved in this interval.
      if(it->Log==0){
	// degenerate linear function. since the Linear coef is never
	// negative, we know that this function must be increasing or
	// numerically constant on this interval.
	// g(x) = Linear*e^x + Constant,
	if(verbose)printf("DEGENERATE LINEAR FUNCTION IN MIN LESS\n");
	if(it->Linear==0){
	  if(verbose){
	    printf("Constant interval\n");
	    it->print();
	  }
	  // store this numerically constant interval.
	  piece_list.emplace_back
	    (it->Linear, it->Log, it->Constant, prev_min_log_mean, it->max_log_mean,
	     it->data_i, true); // equality constraint active on convex piece.
	  prev_min_log_mean = it->max_log_mean;
	}else{
	  if(verbose){
	    printf("Increasing interval\n");
	    it->print();
	  }
	  // don't store this interval, but store its min cost as a
	  // constant.
	  prev_min_cost = it->getCost(it->min_log_mean);
	  prev_data_i = it->data_i;
	}
      }else{
	double mu = it->argmin();
	if(verbose)printf("argmin=%a\n", mu);
	if(mu <= it->min_log_mean){
	  /* The minimum is achieved on the left or before this
	     interval, so this function is always increasing in this
	     interval. We don't need to store it, but we do need to keep
	     track of the minimum cost, which occurs at the min mean
	     value in this interval. */
	  prev_min_cost = it->getCost(it->min_log_mean);
	  prev_data_i = it->data_i;
	}else if(mu < it->max_log_mean){
	  // Minimum in this interval, so add a convex piece up to the
	  // min, and keep track of the min cost to create a constant
	  // piece later. NB it is possible that prev_min_log_mean==mu in
	  // which case we do not need to store the convex piece.
	  if(prev_min_log_mean < mu){
	    piece_list.emplace_back
	      (it->Linear, it->Log, it->Constant, prev_min_log_mean, mu,
	       it->data_i, true); // equality constraint active on convex piece.
	  }
	  prev_min_log_mean = mu;
	  prev_min_cost = it->getCost(mu);
	  if(verbose)printf("prev_min_cost=%a\n", prev_min_cost);
	  prev_data_i = it->data_i;
	}else{
	  // Minimum after this interval, so this function is
	  // decreasing on this entire interval, and so we can just
	  // store it as is.
	  piece_list.emplace_back
	    (it->Linear, it->Log, it->Constant, prev_min_log_mean, it->max_log_mean,
	     it->data_i, true); // equality constraint active on convex piece.
	  prev_min_log_mean = it->max_log_mean;
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
	// points between the constant function prev_min_log_mean and a
	// non-degenerate Poisson loss function piece. 
	if(it->has_two_roots(prev_min_cost)){
	  // There are two mean values where the Poisson loss piece
	  // intersects the constant function prev_min_log_mean, but we
	  // are only concerned with the first mean value (the
	  // lesser of the two).
	  double mu = it->get_smaller_root(prev_min_cost);
	  if(it->min_log_mean < mu && mu < it->max_log_mean){
	    // The smaller intersection point occurs within the
	    // interval, so the constant interval ends here, and we
	    // can store it immediately.
	    piece_list.emplace_back
	      (0, 0, prev_min_cost, prev_min_log_mean, mu, prev_data_i,
	       false);// equality constraint inactive on constant piece.
	    prev_min_cost = INFINITY;
	    prev_min_log_mean = mu;
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
      (0, 0, prev_min_cost, prev_min_log_mean, it->max_log_mean, prev_data_i,
       false);//equality constraint inactive on constant piece.
  }
}

void PiecewisePoissonLossLog::set_to_min_more_of(PiecewisePoissonLossLog *input){
  int verbose=false;
  int prev_data_i = -2;
  piece_list.clear();
  PoissonLossPieceListLog::iterator it = input->piece_list.end();
  it--;
  double prev_min_cost = INFINITY, prev_max_log_mean = it->max_log_mean;
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
	  (it->Linear, it->Log, it->Constant, it->min_log_mean, prev_max_log_mean,
	   it->data_i, true); // equality constraint active on convex piece.
	prev_max_log_mean = it->min_log_mean;
      }else{
	double mu = it->argmin();
	if(it->max_log_mean <= mu){
	  /* The minimum is achieved after this interval, so this
	     function is always decreasing in this interval. We don't
	     need to store it. */
	  prev_min_cost = it->getCost(it->max_log_mean);
	  prev_data_i = it->data_i;
	}else if(it->min_log_mean < mu){
	  // Minimum in this interval, so add a convex piece up to the
	  // min, and keep track of the min cost to create a constant
	  // piece later. NB it is possible that mu==prev_max_log_mean, in
	  // which case we do not need to save the convex piece.
	  if(mu < prev_max_log_mean){ 
	    piece_list.emplace_front
	      (it->Linear, it->Log, it->Constant, mu, prev_max_log_mean, 
	       it->data_i, true); // equality constraint active on convex piece.
	  }
	  prev_max_log_mean = mu;
	  prev_min_cost = it->getCost(mu);
	  prev_data_i = it->data_i;
	}else{
	  // Minimum before this interval, so this function is
	  // increasing on this entire interval, and so we can just
	  // store it as is.
	  piece_list.emplace_front
	    (it->Linear, it->Log, it->Constant, it->min_log_mean, prev_max_log_mean,
	     it->data_i, true); // equality constraint active on convex piece.
	  prev_max_log_mean = it->min_log_mean;
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
	// points between the constant function prev_min_log_mean and a
	// non-degenerate Poisson loss function piece. 
	if(it->has_two_roots(prev_min_cost)){
	  // There are two mean values where the Poisson loss piece
	  // intersects the constant function prev_min_log_mean, but we
	  // are only concerned with the second mean value (the
	  // greater of the two). 
	  mu = it->get_larger_root(prev_min_cost);
	}//if(there are two roots
      }//if(Log is zero
      if(it->min_log_mean < mu && mu < it->max_log_mean){
	// The smaller intersection point occurs within the
	// interval, so the constant interval ends here, and we
	// can store it immediately.
	piece_list.emplace_front
	  (0, 0, prev_min_cost,
	   mu, prev_max_log_mean,
	   prev_data_i,
	   false);// equality constraint inactive on constant piece.
	prev_min_cost = INFINITY;
	prev_max_log_mean = mu;
	it++;
      }//if(mu occurs in interval
    }//if(prev_min_cost is finite
  }//while(it
  if(prev_data_i != -2){
    // ending on a constant piece -- we never have a convex piece at
    // the start, because the end is the min observed data.
    piece_list.emplace_front
      (0, 0, prev_min_cost,
       it->min_log_mean, prev_max_log_mean,
       prev_data_i,
       false);//equality constraint inactive on constant piece.
  }
}

void PiecewisePoissonLossLog::add(int Linear, int Log, double Constant){
  PoissonLossPieceListLog::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->Linear += Linear;
    it->Log += Log;
    it->Constant += Constant;
  }
}

void PiecewisePoissonLossLog::set_prev_seg_end(int prev_seg_end){
  PoissonLossPieceListLog::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->data_i = prev_seg_end;
  }
}

void PiecewisePoissonLossLog::findMean
(double log_mean, int *seg_end, bool *equality_constraint_active){
  PoissonLossPieceListLog::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    if(it->min_log_mean <= log_mean && log_mean <= it->max_log_mean){
      *seg_end = it->data_i;
      *equality_constraint_active = it->equality_constraint_active;
      return;
    }
  }
}

double PiecewisePoissonLossLog::findCost(double mean){
  PoissonLossPieceListLog::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    if(it->min_log_mean <= mean && mean <= it->max_log_mean){
      int verbose = 0;
      return it->getCost(mean);
    }
  }
}

void PiecewisePoissonLossLog::print(){
  PoissonLossPieceListLog::iterator it;
  printf("%10s %10s %10s %10s %10s %10s\n",
	 "Linear", "Log", "Constant", "min_log_mean", "max_log_mean", "data_i");
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->print();
  }
}

void PoissonLossPieceLog::print(){
  printf("%10d %10d %20f %60.55f %60.55f %d\n",
	 Linear, Log, Constant,
	 min_log_mean, max_log_mean, data_i);
}

void PiecewisePoissonLossLog::Minimize(double *best_cost,
	      double *best_log_mean,
	      int *data_i,
	      bool *equality_constraint_active){
  double candidate_cost, candidate_log_mean;
  int verbose=false;
  PoissonLossPieceListLog::iterator it;
  *best_cost = INFINITY;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    candidate_log_mean = it->argmin();
    if(it->min_log_mean <= candidate_log_mean &&
       candidate_log_mean <= it->max_log_mean){
      candidate_cost = it->getCost(candidate_log_mean);
      if(candidate_cost < *best_cost){
	*best_cost = candidate_cost;
	*best_log_mean = candidate_log_mean;
	*data_i = it->data_i;
	*equality_constraint_active = it->equality_constraint_active;
      }
    }
  }
}

// check that this function is the minimum on all pieces.
int PiecewisePoissonLossLog::check_min_of
(PiecewisePoissonLossLog *prev, PiecewisePoissonLossLog *model){
  PoissonLossPieceListLog::iterator it;
  int verbose = 0;
  for(it = piece_list.begin(); it != piece_list.end(); it++){
    if(it != piece_list.begin()){
      PoissonLossPieceListLog::iterator pit = it;
      pit--;
      if(pit->max_log_mean != it->min_log_mean){
	printf("prev->max_log_mean != it->min_log_mean min\n");
	return 3;
      }
    }
    if(it->max_log_mean <= it->min_log_mean){
      printf("max_log_mean<=min_log_mean=%15.10f min\n", it->min_log_mean);
      return 2;
    }
    double mid_mean = (it->min_log_mean + it->max_log_mean)/2;
    double cost_min = it->getCost(mid_mean);
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
      PoissonLossPieceListLog::iterator pit = it;
      pit--;
      if(pit->max_log_mean != it->min_log_mean){
	printf("prev->max_log_mean != it->min_log_mean prev\n");
	return 3;
      }
    }
    if(it->max_log_mean <= it->min_log_mean){
      printf("max_log_mean<=min_log_mean=%15.10f prev\n", it->min_log_mean);
      return 2;
    }
    double mid_mean = (it->min_log_mean + it->max_log_mean)/2;
    double cost_prev = it->getCost(mid_mean);
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
      PoissonLossPieceListLog::iterator pit = it;
      pit--;
      if(pit->max_log_mean != it->min_log_mean){
	printf("prev->max_log_mean != it->min_log_mean model\n");
	return 3;
      }
    }
    if(it->max_log_mean <= it->min_log_mean){
      printf("max_log_mean<=min_log_mean=%15.10f model\n", it->min_log_mean);
      return 2;
    }
    double mid_mean = (it->min_log_mean + it->max_log_mean)/2;
    double cost_model = it->getCost(mid_mean);
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

void PiecewisePoissonLossLog::set_to_min_env_of
(PiecewisePoissonLossLog *fun1, PiecewisePoissonLossLog *fun2, int verbose){
  PoissonLossPieceListLog::iterator
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
    double last_max_log_mean = piece_list.back().max_log_mean;
    if(it1->max_log_mean == last_max_log_mean){
      it1++;
    }
    if(it2->max_log_mean == last_max_log_mean){
      it2++;
    }
  }
}

bool sameFuns
(PoissonLossPieceListLog::iterator it1,
 PoissonLossPieceListLog::iterator it2){
  return it1->Linear == it2->Linear &&
    it1->Log == it2->Log &&
    it1->Constant == it2->Constant;
}

void PiecewisePoissonLossLog::push_min_pieces
(PiecewisePoissonLossLog *fun1,
 PiecewisePoissonLossLog *fun2,
 PoissonLossPieceListLog::iterator it1,
 PoissonLossPieceListLog::iterator it2,
 int verbose){
  bool same_at_left;
  double last_min_log_mean;
  PoissonLossPieceListLog::iterator prev2 = it2;
  prev2--;
  PoissonLossPieceListLog::iterator prev1 = it1;
  prev1--;
  if(it1->min_log_mean < it2->min_log_mean){
    //it1 function piece starts to the left of it2.
    same_at_left = sameFuns(prev2, it1);
    last_min_log_mean = it2->min_log_mean;
  }else{
    //it1 function piece DOES NOT start to the left of it2.
    last_min_log_mean = it1->min_log_mean;
    if(it2->min_log_mean < it1->min_log_mean){
      //it2 function piece starts to the left of it1.
      same_at_left = sameFuns(prev1, it2);
    }else{
      //it1 and it2 start at the same min_log_mean value.
      if(it1==fun1->piece_list.begin() &&
	 it2==fun2->piece_list.begin()){
	same_at_left = false;
      }else{
	same_at_left = sameFuns(prev1, prev2);
      }
    }
  }
  PoissonLossPieceListLog::iterator next2 = it2;
  next2++;
  PoissonLossPieceListLog::iterator next1 = it1;
  next1++;
  bool same_at_right;
  double first_max_log_mean;
  if(it1->max_log_mean < it2->max_log_mean){
    // it2 function piece continues to the right of it1.
    same_at_right = sameFuns(next1, it2);
    first_max_log_mean = it1->max_log_mean;
  }else{
    first_max_log_mean = it2->max_log_mean;
    if(it2->max_log_mean < it1->max_log_mean){
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
  if(last_min_log_mean == first_max_log_mean){
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
    push_piece(it1, last_min_log_mean, first_max_log_mean);
    if(verbose)printf("exactly equal over entire interval\n");
    return;
  }
  PoissonLossPieceLog diff_piece
    (it1->Linear - it2->Linear,
     it1->Log - it2->Log,
     it1->Constant - it2->Constant,
     last_min_log_mean, first_max_log_mean,
     -5, false);
  // printf("it1->Constant=%a\nit2->Constant=%a\n",
  // 	 it1->Constant, it2->Constant);
  double mid_mean = (first_max_log_mean + last_min_log_mean)/2;
  double cost_diff_mid = diff_piece.getCost(mid_mean);
  if(diff_piece.Log == 0){
    // g(x) = Linear*e^x + Constant = 0,
    // x = log(-Constant/Linear).
    if(diff_piece.Linear == 0){
      // They are offset by a Constant.
      if(diff_piece.Constant < 0){
	push_piece(it1, last_min_log_mean, first_max_log_mean);
      }else{
	push_piece(it2, last_min_log_mean, first_max_log_mean);
      }
      if(verbose)printf("offset by a constant\n");
      return;
    }
    if(diff_piece.Constant == 0){
      // The only difference is the Linear coef.
      if(diff_piece.Linear < 0){
	push_piece(it1, last_min_log_mean, first_max_log_mean);
      }else{
	push_piece(it2, last_min_log_mean, first_max_log_mean);
      }
      if(verbose)printf("only diff is linear coef\n");
      return;
    }
    double log_mean_at_equal_cost = log(-diff_piece.Constant/(double)diff_piece.Linear);
    if(last_min_log_mean < log_mean_at_equal_cost &&
       log_mean_at_equal_cost < first_max_log_mean){
      // the root is in the interval, so we need to add two intervals.
      if(0 < diff_piece.Linear){
	push_piece(it1, last_min_log_mean, log_mean_at_equal_cost);
	push_piece(it2, log_mean_at_equal_cost, first_max_log_mean);
      }else{
	push_piece(it2, last_min_log_mean, log_mean_at_equal_cost);
	push_piece(it1, log_mean_at_equal_cost, first_max_log_mean);
      }
      if(verbose)printf("Log zero with one root in interval\n");
      return;
    }
    // the root is outside the interval, so one is completely above
    // the other over this entire interval.
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_max_log_mean);
    }
    if(verbose)printf("Log zero with no roots in interval\n");
    return;
  }//if(diff->Log == 0
  double cost_diff_left = diff_piece.getCost(last_min_log_mean);
  double cost_diff_right = diff_piece.getCost(first_max_log_mean);
  bool two_roots = diff_piece.has_two_roots(0.0);
  double smaller_mean, larger_mean;
  if(two_roots){
    smaller_mean = diff_piece.get_smaller_root(0.0);
    larger_mean = diff_piece.get_larger_root(0.0);
  }
  if(same_at_right){
    // they are equal on the right, but we don't know if there is
    // another crossing point somewhere to the left.
    if(two_roots){
      // there could be a crossing point to the left.
      double mean_at_equal_cost = smaller_mean;
      if(verbose)printf("smaller_mean=%a\n", mean_at_equal_cost);
      if(last_min_log_mean < mean_at_equal_cost &&
      	 mean_at_equal_cost < first_max_log_mean){
	//the cross point is in the interval.
	if(cost_diff_left < 0){
	  push_piece(it1, last_min_log_mean, mean_at_equal_cost);
	  push_piece(it2, mean_at_equal_cost, first_max_log_mean);
	}else{
	  push_piece(it2, last_min_log_mean, mean_at_equal_cost);
	  push_piece(it1, mean_at_equal_cost, first_max_log_mean);
	}
	if(verbose)printf("equal on the right with one crossing in interval\n");
	return;
      }
    }//if(two_roots
    // Test the cost at the midpoint, since the cost may be equal on
    // both the left and the right.
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_max_log_mean);
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
      if(last_min_log_mean < mean_at_equal_cost &&
      	 mean_at_equal_cost < first_max_log_mean){
	// the crossing point is in this interval.
	if(cost_diff_right < 0){
	  push_piece(it2, last_min_log_mean, mean_at_equal_cost);
	  push_piece(it1, mean_at_equal_cost, first_max_log_mean);
	}else{
	  push_piece(it1, last_min_log_mean, mean_at_equal_cost);
	  push_piece(it2, mean_at_equal_cost, first_max_log_mean);
	}
	if(verbose)printf("equal on the left with crossing in interval\n");
	return;
      }
    }//if(there may be crossing
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_max_log_mean);
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
      last_min_log_mean < larger_mean && larger_mean < first_max_log_mean;
    if(verbose)printf("smaller_mean=%f %a\nlarger_mean=%f %a\n",
		      smaller_mean, smaller_mean,
		      larger_mean, larger_mean);
    bool smaller_inside =
      last_min_log_mean < smaller_mean && smaller_mean < first_max_log_mean;
    if(larger_inside){
      if(smaller_inside && smaller_mean < larger_mean){
	// both are in the interval.
	first_mean = smaller_mean;
	second_mean = larger_mean;
	if(verbose){
	  diff_piece.print();
	  printf("%f and %f in [%f,%f]\n",
		 smaller_mean, larger_mean,
		 last_min_log_mean, first_max_log_mean);
	}
      }else{
	// smaller mean is not in the interval, but the larger is.
	first_mean = larger_mean;
	if(verbose){
	  printf("%f in [%f,%f]\n",
		 first_mean,
		 last_min_log_mean, first_max_log_mean);
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
		 last_min_log_mean, first_max_log_mean);
	}
      }
    }
  }//if(two_roots
  if(second_mean != INFINITY){
    // two crossing points.
    double before_mean = (last_min_log_mean + first_mean)/2;
    double cost_diff_before = diff_piece.getCost(before_mean);
    if(cost_diff_before < 0){
      push_piece(it1, last_min_log_mean, first_mean);
      push_piece(it2, first_mean, second_mean);
      push_piece(it1, second_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_mean);
      push_piece(it1, first_mean, second_mean);
      push_piece(it2, second_mean, first_max_log_mean);
    }
    if(verbose)printf("not equal on the sides, 2 crossing points\n");
  }else if(first_mean != INFINITY){
    // "one" crossing point. actually sometimes we have last_min_log_mean
    // < first_mean < first_max_log_mean but cost_diff_before and
    // cost_diff_after have the same sign! In that case we need to
    // just push one piece.
    double before_mean = (last_min_log_mean + first_mean)/2;
    double cost_diff_before = diff_piece.getCost(before_mean);
    if(verbose)printf("cost_diff_before(%.55f)=%f\n", before_mean, cost_diff_before);
    double after_mean = (first_max_log_mean + first_mean)/2;
    double cost_diff_after = diff_piece.getCost(after_mean);
    if(verbose)printf("cost_diff_after(%.55f)=%f\n", after_mean, cost_diff_after);
    if(cost_diff_before < 0){
      if(cost_diff_after < 0){
	// f1-f2<0 meaning f1<f2 on the entire interval, so just push it1.
	push_piece(it1, last_min_log_mean, first_max_log_mean);
      }else{
	push_piece(it1, last_min_log_mean, first_mean);
	push_piece(it2, first_mean, first_max_log_mean);
      }
    }else{//f1(before)-f2(before)>=0 meaning f1(before)>=f2(before)
      if(cost_diff_after < 0){
	//f1(after)-f2(after)<0 meaning f1(after)<f2(after)
	push_piece(it2, last_min_log_mean, first_mean);
	push_piece(it1, first_mean, first_max_log_mean);
      }else{
	//f1(after)-f2(after)>=0 meaning f1(after)>=f2(after)
	push_piece(it2, last_min_log_mean, first_max_log_mean);
      }
    }
    if(verbose)printf("not equal on the sides, 1 crossing point\n");
  }else{
    // "zero" crossing points. actually there may be a crossing point
    // in the interval that is numerically so close as to be identical
    // with last_min_log_mean or first_max_log_mean.
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_max_log_mean);
    }
    if(verbose)printf("not equal on the sides, zero crossing points\n");
  }
}

void PiecewisePoissonLossLog::push_piece
(PoissonLossPieceListLog::iterator it, double min_log_mean, double max_log_mean){
  PoissonLossPieceListLog::iterator last=piece_list.end();
  if(piece_list.size() && sameFuns(--last, it)){
    //it is the same function as last, so just make last extend
    //further to the right.
    last->max_log_mean = max_log_mean;
  }else{
    //it is a different function than last, so push it on to the end
    //of the list.
    piece_list.emplace_back
      (it->Linear, it->Log, it->Constant,
       min_log_mean, max_log_mean,
       it->data_i, it->equality_constraint_active);
  }
}
