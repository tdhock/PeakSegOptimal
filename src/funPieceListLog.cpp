/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "funPieceListLog.h"
#include <list>
#include <math.h>
#include <stdio.h>

#define NEWTON_EPSILON 1e-30
#define NEWTON_STEPS 100
#define PREV_NOT_SET (-3)

#define ABS(x) ((x)<0 ? -(x) : (x))

/*
* Square loss methods
*
* @author Sean Jewell 27 April 2017
*/

SquareLossPiece::SquareLossPiece
  (double a, double b, double c, double m, double M, int i, double prev){
  Square = a;
  Linear = b;
  Constant = c;
  min_mean = m;
  max_mean = M;
  data_i = i;
  prev_mean = prev;
}

SquareLossPiece::SquareLossPiece(){
}

bool SquareLossPiece::has_two_roots(double equals){
  // are there two solutions to the equation 
  // Square * u ^ 2 + 
  // Linear * u + Constant = equals ? 
  double delta = Linear * Linear - 4 * Square * (Constant - equals);
  if (delta > 0) {
    return true;
  } else {
    return false;
  }
}

double SquareLossPiece::SquareLoss(double mean){
  return Square * mean * mean + Linear * mean + Constant;
}

double SquareLossPiece::get_larger_root(double equals){
  double delta = Linear * Linear - 4 * Square * (Constant - equals);
  if (Square > 0) {
    return (-Linear + sqrt(delta)) / (2 * Square);
  } else {
    return (-Linear - sqrt(delta)) / (2 * Square);
  }
}

double SquareLossPiece::get_smaller_root(double equals){
  double delta = Linear * Linear - 4 * Square * (Constant - equals);
  if (Square < 0) {
    return (-Linear + sqrt(delta)) / (2 * Square);
  } else {
    return (-Linear - sqrt(delta)) / (2 * Square);
  }
}

double SquareLossPiece::argmin_mean(){
  // f(u) = Square * u ^ 2 + Linear * u + Constant,
  // f'(u)= 2 * Square * u + Linear = 0 means
  // u = -Linear / (2 * Square)
  
  if (Square > 0 || Square < 0) {
    return - Linear / (2 * Square);  
  } else if(Square == 0 & Linear > 0) {
    return min_mean;
  } else if(Square == 0 & Linear < 0) {
    return max_mean; 
  } else if(Square == 0 & Linear == 0) {
   return min_mean; 
  }
  throw;
}

/// what does this function do differently than argmin_mean()
double SquareLossPiece::argmin(){
  return argmin_mean();
}

double SquareLossPiece::getCost(double mean){
  if (mean < INFINITY) {
    return Square * mean * mean + Linear * mean + Constant;  
  } else if (mean == INFINITY) {
    if (Square > 0) {
      return INFINITY; 
    }
    if (Square < 0) {
      return -INFINITY;
    }
    if (Square == 0) {
      if (Linear > 0) {
        return INFINITY; 
      if (Linear < 0) {
        return -INFINITY;
      } 
      if (Linear == 0) {
        return Constant;
      }
    }
  }
}
  throw;
}

void PiecewiseSquareLoss::set_to_min_less_of
  (PiecewiseSquareLoss *input, double EPS, int verbose){
  piece_list.clear();
  SquareLossPieceList::iterator it = input->piece_list.begin();
  SquareLossPieceList::iterator next_it;
  double prev_min_cost = INFINITY;
  double prev_min_mean = it->min_mean;
  double prev_best_mean;
  while(it != input->piece_list.end()){
    
    if (it -> min_mean == EPS && it -> max_mean == EPS) {
      if (verbose) { printf("hitting EPS interval \n"); }
      
      // return constant function over two segments 
      piece_list.emplace_back
        (0, 0, it -> Constant,
         EPS, EPS, it -> data_i,
         it -> prev_mean);
      prev_min_cost = it -> Constant; 
      prev_best_mean = EPS;
      it++;
    } else {
    
    double left_cost = it->getCost(it->min_mean);
    double right_cost = it->getCost(it->max_mean);
    if(verbose)printf("left_cost=%f right_cost=%f\n", left_cost, right_cost);
    if(prev_min_cost == INFINITY){
      // Look for min achieved in this interval.
      if(verbose){
        printf("Searching for min in\n");
        it->print();
      }
      double next_left_cost;
      next_it = it;
      next_it++;
      double mu = it->argmin();
      double mu_cost = it->getCost(mu);
      bool next_ok;
      if(next_it == input->piece_list.end()){
        next_ok = true;
      }else{
        next_left_cost = next_it->getCost(next_it->min_mean);
        next_ok = NEWTON_EPSILON < next_left_cost-mu_cost; // probably can remove this?
      }
      if(verbose){
        printf("min cost=%f at mean=%f\n", mu_cost, mu);
        printf("next_left_cost-mu_cost=%e right_cost-mu_cost=%e\n", next_left_cost-mu_cost, right_cost-mu_cost);
      }
      bool cost_ok = NEWTON_EPSILON < right_cost-mu_cost && next_ok;
      if(mu <= it->min_mean && cost_ok){
        /* The minimum is achieved on the left or before this
        interval, so this function is always increasing in this
        interval. We don't need to store it, but we do need to keep
        track of the minimum cost, which occurs at the min mean
        value in this interval. */
        if(verbose)printf("min before interval\n");
        prev_min_cost = it->getCost(it->min_mean);
        prev_best_mean = it->min_mean;
      }else if(mu < it->max_mean && cost_ok){
        // Minimum in this interval, so add a convex piece up to the
        // min, and keep track of the min cost to create a constant
        // piece later. NB it is possible that prev_min_log_mean==mu in
        // which case we do not need to store the convex piece.
        if(verbose){
          printf("min in this interval at mean=%f cost=%f\n", mu, mu_cost);
          printf("right_cost=%f right-constant=%e\n", right_cost, right_cost-mu_cost);
          printf("next_left_cost=%f next-constant=%e\n", next_left_cost, next_left_cost-mu_cost);
        }
        if(prev_min_mean < mu){
          piece_list.emplace_back
          (it->Square, it->Linear, it->Constant, prev_min_mean, mu,
           PREV_NOT_SET, INFINITY); // equality constraint active on convex piece.
        }
        prev_min_mean = mu;
        prev_best_mean = mu;
        prev_min_cost = mu_cost;
        if(verbose)printf("prev_min_cost=%f\n", prev_min_cost);
      }else{
        // Minimum after this interval, so this function is
        // decreasing on this entire interval, and so we can just
        // store it as is.
        if(verbose)printf("min after interval\n");
        piece_list.emplace_back
          (it->Square, it->Linear, it->Constant, prev_min_mean, it->max_mean,
           PREV_NOT_SET, INFINITY); // equality constraint active on convex piece.
        prev_min_mean = it->max_mean;
      }//if(non-degenerate mu in interval
    } else{//prev_min_cost is finite
      // Look for a function with prev_min_cost in its interval.
      if(verbose){
        printf("Searching for intersection with %f\n", prev_min_cost);
        printf("cost at limits=[%f,%f] cost-constant=[%e,%e]\n",
               left_cost, right_cost,
               left_cost-prev_min_cost, right_cost-prev_min_cost);
        it->print();
      }
      // SWJ: Update this comment....still should look for 2 possible roots
      // Theoretically there can be zero, one, or two intersection
      // points between the constant function prev_min_log_mean and a
      // non-degenerate Poisson loss function piece.
      if(it->has_two_roots(prev_min_cost)){
        // There are two mean values where the Poisson loss piece
        // intersects the constant function prev_min_log_mean, but we
        // are only concerned with the first mean value (the
        // lesser of the two).
        double mu = it->get_smaller_root(prev_min_cost);
        if(it->min_mean < mu && mu < it->max_mean){
          // The smaller intersection point occurs within the
          // interval, so the constant interval ends here, and we
          // can store it immediately.
          piece_list.emplace_back
          (0, 0, prev_min_cost,
           prev_min_mean, mu, PREV_NOT_SET,
           prev_best_mean);// equality constraint inactive.
          prev_min_cost = INFINITY;
          prev_min_mean = mu;
          it--;
        }//if(mu in interval)
      }//if(has two roots
      if(right_cost <= prev_min_cost + NEWTON_EPSILON &&
         prev_min_cost < INFINITY){
        //ends exactly/numerically on the right.
        if(verbose)printf("constant numerically equal on right\n");
        piece_list.emplace_back
          (0, 0, prev_min_cost,
           prev_min_mean, it->max_mean,
           PREV_NOT_SET,
           prev_best_mean);
        prev_min_cost = INFINITY;
        prev_min_mean = it->max_mean;
      }
    }//if(prev_min_cost is finite
    it++;
    if(verbose){
      printf("current min-less-------------------\n");
      print();
    }
    }
  }//while(it
  if(prev_min_cost < INFINITY){
    // ending on a constant piece.
    it--;
    piece_list.emplace_back
      (0, 0, prev_min_cost,
       prev_min_mean, it->max_mean, PREV_NOT_SET,
       prev_best_mean);//equality constraint inactive on constant piece.
  }
}

void PiecewiseSquareLoss::set_to_clean(PiecewiseSquareLoss *input, double EPS, int verbose) {
  piece_list.clear();
  SquareLossPieceList::iterator it = input->piece_list.begin();
  while(it != input->piece_list.end()){
    
    // three cases
    // 1. in the EPS region
    // 2. scaling moves into EPS region
    // 3. no scaling problems
    
    // in the EPS region -- no scaling
    if (it -> min_mean == EPS && it -> max_mean == EPS) {
      piece_list.emplace_back(0, 0, it -> Constant,
                              EPS, EPS, it -> data_i, it -> prev_mean);
    } else {
      
      double Square =  it -> Square;
      double Linear = it -> Linear;
      double min_mean = it -> min_mean;
      double max_mean = it -> max_mean;
      
      if (max_mean > EPS) {
        piece_list.emplace_back(Square, Linear, it -> Constant,
                                min_mean, max_mean,it -> data_i, it -> prev_mean);              
      }
    }
    it++;
  } //end while
  
}






void PiecewiseSquareLoss::set_to_scaled_of(PiecewiseSquareLoss *input,
                                           double gam, double EPS, int verbose) {
  piece_list.clear();
  SquareLossPieceList::iterator it = input->piece_list.begin();
  while(it != input->piece_list.end()){
    
    // three cases
    // 1. in the EPS region
    // 2. scaling moves into EPS region
    // 3. no scaling problems
    
    // in the EPS region -- no scaling
    if (it -> min_mean == EPS && it -> max_mean == EPS) {
      piece_list.emplace_back(0, 0, it -> Constant,
                              EPS, EPS, it -> data_i, it -> prev_mean);
    } else {
    
    double Square =  (it -> Square) / (gam * gam);
    double Linear = (it -> Linear) / gam;
    double min_mean = (it -> min_mean) * gam;
    double max_mean = (it -> max_mean) * gam;
    
    // scaling moves into the EPS region 
    if (min_mean <= EPS) { 
      // need a new max mean for this segment
      // add a new segment
      if (max_mean > EPS) { 
        piece_list.emplace_back(Square, Linear, it -> Constant,
                                min_mean, EPS,it -> data_i, it -> prev_mean);  
        min_mean = EPS; 
      }
      
    }
    
    piece_list.emplace_back(Square, Linear, it -> Constant,
                              min_mean, max_mean,it -> data_i, it -> prev_mean);      
    }
    it++;
  } //end while
  
}


// output PiecewiseSquareLoss
// that is constant on two intervals
// min over eps interval
// min over (eps, infty) 

void PiecewiseSquareLoss::set_to_eps_min_of
  (PiecewiseSquareLoss *input, double EPS, int verbose) {
  piece_list.clear();
  SquareLossPieceList::iterator it = input->piece_list.begin();
  
  // eps interval costs
  double prev_min_cost_EPS = INFINITY;
  double left_most_EPS = INFINITY;
  double right_most_EPS = -INFINITY;
  double right_mean_EPS, left_mean_EPS, prev_best_mean_EPS, mu_cost_EPS, mu_EPS;
  int data_i;
  
  while(it != input->piece_list.end()){
    if (verbose) {
      printf("start new iter of eps min--------------\n");
      printf("Searching for min cost in \n");
      printf("%10s %10s %15s %15s %15s %15s %s\n",
             "Square", "Linear", "Constant",
             "min_mean", "max_mean",
             "prev_mean", "data_i");
      it -> print();
      printf("current values of prev_min_cost_EPS and data_i %f %d \n", 
             prev_min_cost_EPS, data_i);
    }
    if (it -> min_mean == EPS && it -> max_mean == EPS) {
      if (verbose) { printf("hitting EPS interval \n"); }
      
      mu_EPS = EPS; 
      mu_cost_EPS = it -> Constant;
      
      if (mu_cost_EPS < prev_min_cost_EPS) {
        if (verbose) {printf("setting the min cost at pt %d (new cost, prev cost) \t (%f, %f) \n",
            it -> data_i, mu_cost_EPS, prev_min_cost_EPS);} 
        prev_min_cost_EPS = mu_cost_EPS;
        prev_best_mean_EPS = mu_EPS;
        data_i = it -> data_i;
      }
      
      
      if (verbose) {printf("cost at pt %d \t %f \n", it -> data_i, mu_cost_EPS);}
      
    } else if (it -> max_mean <= EPS) {
      // we are in an interval less than EPS
      if (verbose) { printf("max mean less than EPS \n"); }
      // boundary costs
      left_mean_EPS = it -> min_mean;
      
      if (left_mean_EPS < left_most_EPS) {
        left_most_EPS = left_mean_EPS;
      }
      
      right_mean_EPS = it -> max_mean;
      
      if (right_mean_EPS > right_most_EPS) {
        right_most_EPS = right_mean_EPS;
      }
      
      // determine argmin of this segment
      mu_EPS = it->argmin();
      // argmin occurs in interval
      if (mu_EPS >= left_mean_EPS && mu_EPS <= right_mean_EPS) {
        mu_cost_EPS = it->getCost(mu_EPS);
        // printf("data_i %d (within interval) \n", it -> data_i);
      } else {
        double left_cost_EPS = it -> getCost(left_mean_EPS);
        double right_cost_EPS = it -> getCost(right_mean_EPS);
        
        if (right_cost_EPS < left_cost_EPS) {
          mu_cost_EPS = right_cost_EPS;
          mu_EPS = right_mean_EPS;
          if (verbose) { printf("data_i %d (right interval) \n", it -> data_i);}
        } else {
          mu_cost_EPS = left_cost_EPS;
          mu_EPS = left_mean_EPS;
          if (verbose) {printf("data_i %d (left interval) \n", it -> data_i);}
        }
      }
      
      if (mu_cost_EPS < prev_min_cost_EPS) {
        if (verbose) {printf("setting the min cost at pt %d (new cost, prev cost) \t (%f, %f) \n",
            it -> data_i, mu_cost_EPS, prev_min_cost_EPS);}  
        prev_min_cost_EPS = mu_cost_EPS;
        prev_best_mean_EPS = mu_EPS;
        data_i = it -> data_i;
      }
      
      if (verbose) {printf("cost at pt %d \t %f \n", it -> data_i, mu_cost_EPS);}
      
      
    } 
    it++;
  } // end loop over components
  
  
  // return constant @ EPS interval
    piece_list.emplace_back
    (0, 0, prev_min_cost_EPS,
     EPS, EPS, data_i,
     prev_best_mean_EPS);
  
  // loop over elements and return components that 
  // have max_mean > EPS
  it = input->piece_list.begin();
  while(it != input->piece_list.end()){
  
  if (it -> max_mean > EPS) {
    piece_list.emplace_back
    (it -> Square, it -> Linear, it -> Constant,
     it -> min_mean, it -> max_mean, it -> data_i,
     it -> prev_mean);  
  }
  it++;
  }
  
}



void PiecewiseSquareLoss::set_to_unconstrained_min_of
  (PiecewiseSquareLoss *input, double EPS, int verbose) {
  piece_list.clear();
  SquareLossPieceList::iterator it = input->piece_list.begin();
  // outside of eps interval costs
  double prev_min_cost = INFINITY;
  double left_most = INFINITY;
  double right_most = -INFINITY;
  double right_mean, left_mean, prev_best_mean, mu_cost, mu;
  
  while(it != input->piece_list.end()){
    if (verbose) {
      printf("start new iter of set to unconstrained min of--------------\n");
      printf("Searching for min cost in \n");
      printf("%10s %10s %15s %15s %15s %15s %s\n",
             "Square", "Linear", "Constant",
             "min_mean", "max_mean",
             "prev_mean", "data_i");
      it -> print();

    }
    if (it -> min_mean == EPS && it -> max_mean == EPS) {
      if (verbose) { printf("hitting EPS interval \n"); }
      
      // return constant function over two segments 
      piece_list.emplace_back
        (0, 0, it -> Constant,
         EPS, EPS, it -> data_i,
         it -> prev_mean);
      
    }
      
    // boundary costs
    left_mean = it -> min_mean;
    
    if (left_mean < left_most) {
      left_most = left_mean;
    }
    
    right_mean = it -> max_mean;
    
    if (right_mean > right_most) {
      right_most = right_mean;
    }
    
    // determine argmin of this segment
    mu = it->argmin();
    // argmin occurs in interval
    if (mu >= left_mean && mu <= right_mean) {
      mu_cost = it->getCost(mu);
    } else {
      double left_cost = it -> getCost(left_mean);
      double right_cost = it -> getCost(right_mean);
      
      if (right_cost < left_cost) {
        mu_cost = right_cost;
        mu = right_mean;
      } else {
        mu_cost = left_cost;
        mu = left_mean;
      }
    }
    
    if (mu_cost < prev_min_cost) {
      prev_min_cost = mu_cost;
      prev_best_mean = mu;
    }
    it++;
  } // end loop over components
  
    piece_list.emplace_back
    (0, 0, prev_min_cost,
     left_most, right_most, PREV_NOT_SET,
     prev_best_mean);
  
  if (verbose) {
    printf("interval [%f, %f]\n", left_most, right_most);
    printf("Minimum cost %f \n", prev_min_cost);
    printf("------------------------------------------\n");
  }
}

void PiecewiseSquareLoss::add(double Square, double Linear, double Constant){
  SquareLossPieceList::iterator it;
    for(it=piece_list.begin(); it != piece_list.end(); it++){
      it->Square += Square;
      it->Linear += Linear;
      it->Constant += Constant;
    } 
  }

void PiecewiseSquareLoss::add_protect(double Square, double Linear, double Constant, double EPS, bool eps_segment){
  SquareLossPieceList::iterator it;
  if (eps_segment) {
    for(it=piece_list.begin(); it != piece_list.end(); it++){
      if (it -> min_mean == EPS && it -> max_mean == EPS) {
        it->Constant += Constant;  
      }
    } 
  } else {
    for(it=piece_list.begin(); it != piece_list.end(); it++){
        if (it -> max_mean > EPS) {
          it->Square += Square;
          it->Linear += Linear;
          it->Constant += Constant;   
        }
    } 
  }
}

void PiecewiseSquareLoss::add_penalty(double penalty, double EPS){
  SquareLossPieceList::iterator it;
    for(it=piece_list.begin(); it != piece_list.end(); it++){
      if (it -> max_mean > EPS) {
        it->Constant += penalty;   
      }
    } 
}



void PiecewiseSquareLoss::multiply(double x){
  SquareLossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->Square *= x;
    it->Linear *= x;
    it->Constant *= x;
  }
}

void PiecewiseSquareLoss::scale(double gam){
  SquareLossPieceList::iterator it;
  
  //  double Square =  (it -> Square) / (gam * gam);
  //  double Linear = (it -> Linear) / gam;
  //  double min_mean = (it -> min_mean) * gam;
  //  double max_mean = (it -> max_mean) * gam;
  
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->Square /= (gam * gam);
    it->Linear /= gam;
    it->min_mean *= gam;
    it->max_mean *= gam;
  }
}



void PiecewiseSquareLoss::set_prev_seg_end(int prev_seg_end, double EPS){
  SquareLossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    if (it -> max_mean > EPS) {
      it->data_i = prev_seg_end; 
    }
  }
}

void PiecewiseSquareLoss::findMean
  (double mean, int *seg_end, double *prev_mean){
  SquareLossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    //	  printf("looking for mean %f in interval [%f, %f]\n", mean, it -> min_mean, it -> max_mean);
    if(it->min_mean <= mean && mean <= it->max_mean){
      *seg_end = it->data_i;
      *prev_mean = it->prev_mean;
      //      printf("peel off mean %f \n", &(it -> prev_mean));
      return;
    }
  }
}

double PiecewiseSquareLoss::findCost(double mean){
  SquareLossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    if(it->min_mean <= mean && mean <= it->max_mean){
      int verbose = 0;
      return it->getCost(mean);
    }
  }
  throw;
}

void PiecewiseSquareLoss::print(){
  SquareLossPieceList::iterator it;
  printf("%10s %10s %15s %15s %15s %15s %s\n",
         "Square", "Linear", "Constant",
         "min_mean", "max_mean",
         "prev_mean", "data_i");
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->print();
  }
}

void PiecewiseSquareLoss::checkStable(double MAX){
  SquareLossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    if (it -> Square > MAX) {
      printf("Numerically unstable in interval:\n");
      printf("%10s %10s %15s %15s %15s %15s %s\n",
             "Square", "Linear", "Constant",
             "min_mean", "max_mean",
             "prev_mean", "data_i");
      it -> print(); 
      throw (double) it -> Square; 
    }
  }
}

void SquareLossPiece::print(){
  printf("%.20e %.20e %.20e %.20e %.20e %15f %d\n",
         Square, Linear, Constant,
         min_mean, max_mean,
         prev_mean, data_i);
}

void PiecewiseSquareLoss::Minimize
  (double *best_cost,
   double *best_mean,
   int *data_i,
   double *prev_mean){
  double candidate_cost, candidate_mean;
  int verbose=false;
  SquareLossPieceList::iterator it;
  *best_cost = INFINITY;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    candidate_mean = it->argmin();
    if(candidate_mean < it->min_mean){
      candidate_mean = it->min_mean;
    }else if(it->max_mean < candidate_mean){
      candidate_mean = it->max_mean;
    }
    candidate_cost = it->getCost(candidate_mean);
    if(candidate_cost < *best_cost){
      *best_cost = candidate_cost;
      *best_mean = candidate_mean;
      *data_i = it->data_i;
      *prev_mean = it->prev_mean;
    }
  }
}

double MidMean(double first, double second) {
  if (first == -INFINITY && second == INFINITY) {
    return 0.0;
  } else if (first == -INFINITY && second != INFINITY) {
    return second - 1;
  } else if (first != -INFINITY && second == INFINITY) {
    return first + 1;
  } else {
    return ((first + second)  / 2);
  }
}

// check that this function is the minimum on all pieces.
int PiecewiseSquareLoss::check_min_of
  (PiecewiseSquareLoss *prev, PiecewiseSquareLoss *model){
  SquareLossPieceList::iterator it;
  for(it = piece_list.begin(); it != piece_list.end(); it++){
    if(it != piece_list.begin()){
      SquareLossPieceList::iterator pit = it;
      pit--;
      if(pit->max_mean != it->min_mean){
        printf("prev->max_mean != it->min_mean min\n");
        return 3;
      }
    }
    if(it->max_mean -  it->min_mean <= -NEWTON_EPSILON){
      printf("max_mean (=%15.10f) <=min_mean(=%15.10f) min\n", it -> max_mean, it->min_mean);
      return 2;
    }
    double mid_mean = MidMean(it -> min_mean, it -> max_mean); //(it->min_mean + it->max_mean)/2;
    //    printf("printing out mid mean %10.5f\n", mid_mean);
    if(-INFINITY < mid_mean){
      double cost_min = it->getCost(mid_mean);
      double cost_prev = prev->findCost(mid_mean);
      if(cost_prev+1e-6 < cost_min){
        printf("prev(%f)=%f\n", mid_mean, cost_prev);
        prev->print();
        printf("min(%f)=%f\n", mid_mean, cost_min);
        print();
        return 1;
      }
      double cost_model = model->findCost(mid_mean);
      if(cost_model+1e-6 < cost_min){
        printf("model(%.20e)=%f\n", mid_mean, cost_model);
        model->print();
        printf("min(%.20e)=%f\n", mid_mean, cost_min);
        print();
        return 1;
      }
    }
  }
  for(it = prev->piece_list.begin(); it != prev->piece_list.end(); it++){
    if(it != prev->piece_list.begin()){
      SquareLossPieceList::iterator pit = it;
      pit--;
      if(pit->max_mean != it->min_mean){
        printf("prev->max_mean != it->min_mean prev\n");
        return 3;
      }
    }
    if(it->max_mean - it->min_mean <= -NEWTON_EPSILON){
      printf("max_mean<=min_mean=%15.10f prev\n", it->min_mean);
      return 2;
    }
    double mid_mean = MidMean(it -> min_mean, it -> max_mean);
    if(-INFINITY < mid_mean){
      double cost_prev = it->getCost(mid_mean);
      double cost_min = findCost(mid_mean);
      if(cost_prev+1e-6 < cost_min){
        printf("prev(%f)=%f\n", mid_mean, cost_prev);
        prev->print();
        printf("min(%f)=%f\n", mid_mean, cost_min);
        print();
        return 1;
      }
    }
  }
  for(it = model->piece_list.begin(); it != model->piece_list.end(); it++){
    if(it != model->piece_list.begin()){
      SquareLossPieceList::iterator pit = it;
      pit--;
      if(pit->max_mean != it->min_mean){
        printf("prev->max_mean != it->min_mean model\n");
        return 3;
      }
    }
    if(it->max_mean - it->min_mean <= -NEWTON_EPSILON){
      printf("max_mean=%15.10f<=min_mean=%15.10f model\n", it-> max_mean, it->min_mean);
      return 2;
    }
    double mid_mean = MidMean(it -> min_mean, it -> max_mean);
    if(-INFINITY < mid_mean){
      double cost_model = it->getCost(mid_mean);
      double cost_min = findCost(mid_mean);
      if(cost_model+1e-6 < cost_min){
        printf("model(%f)=%f\n", mid_mean, cost_model);
        model->print();
        printf("min(%f)=%f\n", mid_mean, cost_min);
        print();
        return 1;
      }
    }
  }
  return 0;
}

void PiecewiseSquareLoss::set_to_min_env_of
  (PiecewiseSquareLoss *fun1, PiecewiseSquareLoss *fun2, double EPS, int verbose){
  SquareLossPieceList::iterator
  it1 = fun1->piece_list.begin(),
  it2 = fun2->piece_list.begin();
  if(verbose){
    printf("computing min env of:\n");
    printf("=min-less/more\n");
    fun1->print();
    printf("=cost model\n");
    fun2->print();
  }
  piece_list.clear();
  
  if (it1 -> min_mean == it1 -> max_mean &&
      it2 -> min_mean == it2 -> max_mean) {
  
  // take care of initial eps piece 
  if (it1 -> Constant < it2 -> Constant) {
    piece_list.emplace_back(0, 0, it1->Constant,
     it1 -> min_mean, it1 -> max_mean,
     it1->data_i, it1->prev_mean);
  } else { 
    piece_list.emplace_back(0, 0, it2->Constant,
                            it2 -> min_mean, it2 -> max_mean,
                            it2->data_i, it2->prev_mean);
  }
  
  it1++;
  it2++;
  
  // start at second element of both of these lists
  // fun1 -> piece_list.(it1++);
  // fun2 -> piece_list.erase(it2++);
  
  }
  

  while(it1 != fun1->piece_list.end() &&
        it2 != fun2->piece_list.end()){
    push_min_pieces(fun1, fun2, it1, it2, EPS, verbose);
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

bool sameFunsSquare
  (SquareLossPieceList::iterator it1,
   SquareLossPieceList::iterator it2){
  return it1->Linear == it2->Linear &&
    it1->Square == it2->Square &&
    ABS(it1->Constant - it2->Constant) < NEWTON_EPSILON;
}

void PiecewiseSquareLoss::push_min_pieces
  (PiecewiseSquareLoss *fun1,
   PiecewiseSquareLoss *fun2,
   SquareLossPieceList::iterator it1,
   SquareLossPieceList::iterator it2,
   double EPS, 
   int verbose){
  bool same_at_left;
  double last_min_mean;
  SquareLossPieceList::iterator prev2 = it2;
  prev2--;
  SquareLossPieceList::iterator prev1 = it1;
  prev1--;
  

  if(it1->min_mean < it2->min_mean){
    //it1 function piece starts to the left of it2.
    same_at_left = sameFunsSquare(prev2, it1);
    last_min_mean = it2->min_mean;
  }else{
    //it1 function piece DOES NOT start to the left of it2.
    last_min_mean = it1->min_mean;
    if(it2->min_mean < it1->min_mean){
      //it2 function piece starts to the left of it1.
      same_at_left = sameFunsSquare(prev1, it2);
    }else{
      //it1 and it2 start at the same min_mean value.
      if(it1==fun1->piece_list.begin() &&
         it2==fun2->piece_list.begin()){
        same_at_left = false;
        }else{
         same_at_left = sameFunsSquare(prev1, prev2);
          if (prev1 -> min_mean == EPS && 
              prev1 -> max_mean == EPS &&
              prev2 -> min_mean == EPS &&
              prev2 -> max_mean == EPS)
            same_at_left = false;
       }
    }
  }
  SquareLossPieceList::iterator next2 = it2;
  next2++;
  SquareLossPieceList::iterator next1 = it1;
  next1++;
  bool same_at_right;
  double first_max_mean;
  if(it1->max_mean < it2->max_mean){
    if(verbose)printf("it2 function piece continues to the right of it1.\n");
    same_at_right = sameFunsSquare(next1, it2);
    first_max_mean = it1->max_mean;
  }else{
    first_max_mean = it2->max_mean;
    if(it2->max_mean < it1->max_mean){
      if(verbose)printf("it2 function piece ends before it1.\n");
      same_at_right = sameFunsSquare(it1, next2);
    }else{
      if(verbose)printf("it2 and it1 end at same max_mean.\n");
      if(next1==fun1->piece_list.end() &&
         next2==fun2->piece_list.end()){
        if(verbose)printf("at the end so they can't be equal after this interval.\n");
        same_at_right = false;
      }else{
        if(verbose){
          printf("comparing next function pieces.\n");
          next1->print();
          next2->print();
        }
        same_at_right = sameFunsSquare(next1, next2);
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
    
    // if ((it1 -> Constant) < (it2 -> Constant)) {
    //   push_piece(it1, last_min_mean, first_max_mean);
    // } else {
    //   push_piece(it2, last_min_mean, first_max_mean);
    // }
    return;
  }
  if(sameFunsSquare(it1, it2)){
    // The functions are exactly equal over the entire interval so we
    // can push either of them.
    push_piece(it1, last_min_mean, first_max_mean);
    if(verbose)printf("exactly equal over entire interval\n");
    return;
  }
  SquareLossPiece diff_piece
    (it1->Square - it2->Square,
     it1->Linear - it2->Linear,
     it1->Constant - it2->Constant,
     last_min_mean, first_max_mean,
     -5, false);
  // Evaluate the middle in the original space, to avoid problems when
  // SWJ: need to be careful here in case either first_max_mean of last_min_mean
  // is -Inf or Inf
  // If either is Inf or -Inf then
  
  double mid_mean;
  if (first_max_mean != INFINITY && last_min_mean != -INFINITY) {
    mid_mean = (first_max_mean + last_min_mean) / 2;
  } else if (first_max_mean != INFINITY && last_min_mean == -INFINITY) {
    mid_mean = first_max_mean - 1;
  } else if (first_max_mean == INFINITY && last_min_mean != -INFINITY) {
    mid_mean = last_min_mean + 1;
  } else  {
    mid_mean = 0;
  }
  
  double cost_diff_mid = diff_piece.getCost(mid_mean);
  // Easy case of equality on both left and right.
  if(same_at_left && same_at_right){
    if(verbose)printf("Same on both the left and the right\n");
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    return;
  }
  
  // add back the easy degenerate cases that do not require root finding
  
  // Easy degenerate cases that do not require root finding.
  if(diff_piece.Square == 0){
    // g(x) = Linear * x + Constant = 0,
    // x = - Constant / Linear
    if(diff_piece.Linear == 0){
      // They are offset by a Constant.
      if(diff_piece.Constant < 0){
        push_piece(it1, last_min_mean, first_max_mean);
      }else{
        push_piece(it2, last_min_mean, first_max_mean);
      }
      if(verbose)printf("offset by a constant=%e\n", diff_piece.Constant);
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
      if(verbose)printf("Square zero with one root in interval\n");
      return;
    }
    // the root is outside the interval, so one is completely above
    // the other over this entire interval.
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    if(verbose)printf("Square zero with no roots in interval\n");
    return;
  }
  
  
  double cost_diff_left = diff_piece.getCost(last_min_mean);
  double cost_diff_right = diff_piece.getCost(first_max_mean);
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
      double mean_at_crossing = smaller_mean;
      double mean_between_zeros = (mean_at_crossing + first_max_mean)/2;
      double cost_between_zeros = diff_piece.getCost(mean_between_zeros);
      double mean_at_optimum = diff_piece.argmin();
      if(verbose){
        printf("cost_diff(left:%e)=%e\n", last_min_mean, cost_diff_left);
        printf("cost_diff(cross:%e)=%e\n", mean_at_crossing, diff_piece.getCost(mean_at_crossing));
        printf("cost_diff(between:%e)=%e\n", mean_between_zeros, cost_between_zeros);
        printf("cost_diff(optimum:%e)=%e\n", mean_at_optimum, diff_piece.getCost(mean_at_optimum));
        printf("cost_diff(right:%e)=%e\n", first_max_mean, cost_diff_right);
      }
      if(last_min_mean < mean_at_crossing &&
         mean_at_crossing < mean_at_optimum &&
         mean_at_optimum < first_max_mean){
        //the cross point is in the interval.
        if(cost_diff_left < 0){
          push_piece(it1, last_min_mean, mean_at_crossing);
          push_piece(it2, mean_at_crossing, first_max_mean);
        }else{
          push_piece(it2, last_min_mean, mean_at_crossing);
          push_piece(it1, mean_at_crossing, first_max_mean);
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
      double mean_at_crossing = larger_mean;
      double mean_at_optimum = diff_piece.argmin();
      if(verbose)printf("mean_at_crossing=%f\n", mean_at_crossing);
      if(verbose)printf("last_min_mean=%f\n", last_min_mean);
      if(verbose)printf("mean_at_optimum=%f\n", mean_at_optimum);
      if(verbose)printf("first_max_mean=%f\n", first_max_mean);
      if(last_min_mean < mean_at_optimum &&
         mean_at_optimum < mean_at_crossing &&
         mean_at_crossing < first_max_mean){
        // the crossing point is in this interval.
        if(cost_diff_right < 0){
          push_piece(it2, last_min_mean, mean_at_crossing);
          push_piece(it1, mean_at_crossing, first_max_mean);
        }else{
          push_piece(it1, last_min_mean, mean_at_crossing);
          push_piece(it2, mean_at_crossing, first_max_mean);
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
   // printf("same at left value %d \n", same_at_left);
    if (verbose) printf("cost diff left, mid, right (%f, %f, %f) \n", cost_diff_left, cost_diff_mid, cost_diff_right);
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
      last_min_mean < smaller_mean &&
      // 0 < smaller_mean && // breaks for negative data
      smaller_mean < first_max_mean;
    if(larger_inside){
      if(smaller_inside && smaller_mean < larger_mean){
        // both are in the interval.
        first_mean = smaller_mean;
        second_mean = larger_mean;
        if(verbose){
          diff_piece.print();
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
    // SWJ: Need to be careful for left-most limit as this could be -inf
    // SWJ: Can we just check the sign of the Square coef?
    //    double before_mean;
    //    if (last_min_mean != -INFINITY) {
    //      before_mean = (last_min_mean + first_mean )/2;
    //    } else {
    //      before_mean = first_mean - 1; // should be able to compare anywhere in the interval (-inf, first_mean] ?
    //    }
    //
    //    double cost_diff_before = diff_piece.getCost(before_mean);
    double squareDiff = diff_piece.Square;
    if(squareDiff < 0){
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
    // "one" crossing point. actually sometimes we have last_min_log_mean
    // < first_log_mean < first_max_log_mean but cost_diff_before and
    // cost_diff_after have the same sign! In that case we need to
    // just push one piece.
    
    // need to be careful here, too
    // last_min_mean could be -Inf
    double before_mean;
    if (last_min_mean != -INFINITY) {
      before_mean = (last_min_mean + first_mean) / 2;
    } else {
      before_mean = first_mean - 1;
    }
    
    double cost_diff_before = diff_piece.getCost(before_mean);
    if(verbose){
      printf("cost_diff_before(%.55f)=%f\n", before_mean, cost_diff_before);
    }
    
    // SWJ: need to be careful as first_max_mean could be inf
    double after_mean;
    
    if (first_max_mean != INFINITY) {
      after_mean = (first_max_mean + first_mean)/2;
    } else {
      after_mean = first_mean + 1; // should be able to check anywhere in the interval (first_mean, inf)?
    }
    
    
    double cost_diff_after = diff_piece.getCost(after_mean);
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
    // with last_min_log_mean or first_max_log_mean.
    if(verbose){
      printf("not equal on the sides, zero crossing points\n");
      printf("cost_diff left=%e mid=%e right=%e\n",
             cost_diff_left, cost_diff_mid, cost_diff_right);
    }
    double cost_diff;
    if(first_max_mean == INFINITY){
      // if we are at the last interval and the right limit is
      // infinity, then it should be fine to compare the cost at any
      // point in the interval.
      
      // SWJ: I'm not sure that's true
      // Try something else instead, throw error just in case
      double minimizer = diff_piece.argmin();
      bool min_inside =
        last_min_mean < minimizer && minimizer < first_max_mean;
      if (min_inside) {
        cost_diff = diff_piece.getCost(minimizer + 1);
      } else {
        cost_diff = diff_piece.getCost(last_min_mean + 1);
      }
      
    }else{
      if(ABS(cost_diff_mid) < NEWTON_EPSILON){
        cost_diff = cost_diff_right;
      }else{
        cost_diff = cost_diff_mid;
      }
    }
    if(cost_diff < 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
  }
}

void PiecewiseSquareLoss::push_piece
  (SquareLossPieceList::iterator it, double min_mean, double max_mean){
  if(max_mean <= min_mean){
    return;
  }
  SquareLossPieceList::iterator last=piece_list.end();
  --last;
  if(piece_list.size() &&
  sameFunsSquare(last, it) &&
  it->prev_mean == last->prev_mean &&
  it->data_i == last->data_i){
    //it is the same function as last, so just make last extend
    //farther to the right.
    last->max_mean = max_mean;
  }else{
    //it is a different function than last, so push it on to the end
    //of the list.
    piece_list.emplace_back
    (it->Square, it->Linear, it->Constant,
     min_mean, max_mean,
     it->data_i, it->prev_mean);
  }
}
