/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <R.h>
#include <math.h>
#include <list>
#include <vector>

#include "funPieceListLog.h"
#include "LossPiece.h"
#include "NormalFPOPisotonic.h"

#define EPSILON 1e-12
#define ABS(x) ((x)<0 ? -(x) : (x))

class NormalLossPiece : public LossPiece {
public:
  double Quadratic, Linear, Constant;
  NormalLossPiece();
  NormalLossPiece(double q, double l, double c,
                  double m, double M, int i, double prev);
  double getCost(double mean);
  double argmin();
  double argmin_interval();
  bool has_two_roots(double equals);
  double get_smaller_root(double equals);
  double get_larger_root(double equals);
  void print();
};

typedef std::list<NormalLossPiece> NormalLossPieceList;

class PiecewiseNormalLoss {
public:
  NormalLossPieceList piece_list;
  void set_to_min_less_of(PiecewiseNormalLoss *input);
  void set_to_min_env_of(PiecewiseNormalLoss *fun1,
                         PiecewiseNormalLoss *fun2, int verbose);
  void push_min_pieces(PiecewiseNormalLoss *fun1,
                       PiecewiseNormalLoss *fun2,
                       NormalLossPieceList::iterator it1,
                       NormalLossPieceList::iterator it2, int verbose);
  void push_piece(NormalLossPieceList::iterator it,
                  double min_mean, double max_mean);
  void add(double Quadratic, double Linear, double Constant);
  void multiply(double x);
  void set_prev_seg_end(int prev_seg_end);
  void findMean(double mean, int *seg_end, double *prev_mean);
  double findCost(double mean);
  void Minimize(double *best_cost, double *best_mean,
                int *data_i, double *prev_mean);
  void print();
};

bool sameFuns(NormalLossPieceList::iterator it1,
              NormalLossPieceList::iterator it2);

NormalLossPiece::NormalLossPiece()
  : LossPiece(), Quadratic(0), Linear(0), Constant(0) {}

NormalLossPiece::NormalLossPiece
(double q, double l, double c, double m, double M, int i, double prev)
  : LossPiece(m, M, i, prev), Quadratic(q), Linear(l), Constant(c) {}

double NormalLossPiece::getCost(double mean){
  return Quadratic*mean*mean + Linear*mean + Constant;
}

double NormalLossPiece::argmin(){
  if(Quadratic == 0){
    if(Linear < 0) return INFINITY;
    if(Linear > 0) return -INFINITY;
    return 0;
  }
  return -Linear / (2.0*Quadratic);
}

double NormalLossPiece::argmin_interval(){
  double a = argmin();
  if(a < min_mean) return min_mean;
  if(a > max_mean) return max_mean;
  return a;
}

bool NormalLossPiece::has_two_roots(double equals){
  if(Quadratic == 0) return false;
  double disc = Linear*Linear - 4.0*Quadratic*(Constant - equals);
  return disc > EPSILON;
}

double NormalLossPiece::get_smaller_root(double equals){
  double a = Quadratic, b = Linear, c = Constant - equals;
  if(a == 0){
    if(b == 0) return -INFINITY;
    return -c/b;
  }
  double disc = b*b - 4.0*a*c;
  if(disc < 0) disc = 0;
  double r1 = (-b - sqrt(disc)) / (2.0*a);
  double r2 = (-b + sqrt(disc)) / (2.0*a);
  return (r1 < r2) ? r1 : r2;
}

double NormalLossPiece::get_larger_root(double equals){
  double a = Quadratic, b = Linear, c = Constant - equals;
  if(a == 0){
    if(b == 0) return INFINITY;
    return -c/b;
  }
  double disc = b*b - 4.0*a*c;
  if(disc < 0) disc = 0;
  double r1 = (-b - sqrt(disc)) / (2.0*a);
  double r2 = (-b + sqrt(disc)) / (2.0*a);
  return (r1 > r2) ? r1 : r2;
}

void NormalLossPiece::print(){
  Rprintf("%.10e %.10e %.10e %15f %15f %15f %d\n",
          Quadratic, Linear, Constant,
          min_mean, max_mean, prev_mean, data_i);
}

void PiecewiseNormalLoss::add
(double Quadratic, double Linear, double Constant){
  for(auto it = piece_list.begin(); it != piece_list.end(); it++){
    it->Quadratic += Quadratic;
    it->Linear += Linear;
    it->Constant += Constant;
  }
}

void PiecewiseNormalLoss::multiply(double x){
  for(auto it = piece_list.begin(); it != piece_list.end(); it++){
    it->Quadratic *= x;
    it->Linear *= x;
    it->Constant *= x;
  }
}

void PiecewiseNormalLoss::set_prev_seg_end(int prev_seg_end){
  for(auto it = piece_list.begin(); it != piece_list.end(); it++){
    it->data_i = prev_seg_end;
  }
}

void PiecewiseNormalLoss::findMean
(double mean, int *seg_end, double *prev_mean_out){
  for(auto it = piece_list.begin(); it != piece_list.end(); it++){
    if(it->min_mean <= mean && mean <= it->max_mean){
      *seg_end = it->data_i;
      *prev_mean_out = it->prev_mean;
      return;
    }
  }
  Rf_error("findMean: mean=%e outside all piece intervals", mean);
}

double PiecewiseNormalLoss::findCost(double mean){
  for(auto it = piece_list.begin(); it != piece_list.end(); it++){
    if(it->min_mean <= mean && mean <= it->max_mean){
      return it->getCost(mean);
    }
  }
  return INFINITY;
}

void PiecewiseNormalLoss::Minimize
(double *best_cost, double *best_mean,
 int *data_i, double *prev_mean_out){
  double candidate_cost, candidate_mean;
  *best_cost = INFINITY;
  for(auto it = piece_list.begin(); it != piece_list.end(); it++){
    candidate_mean = it->argmin_interval();
    candidate_cost = it->getCost(candidate_mean);
    if(candidate_cost < *best_cost){
      *best_cost = candidate_cost;
      *best_mean = candidate_mean;
      *data_i = it->data_i;
      *prev_mean_out = it->prev_mean;
    }
  }
}

void PiecewiseNormalLoss::print(){
  Rprintf("%10s %10s %15s %15s %15s %15s %s\n",
          "Quadratic", "Linear", "Constant",
          "min_mean", "max_mean", "prev_mean", "data_i");
  for(auto it = piece_list.begin(); it != piece_list.end(); it++){
    it->print();
  }
}

void PiecewiseNormalLoss::set_to_min_less_of
(PiecewiseNormalLoss *input){
  piece_list.clear();

  auto it = input->piece_list.begin();
  double right_bound = it->min_mean;
  double cur_value   = it->getCost(right_bound);

  piece_list.emplace_back(
    0.0, 0.0, cur_value,
    right_bound, right_bound,
    it->data_i, right_bound);

  bool const_piece;
  {
    double arg0 = it->argmin();
    bool not_const = !(it->Quadratic == 0 && it->Linear == 0);
    const_piece = (arg0 <= right_bound && not_const);
  }

  for(; it != input->piece_list.end(); ++it){
    double argmin_clamp = it->clamp_to_interval(it->argmin());
    double piece_min    = it->getCost(argmin_clamp);

    double decr_a = INFINITY, decr_b = INFINITY;
    bool piece_beats_cur = (piece_min + EPSILON < cur_value);
    bool argmin_ahead    = (right_bound < it->argmin());

    if(piece_beats_cur && argmin_ahead){
      if(const_piece){
        double root = it->get_smaller_root(cur_value);
        if(root < right_bound)   root = right_bound;
        if(root < it->min_mean)  root = it->min_mean;
        if(root > argmin_clamp)  root = right_bound;
        decr_a = root;
      } else {
        decr_a = right_bound;
      }
      decr_b = argmin_clamp;
    } else if(ABS(piece_min - cur_value) < EPSILON && argmin_ahead){
      decr_a = right_bound;
      decr_b = it->max_mean;
    }

    if(decr_a != INFINITY && decr_a < it->min_mean) decr_a = it->min_mean;
    bool has_decr = (decr_a != INFINITY && decr_a < decr_b);

    auto &last = piece_list.back();
    if(!has_decr){
      if(last.Quadratic == 0 && last.Linear == 0){
        last.max_mean = it->max_mean;
      } else {
        piece_list.emplace_back(
          0.0, 0.0, cur_value,
          right_bound, it->max_mean,
          last.data_i, right_bound);
      }
    } else {
      if(decr_a > last.min_mean){
        last.max_mean = decr_a;
      } else {
        piece_list.pop_back();
      }

      if(decr_a < decr_b){
        piece_list.emplace_back(
          it->Quadratic, it->Linear, it->Constant,
          decr_a, decr_b,
          it->data_i, INFINITY);
      }

      if(decr_b < it->max_mean){
        double out_value = it->getCost(decr_b);
        piece_list.emplace_back(
          0.0, 0.0, out_value,
          decr_b, it->max_mean,
          it->data_i, decr_b);
      }
    }

    auto &new_last = piece_list.back();
    right_bound = new_last.max_mean;
    cur_value   = new_last.getCost(right_bound);
    const_piece = (new_last.Quadratic == 0 && new_last.Linear == 0);
  }
}

bool sameFuns
(NormalLossPieceList::iterator it1,
 NormalLossPieceList::iterator it2){
  return it1->Quadratic == it2->Quadratic &&
    it1->Linear == it2->Linear &&
    ABS(it1->Constant - it2->Constant) < EPSILON;
}

void PiecewiseNormalLoss::push_piece
(NormalLossPieceList::iterator it,
 double min_mean, double max_mean){
  if(max_mean <= min_mean) return;
  NormalLossPieceList::iterator last = piece_list.end();
  --last;
  if(piece_list.size() &&
     sameFuns(last, it) &&
     it->prev_mean == last->prev_mean &&
     it->data_i == last->data_i){
    last->max_mean = max_mean;
  }else{
    piece_list.emplace_back(
      it->Quadratic, it->Linear, it->Constant,
      min_mean, max_mean,
      it->data_i, it->prev_mean);
  }
}

void PiecewiseNormalLoss::push_min_pieces
(PiecewiseNormalLoss *fun1,
 PiecewiseNormalLoss *fun2,
 NormalLossPieceList::iterator it1,
 NormalLossPieceList::iterator it2,
 int verbose){
  bool same_at_left;
  double last_min_mean;
  NormalLossPieceList::iterator prev2 = it2;
  prev2--;
  NormalLossPieceList::iterator prev1 = it1;
  prev1--;
  if(it1->min_mean < it2->min_mean){
    same_at_left = sameFuns(prev2, it1);
    last_min_mean = it2->min_mean;
  }else{
    last_min_mean = it1->min_mean;
    if(it2->min_mean < it1->min_mean){
      same_at_left = sameFuns(prev1, it2);
    }else{
      if(it1 == fun1->piece_list.begin() &&
         it2 == fun2->piece_list.begin()){
        same_at_left = false;
      }else{
        same_at_left = sameFuns(prev1, prev2);
      }
    }
  }
  NormalLossPieceList::iterator next2 = it2;
  next2++;
  NormalLossPieceList::iterator next1 = it1;
  next1++;
  bool same_at_right;
  double first_max_mean;
  if(it1->max_mean < it2->max_mean){
    same_at_right = sameFuns(next1, it2);
    first_max_mean = it1->max_mean;
  }else{
    first_max_mean = it2->max_mean;
    if(it2->max_mean < it1->max_mean){
      same_at_right = sameFuns(it1, next2);
    }else{
      if(next1 == fun1->piece_list.end() &&
         next2 == fun2->piece_list.end()){
        same_at_right = false;
      }else{
        same_at_right = sameFuns(next1, next2);
      }
    }
  }
  if(last_min_mean == first_max_mean) return;
  if(sameFuns(it1, it2)){
    push_piece(it1, last_min_mean, first_max_mean);
    return;
  }
  NormalLossPiece diff_piece(
    it1->Quadratic - it2->Quadratic,
    it1->Linear    - it2->Linear,
    it1->Constant  - it2->Constant,
    last_min_mean, first_max_mean, -5, 0.0);
  double mid_mean      = (first_max_mean + last_min_mean) / 2;
  double cost_diff_mid = diff_piece.getCost(mid_mean);
  if(same_at_left && same_at_right){
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    return;
  }
  if(diff_piece.Quadratic == 0){
    if(diff_piece.Linear == 0){
      if(diff_piece.Constant < 0){
        push_piece(it1, last_min_mean, first_max_mean);
      }else{
        push_piece(it2, last_min_mean, first_max_mean);
      }
      return;
    }
    double mean_at_equal_cost = -diff_piece.Constant / diff_piece.Linear;
    if(last_min_mean < mean_at_equal_cost &&
       mean_at_equal_cost < first_max_mean){
      if(0 < diff_piece.Linear){
        push_piece(it1, last_min_mean, mean_at_equal_cost);
        push_piece(it2, mean_at_equal_cost, first_max_mean);
      }else{
        push_piece(it2, last_min_mean, mean_at_equal_cost);
        push_piece(it1, mean_at_equal_cost, first_max_mean);
      }
      return;
    }
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    return;
  }
  double cost_diff_left  = diff_piece.getCost(last_min_mean);
  double cost_diff_right = diff_piece.getCost(first_max_mean);
  bool two_roots = diff_piece.has_two_roots(0.0);
  double smaller_mean, larger_mean;
  if(two_roots){
    smaller_mean = diff_piece.get_smaller_root(0.0);
    larger_mean  = diff_piece.get_larger_root(0.0);
  }
  if(same_at_right){
    if(two_roots){
      double mean_at_crossing = smaller_mean;
      double mean_at_optimum  = diff_piece.argmin();
      if(last_min_mean < mean_at_crossing &&
         mean_at_crossing < mean_at_optimum &&
         mean_at_optimum < first_max_mean){
        if(cost_diff_left < 0){
          push_piece(it1, last_min_mean, mean_at_crossing);
          push_piece(it2, mean_at_crossing, first_max_mean);
        }else{
          push_piece(it2, last_min_mean, mean_at_crossing);
          push_piece(it1, mean_at_crossing, first_max_mean);
        }
        return;
      }
    }
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    return;
  }
  if(same_at_left){
    if(two_roots){
      double mean_at_crossing = larger_mean;
      double mean_at_optimum  = diff_piece.argmin();
      if(last_min_mean < mean_at_optimum &&
         mean_at_optimum < mean_at_crossing &&
         mean_at_crossing < first_max_mean){
        if(cost_diff_right < 0){
          push_piece(it2, last_min_mean, mean_at_crossing);
          push_piece(it1, mean_at_crossing, first_max_mean);
        }else{
          push_piece(it1, last_min_mean, mean_at_crossing);
          push_piece(it2, mean_at_crossing, first_max_mean);
        }
        return;
      }
    }
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    return;
  }
  double first_mean_cross = INFINITY, second_mean_cross = INFINITY;
  if(two_roots){
    bool larger_inside  =
      last_min_mean < larger_mean  && larger_mean  < first_max_mean;
    bool smaller_inside =
      last_min_mean < smaller_mean && smaller_mean < first_max_mean;
    if(larger_inside){
      if(smaller_inside && smaller_mean < larger_mean){
        first_mean_cross  = smaller_mean;
        second_mean_cross = larger_mean;
      }else{
        first_mean_cross = larger_mean;
      }
    }else{
      if(smaller_inside){
        first_mean_cross = smaller_mean;
      }
    }
  }
  if(second_mean_cross != INFINITY){
    double before_mean      = (last_min_mean + first_mean_cross) / 2;
    double cost_diff_before = diff_piece.getCost(before_mean);
    if(cost_diff_before < 0){
      push_piece(it1, last_min_mean, first_mean_cross);
      push_piece(it2, first_mean_cross, second_mean_cross);
      push_piece(it1, second_mean_cross, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_mean_cross);
      push_piece(it1, first_mean_cross, second_mean_cross);
      push_piece(it2, second_mean_cross, first_max_mean);
    }
  }else if(first_mean_cross != INFINITY){
    double before_mean      = (last_min_mean + first_mean_cross) / 2;
    double cost_diff_before = diff_piece.getCost(before_mean);
    double after_mean       = (first_max_mean + first_mean_cross) / 2;
    double cost_diff_after  = diff_piece.getCost(after_mean);
    if(cost_diff_before < 0){
      if(cost_diff_after < 0){
        push_piece(it1, last_min_mean, first_max_mean);
      }else{
        push_piece(it1, last_min_mean, first_mean_cross);
        push_piece(it2, first_mean_cross, first_max_mean);
      }
    }else{
      if(cost_diff_after < 0){
        push_piece(it2, last_min_mean, first_mean_cross);
        push_piece(it1, first_mean_cross, first_max_mean);
      }else{
        push_piece(it2, last_min_mean, first_max_mean);
      }
    }
  }else{
    double cost_diff;
    if(first_max_mean == INFINITY){
      cost_diff = diff_piece.getCost(last_min_mean + 1);
    }else{
      if(ABS(cost_diff_mid) < EPSILON){
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

void PiecewiseNormalLoss::set_to_min_env_of
(PiecewiseNormalLoss *fun1, PiecewiseNormalLoss *fun2, int verbose){
  NormalLossPieceList::iterator
    it1 = fun1->piece_list.begin(),
    it2 = fun2->piece_list.begin();
  piece_list.clear();
  while(it1 != fun1->piece_list.end() &&
        it2 != fun2->piece_list.end()){
    push_min_pieces(fun1, fun2, it1, it2, verbose);
    double last_max_mean = piece_list.back().max_mean;
    if(it1->max_mean == last_max_mean) it1++;
    if(it2->max_mean == last_max_mean) it2++;
  }
}

int NormalFPOPisotonic
(double *data_vec, int data_count,
 double penalty,
 double *cost_vec,
 int *end_vec,
 double *mean_vec,
 int *intervals_vec){
  double min_mean = INFINITY, max_mean = -INFINITY;
  for(int data_i = 0; data_i < data_count; data_i++){
    if(data_vec[data_i] < min_mean) min_mean = data_vec[data_i];
    if(max_mean < data_vec[data_i]) max_mean = data_vec[data_i];
  }
  if(min_mean == max_mean){
    return ERROR_MIN_MAX_SAME;
  }

  std::vector<PiecewiseNormalLoss> cost_model_mat(data_count);
  PiecewiseNormalLoss *cost, *cost_prev;
  PiecewiseNormalLoss min_prev_cost;
  int verbose = 0;

  for(int data_i = 0; data_i < data_count; data_i++){
    cost = &cost_model_mat[data_i];
    double y = data_vec[data_i];

    if(data_i == 0){
      // C_1(m) = (m - y_1)^2 = m^2 - 2*y_1*m + y_1^2
      cost->piece_list.emplace_back(
        1.0, -2.0*y, y*y,
        min_mean, max_mean, -1, INFINITY);
    }else{
      // C_t(m) = (m - y_t)^2 + min{ C_{t-1}(m), C^{<=}_{t-1}(m) + lambda }
      // where C^{<=} enforces the non-decreasing (isotonic) constraint.
      min_prev_cost.set_to_min_less_of(cost_prev);
      min_prev_cost.set_prev_seg_end(data_i - 1);
      min_prev_cost.add(0.0, 0.0, penalty);
      cost->set_to_min_env_of(&min_prev_cost, cost_prev, verbose);
      cost->add(1.0, -2.0*y, y*y);
    }
    cost_prev = cost;
  }

  double best_cost, best_mean, prev_mean_val;
  int prev_seg_end;
  for(int i = 0; i < data_count; i++){
    cost = &cost_model_mat[i];
    intervals_vec[i] = cost->piece_list.size();
    cost->Minimize(cost_vec+i, &best_mean, &prev_seg_end, &prev_mean_val);
  }

  cost_model_mat[data_count-1].Minimize(
    &best_cost, &best_mean, &prev_seg_end, &prev_mean_val);

  for(int i = 0; i < data_count; i++){
    mean_vec[i] = INFINITY;
    end_vec[i] = -2;
  }
  mean_vec[0] = best_mean;
  end_vec[0] = prev_seg_end;
  int out_i = 1;
  while(0 <= prev_seg_end && out_i < data_count){
    cost = &cost_model_mat[prev_seg_end];
    if(prev_mean_val != INFINITY){
      best_mean = prev_mean_val;
    }
    cost->findMean(best_mean, &prev_seg_end, &prev_mean_val);
    mean_vec[out_i] = best_mean;
    end_vec[out_i] = prev_seg_end;
    out_i++;
  }
  return 0;
}
