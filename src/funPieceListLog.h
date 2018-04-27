/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <list>

/* 
* Implement the square loss based on PoissonLoss template
* 
* 
* 
* @author Sean Jewell 27 April 2017
*/

class SquareLossPiece {
public:
  // coefs for square loss
  // Square * u ^ 2 + Linear * u + Constant * 1
  double Square;
  double Linear;
  double Constant;
  // base 10 scale
  double min_mean;
  double max_mean;
  int data_i;
  double prev_mean;
  SquareLossPiece();
  SquareLossPiece
    (double a, double b, double c, double m, double M, int i, double);
  double argmin();
  double argmin_mean();
  void print();
  double get_smaller_root(double);
  double get_larger_root(double);
  bool has_two_roots(double);
  double getCost(double mean);
  double SquareLoss(double);
};

typedef std::list<SquareLossPiece> SquareLossPieceList;

class PiecewiseSquareLoss {
public:
  SquareLossPieceList piece_list;
  void set_to_min_less_of(PiecewiseSquareLoss *, double, int);
  void set_to_unconstrained_min_of(PiecewiseSquareLoss *, double, int);
  void set_to_scaled_of(PiecewiseSquareLoss *, double, double, int);
  void set_to_eps_min_of(PiecewiseSquareLoss *, double, int);
  void set_to_clean(PiecewiseSquareLoss *, double, int);
  void set_to_min_env_of
    (PiecewiseSquareLoss *, PiecewiseSquareLoss *, double, int);
  int check_min_of(PiecewiseSquareLoss *, PiecewiseSquareLoss *);
  void push_min_pieces
    (PiecewiseSquareLoss *, PiecewiseSquareLoss *,
     SquareLossPieceList::iterator, SquareLossPieceList::iterator, double, int);
  void push_piece(SquareLossPieceList::iterator, double, double);
  void add(double Square, double Linear, double Constant);
  void add_penalty(double penalty, double EPS);
  void add_protect(double Square, double Linear, double Constant, double EPS, bool eps_segment);
  void multiply(double);
  void scale(double);
  void print();
  void checkStable(double);
  void set_prev_seg_end(int prev_seg_end, double EPS);
  void findMean(double mean, int *seg_end, double *prev_mean);
  double findCost(double mean);
  void Minimize
    (double *best_cost,
     double *best_mean,
     int *data_i,
     double *prev_mean);
};

bool sameFunsSquare(SquareLossPieceList::iterator, SquareLossPieceList::iterator);


