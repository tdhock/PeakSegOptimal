/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <list>

// NOTE: please only define prototypes in this file (do not define
// methods directly -- instead define them in funPieceList.cpp). This
// is more compatible with R's Makefile, which automatically
// recompiles object files when there are changes to *.cpp but not *.h
// files.

class PoissonLossPieceLog {
public:
  double Linear;
  double Log;
  double Constant;
  double min_log_mean;
  double max_log_mean;
  int data_i;
  double prev_log_mean;
  PoissonLossPieceLog();
  PoissonLossPieceLog
    (double li, double lo, double co, double m, double M, int i, double);
  double argmin();
  double argmin_mean();
  void print();
  double get_smaller_root(double);
  double get_larger_root(double);
  bool has_two_roots(double);
  double getCost(double mean);
  double getDeriv(double);
  double PoissonLoss(double);
  double PoissonDeriv(double);
};

typedef std::list<PoissonLossPieceLog> PoissonLossPieceListLog;

class PiecewisePoissonLossLog {
public:
  PoissonLossPieceListLog piece_list;
  void set_to_min_less_of(PiecewisePoissonLossLog *, int);
  void set_to_min_more_of(PiecewisePoissonLossLog *, int);
  void set_to_min_env_of
    (PiecewisePoissonLossLog *, PiecewisePoissonLossLog *, int);
  int check_min_of(PiecewisePoissonLossLog *, PiecewisePoissonLossLog *);
  void push_min_pieces
    (PiecewisePoissonLossLog *, PiecewisePoissonLossLog *,
     PoissonLossPieceListLog::iterator, PoissonLossPieceListLog::iterator, int);
  void push_piece(PoissonLossPieceListLog::iterator, double, double);
  void add(double Linear, double Log, double Constant);
  void multiply(double);
  void print();
  void set_prev_seg_end(int prev_seg_end);
  void findMean(double mean, int *seg_end, double *prev_log_mean);
  double findCost(double mean);
  void Minimize
    (double *best_cost,
     double *best_mean,
     int *data_i,
     double *prev_log_mean);
};

bool sameFuns(PoissonLossPieceListLog::iterator, PoissonLossPieceListLog::iterator);

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
  void set_to_min_less_of(PiecewiseSquareLoss *, int);
  void set_to_unconstrained_min_of(PiecewiseSquareLoss *, int);
  void set_to_scaled_of(PiecewiseSquareLoss *, double, double, int);
  void set_to_min_env_of
    (PiecewiseSquareLoss *, PiecewiseSquareLoss *, int);
  int check_min_of(PiecewiseSquareLoss *, PiecewiseSquareLoss *);
  void push_min_pieces
    (PiecewiseSquareLoss *, PiecewiseSquareLoss *,
     SquareLossPieceList::iterator, SquareLossPieceList::iterator, int);
  void push_piece(SquareLossPieceList::iterator, double, double);
  void add(double Square, double Linear, double Constant);
  void multiply(double);
  void scale(double);
  void print();
  void checkStable(double);
  void set_prev_seg_end(int prev_seg_end);
  void findMean(double mean, int *seg_end, double *prev_mean);
  double findCost(double mean);
  void Minimize
    (double *best_cost,
     double *best_mean,
     int *data_i,
     double *prev_mean);
};

bool sameFunsSquare(SquareLossPieceList::iterator, SquareLossPieceList::iterator);


