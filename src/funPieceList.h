/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <list>

// NOTE: please only define prototypes in this file (do not define
// methods directly -- instead define them in funPieceList.cpp). This
// is more compatible with R's Makefile, which automatically
// recompiles object files when there are changes to *.cpp but not *.h
// files.

class PoissonLossPiece {
 public:
  double Linear;
  double Log;
  double Constant;
  double min_mean;
  double max_mean;
  int data_i;
  bool equality_constraint_active;
  PoissonLossPiece
    (double li, double lo, double co, double m, double M, int i, bool a);
  double getDiscriminant(double add_constant);
  double discriminant2mean_principal(double discriminant);
  double discriminant2mean_secondary(double discriminant);
  double discriminant2mean_larger(double discriminant);
  double discriminant2mean_smaller(double discriminant);
  double getMinMean();
  double PoissonLoss(double mean);
  double PoissonDeriv(double);
};

typedef std::list<PoissonLossPiece> PoissonLossPieceList;

class PiecewisePoissonLoss {
 public:
  PoissonLossPieceList piece_list;
  void set_to_min_less_of(PiecewisePoissonLoss *);
  void set_to_min_more_of(PiecewisePoissonLoss *);
  void set_to_min_env_of(PiecewisePoissonLoss *, PiecewisePoissonLoss *);
  void push_min_pieces
    (PiecewisePoissonLoss *, PiecewisePoissonLoss *,
     PoissonLossPieceList::iterator, PoissonLossPieceList::iterator);
  void push_piece(PoissonLossPieceList::iterator, double, double);
  void add(double Linear, double Log, double Constant);
  void print();
  void set_prev_seg_end(int prev_seg_end);
  void findMean(double mean, int *seg_end, bool *equality_constraint_active);
  void Minimize(double *best_cost,
		double *best_mean,
		int *data_i,
		bool *equality_constraint_active);
};

bool sameFuns(PoissonLossPieceList::iterator, PoissonLossPieceList::iterator);
