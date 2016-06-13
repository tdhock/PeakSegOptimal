/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <list>
#include <math.h>

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
  void addPiece(PoissonLossPiece new_piece){
    PoissonLossPieceList::iterator it;
    for(it=piece_list.begin(); it != piece_list.end(); it++){
      PoissonLossPiece it_piece = *it;
      it_piece.Linear += new_piece.Linear;
      it_piece.Log += new_piece.Log;
      it_piece.Constant += new_piece.Constant;
    }
  }
  void findMean(double mean, int *seg_end, bool *equality_constraint_active){
    PoissonLossPieceList::iterator it;
    PoissonLossPiece piece;
    for(it=piece_list.begin(); it != piece_list.end(); it++){
      piece = *it;
      if(piece.min_mean <= mean && mean <= piece.max_mean){
	*seg_end = piece.data_i;
	*equality_constraint_active = piece.equality_constraint_active;
      }
    }
  }    
  void Minimize(double *best_cost,
		double *best_mean,
		int *data_i,
		bool *equality_constraint_active){
    double candidate_cost, candidate_mean;
    PoissonLossPieceList::iterator it;
    PoissonLossPiece piece;
    *best_cost = INFINITY;
    for(it=piece_list.begin(); it != piece_list.end(); it++){
      piece = *it;
      candidate_mean = piece.getMinMean();
      if(piece.min_mean <= candidate_mean && candidate_mean <= piece.max_mean){
	candidate_cost = piece.PoissonLoss(candidate_mean);
	if(candidate_cost < *best_cost){
	  *best_cost = candidate_cost;
	  *best_mean = candidate_mean;
	  *data_i = piece.data_i;
	  *equality_constraint_active = piece.equality_constraint_active;
	}
      }
    }
  }    
};
