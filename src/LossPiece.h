#ifndef LOSS_PIECE_H
#define LOSS_PIECE_H

class LossPiece {
public:
  double min_mean, max_mean;
  int data_i;
  double prev_mean;

  LossPiece() : min_mean(0), max_mean(0), data_i(0), prev_mean(0) {}
  LossPiece(double m, double M, int i, double prev)
    : min_mean(m), max_mean(M), data_i(i), prev_mean(prev) {}
  double clamp_to_interval(double x) const {
    if(x < min_mean) return min_mean;
    if(x > max_mean) return max_mean;
    return x;
  }
};

#endif
