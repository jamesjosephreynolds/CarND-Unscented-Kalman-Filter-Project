#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A coordinate system transformation.
  */
  VectorXd Polar2Cartesian(const VectorXd& radar_meas);
  
  /**
   * Generate sigma points from the current state vector
   */
  void GenSigmaPts(MatrixXd& Xsig_aug, const VectorXd& x, const MatrixXd& P, double std_a,
                             double std_yawdd, double lambda, int n_x, int n_aug);
  
  /**
   * Predict sigma points using f(x,v)
   */
  void PredSigmaPts(MatrixXd& Xpred, const MatrixXd& Xsig_aug, double dt, int n_x, int n_aug);
  
  /*
   * Predict mean state vector x
   */
  void PredMean(VectorXd& x, const MatrixXd& Xpred, VectorXd& w);
  
  /*
   * Predict covariance matrix P
   */
  void PredCovariance(MatrixXd& P, const VectorXd& x, const MatrixXd& Xpred, const VectorXd& w, const int angleIdx);
  
  /*
   * Normalize an angle phi between -pi and pi
   */
  double NormAngle(double phi);
  
  /*
   * Convert a cartesian matrix into a polar matrix
   */
  
  void Cartesian2Polar(const MatrixXd& X, MatrixXd& Z, const int n);
  
};

#endif /* TOOLS_H_ */
