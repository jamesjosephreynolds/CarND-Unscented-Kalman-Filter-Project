#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

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
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

  /**
  * A coordinate system transformation.
  */
  Eigen::VectorXd Polar2Cartesian(const Eigen::VectorXd& radar_meas);
  
  /**
   * Generate sigma points from the current state vector
   */
  void GenSigmaPts(Eigen::MatrixXd& Xsig_aug, const Eigen::VectorXd& x, const Eigen::MatrixXd& P, double std_a,
                             double std_yawdd, double lambda, int n_x, int n_aug);
  
  /**
   * Predict sigma points using f(x,v)
   */
  void PredSigmaPts(Eigen::MatrixXd& Xpred, const Eigen::MatrixXd& Xsig_aug, double dt, int n_x, int n_aug);
  
  /*
   * Predict mean state vector x
   */
  void PredMean(Eigen::VectorXd& x, const Eigen::MatrixXd& Xpred, Eigen::VectorXd& w);
  
  /*
   * Predict covariance matrix P
   */
  void PredCovariance(Eigen::MatrixXd& P, const Eigen::VectorXd& x, const Eigen::MatrixXd& Xpred, const Eigen::VectorXd& w);
  
};

#endif /* TOOLS_H_ */
