#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  // Copy from EKF project
  
  VectorXd rmse(4);
  rmse << 0,0,0,0;
    
  // From Lesson 21
  
  if (estimations.size() != ground_truth.size()){
    std::cout << "Array sizes must match!";
    return rmse;
  }
  else if (estimations.size() == 0){
    std::cout << "Arrays must be non-empty!";
    return rmse;
  }
  
  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    VectorXd tmp = estimations[i]-ground_truth[i];
    tmp = tmp.array()*tmp.array();
    rmse += tmp;
  }
  
  //calculate the mean
  rmse = rmse/estimations.size();
  
  //calculate the squared root
  rmse = rmse.array().sqrt();
  
  return rmse;
}

VectorXd Tools::Polar2Cartesian(const VectorXd& radar_meas) {
  /**
   * Convert polar coordinates to Cartesian.
   */
  VectorXd f_x(5,1);
  f_x.fill(0.0);
  
  // From Lesson 14
  double rho = radar_meas(0);
  double phi = radar_meas(1);
  
  f_x << rho*cos(phi),
        -rho*sin(phi),
         0,
         0,
         0;
  
  return f_x;
}

void Tools::GenSigmaPts(MatrixXd& Xsig_aug, const VectorXd& x, const MatrixXd& P,
                    double std_a, double std_yawdd, double lambda, int n_x, int n_aug) {
  /*
   * Augment state vector and covariance matrix, calculate
   *augmented sigma point matrix
   */
  // From lesson 14
  VectorXd x_aug(n_aug);
  MatrixXd P_aug(n_aug, n_aug);
  
  //create augmented mean state
  x_aug.fill(0.0);
  x_aug.head(n_x) = x;
  
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x,n_x) = P;
  P_aug(n_aug-2,n_aug-2) = std_a*std_a;
  P_aug(n_aug-1,n_aug-1) = std_yawdd*std_yawdd;
  
  //std::cout << "P_aug matrix: " << P_aug << "\n\n";
  
  MatrixXd A = P_aug.llt().matrixL();
  //std::cout << "Matrix A (P^-1):\n" << A << "\n\n";
  
  Xsig_aug.col(0) = x_aug;
  for (int idx = 0; idx < n_aug; ++idx){
      Xsig_aug.col(idx+1) = x_aug+sqrt(lambda+n_aug)*A.col(idx);
      Xsig_aug.col(idx+1+n_aug) = x_aug-sqrt(lambda+n_aug)*A.col(idx);
    //Xsig_aug(3,idx) = NormAngle(Xsig_aug(3,idx));
  }
}

void Tools::PredSigmaPts(MatrixXd& Xpred, const MatrixXd& Xsig_aug, double dt, int n_x, int n_aug) {
  /*
   * Use nonlinear f(x,v) to update sigma points
   */
  //From lesson 21
  
  for (int i = 0; i < (2*n_aug+1); ++i){
    // input values
    double px       = Xsig_aug(0,i); // x-position
    double py       = Xsig_aug(1,i); // y-position
    double v        = Xsig_aug(2,i); // velocity magnitude
    double yaw      = Xsig_aug(3,i); // velocity direction
    double yawd     = Xsig_aug(4,i); // rate of change of direction
    double nu_a     = Xsig_aug(5,i); // process noise, longitudinal acceleration
    double nu_yawdd = Xsig_aug(6,i); // process noise, lateral acceleration
    
    // predicted values
    double px_p, py_p, v_p, yaw_p, yawd_p;
    
    // check for very small yawd (straight line driving)
    if (fabs(yawd) < 0.001){
      // straight line driving, use geometry
      px_p = px + v*dt*cos(yaw);
      py_p = py + v*dt*sin(yaw);
    } else {
      // turning, use calculus
      px_p = px + v/yawd * (sin(yaw+yawd*dt)-sin(yaw));
      py_p = py + v/yawd * (-cos(yaw+yawd*dt)+cos(yaw));
    }
    
    v_p    = v; // constant longitudinal velocity
    yaw_p  = yaw + yawd*dt; // integrate rate of change
    //yaw_p = NormAngle(yaw_p);
    yawd_p = yawd; // constant rate of change of direction
    
    // noise effect
    px_p   += 0.5*nu_a*dt*dt*cos(yaw);
    py_p   += 0.5*nu_a*dt*dt*sin(yaw);
    v_p    += nu_a*dt;
    yaw_p  += 0.5*nu_yawdd*dt*dt;
    //yaw_p = NormAngle(yaw_p);
    yawd_p += nu_yawdd*dt;
    
    // place temp variables back into matrix
    Xpred(0,i) = px_p;
    Xpred(1,i) = py_p;
    Xpred(2,i) = v_p;
    Xpred(3,i) = yaw_p;
    Xpred(4,i) = yawd_p;
    
  }
}

void Tools::PredMean(VectorXd& x, const MatrixXd& Xpred, VectorXd& w) {
  int n_sig = Xpred.cols();
  
  // clear previous value in x
  x.fill(0.0);
  
  // accumulate mean
  for (int i = 0; i < n_sig; ++i) {
    x = x+ w(i)*Xpred.col(i); //weighted column i
  }
}

void Tools::PredCovariance(MatrixXd& P, const VectorXd& x, const MatrixXd& Xpred, const VectorXd& w, const int angleIdx) {
  int n_sig = Xpred.cols();
  VectorXd col;
  // clear previous value in P
  P.fill(0.0);
  
  // accumulate covariance
  for (int i = 0; i < n_sig; ++i){
    col = Xpred.col(i)-x;
    //std::cout << "col pre - " << i << ":\n" << col << "\n\n";
    //normalize angular differences between -pi and pi
    while (col(angleIdx)> M_PI) col(angleIdx)-=2.*M_PI;
    while (col(angleIdx)< -M_PI) col(angleIdx)+=2.*M_PI;
    //if (angleIdx >= 0) {
    //    col(angleIdx) = NormAngle(col(angleIdx));
    //}
    //std::cout << "col post - " << i << ":\n" << col << "\n\n";
    MatrixXd tmp(5,5);
    tmp = w(i)*col*col.transpose();
    //std::cout << "tmp matrix - " << i << ":\n" << tmp << "\n\n";
    P = P + w(i)*col*col.transpose();
  }

}

double Tools::NormAngle(double phi) {
  // normalize phi between -pi and pi
  if ((phi >= -M_PI)&&(phi <= M_PI)){
    // already in range
    return phi;
  }
  else if (phi < 0) {
    // phi too negative
    double tmp;
    tmp = phi - M_PI;
    tmp = -remainder(-tmp,(2*M_PI))+M_PI;
    //std::cout<<"phi raw = "<<phi<<"\n";
    //std::cout<<"phi = "<<phi/M_PI<<"\ntmp = "<<tmp/M_PI<<"\n\n";
    return tmp;
  }
      //TODO
  else {
    // phi too positive
    double tmp;
    tmp = phi + M_PI;
    tmp = remainder(tmp,(2*M_PI));
    tmp -= M_PI;
    //std::cout<<"phi raw = "<<phi<<"\n";
    //std::cout<<"phi = "<<phi/M_PI<<"\ntmp = "<<tmp/M_PI<<"\n\n";
    return tmp;
  }
}

void Tools::Cartesian2Polar(const MatrixXd& X, MatrixXd& Z, const int n) {
  // From lesson 26
  for (int i = 0; i < n; ++i) {
    float rho, phi;
    rho = sqrt(X(0,i)*X(0,i)+X(1,i)*X(1,i));
    if (rho < 0.001) rho = 0.001;
    phi = atan2(X(1,i), X(0,i));
    
    float vx, vy;
    vx = cos(X(3,i))*X(2,i);
    vy = sin(X(3,i))*X(2,i);
    
    float rho_dot = (X(0,i)*vx + X(1,i)*vy)/rho;
    Z.col(i) << rho, phi, rho_dot;
  }
}

