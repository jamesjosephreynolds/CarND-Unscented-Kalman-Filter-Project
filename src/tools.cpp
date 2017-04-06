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
  float rho = radar_meas(0);
  float phi = radar_meas(1);
  
  f_x << rho*cos(phi),
        -rho*sin(phi),
         0,
         0,
         0;
  
  return f_x;
}
