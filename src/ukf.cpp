#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // initialization flag
  is_initialized_ = false;
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.0;//3; //std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.0; //std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // Dimensionality of the state vector: x, y, v, phi, phi_dot
  n_x_ = 5;
  
  // Dimensionality of the augmented state: x, y, v, phi, phi_dot, v_dot, phi_dot_dot
  n_aug_ = 7;
  
  // Lambda parameter for predicting sigma points
  lambda_ = 3.0 - double(n_aug_);
  
  // Consistency measures
  NIS_laser_ = 0.0;
  NIS_radar_ = 0.0;
  
  // Mean and covariance weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(1/(2*(lambda_ + n_aug_)));
  weights_(0) = lambda_/(lambda_ + n_aug_);
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  
  
  
  if (!is_initialized_){
    std::cout<<"Initialize\n";
    
    // initialize time
    time_us_ = meas_package.timestamp_;
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /*TODO*/
      // check for rho not too near 0.0
      if (meas_package.raw_measurements_(0) > 0.001){
        x_ = tools.Polar2Cartesian(meas_package.raw_measurements_);
      } else {
      // default values if first measurement is close to rho = 0.0
        x_(0) = 0.1;
        x_(1) = 0.1;
      }
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // check for (x,y) not too near (0.0, 0.0)
      if ((meas_package.raw_measurements_(0) > 0.001) && (meas_package.raw_measurements_(1) > 0.001)){
        x_(0) = meas_package.raw_measurements_(0);
        x_(1) = meas_package.raw_measurements_(1);
      } else {
        // default values if first measurement is close to (x, y) = (0.0, 0.0)
        x_(0) = 0.1;
        x_(1) = 0.1;
      }
    }
    
    // initialize P matrix, use large numbers for unknown initial x values
    P_ << 1, 0,    0,    0,    0,
          0, 1,    0,    0,    0,
          0, 0, 0.1,    0,    0,
          0, 0,    0, 0.1,    0,
          0, 0,    0,    0, 0.1;

    // done initializing, no need to predict or update
    is_initialized_ = true;
  } else {
    // update time information
    double dt = 0.000001*(time_us_-meas_package.timestamp_);
    time_us_ = meas_package.timestamp_;
    
    /*TODO*/
    Prediction(dt);
    
  }
  
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  // Generate and predict sigma points
  int n_sig = 2*n_aug_+1;
  MatrixXd Xsig_aug(n_aug_, n_sig);
  MatrixXd Xpred(n_x_, n_sig);
  
  tools.GenSigmaPts(Xsig_aug, x_, P_, std_a_, std_yawdd_, lambda_, n_x_, n_aug_);
  tools.PredSigmaPts(Xpred, Xsig_aug, delta_t, n_x_, n_aug_);
  
  // Predict mean and covariance
  tools.PredMean(x_, Xpred, weights_);
  tools.PredCovariance(P_, x_, Xpred, weights_);
  std::cout << x_ << "\n\n";
  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
