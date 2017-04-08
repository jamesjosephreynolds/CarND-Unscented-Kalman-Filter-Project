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
  bool is_initialized_ = false;
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3; //std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3; //std_yawdd_ = 30;

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
    time_us_ = measurement_pack.timestamp_;
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /*TODO*/
      // check for rho not too near 0.0
      if (measurement_pack.raw_measurements_(0) > 0.001){
        x_ = tools.Polar2Cartesian(measurement_pack.raw_measurements_);
      } else {
      // default values if first measurement is close to rho = 0.0
        x_(0) = 0.1;
        x_(1) = 0.1;
      }
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // check for (x,y) not too near (0.0, 0.0)
      if ((measurement_pack.raw_measurements_(0) > 0.001) && (measurement_pack.raw_measurements_(1) > 0.001)){
        x_(0) = measurement_pack.raw_measurements_(0);
        x_(1) = measurement_pack.raw_measurements_(1);
      } else {
        // default values if first measurement is close to (x, y) = (0.0, 0.0)
        x_(0) = 0.1;
        x_(1) = 0.1;
      }
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
  } else {
    // update time information
    float dt = time_us-measurement_pack.timestamp_;
    time_us_ = measurement_pack.timestamp_;
    
    /*TODO*/
    
  }
  
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
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
