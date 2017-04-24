#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

#define CTRVPHIIDX int(3)
#define RADPHIIDX int(1)
#define LIDPHIIDX int(-1)

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
  std_a_ = 3;//3; //std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1; //std_yawdd_ = 30;

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
  
  // Dimensionality of the measurement vector for radar: rho, v, phi
  n_radz_ = 3;
  
  // Dimensionality of the measurement vector for lidar: x, y
  n_lidz_ = 2;
  
  // Dimensionality of the augmented state: x, y, v, phi, phi_dot, v_dot, phi_dot_dot
  n_aug_ = 7;
  
  // Number of sigma points
  n_sig_ = 2*n_aug_+1;
  
  
  // predicted sigma points
  Xpred_ = MatrixXd(n_x_, n_sig_);
  
  // Lambda parameter for predicting sigma points
  lambda_ = 3.0 - double(n_aug_);
  
  // Consistency measures
  NIS_laser_ = 0.0;
  NIS_radar_ = 0.0;
  
  // Mean and covariance weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(0.5/(lambda_ + n_aug_));
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
    //std::cout<<"Initialize\n";
    
    // initialize time and CTRV state
    time_us_ = meas_package.timestamp_;
    x_.fill(0.0);
    
    // populate CTRV state initial value with measurement data
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
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
      if ((fabs(meas_package.raw_measurements_(0) > 0.1)) || fabs((meas_package.raw_measurements_(1)) > 0.1)) {
        x_(0) = meas_package.raw_measurements_(0);
        x_(1) = meas_package.raw_measurements_(1);
      } else {
        // default values if first measurement is close to (x, y) = (0.0, 0.0)
        x_(0) = 0.1;
        x_(1) = 0.1;
      }
    }
    
    // initialize P matrix, use identity
    P_ = MatrixXd::Identity(n_x_, n_x_);

    // done initializing, no need to predict or update
    is_initialized_ = true;
    
  } else {
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      if (use_radar_) {
        // update time information only if radar is enabled
        double dt = 0.000001*(meas_package.timestamp_-time_us_);
        time_us_ = meas_package.timestamp_;
        
        // handle large delta t to keep predictions numerically stable
        while (dt > 0.1)
        {
          const double dt_small = 0.05;
          Prediction(dt_small);
          dt -= dt_small;
        }
    
        // predict and update
        Prediction(dt);
        UpdateRadar(meas_package);
      }
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      if (use_laser_) {
        // update time information only if lidar is enabled
        double dt = 0.000001*(meas_package.timestamp_-time_us_);
        time_us_ = meas_package.timestamp_;
        
        // handle large delta t to keep predictions numerically stable
        while (dt > 0.1)
        {
          const double dt_small = 0.05;
          Prediction(dt_small);
          dt -= dt_small;
        }
        
        // predict and update
        Prediction(dt);
        UpdateLidar(meas_package);
      }
    }
    
    
  }
  
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  
  MatrixXd Xsig_aug(n_aug_, n_sig_);
  
  // Clear previous sigma points and prediction
  Xsig_aug.fill(0.0);
  Xpred_.fill(0.0);
  
  // Generate and predict sigma points
  tools.GenSigmaPts(Xsig_aug, x_, P_, std_a_, std_yawdd_, lambda_, n_x_, n_aug_);
  tools.PredSigmaPts(Xpred_, Xsig_aug, delta_t, n_x_, n_aug_);
  
  // Predict mean and covariance
  tools.PredMean(x_, Xpred_, weights_);
  tools.PredCovariance(P_, x_, Xpred_, weights_, CTRVPHIIDX);
  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  
  // Use implementation from EKF (linear function)
  MatrixXd H(n_lidz_, n_x_);
  MatrixXd R(n_lidz_, n_lidz_);
  VectorXd z(n_lidz_);
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  
  // Linear x -> z transformation matrix
  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;
  
  // Measurement noise covariance matrix
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  
  // Get measurement data
  z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);
  
  // Measurement prediction
  VectorXd z_pred = H * x_;
  
  // Prediction error
  VectorXd y = z - z_pred;
  
  // Update state vector andc covariance matrix
  MatrixXd PHt = P_ * H.transpose();
  MatrixXd S = H * PHt + R;
  MatrixXd K = PHt * S.inverse();
  x_ = x_ + (K * y);
  P_ = (I - K * H) * P_;
  
  // Current NIS value
  NIS_laser_ = y.transpose()*S.inverse()*y;
  
  //std::cout << "NIS: " << NIS_lidar_ << ", mean: " << NIS_mean_ << "\n\n";
    
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  MatrixXd Zpred(n_radz_, n_sig_);
  MatrixXd S(n_radz_, n_radz_);
  MatrixXd R(n_radz_, n_radz_);
  VectorXd z(n_radz_);
  VectorXd z_pred(n_radz_);
  MatrixXd Tc = MatrixXd(n_x_, n_radz_);
  VectorXd x_err(n_x_);
  VectorXd z_err(n_radz_);
  MatrixXd K(n_x_,n_radz_);
  Tc.fill(0.0);
  Zpred.fill(0.0);
  S.fill(0.0);
  
  // Get measurement data
  z << meas_package.raw_measurements_(0),
       meas_package.raw_measurements_(1),
       meas_package.raw_measurements_(2);
  
  // Measurement noise covariance matrix
  R << std_radr_*std_radr_, 0 ,0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
  
  // Transform predicted sigma points into measurement space Xpred -> Zpred
  tools.Cartesian2Polar(Xpred_, Zpred, n_sig_);
  
  // Predict mean and covariance
  tools.PredMean(z_pred, Zpred, weights_);
  tools.PredCovariance(S, z_pred, Zpred, weights_, RADPHIIDX);
  S += R;
  
  // From lesson 29
  for (int i = 0; i < n_sig_; ++i) {
    // distance between predicted sigma point and mean
    x_err = Xpred_.col(i)-x_;
    z_err = Zpred.col(i)-z_pred;
    
    //angle normalization (-pi, pi) in x and z domains
    while (x_err(CTRVPHIIDX)> M_PI) x_err(CTRVPHIIDX) -= 2.*M_PI;
    while (x_err(CTRVPHIIDX)<-M_PI) x_err(CTRVPHIIDX) += 2.*M_PI;
    while (z_err(RADPHIIDX)> M_PI) z_err(RADPHIIDX) -= 2.*M_PI;
    while (z_err(RADPHIIDX)<-M_PI) z_err(RADPHIIDX) += 2.*M_PI;
    
    Tc = Tc + weights_(i)*x_err*z_err.transpose();
  }
  
  //calculate Kalman gain K;
  K = Tc*S.inverse();
  
  // Update state mean and covariance matrix (reuse z_err)
  z_err = z - z_pred;
  
  //angle normalization (-pi, pi) in z domain
  while (z_err(RADPHIIDX) >  M_PI) z_err(RADPHIIDX) -= 2.*M_PI;
  while (z_err(RADPHIIDX) < -M_PI) z_err(RADPHIIDX) += 2.*M_PI;
  
  x_ += K*(z - z_pred);
  P_ -= K*S*K.transpose();
  
  // Current NIS value
  NIS_radar_ = z_err.transpose()*S.inverse()*z_err;
  
  //std::cout << "NIS: " << NIS_radar_ << ", mean: " << NIS_mean_ << "\n\n";
  
}
