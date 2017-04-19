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
  std_yawdd_ = 0.5; //std_yawdd_ = 30;

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
  NIS_mean_ = 0.0;
  n_NIS_ = 0;
  
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
    std::cout<<"Initialize\n";
    
    // initialize time
    time_us_ = meas_package.timestamp_;
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /*TODO*/
      // check for rho not too near 0.0
      if (meas_package.raw_measurements_(0) > 0.001){
        std::cout << "Init measurement radar\n\n";
        x_ = tools.Polar2Cartesian(meas_package.raw_measurements_);
      } else {
      // default values if first measurement is close to rho = 0.0
        std::cout << "Init measurement lidar\n\n";
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
          0, 0, 1,    0,    0,
          0, 0,    0, 1,    0,
          0, 0,    0,    0, 1;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    
    std::cout << "Init x vector: \n" << x_ << "\n\n";
    
  } else {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      if (use_radar_) {
        // update time information
        double dt = 0.000001*(time_us_-meas_package.timestamp_);
        time_us_ = meas_package.timestamp_;
    
        // predict and update
        Prediction(dt);
        //std::cout << "Predicted x vector: \n" << x_ << "\n\n";
        //std::cout << "Radar measurement \n" << meas_package.raw_measurements_<< "\n\n";
        //std::cout << "Radar update x vector: \n" << x_ << "\n\n";
        UpdateRadar(meas_package);
      }
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      if (use_laser_) {
        // update time information
        double dt = 0.000001*(time_us_-meas_package.timestamp_);
        time_us_ = meas_package.timestamp_;
        
        // predict and update
        Prediction(dt);
        //std::cout << "Predicted x vector: \n" << x_ << "\n\n";
        //std::cout << "Lidar measurement: \n" << meas_package.raw_measurements_<< "\n\n";
        //std::cout << "Lidar update x vector: \n" << x_ << "\n\n";
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
  // Generate and predict sigma points
  MatrixXd Xsig_aug(n_aug_, n_sig_);
  //MatrixXd Xpred(n_x_, n_sig_);
  Xsig_aug.fill(0.0);
  Xpred_.fill(0.0);
  
  tools.GenSigmaPts(Xsig_aug, x_, P_, std_a_, std_yawdd_, lambda_, n_x_, n_aug_);
  //std::cout << "Sigma points: \n" << Xsig_aug << "\n\n";
  tools.PredSigmaPts(Xpred_, Xsig_aug, delta_t, n_x_, n_aug_);
  //std::cout << "Predict sigma points: \n" << Xpred_ << "\n\n";
  
  // Predict mean and covariance
  tools.PredMean(x_, Xpred_, weights_);
  //std::cout << "Predicted x mean: \n" << x_ << "\n\n";
  tools.PredCovariance(P_, x_, Xpred_, weights_, CTRVPHIIDX);
  //std::cout << "Predicted P mean: \n" << P_ << "\n\n";
  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  
  MatrixXd H(n_lidz_, n_x_);
  MatrixXd R(n_lidz_, n_lidz_);
  VectorXd z(n_lidz_);
  
  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;
  
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  
  z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);
  
  VectorXd z_pred = H * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd K = P_ * Ht * S.inverse();
  
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;
  
  /*MatrixXd Zpred(n_lidz_, n_sig_);
  MatrixXd S(n_lidz_, n_lidz_);
  MatrixXd R(n_lidz_, n_lidz_);
  VectorXd z(n_lidz_);
  z = meas_package.raw_measurements_;
  
  S.fill(0.0);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  
  
  Zpred = Xpred_.block(0, 0, n_lidz_, n_sig_);
  //std::cout << Xpred_ << "\n\n";
  //std::cout << Zpred << "\n\n";
  
  VectorXd z_pred(n_lidz_);
  tools.PredMean(z_pred, Zpred, weights_);
  tools.PredCovariance(S, z_pred, Zpred, weights_, LIDPHIIDX);
  S += R;
  
  // From lesson 29
  MatrixXd Tc = MatrixXd(n_x_, n_lidz_);
  for (int i = 0; i < n_sig_; ++i){
    VectorXd col = Xpred_.col(i)-x_;
    VectorXd row = Zpred.col(i)-z_pred;
    Tc = Tc + weights_(i)*col*row.transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K(n_x_,n_lidz_);
  K = Tc*S.inverse();
  //update state mean and covariance matrix
  x_ += K*(z-z_pred);
  P_ -= K*S*K.transpose();
  std::cout << "Updated P mean: \n" << P_ << "\n\n";
   */
  
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
  Zpred.fill(0.0);
  z = meas_package.raw_measurements_;
  
  S.fill(0.0);
  R << std_radr_*std_radr_, 0 ,0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
  
  
  tools.Cartesian2Polar(Xpred_, Zpred, n_sig_);
  
  VectorXd z_pred(n_radz_);
  tools.PredMean(z_pred, Zpred, weights_);
  tools.PredCovariance(S, z_pred, Zpred, weights_, RADPHIIDX);
  S += R;
  
  // From lesson 29
  MatrixXd Tc = MatrixXd(n_x_, n_radz_);
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; ++i){
    VectorXd col = Xpred_.col(i)-x_;
    col(CTRVPHIIDX) = tools.NormAngle(col(CTRVPHIIDX));
    VectorXd row = Zpred.col(i)-z_pred;
    row(RADPHIIDX) = tools.NormAngle(row(RADPHIIDX));
    Tc = Tc + weights_(i)*col*row.transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K(n_x_,n_radz_);
  K = Tc*S.inverse();
  //update state mean and covariance matrix
  x_ += K*(z-z_pred);
  P_ -= K*S*K.transpose();
  //std::cout << "Updated P mean: \n" << P_ << "\n\n";
  
  //calculate NIS
  MatrixXd zkp1(n_radz_, 1); //z at time k + 1
  tools.Cartesian2Polar(x_, zkp1, 1);
  VectorXd delta(n_radz_);
  z = zkp1.col(0);
  delta = z - z_pred;
  NIS_radar_ = delta.transpose()*S.inverse()*delta;
  n_NIS_ += 1;
  if (n_NIS_ > 1) {
    NIS_mean_ = (NIS_radar_ + (n_NIS_ - 1)*NIS_mean_)/n_NIS_;
  } else {
    NIS_mean_ = NIS_radar_;
  }
  std::cout << "NIS: " << NIS_radar_ << ", mean: " << NIS_mean_ << "\n\n";
  
}
