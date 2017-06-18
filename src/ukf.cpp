#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

double epsilon = 1e-5;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */


  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  

  
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
  // Prediction step

  // Update step
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

  // 0. Generate sigma points - augmented

  //create augmented mean vector
  VectorXd X_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  MatrixXd Q(2, 2);
  Q << std_a_ * std_a_, 0,
       0, std_yawdd_ * std_yawdd_;

  P_aug.fill(0);

  //create augmented mean state
  X_aug.head(5) = x_;
  X_aug(5) = 0;
  X_aug(6) = 0;

  //create augmented covariance matrix

  P_aug.topLeftCorner( 5, 5 ) = P_;
  P_aug.bottomRightCorner( 2, 2 ) = Q;

  // std::cout << "P_aug = " << std::endl << P_aug << std::endl;

  //create square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();

  double sql_na = sqrt(lambda_ + n_aug_);

  //create augmented sigma points

  Xsig_aug.col(0) = X_aug;
  for (int i = 0; i < 7; i++) {
      Xsig_aug.col(i+1)           = X_aug + sql_na * A_aug.col(i);
      Xsig_aug.col(i + 1 + n_aug_) = X_aug - sql_na * A_aug.col(i);
  }


  // 1. Predict sigma points

  for (int i = 0; i < 2*n_aug_+1; i++) {

      VectorXd col_i = Xsig_aug.col(i);

      double vk = col_i(2);
      double phik = col_i(3);
      double phidotk = col_i(4);
      double n_ak = col_i(5);
      double n_phi = col_i(6);

      VectorXd offset = VectorXd(5);
      VectorXd process_noise = VectorXd(5);

      if (phidotk == 0) {
          phidotk = epsilon;
      }

      double v_phidot = vk / phidotk;

      offset << v_phidot * ( sin(phik + phidotk * delta_t) - sin(phik) ),
                v_phidot * ( -cos(phik + phidotk * delta_t) + cos(phik) ),
                0,
                phidotk * delta_t,
                0;

      process_noise << delta_t*delta_t*0.5 * cos(phik) * n_ak,
                       delta_t*delta_t*0.5 * sin(phik) * n_ak,
                       n_ak * delta_t,
                       delta_t*delta_t*0.5 * n_phi,
                       delta_t * n_phi;


      Xsig_pred_.col(i) = col_i.head(5) + offset + process_noise;
  }

  // 2. Predict next state

  //set weights
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  
  for (int i = 1; i < 2*n_aug_+1; i++) {
    weights_(i) = 1/(2 * (lambda_ + n_aug_)); 
  }
  
  //predict state mean
  
  x_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ += (Xsig_pred_.col(i).array() * weights_(i)).matrix();
  }
  
  //predict state covariance matrix
  P_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    MatrixXd error = Xsig_pred_.col(i) - X_aug;
    P_ += weights_(i) * error * error.transpose();
  }
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
  int n_z_ = 3;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);

  MatrixXd z = meas_package.raw_measurements_;

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  for (int i = 0; i < 2*n_aug_+1; i++) {
      
      MatrixXd col = Xsig_pred_.col(i);
      
      double px = col(0);
      double py = col(1);
      double v  = col(2);
      double phi= col(3);
      
      
      VectorXd sig(3);
      
      sig << sqrt(px * px + py * py),
             atan2(py, px),
             (px * cos(phi) * v + py * sin(phi) * v) / sqrt(px * px + py * py);
      
      // TODO: normalize atan2
      
      Zsig.col(i) = sig;
  }

  //calculate mean predicted measurement
  
  for (int i = 0; i < 2*n_aug_+1; i++) {
      z_pred += weights_(i)* Zsig.col(i);
  }
 
  MatrixXd R(3, 3);
  R << std_radr_ * std_radr_,   0,   0,
       0, std_radphi_ * std_radphi_, 0,
       0,  0,  std_radrd_ * std_radrd_;
  
  //calculate measurement covariance matrix S
  for (int i = 0; i < 2*n_aug_+1; i++) {
      MatrixXd error = Zsig.col(i) - z_pred;
      S += weights_(i) * error * error.transpose();
  }
  
  S += R;
    

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  
  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * ((Zsig.col(i) - z_pred).transpose());
  }
  
  // TODO: normalize Zsig - Zpred [1]

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K * S * K.transpose();
}

