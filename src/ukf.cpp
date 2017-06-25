#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

double epsilon = 1e-5;

inline double normalizeAngle(double angle) {
    while (angle > M_PI) {
        angle -= 2. * M_PI;
    }
    while (angle < -M_PI) {
        angle += 2. * M_PI;
    }
    
    return angle;
}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;
    
    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = false;
    
    is_initialized_ = false;
    
    // initial state vector
    x_ = VectorXd(5);
    
    // initial covariance matrix
    P_ = MatrixXd(5, 5);
    
    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 6;
    
    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 5;
    
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
    
    // Number of state dimensions (x, y, long-vel, yaw, yawrate)
    n_x_ = 5;
    // Number of augmented state dimensions (n_x_ + process noise for x, y)
    n_aug_ = 7;
    
    // Lambda values for x_ and Xaug matrices
    lambda_x_ = 3 - n_x_;
    lambda_aug_ = 3 - n_aug_;
    
    // Weight matrices
    weights_ = VectorXd(2 * n_aug_ + 1);
    
    // Predicted sigma points matrices
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    
    // Initial covariance matrix
    P_ << 0.15, 0, 0, 0, 0,
    0, 0.15, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1;
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
    
    if (!is_initialized_) {
        time_us_ = meas_package.timestamp_;
        std::cout << "Is initialized @ " << time_us_ << std::endl;
        
        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1],
            epsilon,
            epsilon,
            epsilon;
        } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            double rho = meas_package.raw_measurements_[0];
            double phi = meas_package.raw_measurements_[1];
            
            x_ << rho * cos(phi),
            rho * sin(phi),
            epsilon,
            epsilon,
            epsilon;
        }
        
        is_initialized_ = true;
        return;
    }
    
    double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.;
    
    time_us_ = meas_package.timestamp_;
    
    // Prediction step
    
    Prediction(delta_t);
    
    // Update step
    
    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        std::cout << "RADAR" << std::endl;
        UpdateRadar(meas_package);
    } else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
        std::cout << "LIDAR" << std::endl;
        UpdateLidar(meas_package);
    }
    
    std::cout << "Update done" << std::endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    VectorXd X_aug = VectorXd(7);
    MatrixXd P_aug = MatrixXd(7, 7);
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    
    
    //
    MatrixXd Q(2, 2);
    Q << std_a_ * std_a_, 0,
    0, std_yawdd_ * std_yawdd_;
    
    X_aug.fill(0);
    P_aug.fill(0);
    Xsig_aug.fill(0);
    
    //create augmented mean state
    X_aug.head(5) = x_;
    X_aug(5) = 0;
    X_aug(6) = 0;
    
    //create augmented covariance matrix
    P_aug.fill(0);
    P_aug.topLeftCorner( 5, 5 ) = P_;
    P_aug.bottomRightCorner( 2, 2 ) = Q;
    
    // std::cout << "P_aug = " << std::endl << P_aug << std::endl;
    
    //create square root matrix
    MatrixXd A_aug = P_aug.llt().matrixL();
    
    double sql_na = sqrt(lambda_x_ + n_aug_);
    
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
        
        
        offset.fill(0);
        process_noise.fill(0);
        
        if (phidotk == 0) {
            phidotk = epsilon;
        }
        
        double v_phidot;
        if (fabs(phidotk) < epsilon) {
            v_phidot = vk;
        } else {
            v_phidot = vk / phidotk;
        }
        
        
        
        
        offset << v_phidot * ( sin(phik + phidotk * delta_t) - sin(phik) ),
        v_phidot * ( -cos(phik + phidotk * delta_t) + cos(phik) ),
        0,
        0, // phidotk * delta_t
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
    weights_(0) = lambda_aug_ / (lambda_aug_ + n_aug_);
    
    for (int i = 1; i < 2*n_aug_+1; i++) {
        weights_(i) = 1/(2 * (lambda_aug_ + n_aug_));
    }
    
    
    //predict state mean
    
    x_.fill(0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        x_ += weights_(i) * Xsig_pred_.col(i);
    }
    std::cout << "Estimated state: " << std::endl << x_ << std::endl;
    
    
    //predict state covariance matrix
    P_.fill(0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd error = Xsig_pred_.col(i) - x_;
        error(3) = normalizeAngle(error(3));
        P_ += weights_(i) * error * error.transpose();
    }
    // std::cout << "P_" << std::endl << P_ << std::endl;
    
    std::cout << "Prediction complete" << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    
    // Take in new measurements (x, y)
    MatrixXd z = meas_package.raw_measurements_;
    int n_z_ = 2;
    
    // Transform predicted sigma points to measurement space (px, py)
    MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
    
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd col = Xsig_pred_.col(i);
        VectorXd z_col = Zsig.col(i);
        
        z_col(0) = col(0);
        z_col(1) = col(1);
        
        Zsig.col(i) = z_col;
    }
    
    
    // Find the mean predicted measurement
    MatrixXd z_pred = VectorXd(n_z_);
    z_pred.fill(0);
    
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred += weights_(i) * Zsig.col(i);
    }
    
    
    // Find the measurement covariance
    MatrixXd S = MatrixXd(n_z_, n_z_);
    
    S.fill(0);
    
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd error = Zsig.col(i) - z_pred;
        S += weights_(i) * error * error.transpose();
    }
    
    
    // Add the sensor's measurement inaccuracy
    MatrixXd R = MatrixXd(n_z_, n_z_);
    
    R << std_laspx_ * std_laspx_, 0,
    0, std_laspy_ * std_laspy_;
    
    S += R;
    
    
    // Calculate the cross-correlation matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z_);
    Tc.fill(0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd error = Zsig.col(i) - z_pred;
        VectorXd x_err = Xsig_pred_.col(i) - x_;
        x_err(3) = normalizeAngle(x_err(3));
        
        Tc += weights_(i) * x_err * error.transpose();
    }
    
    // Find the Kalman gain
    MatrixXd K = Tc * S.inverse();
    
    // Update x and p matrices
    x_ += K * (z - z_pred);
    P_ -= K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    int n_z_ = 3;
    MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
    VectorXd z_pred = VectorXd(n_z_);
    MatrixXd S = MatrixXd(n_z_,n_z_);
    MatrixXd z = meas_package.raw_measurements_;
    
    // Transform the sigma points into measurement space
    //  (x,y,v,phi,phidot) -> (rho, phi, phidot)
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
        
        Zsig.col(i) = sig;
    }
    
    
    // Find the mean predicted measurement
    z_pred.fill(0.);
    for (int i = 0; i < 2*n_aug_+1; i++) {
        z_pred += weights_(i)* Zsig.col(i);
    }
    
    // Radar inaccuracy matrix
    MatrixXd R(3, 3);
    R << std_radr_ * std_radr_,   0,   0,
    0, std_radphi_ * std_radphi_, 0,
    0,  0,  std_radrd_ * std_radrd_;
    
    // Find the measurement covariance matrix
    S.fill(0);
    for (int i = 0; i < 2*n_aug_+1; i++) {
        MatrixXd error = Zsig.col(i) - z_pred;
        
        error(1) = normalizeAngle(error(1));
        
        S += weights_(i) * error * error.transpose();
    }
    
    S += R;
    
    
    // Find cross correlation
    MatrixXd Tc = MatrixXd(n_x_, n_z_);
    Tc.fill(0);
    
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        
        VectorXd x_err = Xsig_pred_.col(i) - x_;
        x_err(3) = normalizeAngle(x_err(3));
        
        // Normalize the z_err
        
        VectorXd z_err = Zsig.col(i) - z_pred;
        z_err(1) = normalizeAngle(z_err(1));
        
        Tc += weights_(i) * x_err * z_err.transpose();
    }
    
    // Find the Kalman gain K
    MatrixXd K = Tc * S.inverse();
    MatrixXd z_err = z - z_pred;
    
    // Normalize phi between -pi, pi
    z_err(1) = normalizeAngle(z_err(1));
    
    // Update x_, P_
    x_ = x_ + K * z_err;
    P_ = P_ - K * S * K.transpose();
    
    std::cout << "Radar update: x_ = " << std::endl << x_ << std::endl;
    
    // TODO: Find Radar NIS
}

