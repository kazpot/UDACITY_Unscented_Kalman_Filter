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

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 3;

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
  
  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  //set predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //set sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  //set current NIS for radar
  NIS_radar_ = 0;

  //set current NIS for laser
  NIS_laser_ = 0;

  //set weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.5/(n_aug_ + lambda_));
  weights_(0) = lambda_/(lambda_ + n_aug_);

  is_initialized_ = false;
  previous_timestamp_ = 0.0;
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
  
    // init covariance matrix
    P_ <<    1, 0, 0, 0, 0,
           0, 1, 0, 0, 0,
           0, 0, 1, 0, 0,
           0, 0, 0, 1, 0,
           0, 0, 0, 0, 1;

    //init timestamp
    previous_timestamp_ = meas_package.timestamp_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];

      double px = rho * cos(phi); 
      double py = rho * sin(phi);

      x_ << px,py,0,0,0; 
    }else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];

      x_ << px,py,0,0,0;
    }


    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  //calculate dt - expressed in seconds
  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; 
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
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

  1. generate augmented sigma points
  2. predict sigma points
  3. predict mean and covariance
  */
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;
  
  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  
  for (int i=0; i < n_aug_; i++){
     Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
     Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  for (int i=0; i < 2 * n_aug_ + 1; i++){
      double p_x      = Xsig_aug(0,i);
      double p_y      = Xsig_aug(1,i);
      double v        = Xsig_aug(2,i);
      double yaw      = Xsig_aug(3,i);
      double yawd     = Xsig_aug(4,i);
      double nu_a     = Xsig_aug(5,i);
      double nu_yawdd = Xsig_aug(6,i);
      
      double px_p;
      double py_p;

      //avoid zero devision error
      if(fabs(yawd) > 0.0001){
          px_p = p_x + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
          py_p = p_y + v/yawd * (cos(yaw)-cos(yaw + yawd * delta_t));
      }else{
          px_p = p_x + v * cos(yaw) * delta_t;
          py_p = p_y + v * sin(yaw) * delta_t;
      }
      
      double v_p = v;
      double yaw_p = yaw + yawd * delta_t;
      double yawd_p = yawd;
      
      //add noise
      px_p = px_p + 0.5 * delta_t * delta_t * cos(yaw) * nu_a;
      py_p = py_p + 0.5 * delta_t * delta_t * sin(yaw) * nu_a;
      v_p = v_p + delta_t * nu_a;
      yaw_p = yaw_p + 0.5 * delta_t * delta_t * nu_yawdd;
      yawd_p = yawd_p + delta_t * nu_yawdd;
      
      Xsig_pred_(0,i) = px_p;
      Xsig_pred_(1,i) = py_p;
      Xsig_pred_(2,i) = v_p;
      Xsig_pred_(3,i) = yaw_p;
      Xsig_pred_(4,i) = yawd_p;
  }

  x_.fill(0.0);
  for(int i=0; i < 2 * n_aug_ + 1; i++){
      x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  
  P_.fill(0.0);
  for(int i=0; i < 2 * n_aug_ + 1; i++){
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      
      //angle normalization
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3) += 2. * M_PI;
      
      P_ = P_ + weights_(i) * x_diff * x_diff.transpose(); 
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
  VectorXd z = meas_package.raw_measurements_;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  for(int i=0; i < 2 * n_aug_ + 1; i++){
      double px  = Xsig_pred_(0,i);
      double py  = Xsig_pred_(1,i);
      
      Zsig(0,i) = px;
      Zsig(1,i) = py;
  }
  
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for(int i=0; i < 2 * n_aug_ + 1; i++){
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  
  //calculate measurement covariance matrix
  S.fill(0.0);
  for(int i=0; i < 2 * n_aug_ + 1; i++){
      VectorXd z_diff = Zsig.col(i) - z_pred;
      S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  
  //calculate measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);;
  R << std_laspx_ * std_laspx_, 0,
       0, std_laspy_ * std_laspy_;
  S = S + R;

  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for(int i=0; i < 2 * n_aug_ + 1; i++){
      VectorXd z_diff = Zsig.col(i) - z_pred;
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      
      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  //Kalman gain
  MatrixXd K = Tc * S.inverse();

  //calculate NIS
  NIS_laser_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
  
  //update state mean and covariance matrix
  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K * S * K.transpose();

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
  VectorXd z = meas_package.raw_measurements_;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  for(int i=0; i < 2 * n_aug_ + 1; i++){
      double px  = Xsig_pred_(0,i);
      double py  = Xsig_pred_(1,i);
      double v   = Xsig_pred_(2,i);
      double yaw = Xsig_pred_(3,i);
      
      double rho  = sqrt(px * px + py * py);
      double phi  = atan2(py,px);
      double rhod = (px * cos(yaw) * v + py * sin(yaw) * v)/sqrt(px * px + py * py);
      
      Zsig(0,i) = rho;
      Zsig(1,i) = phi;
      Zsig(2,i) = rhod;
  }
  
  //calculate mean prdicted measurement
  z_pred.fill(0.0);
  for(int i=0; i < 2 * n_aug_ + 1; i++){
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  
  //calculate measurement covariance matrix
  S.fill(0.0);
  for(int i=0; i < 2 * n_aug_ + 1; i++){
      VectorXd z_diff = Zsig.col(i) - z_pred;
      
      //angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
      
      S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  
  //calculate measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);;
  R << std_radr_ * std_radr_, 0, 0,
       0, std_radphi_ * std_radphi_,0,
       0,0,std_radrd_ * std_radrd_;
  S = S + R;

  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for(int i=0; i < 2 * n_aug_ + 1; i++){
      VectorXd z_diff = Zsig.col(i) - z_pred;
      
      //angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
      
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      
      //angle normalization
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
      
      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  //calculate Kalman gain
  MatrixXd K = Tc * S.inverse();
  
  //calculate NIS
  NIS_radar_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);

  //update state mean and covariance matrix
  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K * S * K.transpose();

}
