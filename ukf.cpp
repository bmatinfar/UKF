#include "ukf.h"
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
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.28; // trail and error

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.8; //trial and error

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
    
    is_initialized_ = false;
    time_us_ = 0;
    
    n_x_ = 5;
    
    n_aug_ = n_x_ + 2;
    
    lambda_ = 3 - n_x_;
    
    int n_z = 3;
    
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    Xsig_pred_.fill(0.0);
    
    H = MatrixXd(2, n_x_);
    R_Lidar = MatrixXd(2, 2);
    
    R_Radar= MatrixXd(n_z,n_z);
}

UKF::~UKF() {
    myfile.close();
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

    if(!is_initialized_)
    {
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            x_(0) = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
            x_(1) = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
            //double vx = meas_package.raw_measurements_[2] * cos(meas_package.raw_measurements_[1]);
            //double vy = meas_package.raw_measurements_[2] * sin(meas_package.raw_measurements_[1]);
            x_(2) = 0.0; //sqrt(vx*vx + vy*vy);
            x_(3) = 0.0;
            x_(4) = 0.0;
            
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            x_(0) = meas_package.raw_measurements_[0];
            x_(1) = meas_package.raw_measurements_[1];
            x_(2) = 0.0;
            x_(3) = 0.0;
            x_(4) = 0.0;
        }


        //Initial state covariance matrix P
        P_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0,
              0, 0, 1 , 0, 0,
              0, 0, 0, 1 , 0,
              0, 0, 0, 0, 1 ;
        
        H<< 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;

        R_Lidar << std_laspx_ * std_laspx_, 0,
        0, std_laspy_ * std_laspy_;

        R_Radar <<    std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0,std_radrd_*std_radrd_;
        
        time_us_ = meas_package.timestamp_;

        is_initialized_ = true;
        myfile.open("RadarNIS");
        return;
    }
    
    //compute the time elapsed between the current and previous measurements
    double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
    time_us_ = meas_package.timestamp_;
    Prediction(dt);
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        UpdateRadar(meas_package);
    } else {
        
        UpdateLidar(meas_package);
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

    //cout<<"State before prediction: " <<x_<<endl;
    
    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);
    
    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    
    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    
    //create augmented mean state
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;
    
    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;
    
    //MatrixXd Q = MatrixXd(2,2);
    //Q << std_a_*std_a_, 0 , 0, std_yawdd_*std_yawdd_;
    
    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();
    
    //create augmented sigma points
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }
    
    //predict sigma points
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        const double p_x = Xsig_aug(0,i);
        const double p_y = Xsig_aug(1,i);
        const double v = Xsig_aug(2,i);
        const double yaw = Xsig_aug(3,i);
        const double yawd = Xsig_aug(4,i);
        const double nu_a = Xsig_aug(5,i);
        const double nu_yawdd = Xsig_aug(6,i);
        
        //predicted state values
        double px_p, py_p;
        
        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;
        
        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;
        
        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;
        
        //write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
    }
    
    //create vector for weights
    VectorXd weights = CreateWeights(lambda_);
    
    //predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x_ = x_+ weights(i) * Xsig_pred_.col(i);
    }
    
    //cout<<"predicted state:" << x_ << endl;

    //predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        
        //angle normalization
        NormalizeAngle(x_diff(3));
        
        P_ = P_ + weights(i) * x_diff * x_diff.transpose() ;
    }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

    MatrixXd I = MatrixXd::Identity(5, 5);
    VectorXd y = meas_package.raw_measurements_ - H * x_;
    MatrixXd Ht = H.transpose();
    MatrixXd PHt = P_ * Ht;
    MatrixXd S_Lidar = H * PHt + R_Lidar;
    MatrixXd Si = S_Lidar.inverse();
    MatrixXd K =  PHt * Si;
    
    x_ = x_ + K * y;
    P_ = (I - K*H)* P_.transpose();
    
    //cout<<"updated state from Lidar measurement;" << x_ <<endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
 
    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;
    
    //define spreading parameter
    double lambda = 3 - n_aug_;
    
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig.fill(0);

    
    //create vector for weights
    VectorXd weights = CreateWeights(lambda);
    
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
        
        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
        if (p_x < 0.0001 & p_y < 0.0001)
            Zsig(1,i) = 0;
        else
            Zsig(1,i) = atan2(p_y,p_x);                                 //phi
        if(Zsig(0,i) < 0.001)
            Zsig(2,i) = 0;
        else
            Zsig(2,i) = (p_x*v1 + p_y*v2 ) / Zsig(0,i);   //r_dot
    }
    
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights(i) * Zsig.col(i);
    }

    //cout<<"predicted state for Radar measurement:" <<z_pred <<endl;
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        //angle normalization
        NormalizeAngle(z_diff(1));
        
        S = S + weights(i) * z_diff * z_diff.transpose();
    }
    
    //add measurement noise covariance matrix
    S = S + R_Radar;
    
    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    //calculate cross correlation matrix
    Tc.fill(0.0);

    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        //angle normalization
        NormalizeAngle(z_diff(1));
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        
        //angle normalization
        NormalizeAngle(x_diff(3));
        
        Tc = Tc + weights(i) * x_diff * z_diff.transpose();
    }
    
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    //cout<< "raw measurement in radar:" <<meas_package.raw_measurements_ << endl;
    
    //residual
    VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
    
    //angle normalization
    NormalizeAngle(z_diff(1));
    
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    
    
    double nis =  z_diff.transpose() * S.inverse() * z_diff;
    
    //cout<<"updated state from Radar measurement:" << x_ <<endl;
    
    //cout<< "nis for Radar:" << nis<<endl;

    myfile << nis <<endl;

}

VectorXd UKF::CreateWeights(double lambda)
{
    //create vector for weights
    VectorXd weights = VectorXd(2*n_aug_+1);

    // set weights
    double weight_0 = lambda/(lambda+n_aug_);
    weights(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
        double weight = 0.5/(n_aug_+lambda);
        weights(i) = weight;
    }    
    return weights;
}

void UKF::NormalizeAngle(double& phi)
{
    phi = atan2(sin(phi), cos(phi));
}

