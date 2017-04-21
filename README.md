# Unscented Kalman Filter Project #

This project is originally forked from https://github.com/udacity/CarND-Unscented-Kalman-Filter-Project.  This repository includes starter code, that is used herein.

## Initialization ##

As with the EKF case, there are two possibilities for initializing the state estimate, depening on which sensor provides the initial measurement.

If the first available measurement comes from the lidar sensor, then the initialization is simple.  We simply take the measured `x` and `y` positions as the initial values for `px` and `py`, and assume `v`, `phi` and `phi_dot`  are initially `0.0`.

If the first available measurement comes from the radar sensor, then we need to convert the polar coordinates `rho` and `phi` into Cartesian coordinates.  Unfortunately, we don't know the direction of `rho_dot`, so we can't use this information to initialize `vx` or `vy`.  Like the lidar case, assume `v`, `phi` and `phi_dot` are initially `0.0`.

*From ukf.cpp*
```C++
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
  if ((meas_package.raw_measurements_(0) > 0.001) && (meas_package.raw_measurements_(1) > 0.001)) {
    x_(0) = meas_package.raw_measurements_(0);
    x_(1) = meas_package.raw_measurements_(1);
  } else {
    // default values if first measurement is close to (x, y) = (0.0, 0.0)
    x_(0) = 0.1;
    x_(1) = 0.1;
  }
}
```

*From tools.cpp*
```C++
VectorXd Tools::Polar2Cartesian(const VectorXd& radar_meas) {
  /**
   * Convert polar coordinates to Cartesian.
   */
  VectorXd f_x(5);
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
```

My UKF implementation satisfies the rubric accuracy criteria for both datasets.  For the case of *obj_pose-laser-radar-synthetic-input.txt* the results are greatly improved with the UKF when compared to the EKF.

#### RMSE for sample-laser-radar-measurement-data-1.txt ####
|State |RMSE UKF |RMSE EKF |RMSE Limit                     
|:-----|:--------|:--------|:---------
|`px` |0.053 |0.065 |0.090 
|`py` |0.064 |0.062 |0.090 
|`vx` |0.529 |0.544 |0.650 
|`vy` |0.547 |0.544 |0.650 

#### RMSE for obj_pose-laser-radar-synthetic-input.txt ####
|State |RMSE UKF |RMSE EKF |RMSE Limit                     
|:-----|:--------|:--------|:---------
|`px` |0.073 |0.140 |0.090 
|`py` |0.085 |0.665 |0.100 
|`vx` |0.357 |0.579 |0.400 
|`vy` |0.244 |1.634 |0.300 
