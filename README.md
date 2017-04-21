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
