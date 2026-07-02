"""Relative navigation and sensor models.

Modules
-------
ekf      Extended Kalman filter on CW/TH relative dynamics; fixed-dimension,
         bounded-iteration implementation targeted at the static subset.
sensors  Measurement models with error budgets: relative GNSS (carrier-phase
         differential), vision-based bearing/range (docking camera + fiducial
         markers), lidar/rangefinder, and star-tracker/gyro attitude.
"""
