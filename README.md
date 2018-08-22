# CalibrateCamera #

Welcome to my project on calibrating a camera for computer vision using Matlab. The aim is to find optimal intrinsic parameters 
of the K matrix. 

This is accomplished by first simulating a scene and finding point correspondences to create a Homography matrix. 

Then RANSAC is used to minimise outliers and finally optimisation scripts are written to fit the optimal parameters
to the correspondence points.

This project is outlined in the report tinyurl.com/b1-nihaar
