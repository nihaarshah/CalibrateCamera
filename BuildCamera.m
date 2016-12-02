function [KMatrix, CameraHeight, CameraWidth] = BuildCamera

Parameters = DefineParameters;
KMatrix = SingleVectorCameraModel( Parameters );
CameraHeight = Parameters(1);
CameraWidth = Parameters(2);