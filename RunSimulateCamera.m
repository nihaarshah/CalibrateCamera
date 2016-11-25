% Simulates viewing a 3D object (a cube) to demonstrate the 
% structures built up to calibrate a camera

%Construct a camera
[KMatrix, CameraHeight, CameraWidth] = BuildCamera;

%Construct an object in its own frame
Cube = BuildCube;

%Position the object in space
T_ow = PositionObject;

%Position the camera so that it is likely that the object can be seen
T_cw = PositionCamera(T_ow);

%Look at what is produced
ViewCamera(Cube, T_ow, KMatrix, CameraHeight, CameraWidth, T_cw)

