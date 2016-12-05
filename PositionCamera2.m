function[T_cw] =  PositionCamera2(T_ow, Distance)
% PositionCamera
% Generates a random camera frame that almost points at the origin of the
% calibration object back along the object's z axis (which is perturbed).
% Aim is to 'fill' the camera image with calibration grid.

% T_ow is the 4x4 frame of the calibration grid in world coordinates.
% Distance is the base distance to the grid from the camera along the
% normal and a random factor of 0.1 times the distance is added.

% T_cw is the 4x4 camera frame in world coordinates

% Assign space for the camera frame
T_cw = zeros(4);

% Set the homogenous multiplier to 1
T_cw(4,4) = 1;

% extract the object origin
ObjectOrigin = T_ow(1:3,4);

% View the camera from the 'Distance' with a bit of randomness along the
% negative z-axis of the grid. This vector will be almost parallel to the
% camera z axis.

InitialViewVector = -Distance*T_ow(1:3,3) + 0.1*Distance*rand(3,1);

% Define the origin of the camera frame in world coordinates
T_cw(1:3,4) = ObjectOrigin - InitialViewVector;

% Define the camera z-axis as pointing at the object origin
Normz = norm(InitialViewVector);

if Normz < eps
    error('Unable to normalize the camera z axis');
end
% Define a unit vector pointing at the object
InitialCameraz = InitialViewVector/Normz;

% Pertub the initial z axis a bit for more luck
CameraZ = InitialCameraz -  0.01*(rand(3,1) - 0.5);

% ...and normalize again (no need to check norm)
CameraZ = CameraZ / norm(CameraZ);

% Define a random camera x-axis
CameraX = rand(3,1);

% Project out the z axis
CameraX = CameraX - (CameraZ' * CameraX)*CameraZ;
% normalize the x axis
Normx = norm(CameraX);

if Normx < eps
    error ('Unable to normalize the camera x axis')
end

CameraX = CameraX/Normx;

% Define the y axis
CameraY = cross(CameraZ, CameraX);

% Complete the transformation matrix
T_cw(1:3,1:3) = [CameraX CameraY CameraZ];

end
