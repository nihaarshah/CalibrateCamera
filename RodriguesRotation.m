function [ R ] = RodriguesRotation(RotationAxis,RotationAngle ) %RodriguesRotation
% Generates a rotation matrix given the axis of rotation and the
% angle of rotation in radians. This is not an angle-axis
% representation (where angle is encoded as the norm of the axis). %
% RotationAxis must be a vector with 3 elements and must be non
% zero. The code normalises the vector and thus it does not matter
% if the vector is not of unit lenth.
% RotationAngle is in radians. %
% The algorithm used is Rodrigues Rotation Formula


% Peform tests on correctness of input if length(RotationAxis) ~=3

if length(RotationAxis) ~=3
    
   error('The rotation axis must only have 3 elements')
end

NormAxis = norm(RotationAxis);
if NormAxis < eps
error('Cannot normalise the axis reliably');
end

% This step may not be necessary if the axis is already unit norm RotationAxis = RotationAxis / NormAxis;
% Formulate the cross product matrix
K = [ 0 -RotationAxis(3) RotationAxis(2); ...
RotationAxis(3) 0 -RotationAxis(1) ; ... 
-RotationAxis(2) RotationAxis(1) 0];
% Compute the rotation matrix using Rodrigues Rotation Formula
R = eye(3) + sin(RotationAngle)*K + (1 - cos(RotationAngle))*K^2;
end
