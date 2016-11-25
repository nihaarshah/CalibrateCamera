function [T_ow] = PositionObject()
%PositionObject
%Defines a 4x4 homogenous transformation that places an object 
%in space at a random location in world coorinates and at a 
%random orientation. T_ow means T:O->W i.e. the transformation 
% of Points in object coordinates into Points in World Coordinates.

T_ow = zeros(4,4);

%The bottom row is [0 0 0 1]
T_ow(4,:) = [0 0 0 1];

%Put some random location
T_ow(1:3, 4) = 2000*rand(3,1);

% With a random orientation
T_ow(1:3, 1:3) = RandomRotationMatrix;

end


