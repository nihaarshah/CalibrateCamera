function [ T_cw ] = FillImage(T_ow,KMatrix,GridWidth,CameraHeight,CameraWidth)

% FillImage Takes in a Calibration grid definition and positions the
% camera so that the image of the grid fills the camera image.

% T_ow is the Calibration Grid frame KMatrix is the intrinic camera model
% Gridwidth is the length of a size of the whole grid im mm
% CameraHeight and CameraWidth are the sizes of the image in pixels

% Define the Left Top, Left Bottom, Right Top and Right Bottom of the grid .
GridCorners = [-GridWidth/2 -GridWidth/2 GridWidth/2 GridWidth/2;
    GridWidth/2 -GridWidth/2 GridWidth/2 -GridWidth/2;
    0 0 0 0;
    1 1 1 1];

% Compute the positions of the GridCorners in the world
GridCorners = T_ow*GridCorners;

% We have a 1m by 1m grid somewhere in space and we need to view
% the grid from the camera. We view from a random location based
% on a 'distance' called CameraBaseDistance. this is an initial estimate.
CameraBaseDistance = 2000;

% Keep reducing the distance until all 4 corners are outside the image
% and InsideImage is a flag that records the failure of our objective
% thus triggering a move towards the object and a retry.
InsideImage = 1;
while InsideImage == 1
    InsideImage = 0;
    % PositionCamera makes the camera 'almost' point at the object, but
    % introduces randomness so we have a range of angles betweeen the
    % grid's normal and the camera z-axis.
    
    T_cw = PositionCamera2(T_ow,CameraBaseDistance);
    
    % Compute where corners are in the unit camera frame
    UnitCorners = (T_cw \ GridCorners);
    % Project into the unit camera
    UnitCorners = UnitCorners(1:3,:);
    % And convert to camera pixels and label them as homogeneous.
    HomogeneousCorners = KMatrix * UnitCorners;
    % Dehomogenise the image of the corners of the grid
    Corners = zeros(2,4);
    for j = 1:4
        
        Corners(1:2,j) = ...
            HomogeneousCorners(1:2,j) / HomogeneousCorners(3,j);
        if Corners(1,j) > 0 && Corners(1,j) < CameraWidth-1
            % The u component is not outside the image - keep searching
            InsideImage = 1;
        end
        if Corners(2,j) > 0 && Corners(2,j) < CameraHeight-1
            InsideImage = 1;
        end
    end
    
    % Move the camera nearer to the object if any of the corners are
    % inside the image. The next time around we have a better chance
    % of achieving our objective.
    if InsideImage == 1
        CameraBaseDistance = CameraBaseDistance * 0.9;
    end
    
end
end