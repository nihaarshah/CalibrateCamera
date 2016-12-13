%Example script for estimation task 1

%This script performs the following actions:
% 1. Constructs a camera model loosely based on an iPhone 6
% 2. Constructs a calibration grid of 1m on a side with 10mm grid spacing.
% 3. Positions the grid somewhere in space
% 4. Places the camera in a 'random' location so that the grid fills the
% camera image. Check this location to see if the grid outline is outside
% the image; if not, move the camera nearer to the grid. This is NOT a smart
% way of positioning things, but it models the likely begavior of a human
% user.
% 5. Compute the image of the grid and construct a set of grid-image
% correspondences and add 'noise' to the image.
% 6. Convert some of the correspondences to outliers.
% 7. Carry out a RANSAC estimation of the homography

% 1. COnstruct the camera model
[KMatrix, CameraHeight, CameraWidth] = BuildCamera;

%2. Construct a 1m by 1m grid with 10mm tiles in the grid frame
% The grid is a set of 4-element vectors [x y 0 1]'
GridWidth = 100;
GridIncrement = 10;
CalibrationGrid = BuildGrid(GridIncrement, GridWidth);

%3. Choose somehwere in space for the grid
% T_ow is the 4x4 transformation matrix from the grid to the world
T_ow = PositionGrid();


% 4. Choose a 'random' location for the camera that fills the image.
% T_cw is the 4x4 ransformation matrix from camera to the world

T_cw = FillImage (T_ow, KMatrix,GridWidth, CameraHeight, CameraWidth);

% 5. Now fill the camera with a noisy image of the grid and generate the
% point correspondences 
% Correspond is a set of pairs of vectors of the form [[u v]' [x y]']
%for each grid corner that lies inside the image.

Correspond = BuildNoisyCorrespondences(T_ow, T_cw, CalibrationGrid, KMatrix ...
,CameraHeight, CameraWidth);%Not sure how to use CameraWidth and Height. Not added noise yet

% 6. Add in some 'outliers' by replacing [u v]' with a point
% somewhere in the image.
% Define the outlier probability
pOutlier = 0.05;

for j=1:length(Correspond)
    if rand < pOutlier
        % This is an outlier - so put the point anywhere in the image.
        Correspond(1, j) = rand * (CameraWidth-1);
        Correspond(2, j) = rand * (CameraHeight-1);
    end
end

figure(1)
plot(Correspond(1,:), Correspond(2,:), '.')
title('The noisy measurements of the title corners')
axis ij

% 7 Perform the Ransac estimation - output the result for the inspection
% If the Ransac fails, it returns a zero Homography
Maxerror = 3;   %The maximum error allowed before rejecting a point
RansacRuns = 50;    %The number of runs when creating the consensus set
[Homog, BestConsensus] = RansacHomog(Correspond, Maxerror, RansacRuns);

% If you want to test the result, we can construct the homography for the
% system from its definition

% First find the object frame in the camera frame
T_oc = T_cw \ T_ow;
%Construct the non-normalized homography from K*[x y t]
OrigHomog = KMatrix * [T_oc(1:3,1) T_oc(1:3,2) T_oc(1:3,4)];
%And normalize so that (3,3) is 1.0-output for inspection
OrigHomog = OrigHomog / OrigHomog(3,3);