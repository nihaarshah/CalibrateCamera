function [Correspond] = BuildNoisyCorrespondences(T_ow, T_cw,...
    CalibrationGrid, KMatrix,CameraHeight, CameraWidth,mu,std)
% This function takes points from our simulated calibration grid in xy
% coordinates and transforms this into camera-sensor uv coordinates using
% the Homogenous transformation T_oc = T_ow\T_cw followed by multiplying
% with our original K-Matrix. We shall disregard any pixel values beyond
% the camera's dimensions. Finally we add Gaussian noise whose variance is
% passed as a flag which can be changed from the script.

xy = CalibrationGrid(1:2,:);

% Construction of grid point to image point
CalibrationGrid = T_ow * CalibrationGrid;
CalibrationGrid = T_cw \ CalibrationGrid; % output: [x y z(=0) 1]
% Transform the x y coordinate to u v coordinates
CalibrationImage = KMatrix * CalibrationGrid(1:3,:);

% Normalize the uv coordinates with respect to the z coordinate
CalibrationImage(1,:) = CalibrationImage(1,:) ./ CalibrationImage(3,:);
CalibrationImage(2,:) = CalibrationImage(2,:) ./ CalibrationImage(3,:);
CalibrationImage = CalibrationImage(1:2,:);

% Size of CalibrationGridImage so we can add noise
[m,n] = size(CalibrationImage);

% Adding Gaussian noise with a changeable variance to the image
CalibrationImage = CalibrationImage + random('norm',mu,std,m,n);

Correspond = zeros(size(CalibrationGrid));
TempCorrespond = zeros(size(Correspond));

Correspond(1:2,:) = CalibrationImage(1:2,:);
Correspond(3:4,:) = xy(1:2,:);

n = 0;
for i = 1:length(Correspond)
    
    if Correspond(1,i)>=0 && Correspond(1,i) <CameraWidth...
            && Correspond(2,i)>=0 && Correspond(2,i)<CameraHeight
        n = n+1;
        TempCorrespond(:,n) = Correspond(:,i);
    end
end
Correspond = TempCorrespond(:,1:n);
end