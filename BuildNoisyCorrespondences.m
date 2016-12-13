function [Correspond] = BuildNoisyCorrespondences(T_ow, T_cw, CalibrationGrid, KMatrix ...
    ,CameraHeight, CameraWidth)

T_oc = T_cw \ T_ow;
Perspectivity = horzcat(T_oc(1:3,1:2), T_oc(1:3,4));  % [r1 r2 t]
Homography = KMatrix * Perspectivity;

% Calibration grid is [x y 0 1] for each intersection grid point
[m n] = size(CalibrationGrid);%Will need to run for loop number of row times
Correspond =[]; %Initializing

for j=1:m
    qo=horzcat(CalibrationGrid(j,1:2),CalibrationGrid(j,4)); % [x y 1]
    qc =Homography * qo'; % [u v 1]'*S?
    if (0<qc(1)<CameraWidth && 0<qc(2)<CameraHeight)
        
        noise = random('norm',0,0.0,3,1);
        qcnoisy = qc+noise;
        newxy = qo(1,1:2);  % [x y]
        newuv = qcnoisy(1:2,1); % [u v]'
        newuvxy = vertcat(newuv,newxy');
        Correspond = horzcat(newuvxy,Correspond);
    end
    
    % output is Correspond = [ u1 x1 u2 x2 u3 x3
    %                         v1 y1 v2 y2 v3 y3 ]
    
    
    % 4 row vector with noise. Multiply each term of the uvxy matrix with a
    % noisy term from x. But we will multiply each column with same noise
    % (maybe not what we want?)
    %  x = random('norm', 0, 0.1, 4, 1);
    %  Correspond = Correspond * diag(x);
    
end

%Correspond = horzcat(GridCornersC(1:2,1:4),GridCorners(1:2,1:4));


% Correspond = [GridCornersC(1,1:4) GridCornersC(2,1:4)
%                GridCorners(1,1:4) GridCorners(2,1:4)];
% Not sure what arrangement is expected from matrix currently we have
% [u1 u2 u3 u4 v1 v2 v3 v4; x1 x2 x3 x4 y1 y2 y3 y4] where u,v camera
% cordinates and xy are object(grid) coorinates

%Also not sure how to use cameraheight and camerawidth. is it because we
%are concerned with grid corners inside the image and not necessarily the
%actual grid corners?


