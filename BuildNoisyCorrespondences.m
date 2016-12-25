function [Correspond] = BuildNoisyCorrespondences(T_ow, T_cw, CalibrationGrid, KMatrix ...
    ,CameraHeight, CameraWidth,mu,std)
 
% CalibrationGrid = CalibrationGrid';
% To generate a random matrix which has the same size as the
% calibrationgrid and add noise to the calibration grid
 
s1 = size(CalibrationGrid);
 xy = CalibrationGrid(1:2,:);
% Introduces noise with a random set variance
 noise = rand();
 
 
% Construction of grid point to image point 
CalibrationGrid = T_ow * CalibrationGrid;
CalibrationGrid = T_cw \ CalibrationGrid; % [x y z(=0) 1]
 
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
 
 
 
Correspond = zeros(s1);
 
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
 
% for j=1:m
%     qo=horzcat(CalibrationGrid(j,1:2),CalibrationGrid(j,4)); % [x y 1]
%     qc =Homography * qo'; % [u v 1]'*S?
%     qc(1) = qc(1)/qc(3); % normalizing
%     qc(2) = qc(2)/qc(3); % normalizing
% %      if (0<qc(1)<CameraWidth && 0<qc(2)<CameraHeight)
%         
%         %noise = random('norm',0,0.0,3,1);
%         %qcnoisy = qc+noise;
%         newxy = qo(1,1:2);  % [x y]
%         newuv = qc(1:2,1); % [u v]'
%         newuvxy = vertcat(newuv,newxy');
%         Correspond = horzcat(newuvxy,Correspond);
% %      end
%     
% end
 
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
 
 
 

