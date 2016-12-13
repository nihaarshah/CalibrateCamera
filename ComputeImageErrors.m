function [Error] = ComputeImageErrors(KMatrix,RotAxis,...
    Translation,Correspond, BestConsensusIndex)

Theta = norm(RotAxis); % Maybe the axis of rotation is theta * unit 
    % length axis of rotation so taking the norm will give us the magnitude
    % of this axis which is theta?
n = 4;
x = RotAxis(1)/Theta;
y = RotAxis(2)/Theta;
z = RotAxis(3)/Theta;

K =[0 -z y; z 0 -x; -y x 0]; %Different K matrix 

% Theta = Phi/RotAxis - n*pi; What i thought earlier

R = [1 0 0;0 1 0; 0 0 1] + sin(Theta)*K + (1-cos(Theta))*[K]^2;

Perspectivity = [R Translation; [0 0 0] 1];

Homography = KMatrix * Perspectivity;

Estimates = Homography * Correspond(3:4,BestConsensusIndex); % [u' v'] or the estimated u v coordinates

Error = (Estimates - Correspond(1:2,BestConsensusIndex))^2; % error in pixels



