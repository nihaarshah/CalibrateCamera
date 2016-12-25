function [Homography] = CalculateHomog(RotAxis,KMatrix,Translation)

Translation = Translation';
Theta = norm(RotAxis); % Maybe the axis of rotation is theta * unit 
    % length axis of rotation so taking the norm will give us the magnitude
    % of this axis which is theta?
% n = 4;
% x = RotAxis(1)/Theta;
% y = RotAxis(2)/Theta;
% z = RotAxis(3)/Theta;
% 
% K =[0 -z y; z 0 -x; -y x 0]; %Different K matrix 

% Theta = Phi/RotAxis - n*pi; What i thought earlier

% R = [1 0 0;0 1 0;0 0 1] + sin(Theta)*K + (1-cos(Theta))*K^2;

[ R ] = RodriguesRotation(RotAxis,Theta);

P(1:3,1:3) = R;
P(1:3,4) = Translation;
% P(4,1:3) = [0 0 0];
% P(4,4) = 1;

% P = [R Translation; [0 0 0] 1];
Perspectivity = [P(1:3,1:2),P(1:3,4)];
% One = ones(size(Correspond));
% Correspond = vertcat(Correspond,One);
Homography = KMatrix * Perspectivity;
end