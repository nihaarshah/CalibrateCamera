
%Parameters come from KMatrix, Rotation axis and Translation vector
K = sym('k',[3,3]);
RotAxis = sym('rot',[1,3]);
t = sym('t',[3,1]);
Theta = norm(RotAxis);
RotAxis = RotAxis/Theta;

n = lenghth(BestConsensus);

% Grid point coordinate
syms x y

% L is a matrix made of rotaxis parameters. It is used to find the rotation
% matrix
L = [0, -rot(1,3), rot(1,2);rot(1,3),0,-rot(1,1);-rot(1,2),rot(1,1),0];
I = eye(3); % identity matrix 3x3


% Assembling the rotation matrix from the 3 rotation axis parameters
% packed in the L matrix
R = I + sin(Theta)*L + (1-cos(Theta))*L^2;

% Forming the Perspectivity from the Rotation axis parameters packed inside
% r1 and r2 and the Translation vector. Thus P contains all 6 extrinsic
% parameters-Rotation and translation
r1 = R(1:3,1);
r2 = R(1:3,2);
P = [r1,r2,t];

for i = 1:n
    [x,y]' = Correspond(3:4,BestConsensus(n));
    f = K * P * [x y 1]';
    J = jacobian(f,[K(1,1),K(1,2),K(1,3),K(2,2),K(2,3),rot(1,1),rot(2,1),rot(3,1),t(1,1),t(2,1),t(3,1)]);
end


KMatrix = [1000,0,500;0,1000,500;0,0,1];
RotAxis = [1 2 1]';
Translation = [1 2 1]';


J1 = subs(J,K,KMatrix);
J2 = subs(J1,r1,RotAxis);
J3 = subs(J2,t,Translation);
