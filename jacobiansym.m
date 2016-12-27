function [J4] = JacobianSym(KMatrixValues,RotValues,...
    TranslationValues,CorrValues,BestConsensus)
 
% BestConsensus = [1,2];
n = length(BestConsensus);

%Parameters come from KMatrix, Rotation axis and Translation vector
K = sym('k',[3,3]);
RotAxis = sym('rot',[1,3]);
% RotAxisNormalized = sym('rotN',[1,3]);
t = sym('t',[3,1]);
Theta = norm(RotAxis);
% RotAxisNormalized = RotAxis/Theta;
Correspond=sym('cor',[4,n]);
% RotAxisNormalized = RotAxis/Theta;
fin = sym('fin',[2*n,1]);

% Grid point coordinate
syms x y

% L is a matrix made of rotaxis parameters. It is used to find the rotation
% matrix
L = [0, -RotAxis(1,3)/Theta, RotAxis(1,2)/Theta;RotAxis(1,3)/Theta,0,-RotAxis(1,1)/Theta;-RotAxis(1,2)/Theta,RotAxis(1,1)/Theta,0];
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
    XY = Correspond(3:4,BestConsensus(n));
    f = K * P * [XY' 1]';
%     J = jacobian(f,[k1_1,k1_2,k1_3,k2_2,k2_3,rotN1_1,rotN1_2,rotN1_3,t1_1,t2_1,t3_1]);
    f = f/f(3);
    fin = [f(1:2,1);fin];
    
end
Variables = [K(1,1),K(1,2),K(1,3),K(2,2),K(2,3),RotAxis(1,1),RotAxis(1,2),RotAxis(1,3),t(1,1),t(2,1),t(3,1)];
J = jacobian(fin,Variables);


% Assigning values to the symbols
% KMatrixValues = [1000,0,500;0,1000,500;0,0,1];
% RotValues = [1 2 1];
% TranslationValues = [1 2 1]';
% CorrValues = [1 2 3 1;1 1 1 1]; 

J1 = subs(J,K,KMatrixValues);
J2 = subs(J1,RotAxis,RotValues');
J3 = subs(J2,t,TranslationValues);
J4 = subs(J3,Correspond,CorrValues);


% Ls = subs(f,{RotAxisNormalized,t,Correspond,K},{RotValues,TranslationValues,CorrValues,KMatrixValues});

