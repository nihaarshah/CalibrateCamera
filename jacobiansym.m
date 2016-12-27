function [J4,J8] = JacobianSym(KMatrixValues,RotValues,...
    TranslationValues,CorrValues,BestConsensus)
 
% BestConsensus = [1,2];
n = length(BestConsensus);

%Parameters come from KMatrix, Rotation axis and Translation vector
K = sym('k',[3,3]);
RotAxis = sym('rot',[1,3]);
t = sym('t',[3,1]);
Theta = norm(RotAxis);
Correspond=sym('cor',[4,n]);
fin = sym('fin',[2*n,1]);


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
fin =[];
for i = 1:n
    XY = Correspond(3:4,BestConsensus(i));
    f = K * P * [XY' 1]';
%     J = jacobian(f,[k1_1,k1_2,k1_3,k2_2,k2_3,rotN1_1,rotN1_2,rotN1_3,t1_1,t2_1,t3_1]);
    f = f/f(3);
    fin = [fin;f(1:2,1)];
    
end
Variables = [K(1,1),K(1,2),K(1,3),K(2,2),K(2,3),RotAxis(1,1),RotAxis(1,2),RotAxis(1,3),t(1,1),t(2,1),t(3,1)];

JKMat = jacobian(fin,Variables(1,1:5));
JFram = jacobian(fin,Variables(1,6:11));

% Subbing in the K parameter values
J1 = subs(JKMat,K,KMatrixValues);
J2 = subs(J1,RotAxis,RotValues');
J3 = subs(J2,t,TranslationValues);
J4 = subs(J3,Correspond,CorrValues);

J5 = subs(JFram,K,KMatrixValues);
J6 = subs(J5,RotAxis,RotValues');
J7 = subs(J6,t,TranslationValues);
J8 = subs(J7,Correspond,CorrValues);


