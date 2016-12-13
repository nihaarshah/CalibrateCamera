function [KMatrixJacobBlock, FrameJacobBlock] = SingleImageJacobian(KMatrix,anglevector,translationvector, Correspond, BestConsensusIndices)

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
 

KParameter(1) = KMatrix(1,1);
KParameter(2) = KMatrix(1,2);
KParameter(3) = KMatrix(1,3);
KParameter(4) = KMatrix(2,2);
KParameter(5) = KMatrix(2,3);

nBestConsensus = size(BestConsensusIndices);

for block1row = 1:nBestConsensus
    %Estimate image location for that particular index 
    Estimates = Homography * Correspond(3:4,BestConsensusIndices(block1row));
     for block1col = 1:5             
        B1(block1row,block1col)=ForwardDifference(Estimates,KParameter(block1col));
    end 
end

RotParameter(1) = anglevector(1);
RotParameter(2) = anglevector(2);
RotParameter(3) = anglevector(3);
RotParameter(4) = translationvector(1);
RotParameter(5) = translationvector(2);
RotParameter(6) = translationvector(3);

    
for block2row = 1:nBestConsensus
    Estimates = Homography * Correspond(3:4,BestConsensusIndices(block2row));
    for block2col = 1:6
        B2(block2row,block2col)=ForwardDifference(Estimates,RotParameter(block2col));
    end
end

KMatrixJacobBlock = B1;
FrameJacobBlock = B2;

end