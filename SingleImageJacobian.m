function [KMatrixJacobBlock, FrameJacobBlock] = SingleImageJacobian(KMatrix,RotAxis,Translation, Correspond, BestConsensusIndices)
% Output two blocks of matrices. One block is a first derivative of the
% Estimated image coordinates u and v with respect each of the K matrix
% parameters. There are n*2 rows where n is the nymber of elements in the
% BestConsenus vetcor. There are 5 columns each representing the partial
% derivative with respect to one of tge 5 Kmatrix parameters.
% The second block is similar except the parameters are now the 6 parametsr
% of the Rotation (3 axes) and translation(3 axes).


% Theta = norm(RotAxis); % Maybe the axis of rotation is theta * unit 
%     % length axis of rotation so taking the norm will give us the magnitude
%     % of this axis which is theta?
% n = 4;
% x = RotAxis(1)/Theta;
% y = RotAxis(2)/Theta;
% z = RotAxis(3)/Theta;
%  
% K =[0 -z y; z 0 -x; -y x 0]; %Different K matrix 
 
% Theta = Phi/RotAxis - n*pi; What i thought earlier
 
% R = [1 0 0;0 1 0; 0 0 1] + sin(Theta)*K + (1-cos(Theta))*[K]^2;
%  
% Perspectivity = [R Translation; [0 0 0] 1];
%  
% Homography = KMatrix * Perspectivity;
 
Parameter(1) = KMatrix(1,1);
Parameter(2) = KMatrix(1,2);
Parameter(3) = KMatrix(1,3);
Parameter(4) = KMatrix(2,2);
Parameter(5) = KMatrix(2,3);
Parameter(6) = RotAxis(1);
Parameter(7) = RotAxis(2);
Parameter(8) = RotAxis(3);
Parameter(9) = Translation(1);
Parameter(10) = Translation(2);
Parameter(11) = Translation(3);

BestConsensusIndices = BestConsensusIndices(BestConsensusIndices~=0); %find non zero

nBestConsensus = length(BestConsensusIndices);
One = ones(length(Correspond));
Correspond = vertcat(Correspond,One);
[Homography] = CalculateHomog(RotAxis,KMatrix,Translation); % Calculating the Homography for finding initial estimate

B1 = [];
for row = 1:2:((2*nBestConsensus)-1)
    consindex = (row+1)/2;
    for i = 1:5
        % Calculate initial estimate fk(p)
        InitialEstimate = Homography * Correspond(3:5,BestConsensusIndices(consindex));% [u0' v0' s0]
        % Need to normalize
        InitialEstimate(1) = InitialEstimate(1)/InitialEstimate(3);
        InitialEstimate(2)= InitialEstimate(2)/InitialEstimate(3);
        InitialEstimate = InitialEstimate(1:2,1); % [u0' v0']
        
        % Perturb K parameters
        dp = zeros(11);
        dp(i) = dp(i) + 1.0e-1;
        NewParameters = Parameter + dp;
        
        % Find new estimates after perturbation
        NewKMatrix=CalculateNewKMatrix(NewParameters);
        [NewHomography] = CalculateHomog(RotAxis,NewKMatrix,Translation);
        NewEstimate = NewHomography * Correspond(3:5,BestConsensusIndices(consindex)); % [u1' v1' s1]
        NewEstimate(1) = NewEstimate(1)/NewEstimate(3);
        NewEstimate(2)= NewEstimate(2)/NewEstimate(3);
        NewEstimate = NewEstimate(1:2,1); % [u1' v1']
        
        % Find derivative
        dfdp = (NewEstimate-InitialEstimate)/dp(i); % [du dv]
        B1(row:row+1,i) = dfdp;
    end
    
end


B2 = [];
for row = 1:2:((2*nBestConsensus)-1)
    consindex = (row+1)/2;
    for j = 6:11
        % Calculate initial estimate fk(p)
        InitialEstimate = Homography * Correspond(3:5,BestConsensusIndices(consindex));% [u0' v0' s0]
        InitialEstimate(1) = InitialEstimate(1)/InitialEstimate(3);
        InitialEstimate(2)= InitialEstimate(2)/InitialEstimate(3);
        InitialEstimate = InitialEstimate(1:2,1); % [u0' v0']
        
        % Perturb K parameters
        dp = zeros(11);
        dp(j) = dp(j) + 1.0e-1;
        NewParameters = Parameter + dp;
        
        % Find new estimates after perturbation
        [NewRotAxis,NewTranslation] = CalculateNewRotTrans(NewParameters);
        [NewHomography] = CalculateHomog(NewRotAxis,KMatrix,NewTranslation);
        NewEstimate = NewHomography * Correspond(3:5,BestConsensusIndices(consindex)); % [u1' v1' s1]
        NewEstimate(1) = NewEstimate(1)/NewEstimate(3);
        NewEstimate(2)= NewEstimate(2)/NewEstimate(3);
        NewEstimate = NewEstimate(1:2,1); % [u1' v1']
        
        % Find derivative
        dfdp = (NewEstimate-InitialEstimate)/dp(j); % [du dv]
        B2(row:row+1,j-5) = dfdp;
    end
end

% for block1row = 1:nBestConsensus
%     %Estimate image location for that particular index 
%     Estimates = Homography * Correspond(3:4,BestConsensusIndices(block1row));
%      for block1col = 1:5             
%         B1(block1row,block1col)=ForwardDifference(Estimates,KParameter(block1col));
%     end 
% end
% 
% 
%     
% for block2row = 1:nBestConsensus
%     Estimates = Homography * Correspond(3:4,BestConsensusIndices(block2row));
%     for block2col = 1:6
%         B2(block2row,block2col)=ForwardDifference(Estimates,RotParameter(block2col));
%     end
% end

KMatrixJacobBlock = B1;
FrameJacobBlock = B2;

end