function [KMatrixJacobBlock, FrameJacobBlock] = SingleImageJacobian(KMatrix,RotAxis,Translation, Correspond, BestConsensusIndices)
% Output two blocks of matrices. One block is a first derivative of the
% Estimated image coordinates u and v with respect each of the K matrix
% parameters. There are n*2 rows where n is the nymber of elements in the
% BestConsenus vetcor. There are 5 columns each representing the partial
% derivative with respect to one of tge 5 Kmatrix parameters.
% The second block is similar except the parameters are now the 6 parametsr
% of the Rotation (3 axes) and translation(3 axes).
 
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
        dp = ones(11);
        dp(i) = 1+(1e-2);
        NewParameters = Parameter * dp;
        
        % Find new estimates after perturbation
        NewKMatrix=CalculateNewKMatrix(NewParameters);
        [NewHomography] = CalculateHomog(RotAxis,NewKMatrix,Translation);
        NewEstimate = NewHomography * Correspond(3:5,BestConsensusIndices(consindex)); % [u1' v1' s1]
        NewEstimate(1) = NewEstimate(1)/NewEstimate(3);
        NewEstimate(2)= NewEstimate(2)/NewEstimate(3);
        NewEstimate = NewEstimate(1:2,1); % [u1' v1']
        
        % Find derivative
        dfdp = (NewEstimate-InitialEstimate)/(Parameter(i)*(dp(i)-1)); % [du dv]
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
        dp = ones(11);
        dp(j) = 1+(1e-2);
        NewParameters = Parameter * dp;
        
        % Find new estimates after perturbation
        [NewRotAxis,NewTranslation] = CalculateNewRotTrans(NewParameters);
        [NewHomography] = CalculateHomog(NewRotAxis,KMatrix,NewTranslation);
        NewEstimate = NewHomography * Correspond(3:5,BestConsensusIndices(consindex)); % [u1' v1' s1]
        NewEstimate(1) = NewEstimate(1)/NewEstimate(3);
        NewEstimate(2)= NewEstimate(2)/NewEstimate(3);
        NewEstimate = NewEstimate(1:2,1); % [u1' v1']
        
        % Find derivative
        dfdp = (NewEstimate-InitialEstimate)/(Parameter(j)*(dp(j)-1)); % [du dv]
        B2(row:row+1,j-5) = dfdp;
    end
end

KMatrixJacobBlock = B1;
FrameJacobBlock = B2;

end