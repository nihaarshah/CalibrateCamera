function [KMatJacob, ImgJacob] = SingleImageJacobianIlker(KMatrix,Phi,...
    Trans,Correspond,Consensus)


%Calculate the initial error vector

InitialErrorVector = ComputeImageErrors(KMatrix,Phi,Trans,Correspond,Consensus);

Jacob = zeros(2*nnz(Consensus),11);

p = [KMatrix(1,1:3)';KMatrix(2,2:3)';Phi;Trans];

dp = zeros(11,1);

dp(1) = 1.0e-4;
dp(2) = 1.0e-4;
dp(3) = 1.0e-4;
dp(4) = 1.0e-4;
dp(5) = 1.0e-4;

dp(6:11) = [Phi/100;Trans/100];

for j = 1:length(dp)
    p(j) = p(j) + dp(j);
    
    NewKMatrix = zeros(3);
    NewKMatrix(1,1:3) = p(1:3)';
    NewKMatrix(2,2:3) = p(4:5)';
    NewKMatrix(3,3) = 1;
    
    NewPhi = p(6:8);
    
    NewTrans = p(9:11);
    
    % Re-compute image errors
    ErrorVector = ComputeImageErrors(NewKMatrix,NewPhi,NewTrans,Correspond,Consensus);
    
    Jacob(:,j) = (ErrorVector - InitialErrorVector)./dp(j);
    
    %Undo the perturbation
    p(j) = p(j) - dp(j);
end
KMatJacob = Jacob(:,1:5);
ImgJacob = Jacob(:,6:11);
end
