function [ErrorVec] = ComputeImageErrors(KMatrix,RotAxis,...
    Translation,Correspond,ConsensusVector)
% This function calculates the difference between [u' v'] estimates and [u
% v] measured. The error is a 2 x 1 column vector. The loop runs n
% times if there are number n non-zero entries in the BestConsensus
% vector. The error calculated in each iteration is stacked in a column
% vector of size 2n x 1.


ConsensusVector = ConsensusVector(ConsensusVector~=0); %find non zero
sz = length(ConsensusVector);

ErrorVec=[];
One = ones(size(Correspond));
Correspond = vertcat(Correspond,One); % Adding a row of 1s in the bottom of 
% Correspond so now [u v x y 1] is one column of 5 x n matrix correspond

[Homography]=CalculateHomog(RotAxis,KMatrix,Translation);



for i = 1:2:(2*sz-1)
    fk=Homography * Correspond(3:5,ConsensusVector((i+1)/2)); % [u' v' s]' or the estimated u v coordinates
    % fk = fk(1:2,(i+1)/2); % [u' v' s]'
    % we normalize
    fk(1) = fk(1)/fk(3);
    fk(2) = fk(2)/fk(3);
        
    % fk = norm(Estimates(i));
    mk = Correspond(1:2,ConsensusVector((i+1)/2)); % [u v]'
    Error(i:i+1,1) = fk(1:2,1) - mk(1:2,1); % error in pixels [du dv]'
    ErrorVec = [ErrorVec; Error(i:i+1,1)]; % stacking [du dv]'
end

end


