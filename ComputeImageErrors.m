function [ErrorVec] = ComputeImageErrors(KMatrix,RotAxis,...
    Translation,Correspond,ConsensusVector)

ConsensusVector = ConsensusVector(ConsensusVector~=0); %find non zero
% Theta = norm(RotAxis); % Maybe the axis of rotation is theta * unit
%     % length axis of rotation so taking the norm will give us the magnitude
%     % of this axis which is theta?
% n = 4;
% x = RotAxis(1)/Theta;
% y = RotAxis(2)/Theta;
% z = RotAxis(3)/Theta;
%
% K =[0 -z y; z 0 -x; -y x 0]; %Different K matrix
%
% % Theta = Phi/RotAxis - n*pi; What i thought earlier
%
% R = [1 0 0;0 1 0; 0 0 1] + sin(Theta)*K + (1-cos(Theta))*[K]^2;
%
% P = [R Translation; [0 0 0] 1];
% Perspectivity = [P(1:3,1:2),P(1:3,4)];
% One = ones(size(Correspond));
% Correspond = vertcat(Correspond,One);
% Homography = KMatrix * Perspectivity;

ErrorVec=[];
One = ones(size(Correspond));
Correspond = vertcat(Correspond,One);

[Homography]=CalculateHomog(RotAxis,KMatrix,Translation);
sz = length(ConsensusVector);
% fk =zeros(1:3,1);
% mk = zeros(1:2,1);

for i = 1:2:(2*sz-1)
    fk=Homography * Correspond(3:5,ConsensusVector((i+1)/2)); % [u' v' s]' or the estimated u v coordinates
%     fk = fk(1:2,(i+1)/2); % [u' v' s]'
    % shall we normalize???
    fk(1) = fk(1)/fk(3);
    fk(2) = fk(2)/fk(3);
        
    % fk = norm(Estimates(i));
    mk = Correspond(1:2,ConsensusVector((i+1)/2)); % [u v]'
    Error(i:i+1,1) = fk(1:2,1) - mk(1:2,1); % error in pixels [du dv]'
    ErrorVec = [ErrorVec; Error(i:i+1,1)]; % stacking [du dv]'
end

end


