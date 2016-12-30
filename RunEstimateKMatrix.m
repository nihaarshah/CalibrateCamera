% This script runs the estimation of the KMatrix %
%
% 1. Constructs a camera model loosely based on an iPhone6
% 2. Constructs a calibration grid 1m on a side with 10mm grid spacing.
% 3. Positions the grid somewhere in space.
%
% Perform the following actions on each image, repeating an image
% if the homography estimation failed.
% 4. Places the camera somewhere in space to generate a full image of tile
% square corner locations.
% 5. Generate the noisy image of the grid.
% 6. Add in some outliers.
% 7. Perform a Ransac estimate of the homography

% Once the homographies have been estimated
% 8. Build the regressor for estimating the K-matrix
% 9. Carry out the Cholesky factorization and invert.

% The number of calibration images to use
nImages = 6;

% 1. Construct the Camera model
[KMatrix , CameraHeight , CameraWidth] = BuildCamera;

% 2. Construct a 1m by 1m grid with 10mm tiles in the grid frame
% The grid is a set of 4-element vectors [x y 0 1]'.
GridWidth = 1000;
GridIncrement = 10;
CalibrationGrid = BuildGrid(GridIncrement ,GridWidth);

% 3. Choose somewhere in space for the grid
% T_ow is the 4x4 tranformation matrix from grid to world.
T_ow = PositionGrid();

% Define the scaling to use
if CameraHeight > CameraWidth
    CameraScale = 2.0 / CameraHeight;
else
    CameraScale = 2.0 / CameraWidth;
end
GridScale = 2.0/GridWidth;
ErrorInEstimationVector = [];
RansacRunsVector = [];
KSum = zeros(3,3);
%This loop if to change the Ransac runs in each iteration and store the
% error between the estimated and real K matrix each time
for RansacRuns = 5:15
% for KMatIterations = 1:10
% Generate the calibration images and the homographies
% Store Homographies and consensus sets in a Matlab Cell Array
% called HomogData.
HomogData = cell(nImages,3);

for CalImage = 1: nImages
    % Keep looking for homographies until we get a non-zero result.
    % 'Estimating ' is a toggle.
    Estimating = 1;
    while Estimating == 1
        % The default is 'Success', i.e. Estimating = 0
        Estimating = 0;
        
        % 4 Choose a 'random' location for the camera that fills the image .
        % T_cw is the 4x4 transformation matrix from camera to world
        T_cw = FillImage(T_ow,KMatrix,GridWidth,CameraHeight,CameraWidth);
        % 5 We now fill the camera with a noisy image of the grid and generate
        % the point correspondences.
        
        % Correspond is a set of pairs of vectors of the form [[u v]' [x y ]']
        % for each grid corner that lies inside the image.
        mu = 0;
        std = 0.3;
        Correspond = BuildNoisyCorrespondences(T_ow,T_cw,CalibrationGrid,...
            KMatrix,CameraHeight,CameraWidth,mu,std);
        
        % 6. Add in some 'outliers' by replacing [u v]' with a point
        % somewhere in the image.
        % Define the Outlier probability
        pOutlier = 0.10;
        for j = 1:length(Correspond)
            r = rand;
            
            if r < pOutlier
                
                Correspond(1,j) = rand * (CameraWidth-1);
                Correspond(2,j) = rand * (CameraHeight-1);
            end
        end
        
        % Now scale the grid and camera to [-1, 1] to improve
        % the conditioning of the Homography estimation.
        Correspond(1:2,:) = Correspond(1:2,:)*CameraScale - 1.0;
        Correspond(3:4,:) = Correspond(3:4,:)*GridScale;
        
        
        
        % 7. Perform the Ransac estimation - output the result for inspection
        % If the Ransac fails it retuns a zero Homography
        Maxerror = 3; % The maximum error allowed before rejecting a point.
        
        % I am using a variance of 0.5 pixels so sigma is sqrt(0.5)
        % 3 pixels in the *NORM* is 3 sigma as there are 2 errors involved
        % in the norm (u and v).
        % Note: The above is in pixels - so scale before Ransac! 
        
%         RansacRuns = 50;
        % The number of runs when creating the consensus set.
        
        [Homog,BestConsensus]=RansacHomog(Correspond,Maxerror*CameraScale,RansacRuns);
        
        
        if Homog(3,3) > 0
            % This image worked. So record the homography and the
            % consensus set
            HomogData{CalImage,1} =  Homog;
            HomogData{CalImage,2} = Correspond;
            HomogData{CalImage,3} = BestConsensus;
            
        else
            % estimation failed so go around again
            Estimating = 1;
        end
        
        
        
    end % end of the whole Estimating ==1 loop
end % end of the nImages loop


% 8. Build the regressor for estimating the Cholesky product
Regressor = zeros(2*nImages,6);

for CalImage = 1:nImages
    r1 = 2* CalImage -1;
    r2 = 2*CalImage;
    Regressor(r1:r2,:) = KMatrixRowPair(HomogData{CalImage,1});
end

% Find the kernel
[U, D, V] = svd(Regressor,'econ');
D = diag(D);
[M, I] = min(D);
% K is the estimate of the kernel
K = V(:,I);

% The matrix to be constructed needs to be positive definite
% It is necessary that K(1) be positive.
if K(1) < 0
    K = -K;
end

% Construct the matrix Phi from the kernel
Phi = zeros(3);

Phi(1,1) = K(1);
Phi(1,2) = K(2);
Phi(1,3) = K(3);
Phi(2,2) = K(4);
Phi(2,3) = K(5);
Phi(3,3) = K(6);

% Add in the symmetric components
Phi(2,1) = Phi(1,2);
Phi(3,1) = Phi(1,3);
Phi(3,2) = Phi(2,3);

% Check if the matrix is positive definite
e = eig(Phi);
for j = 1:3
    if e(j)<= 0
        error('The Cholesky product is not positive definite')
    end
end

% 9. Carry out the Cholesky factorization
KMatEstimated = chol(Phi);
% Invert the factor
KMatEstimated = KMatEstimated \ eye(3);

% The scaling of the grid has no impact on the scaling of
% the K-matrix as the vector 't' takes no part in the estimate
% of Phi. Only the image scaling has an impact.

% first normalize the KMatrix
KMatEstimated = KMatEstimated / KMatEstimated(3,3);

% Add 1.0 to the translation part of the image
KMatEstimated(1,3) = KMatEstimated(1,3) + 1;
KMatEstimated(2,3) = KMatEstimated(2,3) + 1;

% Rescale back to pixels
KMatEstimated(1:2,1:3) = KMatEstimated(1:2,1:3) / CameraScale;



% KSum = KMatEstimated+KSum;
% 
% end % end the loop for K iterations
% KEstimatedAverage = KSum/10;

% error between the Kestimated and actual K
ErrorInEstimation = norm(KMatrix-KMatEstimated,1);

ErrorInEstimationVector = [ErrorInEstimationVector,ErrorInEstimation];
RansacRunsVector = [RansacRunsVector,RansacRuns];
end % end of loop that changes ransac runs

plot(RansacRunsVector,ErrorInEstimationVector,'black*');
ylabel('||K estimated - Original K.||');
xlabel('number of RanSac runs')


BestConsensus = BestConsensus(BestConsensus~=0);
X=['The best consensus for Maxerror of ',num2str(Maxerror),' has ',num2str(length(BestConsensus)),' elements in it.'];
BestConsensus;
disp(X);
