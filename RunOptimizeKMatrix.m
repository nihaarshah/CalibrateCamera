% This script performs a full estimation and optimisation of
% a camera K-matrix. %
% 1. Construct a camera model loosely based on an iPhone6
% 2. Construct a calibration grid 1m on a side with 10mm grid spacing.
% 3. Position the grid somewhere in space. %
% Perform the following actions on each image, repeating an image
% if the homography estimation failed.
%
% 4. Place the camera somewhere in space to generate a full image of tile
% square corner locations.
% 5. Generate the noisy image of the grid.
% 6. Add in some outliers.
% 7. Perform a Ransac estimate of the homography
%
% Once the homographies have been estimated
% 8. Build the regressor for estimating the K-matrix
% 9. Carry out the Cholesky factorization and invert. This generates an
% initial model of the K-Matrix


% The number of calibration images to use
nImages = 6;

% 1. Construct the Camera model
[KMatrix , CameraHeight , CameraWidth] = BuildCamera;

% 2. Construct a 1m by 1m grid with 10mm tiles in the grid frame
% The grid is a set of 4-element vectors [x y 0 1]'.
GridWidth = 100;
GridIncrement = 10;
CalibrationGrid = BuildGrid(GridIncrement,GridWidth);


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
% Generate the calibration images and the homographies
% Store Homographies and consensus sets in a Matlab Cell Array
% called HomogData.

HomogData = cell(nImages ,3);

% Define the numbers to access the DataCell
NHOMOGRAPHY = 1;
NCORRESPOND = 2;
NCONSENSUS = 3;

for CalImage = 1: nImages
    % Keep looking for homographies until we get a non-zero result.
    % 'Estimating ' is a toggle.
    Estimating = 1;
    
    
    while Estimating == 1
        % The default is 'Success', i.e. Estimating = 0
        Estimating = 0;
        % 4 Choose a 'random' location for the camera that fills the
        % image.
        % T_cw is the 4x4 transformation matrix from camera to world
        T_cw = FillImage(T_ow,KMatrix,GridWidth,CameraHeight,CameraWidth);
        
        
        % 5 We now fill the camera with a noisy image of the grid and
        % generate the point correpondences.
        % Correspond is a set of pairs of vectors of the form
        % [[u v]' [x y]'] for each grid corner that lies inside the
        % image.
        mu = 0;
        std = 0;
        Correspond = BuildNoisyCorrespondences(T_ow,T_cw,...
            CalibrationGrid ,KMatrix ,CameraHeight ,CameraWidth,mu,std);
        
        
        % 6. Add in some 'outliers' by replacing [u v]' with a point
        % somewhere in the image.
        % Define the Outlier probability
        pOutlier = 0.00;
        
        for j = 1:length(Correspond)
            r = rand;
            if r < pOutlier
                Correspond(1,j) = rand * (CameraWidth - 1);
                Correspond(2,j) = rand * (CameraHeight - 1);
            end
        end
        
        
        % Now scale the grid and camera to [-1, 1] to improve
        % the conditioning of the Homography estimation.
        Correspond(1:2,:) = Correspond(1:2,:)*CameraScale - 1.0;
        Correspond(3:4,:) = Correspond(3:4,:)*GridScale;
        
        % 7. Perform the Ransac estimation
        % If the Ransac fails it retuns a zero Homography
        
        % The maximum error allowed before rejecting a point.
        Maxerror = 3.0;
        % I am using a variance of 0.5 pixels so sigma is sqrt(0.5)
        % 3 pixels in the *NORM* is 3 sigma as there are 2 errors involved
        % in the norm (u and v).
        % Note: The above is in pixels - so scale before Ransac!
        % The number of runs when creating the consensus set.
        RansacRuns = 50;
        
        [Homog, BestConsensus] = ...
            RansacHomog(Correspond,Maxerror*CameraScale,RansacRuns);
        
        BestConsensus=BestConsensus(BestConsensus~=0);
        
        if Homog (3 ,3) > 0
            % This image worked. So record the homography and the
            % consensus set
            HomogData{CalImage ,NHOMOGRAPHY} = Homog;
            HomogData{CalImage ,NCORRESPOND} = Correspond;
            HomogData{CalImage ,NCONSENSUS} = BestConsensus;
        else
            % The estimate failed. So go around again.
            Estimating = 1;
        end
        
    end %end of the while Estimating ==1 loop
end % end of the nImages loop

% 8. Build the regressor for the estimating the Cholesky product
Regressor = zeros(2*nImages,6);
for CalImage = 1:nImages
    r1 = 2*CalImage-1;
    r2 = 2*CalImage;
    Regressor(r1:r2,:) = KMatrixRowPair(HomogData{CalImage,1});
end

% Find the kernel
[U, D, V] = svd(Regressor, 'econ');
D = diag(D);
[M, I] = min(D);
% K is the estimate of the kernel
K = V(:,I);

% The matrix to be constriucted needs to be positive definite
% It is necessarty that K(1) be positive

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
    if e(j) <= 0
        error('The cholesky product is not positive definite')
    end
end

% 9. Carry out the Cholesky factorization
KMatEstimated = chol(Phi);

% Invert the factor
KMatEstimated = KMatEstimated\eye(3);

% The scaling of the grid has no impact on the scaling of the KMatrix as
% the vector 't' takes no part in the estimate of Phi. Only the image
% scaling has an impact.

% first normalize the KMatrix
KMatEstimated = KMatEstimated/KMatEstimated(3,3);

% Optimize the K matrix
OptKMatrix = OptimiseKMatrix(KMatEstimated, HomogData);

% Add 1.0 to trhe translation part of the image
KMatEstimated(1,3) = KMatEstimated(1,3) + 1;
KMatEstimated(2,3) = KMatEstimated(2,3) + 1;

% Rescale back to pixels
KMatEstimated(1:2,1:3) = KMatEstimated(1:2,1:3)/CameraScale;

% Add 1.0 to trhe translation part of the image
OptKMatrix(1,3) = OptKMatrix(1,3) + 1;
OptKMatrix(2,3) = OptKMatrix(2,3) + 1;

% Rescale back to pixels
OptKMatrix(1:2,1:3) = OptKMatrix(1:2,1:3)/CameraScale;

