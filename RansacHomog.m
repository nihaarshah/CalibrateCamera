function [Homog , BestConsensus] = RansacHomog(Correspond, MaxError, NRuns)
%RansacHomog runs a ransac estimation of the homography
% passed as pairs of points in Correspond where Correspond
% is a set of 4-vetcors of the form [[u v]' ; [x y]']

% MaxError is the maximum error for a point to be rejected in the consensus
% set
% NRuns is the nymber of times to run the estimator 

% Homof is the homography that has been identified. If no homography can be
% computed then a 3x3 zero matrix is returned

% 1. For each of NRuns, choose as set of 4 points.
% 2. Use the 4 points to generate a regression matric and a data vector
% 3. If the regressor is a full rank, estimate this homography
% 4. Go through the points and accept those whose error norm is less than
% that set by MaxError.
% 5. Record the set of points in the largest consensus set.
% 6. If the consensus set is not empty, carry out the least squares
% estimate of the homography using svd

% The number of points available
n = length(Correspond);

% Allocate space for the homography
Homog = zeros(3);

% The number of points in the best consenus 
nBest = 0;

for Runs = 1:NRuns
    
    RankTest = 1;
    while RankTest == 1
        % RankTest is set to 1 if the 4 points do not give 
        % a full rank matrix in the estimator
        
        RankTest = 0;
        
        % 1 Choose a sample set
        SamplePoints = zeros(1,4);
        
        % The first point is a random integer between 1 and inclusive
        SamplePoints(1) = 1 + fix(n*rand);
        if SamplePoints(1) > n
            SamplePoints(1) = 1;
        end
        
        % Choose the next three points ensuring that there are no repeats
        for j = 2:4
            % searching is a toggle that is triggered if there is a repeat
            searching = 1;
            
            while searching == 1
                searching = 0;
                % Inititial sample point 
                SamplePoints(j) = 1 + fix(n*rand);
                if SamplePoints(j) > n
                    SamplePoints(j) = n;
                end
                
                % Is the new point a repeat? If so go again
                for k = 1:j-1
                    if SamplePoints(k) == SamplePoints(j)
                        searching = 1;
                    end
                end
            end
        end
        
% 6. BestConsenus now contains the largest set of consistent estimates.
% Use this set to estimate the homography using a robust inverse
if nBest > 0
    % The number of measurements in the consenus set 
    Regressor = zeros(2*nBest,1);
    DataVac = zeros(2*nBest,1);
    
    % Construct the regressor
    for j = 1:nBest
        r1 = j*2-1;
        r2 = j*2;
        % HomogRowPair generates 2 rows of the 8 column matrix that
        % multiplies the unknown vector of the homography elements
        % to get the vector of the measurements.
        
        [Regressor(r1:r2,:), DataVec(r1:r2)] =...
            HomogRowPair(Correspond(:,BestConsensus(j)));
        
        % Find the singular value decomposition in order to compute the
        % robust pseudo inverse
        [U, D, V] = svd(Regressor, 'econ');
        
        % The condition number of the computation- a measure of how reiable
        % the inversion is 
        if D(8,8) < eps
            Condition = D(1,1)/D(8,8);
        end
        
        