function [NERRORVECTOR]= ComputeImageErrorsMax(KMatrix, NANGLE, NTRANSLATION, NCORRESPOND, NCONSENSUS)
%This function will generate a 2-coloumn matrix of errors of u and v of 
%each individual point of an actual image.

%The following steps will be exeuted in this function:
%1. Convert the current set of `best' positional parameters to a set of n 3x3
%Perspectivites.
%2. Predict the positions of the points in the image using these transformation
%matrices and the current best camera K-matrix.
%3. Subtract the actual image positions from the predicted image positions to
%generate a set of image errors in pixels.
%4. Store the error into NERRORVECTOR a 2 x nBest

%Counting the number of terms in NCONSENSUS
n = 0;
for j = 1 : length(NCONSENSUS)
    if NCONSENSUS(j) >= 1
       n = n + 1;
    end
end

%Generate a NERRORVECTOR, a 2 x n matrix
NERRORVECTOR = zeros(2*n,1);
for a = 1 : n
    %Extract corresponding point from NCORRESPOND
    b = NCONSENSUS(a);
    XY = NCORRESPOND(3:4, b);
    UV = NCORRESPOND(1:2, b);
    
    %Calculate perspectivies by finding rotaional matrix from Angle-axis
    ROTANG = norm(NANGLE);
    ROTAXIS = NANGLE/norm(NANGLE);
    R = RodriguesRotation(ROTAXIS, ROTANG);
    Perspectivities = [ R(1:3,1:2) NTRANSLATION];
   
    %Calculate predicted u and v and noramlise by value of the third term
    PreUV = KMatrix * Perspectivities * [XY; 1];
    PreUV(1) = PreUV(1) / PreUV(3);
    PreUV(2) = PreUV(2) / PreUV(3);
   
    %Store the difference in NERRORVECTOR
    NERRORVECTOR(2*a-1:2*a) = PreUV(1:2)- UV(1:2);
end

end
    
    