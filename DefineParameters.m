function [Parameters]= DefineParameters
% There must be 8 parameters passed ordered as
% (1) ChipWidth - An integer describing the number of horizontal pixels. 
% (2) ChipHeight - An integer describing the number of vertical pixels. 
% (3) FocalLength - The camera focal length (between 1.0 and 100.0 mm) 
% (4) PixelWidth - The pixel width (between 0.001 and 0.1 mm)
% (5) PixelHeight - The pixel height (between 0.001 and 0.1 mm)
% (6) Skewness - The skewness in u-pixels (between -0.1 and 0.1)
% (7) P_u - The offset to the principal point as a fraction of the width 
% (8) P_v - The offset to the principal point as a fraction of the height

Parameters = zeros(8,1);
Parameters(1)=1000;
Parameters(2)=1000;
Parameters(3)=10;
Parameters(4)=.01;
Parameters(5)=.01;
Parameters(6)=0;
Parameters(7)=.5;
Parameters(8)=.5;