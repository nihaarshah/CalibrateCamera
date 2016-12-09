function [psi] = KMatrixRowPair(Homog)

h11 = Homog(1,1);
h12 = Homog(1,2);
h13 = Homog(1,3);
h21 = Homog(2,1);
h22 = Homog(2,2);
h23 = Homog(2,3);
h31 = Homog(3,1);
h32 = Homog(3,2);
h33 = Homog(3,3);

% this is also the regressor matrix for the jth image
psi = ...
[(h11^2-h12^2) 2*(h11*h21-h12*h22) 2*(h11*h31-h12*h32) (h21^2-h22^2)...
 2*(h21*h31-h22*h32) (h31^2-h32^2);
 (h11*h12) (h11*h22+h21*h12) (h11*h32+h31*h12) (h21*h22) (h21*h32+h31*h22)...
 (h31*h32)];

