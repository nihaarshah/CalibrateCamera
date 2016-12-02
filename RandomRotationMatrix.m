function [ R ] = RandomRotationMatrix()
%RandomRotationMatrix
% Generates a random rotation matrix by starting with
% two random vectors and then normalising, projecting
% out components and then generating the final
% column by using the cross-product.

% Make sure we have a random matrix larger than 'eps'
xnorm = 0;
while xnorm < eps
    x = rand(3,1);
    xnorm = norm(x);
 end
% The 'random ' unit vector
xhat = x / xnorm;

ynorm =0;
while ynorm < eps
 y = rand(3,1);
 % Project out xhat by finding the dot product and removing from y
y = y - (y'*xhat)*xhat;
ynorm = norm(y);
 end
 % The 'random' unit y vector
 yhat = y / ynorm;

 % Find the third vector as the cross product of x and
 zhat = cross(xhat,yhat);

 % Construct the rotation matrix for return
 R = [xhat yhat zhat];
end