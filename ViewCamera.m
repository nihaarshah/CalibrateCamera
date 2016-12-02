function [] = ViewCamera(ObjectLines, T_ow, KMatrix, CameraHeight, CameraWidth, T_cw)
    %ViewCamera
%Takes an object described by a set of lines passed in ObjectLines and
%draws a picture of the camera's view of the object.

%ObjectLines is a 4x2n were each column is a homogenous Point in the 
%objective frame. The object is defined as pairs of Points that should have
%a line drawn between them.
%T_ow is a 4x4 homogenous transformation matrix describing 
%the object's frame.

%KMatrix is the K-Matrix of the camera in pixels
%CameraHeight is the number of vertical pixels
%CameraWidth is the number of horizontal pixels 
%T_cw is the 4x4 Camera frame in world coordinates

%Check sizes
s = size(ObjectLines);
if s(1) ~=4 || mod (s(2), 2) ~= 0;
    error('ObjectLines has an invalid size')
end

s = size(T_ow);

if s(1) ~= 4|| s(2) ~= 4
    error('T_ow has an invalid size')
end

s= size(KMatrix);
if s(1) ~=3 || s(2) ~= 3
    error('KMatrix has an invalid size')
end

s = size(T_cw);
if s(1) ~= 4 || s(2) ~= 4
    error('T_oc has an invalid size')
end

% We could perform other tests to make the code bomb proof

% Transform the object into world co-ordinates
ObjectLines = T_ow * ObjectLines;

%Transform the object into camera coordinates using the backslash operator
ObjectLines = T_cw \ ObjectLines;

%Project out the 4th coordinate and multiply by the KMatrix
ObjectLines = KMatrix * ObjectLines(1:3, :);

% We now have a set of homogenous points representing 2D points.
% We need to normalize these points to get 2D points

s= size(ObjectLines);
for j = 1:s(2)
    ObjectLines(1:2,j) = ObjectLines(1:2, j)/ObjectLines(3,j);
end

%Throw away the normalizing components
ObjectLines = ObjectLines(1:2,:);

%Generating the image
figure(1)
clf

axis([0 CameraWidth 0 CameraHeight])
hold on 

%Count through the point pair
for j =1 : 2 : s(2)
    plot([ObjectLines(1,j) ObjectLines(1, j+1)],[ObjectLines(2,j) ObjectLines(2, j+1)])
end

axis ij

end

    
    