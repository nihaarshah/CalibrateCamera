function [grid] = BuildGrid (GridIncrement, GridWidth)
x =[];
y=[];
for i=-(GridWidth/2):GridIncrement:(GridWidth/2)
    for j=-(GridWidth/2):GridIncrement:(GridWidth/2)
        x=[x i];
        y=[y j];
    end
end

z = zeros(1,length(x));
o = ones(1,length(x));
grid = [x; y; z; o;]; 
% output [ x y 0 1
%          x y 0 1
%          x y 0 1]