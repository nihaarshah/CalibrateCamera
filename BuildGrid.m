function [grid] = BuildGrid (GridIncrement, GridWidth)
x = 0:GridIncrement:GridWidth;
y = x;
z = zeros(1,(GridWidth/GridIncrement) + 1);
o = ones(1,(GridWidth/GridIncrement) + 1);
grid = [x; y; z; o;]';