function [grid] = BuildGrid (GridIncrement, GridWidth)
x =[];
y=[];
for i=0:GridIncrement:GridWidth
    for j=0:GridIncrement:GridWidth
        x=[x i];
        y=[y j];
    end
end

z = zeros(1,length(x));
o = ones(1,length(x));
grid = [x; y; z; o;]';