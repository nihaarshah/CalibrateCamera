function [Regressor, DataVec] = HomogRowPair(Correspond)



 u = Correspond(1,1);
 v = Correspond(2,1);
 x = Correspond(3,1);
 y = Correspond(4,1);
 
 Regressor = [x y 1 0 0 0 -u*x -u*y; 0 0 0 x y 1 -v*x -v*y];
%  Regressor=vertcat(Regressor,newRegressor);
 
 
 DataVec = [u v]';
%  DataVec = vertcat(DataVec,newDatavec);
            

