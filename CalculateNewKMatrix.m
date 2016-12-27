function [NewKMatrix] = CalculateNewKMatrix(Parameters)

NewKMatrix(1,1)=Parameters(1);
NewKMatrix(1,2)=Parameters(2);
NewKMatrix(1,3)=Parameters(3);
NewKMatrix(2,1)=0;
NewKMatrix(2,2)=Parameters(4);
NewKMatrix(2,3)=Parameters(5);
NewKMatrix(3,1)=0;
NewKMatrix(3,2)=0;
NewKMatrix(3,3)=1;