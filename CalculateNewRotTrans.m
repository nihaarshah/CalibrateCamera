function [NewRotAxis,NewTranslation] = CalculateNewRotTrans(NewParameters)
% These parameters are perturbed one at a time and the new translation and
% rotation vectors use the perturbed values as output
NewRotAxis(1)=NewParameters(6);
NewRotAxis(2)=NewParameters(7);
NewRotAxis(3)=NewParameters(8);
NewTranslation(1)=NewParameters(9);
NewTranslation(2)=NewParameters(10);
NewTranslation(3)=NewParameters(11);
end