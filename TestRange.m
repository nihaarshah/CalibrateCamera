function [] = TestRange(Parameter , MinVal , MaxVal , Name)
% TestRange
% Tests the range of Parameter agains MinVal and MaxVal.
% If the test fails an error message results quoting the Parameter % as an identifier and Matlab exits.
if Parameter < MinVal || Parameter > MaxVal
error('Input parameter %s, value %d, was out of range', ...
Name,Parameter)
end
end