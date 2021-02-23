%Temperature function for X1
%
%Inputs:
%   X1: x-coordinate of points to evaluate
%Ouputs:
%   U: temperature at given points
%   UDX1: x-component of gradient of temperature function
%   UDX2: y-component of gradient of temperature function
function [U,DX1,DX2] = temp_int(X1,X2)
    %Calculate temperature
    U = X1;
    %Calculate gradients
    DX1 = 1;
    DX2 = 0*X2;
end%function