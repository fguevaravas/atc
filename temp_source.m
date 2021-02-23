%Temperature function for a point source in an unbounded domain
%
%Inputs:
%   X1: x-coordinate of points to evaluate
%   X2: y-coordinate of points to evaluate
%   t: time to evaluate
%   Y1: x-coordinate of point source
%   Y2: y-coordinate of point source
%   stregnth: constant to multiply temperature by 
%   k: thermal diffusivity
%Ouputs:
%   U: temperature at given points
%   UDX1: x-component of gradient of temperature function
%   UDX2: y-component of gradient of temperature function
function [U,DX1,DX2] = temp_source(X1,X2,t,source_time, Y1, Y2, strength, k)
    %Calculate temperature
    if(t>source_time)
        U = strength/(4*pi*k*(t-source_time))*exp(-((X1-Y1).^2+(X2-Y2).^2)/(4*k*(t-source_time)));
    else
        U = X1*0;
    end
    %Calculate gradients
    if(t>source_time)
        DX1 = -strength*(X1-Y1)./(8*pi*k^2*(t-source_time)^2).*exp(-((X1-Y1).^2+(X2-Y2).^2)/(4*k*(t-source_time)));
    else
        DX1 = 0*X1;
    end
    if(t>source_time)
        DX2 = -strength*(X2-Y2)./(8*pi*k^2*(t-source_time)^2).*exp(-((X1-Y1).^2+(X2-Y2).^2)/(4*k*(t-source_time)));
    else
        DX2 = 0*X2;
    end
end%function

