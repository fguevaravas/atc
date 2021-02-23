%Approximate L2 norm of temperature
%
%Input:
%   U: temperature field
%   N1: number of x grid points
%   N2: number of y grid points
%Output:
%   x: approximate L2 norm of data
function[x] = nap(U,N1,N2) 
if N1 == 1 |N2 ==1
    x = sqrt(sum(U.^2)/N1/N2);
else
    x = sqrt(sum(sum(U.^2))/N1/N2);
end
end %function
