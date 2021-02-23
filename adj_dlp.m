% Derivative of the SLP
%
%Input:
%   f: neumann data
%   pmc1: x-coordinates of difference between points on surface and grid
%   points
%   pmc2: y-coordinates of difference between points on surface and grid
%   points
%   h: difference between points on surface
%   t: current time
%   k: thermal diffusivity
%Ouput:
%   U: derivative of single-layer component of recreated temperature field at t 

function[U] =adj_dlp(f, nv, pmc1, pmc2, h, t, k)
if (t>0)
    DG1 = -pmc1./(8*k^2*pi*t^2).*exp(-(pmc1.^2+pmc2.^2)/(4*k*t));
    DG2 = -pmc2./(8*k^2*pi*t^2).*exp(-(pmc1.^2+pmc2.^2)/(4*k*t));
else
    DG1 = 0*pmc1;
    DG2 = 0*pmc1;
end

% calculate integral
U = (diag(nv(:,1))*DG1+diag(nv(:,2))*DG2)*(f.*h);
end%function