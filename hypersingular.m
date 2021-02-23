% Derivative of double layer potential, hypersingular potential
%
%Input:
%   f: dirichlet data data
%   nv: normal vectors of original object
%   nv_c: normal vectors of where the derivative of the DLP is calculated
%   pmc1: x-coordinates of difference between points on surface and grid
%   points
%   pmc2: y-coordinates of difference between points on surface and grid
%   points
%   h: difference between points on surface
%   t: current time
%   k: thermal diffusivity
%Ouput:
%   U: derivative of double-layer component of recreated temperature field at t 

function[U] = hypersingular(f, nv, nv_c, pmc1, pmc2, h, t, k)
if (t>0)
    D11 = (1-pmc1.^2/(2*k*t))/(8*k^2*pi*t^2).*exp(-(pmc1.^2+pmc2.^2)/(4*k*t));
    D12 = -pmc1.*pmc2./(16*k^3*pi*t^3).*exp(-(pmc1.^2+pmc2.^2)/(4*k*t));
    D21 = D12;
    D22 = (1-pmc2.^2/(2*k*t))/(8*k^2*pi*t^2).*exp(-(pmc1.^2+pmc2.^2)/(4*k*t));
    comp1 = diag(nv_c(:,1))*D11*diag(nv(:,1));
    comp2 = diag(nv_c(:,1))*D21*diag(nv(:,2));
    comp3 = diag(nv_c(:,2))*D12*diag(nv(:,1));
    comp4 = diag(nv_c(:,2))*D22*diag(nv(:,2));
    U = (comp1+comp2+comp3+comp4)*(f.*h);
else
    U=0*pmc1*(f.*h);
end
end%function