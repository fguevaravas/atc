%Double layer potential
%
%Input:
%   f: dirichlet data
%   nv: normal vectors to surface
%   pmc1: x-coordinates of difference between points on surface and grid
%   points
%   pmc2: y-coordinates of difference between points on surface and grid
%   points
%   mode: switch for Green's function for heat kernel or 2d initial
%   condition 
%   h: difference between points on surface
%   t: current time
%   k: thermal diffusivity
%Ouput:
%   U: double-layer component of recreated temperature field at t 

function[U] =dlp(f, nv, pmc1, pmc2, mode, h, t, k)
    if (mode == 'reg')
        %gradient components for heat kernel
        if (t>0)
            DG1 = pmc1./(8*k^2*pi*t^2).*exp(-(pmc1.^2+pmc2.^2)/(4*k*t));
            DG2 = pmc2./(8*k^2*pi*t^2).*exp(-(pmc1.^2+pmc2.^2)/(4*k*t));
        else
            DG1 = 0*pmc1;
            DG2 = 0*pmc1;
        end
    elseif (mode == 'int')
        %gradient components for 2d initial condition
        DG1 = 1/(4*pi)*((2*pmc1./(pmc1.^2+pmc2.^2))-(2*pmc1./(pmc1.^2+pmc2.^2)).*exp(-(pmc1.^2+pmc2.^2)/(4*k*t)));
        DG2 = 1/(4*pi)*((2*pmc2./(pmc1.^2+pmc2.^2))-(2*pmc2./(pmc1.^2+pmc2.^2)).*exp(-(pmc1.^2+pmc2.^2)/(4*k*t)));
    end
    %approximate integral with trapezoidal rule
    U = (DG1.*nv(:,1)'+DG2.*nv(:,2)')*(f.*h);

end%function

