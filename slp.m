%Single layer potential
%
%Input:
%   f: neumann data
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
%   U: single-layer component of recreated temperature field at t 
function [U] = slp(f, pmc1, pmc2, mode, h , t,  k)
    %Heat kernel
    if (mode == 'reg')
        if (t>0)
          G = 1/(4*k*pi*t).*exp(-(pmc1.^2+pmc2.^2)/(4*k*t));
        else
          G = 0*pmc1;
        end
    %Function for 2d initial conditions
    elseif(mode == 'int')
        G = -1/(4*pi)*(expint((pmc1.^2 + pmc2.^2)/(4*k*t)) +log((pmc1.^2+pmc2.^2)/(4*k*t)) + 0.577215664901532);
    end

    %approximate integral with trapezoidal rule
    U =G*(f.*h);

end %function
