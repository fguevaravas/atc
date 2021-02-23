%Applies the inverse of the single layer potential evaluated at the first
%time step. 
%
%Inputs:
%   f: one time step of the RHS of the Calderon system 
%   pmc1: x-coordinates of differences of ctr points on boundary of
%   inclusion
%   pmc2: y-coordinates of differences of ctr points on boundary of
%   inclusion
%   h: distance between centers on boundary of inclusion
%   k: thermal diffustivity
%   dt: time step size
function [U] = slp_inverse(f, pmc1, pmc2, h, k, dt)
t = dt/2;
G = (1/(4*k*pi*t)).*exp(-(pmc1.^2+pmc2.^2)/(4*k*t));
U =(G\f).*(1./h);
end %function
