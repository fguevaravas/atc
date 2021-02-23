% Applies a causal convolution in time i.e.
%   y(t) = \int_0^t A(t-s)*x(s) ds
% where A(t-s) is a linear operator and x(s) is a vector. We use Riemann sum:
%
%   y(i\delta t) = \sum_{j=0}^i A((i-j-1/2) \delta t) x((j+1/2)\delta t)
%
% inputs
%   apply_A    is a function such that apply_A(t,x) = A(t)*x
%
%   X Nx(M+1)  matrix of M+1 time snapshots, 
%              i.e. X = [x(0), x(dt), ... , x(M*dt)]
%
%   dt         time step = T / M
%
% outputs
%   Y Nx(M+1)  matrix of M+1 time snapshots of y:
%              i.e. Y = [y(0), y(dt), ... , y(M*dt)]

function Y=apply_time_convolution(apply_A,X,dt)

% get dimensions and allocate
N = size(X,1); M = size(X,2)-1;
size_ouput = length(apply_A(dt/2, X(:,1)));
Y = zeros(size_ouput, M+1);

for i=1:(M+1)
 % to compute y((i-1)*dt) we need the kernel evaluated at (i-j)*dt
 for j=1:i
   Y(:,i) = Y(:,i) +  dt*apply_A((i-j+1/2)*dt,X(:,j));
 end
end
end%function
