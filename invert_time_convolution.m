% Inverts a causal convolution in time i.e.
%   y(t) = \int_0^t A(t-s)*x(s) ds
% where A(t-s) is a linear operator and x(s) is a vector.
%
% We work on the approximation by the Riemann sum
%
%   y(i\delta t) = \sum_{j=0}^i A((i-j) \delta t) x(j\delta t)
%
% inputs
%   apply_A    is a function s.t. apply_A(t,x) = A(t)*x
%   apply_Ainv is a function s.t. apply_Ainv(t,x) = inv(A(t))*x
%
%   Y Nx(M+1)  matrix of M+1 time snapshots, 
%              i.e. Y = [y(dt/2), y(3dt/2), ... , y(M*dt)]
%
%   dt         time step = T / M
%
% outputs
%   X Nx(M+1)  matrix of M+1 time snapshots of x:
%              i.e. X = [x(dt/2), x(3dt/2), ... , x(M*dt)]
function X=invert_time_convolution(apply_A,apply_Ainv,Y,dt)

% get dimensions and allocate
N = size(Y,1); M = size(Y,2); X = 0*Y;

% we do forward substitution in time
for i=1:M
 tmp = Y(:,i);
 for j=1:(i-1)
  tmp = tmp - dt*apply_A((i-j+1/2)*dt,X(:,j));
 end
 X(:,i) = apply_Ainv(dt/2,tmp)/dt;
end
