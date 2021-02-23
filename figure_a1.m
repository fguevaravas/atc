% Appendix - Figure 1 Active Thermal Cloaking and Mimicking
%
% Demonstrates the interior reproduction problem in 1D

%   config.t_eval: time of evaluation
%   config.N: number of time steps
%   config.Y1: location of point source
%   config.cp: location of cloaking devices
%   config.N1: number of spatial discretizations
%   config.k: thermal diffusvity

config.t_eval = .1;
config.N = 3000;
config.Y1 = 0.1;
config.cp = [.4,.8];
config.N1 = 4000;
config.k = 1;
thickLines;
therm_rep_1D(config)