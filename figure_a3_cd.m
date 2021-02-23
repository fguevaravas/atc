% Appendix Figure 3 - Active Thermal Cloaking and Mimicking
%  computes and plots errors for the exterior reproduction problem
%  with different number of time steps

%   config.t_eval: time of evaluation
%   config.N: number of time steps
%   config.Y1: location of point source
%   config.cp: location of cloaking devices
%   config.N1: number of spatial discretizations
%   config.k: thermal diffusvity

config.t_eval = 1;
config.N = [300, 3000];
config.Y1 = 0.6;
config.cp = [.4,.8];
config.N1 = 1000;
config.k = 1;
thickLines;
therm_rep_1D_Error_study(config)