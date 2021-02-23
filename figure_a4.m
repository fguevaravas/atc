% Appendix Figure 4 - Active Thermal Cloaking and Mimicking
% demonstrates roles of different terms in boundary integral representatio
% when the initial condition is non-zero.
config.k = .3;
config.N1 = 200; config.N2 = 200;
config.t_eval = 0.15;
config.t_pics = [0.01, 0.05,0.1];%these must be at exact time steps
config.cp = 128;
config.N = 300;
thickLines; %settings for plots 
thermal_reproductions_supp_IC(config);
