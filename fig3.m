%Active Thermal Cloaking supplementary images 
config.k = .3;
config.N1 = 200; config.N2 = 200;
config.Y1 = .25; config.Y2 = .25;
config.t_eval = 0.3; %final time for max interval
config.t_min = 0.2; %Time that field will actually be recreated at
config.cp = 128;
config.N = 300;
config.source_time = 0;
config.errors = 0; %1 yes
config.perc = 0.03; 

if floor(config.N*config.t_min/config.t_eval)==config.N*config.t_min/config.t_eval
    thickLines; %settings for plots 
    thermal_reproductions_supp(config);
else
    sprintf('t_min must be a time step')
end