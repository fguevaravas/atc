%Figure 7 interior reproduction problem Active Thermal Cloak
clear all;
config.k = .2;
config.N1= 100; config.N2=100;
config.N = 1000;
config.Y1 = 0;
config.Y2 = 0;
config.t_eval = 1;

thickLines; %Settings for plotting
sd_study(config);