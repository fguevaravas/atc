%Figure 9 interior reproduction problem Active Thermal Cloak
clear all;
config.k=.2; 
config.t_eval = 1; 
config.Y1= 0; config.Y2= 0; 
config.N1= 100; config.N2=100;
config.cp = 100;
config.N = 1000;
config.perc = .03;

thickLines; % Plotting Settings
data_study(config);