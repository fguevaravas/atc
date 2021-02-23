%Code to generate figure 7 exterior reproduction problem Active Thermal Cloak
clear all;
config.k = .2;
config.N1= 100; config.N2=100;
config.N = 1000;
config.Y1 = 0.5;
config.Y2 = 0.55;
config.t_eval = 1;

thickLines;
sd_study(config);