%1D thermal reproductions for a point source
%
%Inputs:
%   config.t_eval: time of evaluation
%   config.N: number of time steps
%   config.Y1: location of point source
%   config.cp: location of cloaking devices
%   config.N1: number of spatial discretizations
%   config.k: thermal diffusvity
%Ouputs:
%   figure(1): original temperature & recreated temperature
%   figure(2): log10 error 
function[] = therm_rep_1D(config)
L1 = 1; 
dt = config.t_eval/config.N;

nv = [1,-1;-1,1]; %interior rep

%Spatial discretization
X1 = linspace(0,L1,config.N1);

%geometry
geo.cp = config.cp;
geo.nv = nv;

%operators
temp = @(X,t) tempsrc1d(config.k,t,X,config.Y1);
dtemp = @(X,t) dertempsrc1d(config.k,t,X,config.Y1);
single_layer = @(t, x) slp_1d(x,t,X1,config.k,geo);
double_layer = @(t, x) dlp_1d(x, t,X1,config.k,geo);

%Get Data for Convolutions
lambda = zeros(2,config.N);
phi = zeros(2,config.N);
for i = 1:config.N
   lambda(1,i) = temp(geo.cp(1),(i-1/2)*dt);
   lambda(2,i) = temp(geo.cp(2),(i-1/2)*dt);
   phi(1,i) = dtemp(geo.cp(1),(i-1/2)*dt);
   phi(2,i) = dtemp(geo.cp(2),(i-1/2)*dt);
end

%Apply SLP & Dlp
U1 = apply_time_convolution(single_layer, phi, dt);
U2 = apply_time_convolution(double_layer, lambda, dt);
if config.Y1 - config.cp(1) <0 ||config.Y1-config.cp(2)>0
    U_rec_all_times = U1-U2; %interior rep 
else
    U_rec_all_times = U2-U1; %exterior rep 
end
U_rec = U_rec_all_times(:,config.N);
U_rec = reshape(U_rec,size(X1));

%Calculate errors
if config.Y1-config.cp(1) <0||config.Y1-config.cp(2)>0
    ext_err = ((X1<geo.cp(1))+(X1>geo.cp(2))).*log10(abs(U_rec));
    int_err = ((X1>geo.cp(1)).*(X1<geo.cp(2))).*log10(abs(U_rec-temp(X1,config.t_eval)));
else
    ext_err = ((X1<geo.cp(1))+(X1>geo.cp(2))).*log10(abs(U_rec-temp(X1,config.t_eval)));
    int_err = ((X1>geo.cp(1)).*(X1<geo.cp(2))).*log10(abs(U_rec));
end
tot_err = int_err + ext_err;
%Plot
figure(1); clf;
hold on;
plot(X1,temp(X1,config.t_eval))
plot(X1,U_rec)
a = ylim;
plot(config.Y1,a(1)+.001,'rp', 'MarkerSize', 15,'MarkerFaceColor', 'r')
set(gca,'XTick', unique([[config.Y1,geo.cp(1),geo.cp(2)],get(gca,'XTick')]));
xline(geo.cp(1),'--');
xline(geo.cp(2),'--');
hold off

figure(2); clf;
hold on;
plot(X1,tot_err)
a = ylim;
plot(config.Y1,a(1)+.001,'rp', 'MarkerSize', 15,'MarkerFaceColor', 'r')
set(gca,'XTick', unique([config.Y1,geo.cp(1),geo.cp(2),get(gca,'XTick')]));
xline(geo.cp(1),'--');
xline(geo.cp(2),'--');
hold off

%%%%%%%%%%%%%%%%%%%%%%%1D Functions 
%1D Single Layer Potential
%  
%Inputs:
%   f: Neumann Data
%   t: time
%   X1: points on line
%   k: thermal diffusvity
%   geo: geometry of cloaking points
%Ouput:
%   U: single layer component of recreated field 
function[U] = slp_1d(f,t,X1,k,geo)
    pmc1 = X1-geo.cp(1);
    pmc2 = X1-geo.cp(2);
    G(:,1) = tempsrc1d(k,t,pmc1,0);
    G(:,2) = tempsrc1d(k,t,pmc2,0);
    U = (G)*(f.*geo.nv(:,1));
end%function

%1D Double Layer Potential
%  
%Inputs:
%   f: Dirichlet Data
%   t: time
%   X1: points on line
%   k: thermal diffusvity
%   geo: geometry of cloaking points
%Ouputs:
%   U: double layer component of recreated field 
function[U] = dlp_1d(f, t,X1,k,geo)
    pmc1 = X1-geo.cp(1);
    pmc2 = X1-geo.cp(2);
    G(:,1) = dertempsrc1d(k,t,pmc1,0);
    G(:,2) = dertempsrc1d(k,t,pmc2,0);
    U = (G.*geo.nv(:,2)')*(f);
end%function
%1D point source temperature
%   
%Inputs: 
%   k: thermal diffusivity
%   t: time
%   X: points on line
%   Y: location of source
%Ouputs:
%   U: temperature
function[U] = tempsrc1d(k,t,X,Y)
    if t ==0
        U = X*0;
    else
        U = 1/sqrt(4*pi*k*t)*exp(-((X-Y).^2)/(4*k*t));
    end
end%function

%Derivative of 1D temperature from point source
%Inputs: 
%   k: thermal diffusivity
%   t: time
%   X: points on line
%   Y: location of source
%Ouput:
%   U: derivative
function[U] = dertempsrc1d(k,t,X,Y)
    if t ==0
        U = X*0;
    else
        U = (X-Y)/2/k/t/sqrt(4*pi*k*t).*exp(-((X-Y).^2)/(4*k*t));
    end
end%function
end%function
