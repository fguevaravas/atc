% Figure 11 Active thermal Cloak
%
% Make a point source inside of a circular region appear in a different
% position from the perspective outside of that region

clear all;

L1 = 1; L2 = 1; % dimensions of domain
k=.2; % thermal diffusivity

t_eval = .2; % where to evaluate the surface integral time
N=250; % number discretization steps in time
dt = t_eval / N; 


Y1= .6; Y2= .4; % position of point source to cloak
Ymove1 = .39; Ymove2 = .6; %position to "move" source to
%Y1 = .5; Y2=.4;
%Ymove1 = .5; Ymove2=.6;

N1 = 200; N2 = 200; % number of grid points in x and y

% Points at which to evaluate
[X1,X2] = ndgrid(linspace(0,L1,N1),linspace(0,L2,N2));
p = [X1(:),X2(:)];

cp = 128; % number of points on boundary of circ region
Rcenter = [L1/2,L2/2]; radius = min (L1,L2)/4; % center and radius of region to be cloaked

geo = circ(cp, Rcenter, radius); % get circle geometry

% Differences between grid points and center of cloaking region
[pmc1, pmc2] = point_diff(p, geo.ctr);

temp = @(X1,X2,T,mode) temp_source(X1,X2,T,0, Y1, Y2, 1, k);  % temperature of original source
temp2 = @(X1,X2,T,mode) temp_source(X1,X2,T,0, Ymove1, Ymove2, 1, k); % temperature of moved source

% Single and Double Layer operators
double_layer = @(t,x ) dlp(x, geo.nv, pmc1, pmc2, 'reg', geo.h, t, k)*k;
single_layer = @(t, x) slp(x, pmc1, pmc2, 'reg', geo.h , t,  k)*k;

% Get data for original source
lambda = zeros(cp,N); % Dirichlet
phi = zeros(cp,N); % Neumann
for i=1:N
    [U, UDX1, UDX2] = temp(geo.ctr(:,1),geo.ctr(:,2),(i-1/2)*dt,0);
    F = UDX1.*geo.nv(:,1) +UDX2.*geo.nv(:,2); 
    lambda(:,i) = U;
    phi(:,i) = F; 
end

% Apply double and single layer potentials
U1 = apply_time_convolution(single_layer, phi, dt);
U2 = apply_time_convolution(double_layer, lambda, dt);

U_t_eval = reshape(U2(:,N)-U1(:,N),size(X1)); % original field in the exterior

% Get data for source to mimic
for i=1:N
    [U, UDX1, UDX2] = temp2(geo.ctr(:,1),geo.ctr(:,2),(i-1/2)*dt,0);
    F = UDX1.*geo.nv(:,1) +UDX2.*geo.nv(:,2); 
    lambda(:,i) = U;  
    phi(:,i) = F; 
end

% Apply double and single layer potentials
U1 = apply_time_convolution(single_layer, phi, dt);
U2 = apply_time_convolution(double_layer, lambda, dt);

U_t_evalmove = reshape(U1(:,N)-U2(:,N),size(X1)); % field to mimic at final time

U_rec = U_t_eval+U_t_evalmove; % total cloaking field

% Calculate error in the exterior 
eps = 0;
ext_mask = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 > (radius+eps)^2;
log_error = log10(abs(temp2(X1,X2,t_eval,0)-temp(X1,X2,t_eval,0)+U_rec).*ext_mask);

% Plot at final time 
cmin = min(min(temp(X1,X2,t_eval,0)));
cmax = max(max(temp(X1,X2,t_eval,0)));
figure(1); clf(); hold on;
[~,h]=contourf(linspace(0,L1,N1),linspace(0,L2,N2),(temp(X1,X2,t_eval,0)-U_rec).',[linspace(cmin,cmax,16)]);
set(h,'edgecolor','none');
plot(geo.ctr(:,1),geo.ctr(:,2),'kx');
axis xy
axis tight
axis square;
hold off;
axis off;
xlim([0 L1]);
ylim([0 L2]);
colorbar;

figure(2); clf(); hold on;
[~,h]=contourf(linspace(0,L1,N1),linspace(0,L2,N2),temp(X1,X2,t_eval,0).',[linspace(cmin,cmax,16)]);
set(h,'edgecolor','none');
axis xy
axis tight
axis off;
axis square;
hold off;
xlim([0 L1]);
ylim([0 L2]);
colorbar;

figure(3); clf(); hold on;
[~,h]=contourf(linspace(0,L1,N1),linspace(0,L2,N2),temp2(X1,X2,t_eval,0).',[linspace(cmin,cmax,16)]);
set(h,'edgecolor','none');
axis xy
axis tight
axis off;
axis square;
hold off;
xlim([0 L1]);
ylim([0 L2]);
colorbar;

figure(4); clf(); hold on;
imagesc(linspace(0,L1,N1),linspace(0,L2,N2),log_error.');
plot(geo.ctr(:,1),geo.ctr(:,2),'kx');
axis xy;
axis tight;
axis off;
axis square;
xlim([0 L1]);
ylim([0 L2]);
colorbar;
hold off;