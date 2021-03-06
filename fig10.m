% Figure 10 Active Thermal Cloak
%
% Cloaking of a kite object from an incident field generated by a point
% source at (0,0) 

L1 = 1; L2 = 1; % dimensions of domain
k=.2; % thermal diffusivity

t_eval = .5; % final time of evaluation
N = 600; % number discretization steps in time

dt = t_eval / N; 

%Y1= 0.38; Y2= 0.28; % position of point source
Y1= .9; Y2 = .3;

N1 = 200; N2 = 200; %dimensions of domain 

% Points at which to evaluate
[X1,X2] = ndgrid(linspace(0,L1,N1),linspace(0,L2,N2));
p = [X1(:),X2(:)];

% Cloaking region geometry
cp_c = 128; % number of points on boundary of cloaking region

t_c = (0:cp_c-1)'*2*pi/cp_c;
geo_c = modified_kite(t_c);
%geo_c = circ(cp_c, [0.5, 0.5], 1/3);
% Inclusion geometry
cp = 128; % number of points on the boundary of the inclusion

t = (0:cp-1)'*2*pi/cp;
geo = kite(t);

% Difference of points for operators
[pmc1, pmc2] = point_diff(p,geo.ctr); % grid points minus inclusion centers 
[y1, y2] = point_diff(geo.ctr,geo.ctr); % inclusion centers minus inclusion centers
[pmc1_c, pmc2_c] = point_diff(p,geo_c.ctr); % grid points minus cloaking centers
[y1_r_c,y2_r_c] = point_diff(geo.ctr, geo_c.ctr); % cloaking centers minus inclusion centers

% Temperature function
temp = @(X1,X2,t,source_time) temp_source(X1,X2,t,source_time, Y1, Y2, 1, k);

% SLP and DLP for both the inclusion and cloaking region
double_layer = @(t,x ) dlp(x, geo.nv, pmc1, pmc2, 'reg', geo.h, t, k)*k;
single_layer = @(t, x) slp(x, pmc1, pmc2, 'reg', geo.h , t,  k)*k;
double_layer_c = @(t,x ) dlp(x, geo_c.nv, pmc1_c, pmc2_c, 'reg', geo_c.h, t, k)*k;
single_layer_c = @(t, x) slp(x, pmc1_c, pmc2_c, 'reg', geo_c.h , t,  k)*k;

% SLP and DLP on the boundary of the inclusion 
double_layer_boundary = @(t,x) dlp(x, geo.nv, y1, y2, 'reg', geo.h, t, k)*k;
single_layer_boundary = @(t,x) slp(x, y1, y2, 'reg', geo.h, t, k)*k;

% SLP inverse on the boundary of the inclusion
slp_inv = @(t,x) slp_inverse(x, y1, y2, geo.h, k, dt)/k;

% SLP and DLP evaluated at the difference between inclusion centers and
% cloaking centers 
single_layer_on_d_omega = @(t,x) slp(x, y1_r_c, y2_r_c, 'reg', geo_c.h, t, k)*k;
double_layer_on_d_omega = @(t,x) dlp(x, geo_c.nv, y1_r_c, y2_r_c, 'reg', geo_c.h, t, k)*k;

% Get data for the scattered field
lambda = zeros(cp,N); % Dirichlet
phi = zeros(cp,N); % Neumann
for i=1:N
    [U, UDX1, UDX2] = temp(geo.ctr(:,1),geo.ctr(:,2),(i-1/2)*dt,0);
    F = UDX1.*geo.nv(:,1) +UDX2.*geo.nv(:,2); 
    lambda(:,i) = U; 
    phi(:,i) = F; 
end

% Find the Neumann data for the scattered field
% V phi = -1/2 g + K g
% g = -lambda
U2_bd = apply_time_convolution(double_layer_boundary,-lambda,dt); 
% Construct the entire rhs of the calderon system
Y = 1/2*lambda+U2_bd;
% Invert the SLP to find the Neumann data
nman_data = zeros(cp,N); % initialize Neumann data for the scattered field
nman_data = invert_time_convolution(single_layer_boundary,slp_inv,Y,dt);

% Find scattered field
U2_ext = apply_time_convolution(double_layer,-lambda,dt);
U1_scat = apply_time_convolution(single_layer,nman_data,dt);
Scat_all_times = (U2_ext-U1_scat);

% Get the data for the cloaking field
lambda_c = zeros(cp_c,N); % Dirichlet 
phi_c = zeros(cp_c,N); % Neumann
for i=1:N
    [U_c, UDX1_c, UDX2_c] = temp(geo_c.ctr(:,1),geo_c.ctr(:,2),(i-1/2)*dt,0);
    F_c = UDX1_c.*geo_c.nv(:,1) +UDX2_c.*geo_c.nv(:,2); 
    lambda_c(:,i) = U_c; 
    phi_c(:,i) = F_c;
end

% Find cloaking field
U1_rec_inc = apply_time_convolution(single_layer_c, phi_c, dt); 
U2_rec_inc = apply_time_convolution(double_layer_c, lambda_c,  dt);
U_cl_all_times = U1_rec_inc - U2_rec_inc;

% Find cloaking field exactly at the boundary of the inclusion
U1_rec_inc = apply_time_convolution(single_layer_on_d_omega, phi_c, dt); 
U2_rec_inc = apply_time_convolution(double_layer_on_d_omega, lambda_c,  dt);
U_cl_on_omega = U1_rec_inc - U2_rec_inc; 

% Find scatter field from -U_cl
% V phi = -1/2 g + K g
% g = -lambda
lambda = -U_cl_on_omega;
U2_bd = apply_time_convolution(double_layer_boundary,-lambda,dt); 
% Construct the entire rhs of the calderon system
Y = 1/2*lambda+U2_bd;
% Find the field scattered by the cloaking field
nman_data = invert_time_convolution(single_layer_boundary,slp_inv,Y,dt);
U2_ext = apply_time_convolution(double_layer,-lambda,dt);
U1_scat = apply_time_convolution(single_layer,nman_data,dt);
Scat_all_times_cl = (U2_ext-U1_scat);

%Plot
t = .5;
t_eval = .5;
N = t_eval/dt;
U_cl = reshape(U_cl_all_times(:,N),size(X1));
U_scat_cl = reshape(Scat_all_times_cl(:,N),size(X1));
U_scat = reshape(Scat_all_times(:,N),size(X1));
cmin = 0;
cmax = max(max(temp(X1,X2,t_eval,0)));
log_error=log10(abs(U_scat+U_scat_cl));

figure(1); clf; 
hold on;
[~,h]=contourf(X1,X2,(temp(X1,X2,t_eval,0)-U_cl+U_scat+U_scat_cl),[-2e300,linspace(cmin,cmax,22)]);
set(h,'edgecolor','none');
h= patch(geo.ctr(:,1),geo.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
plot(geo_c.ctr(:,1),geo_c.ctr(:,2),'rx');
axis square;
colorbar;
caxis([cmin,cmax]);
axis tight;
axis xy;
axis off;
hold off;

figure(2); clf; 
hold on;
[~,h]=contourf(X1,X2,((temp(X1,X2,t_eval,0)+U_scat)),[-2e300,linspace(cmin,cmax,22)]);
set(h,'edgecolor','none');
h= patch(geo.ctr(:,1),geo.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
axis square;
axis tight;
axis xy;
axis off;
colorbar;
caxis([cmin,cmax]);
hold off;

figure(3); clf(); hold on;
imagesc(linspace(0,L1,N1),linspace(0,L2,N2),log_error.');
h= patch(geo.ctr(:,1),geo.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
axis xy;
axis tight;
hold off;
axis square;
axis off;
xlim([0 L1]);
caxis([-10,0])
ylim([0 L2]);
colorbar;

% Evaluate the fields at t=.25
t_eval = .25;
N = t_eval/dt;
U_cl = reshape(U_cl_all_times(:,N),size(X1));
U_scat_cl = reshape(Scat_all_times_cl(:,N),size(X1));
U_scat = reshape(Scat_all_times(:,N),size(X1));
cmax = max(max(temp(X1,X2,t_eval,0)));
log_error=log10(abs(U_scat+U_scat_cl));

figure(4); clf; 
hold on;
[~,h]=contourf(X1,X2,(temp(X1,X2,t_eval,0)-U_cl+U_scat+U_scat_cl),[-2e300,linspace(cmin,cmax,22)]);
set(h,'edgecolor','none');
h= patch(geo.ctr(:,1),geo.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
plot(geo_c.ctr(:,1),geo_c.ctr(:,2),'rx');
axis square;
colorbar;
caxis([cmin,cmax]);
axis tight;
axis xy;
axis off;
hold off;

figure(5); clf; 
hold on;
[~,h]=contourf(X1,X2,(temp(X1,X2,t_eval,0)+U_scat),[-2e300,linspace(cmin,cmax,22)]);
set(h,'edgecolor','none');
h= patch(geo.ctr(:,1),geo.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
axis square;
axis tight;
axis xy;
axis off;
colorbar;
caxis([cmin,cmax]);
hold off;

figure(6); clf(); hold on;
imagesc(linspace(0,L1,N1),linspace(0,L2,N2),log_error.');
h= patch(geo.ctr(:,1),geo.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
axis xy;
axis tight;
hold off;
axis square;
axis off;
xlim([0 L1]);
caxis([-10,0])
ylim([0 L2]);
colorbar;

% Evaluate the fields at t=.05
t_eval = .05;
cmin = 0;
N = t_eval/dt;
U_cl = reshape(U_cl_all_times(:,N),size(X1));
U_scat_cl = reshape(Scat_all_times_cl(:,N),size(X1));
U_scat = reshape(Scat_all_times(:,N),size(X1));
cmax = max(max(temp(X1,X2,t_eval,0)));
log_error=log10(abs(U_scat+U_scat_cl));

figure(7); clf; 
hold on;
[~,h]=contourf(X1,X2,(temp(X1,X2,t_eval,0)-U_cl+U_scat+U_scat_cl),[-2e300,linspace(cmin,cmax,22)]);
set(h,'edgecolor','none');
h= patch(geo.ctr(:,1),geo.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
plot(geo_c.ctr(:,1),geo_c.ctr(:,2),'rx');
axis square;
colorbar;
caxis([cmin,cmax]);
axis tight;
axis xy;
axis off;
hold off;

figure(8); clf; 
hold on;
[~,h]=contourf(X1,X2,(temp(X1,X2,t_eval,0)+U_scat),[-2e300,linspace(cmin,cmax,22)]);
set(h,'edgecolor','none');
h= patch(geo.ctr(:,1),geo.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
axis square;
axis tight;
axis xy;
axis off;
colorbar;
caxis([cmin,cmax]);
hold off;

figure(9); clf(); hold on;
imagesc(linspace(0,L1,N1),linspace(0,L2,N2),log_error.');
h= patch(geo.ctr(:,1),geo.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
axis xy;
axis tight;
hold off;
axis square;
axis off;
xlim([0 L1]);
caxis([-10,0])
ylim([0 L2]);
colorbar;


