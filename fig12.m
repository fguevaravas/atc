% Figure 12 Active Thermal Cloak
%
% Make a circle with Dirichlet boundary conditions appear as a kite with
% Dirichlet Boundary Conditioms

clear all;

L1 = 1; L2 = 1; % dimensions of domain
k=.2; % thermal diffusivity

t_eval = .05; % final time of evaluation
N=180; % number discretization steps in time
dt = t_eval / N; 


Y1= .25; Y2= .5; % position of point source

N1 = 200; N2 = 200; % number of grid points in x and y

% Points at which to evaluate
[X1,X2] = ndgrid(linspace(0,L1,N1),linspace(0,L2,N2));
p = [X1(:),X2(:)];


cp_c = 128; % number of points on boundary of cloaking region
cp_i = 128; % number of points on inclusion to fake
cp_a = 128; % number of points on actual inclusion
Rcenter = [L1/2,L2/2]; radius = min (L1,L2)/4.5; % center and radius of cloak

% Geometry for inclusions and cloak
t_i = (0:cp_i-1)'*2*pi/cp_i;
t_a = (0:cp_a-1)'*2*pi/cp_a;
g_i = flower(t_i); % fake inclusion
g_a = kite(t_a); % actual inclusion
g_c = circ(cp_c, Rcenter,radius); % cloaking region

% Differences between points
[pmc1, pmc2] = point_diff(p,g_i.ctr); % grid points - fake inclusion
[y1_i, y2_i] = point_diff(g_i.ctr,g_i.ctr); % fake inclusion - fake inclusion
[pmc1_c, pmc2_c] = point_diff(p,g_c.ctr); % grid points - cloaking region
[y1_r_c, y2_r_c] = point_diff(g_c.ctr, g_i.ctr); % cloaking region - fake inclusion
[pmc1_a, pmc2_a] = point_diff(p,g_a.ctr); % grid ponts - actual inclusion
[y1_a, y2_a] = point_diff(g_a.ctr, g_a.ctr); % actual inclusion - actual inclusion
[y1_r_a, y2_r_a] = point_diff(g_a.ctr, g_c.ctr); % actual inclusion - cloaking region

% Temperature function
temp = @(X1,X2,t,source_time) temp_source(X1,X2,t,source_time, Y1, Y2, 1, k);

% External SLPs and DLPs
double_layer = @(t,x ) dlp(x, g_i.nv, pmc1, pmc2, 'reg', g_i.h, t, k)*k;
single_layer = @(t, x) slp(x, pmc1, pmc2, 'reg', g_i.h , t,  k)*k;
double_layer_c = @(t,x ) dlp(x, g_c.nv, pmc1_c, pmc2_c, 'reg', g_c.h, t, k)*k;
single_layer_c = @(t, x) slp(x, pmc1_c, pmc2_c, 'reg', g_c.h , t,  k)*k;
double_layer_a = @(t,x ) dlp(x, g_a.nv, pmc1_a, pmc2_a, 'reg', g_a.h, t, k)*k;
single_layer_a = @(t, x) slp(x, pmc1_a, pmc2_a, 'reg', g_a.h , t,  k)*k;

% SLP, DLP, and inverses on the objects
double_layer_boundary_i = @(t,x) dlp(x, g_i.nv, y1_i, y2_i, 'reg', g_i.h, t, k)*k;
single_layer_boundary_i = @(t,x) slp(x, y1_i, y2_i, 'reg', g_i.h, t, k)*k;
slp_inv_i = @(t,x) slp_inverse(x, y1_i, y2_i, g_i.h, k, dt)/k;
double_layer_boundary_a = @(t,x) dlp(x, g_a.nv, y1_a, y2_a, 'reg', g_a.h, t, k)*k;
single_layer_boundary_a = @(t,x) slp(x, y1_a, y2_a, 'reg', g_a.h, t, k)*k;
slp_inv_a = @(t,x) slp_inverse(x, y1_a, y2_a, g_a.h, k, dt)/k;

% SLP, DLP, and derivatives to find scattered field and derivative exactly on boundary 
single_layer_on_d_r = @(t,x) slp(x, y1_r_c, y2_r_c, 'reg', g_i.h, t, k)*k;
double_layer_on_d_r = @(t,x) dlp(x, g_i.nv, y1_r_c, y2_r_c, 'reg', g_i.h, t, k)*k;
adj_double_layer_on_d_r = @(t,x ) adj_dlp(x, g_c.nv, y1_r_c, y2_r_c, g_i.h, t, k)*k;
hyp_sing_on_d_r = @(t,x) hypersingular(x, g_i.nv, g_c.nv, y1_r_c, y2_r_c, g_i.h, t, k)*k;

% SLP and DLP to find scattered field from cloaking field 
single_layer_on_d_omega = @(t,x) slp(x, y1_r_a, y2_r_a, 'reg', g_c.h, t, k)*k;
double_layer_on_d_omega = @(t,x) dlp(x, g_c.nv, y1_r_a, y2_r_a, 'reg', g_c.h, t, k)*k;

% Find scattered field from actual inclusion

lambda = zeros(cp_a,N); % Dirichlet data (Actual inclusion)
phi = zeros(cp_a,N); % Neumann data (Actual inclusion)
nman_data = zeros(cp_a,N); % Unknown Neumann data for scattered field

% Get data
for i=1:N
    [U, UDX1, UDX2] = temp(g_a.ctr(:,1),g_a.ctr(:,2),(i-1/2)*dt,0);
    F = UDX1.*g_a.nv(:,1) +UDX2.*g_a.nv(:,2); 
    lambda(:,i) = U; 
    phi(:,i) = F; 
end

% Find the Neumann data for the scattered field
% V phi = -1/2 g + K g
% g = -lambda
U2_bd = apply_time_convolution(double_layer_boundary_a,-lambda,dt); 
% Construct the entire rhs of the calderon system
Y = 1/2*lambda+U2_bd;
% Find the Neumann data by inverting the SLP
nman_data = invert_time_convolution(single_layer_boundary_a,slp_inv_a,Y,dt);

% Apply SLP and DLP
U2_ext = apply_time_convolution(double_layer_a,-lambda,dt);
U1_scat = apply_time_convolution(single_layer_a,nman_data,dt);

Scat_all_times = (U2_ext-U1_scat); % Scattered field at all times
U_scat = reshape(Scat_all_times(:,N),size(X1)); % Scattered field at final time

% Find cloaking field 

lambda_c = zeros(cp_c,N); % Dirichlet Data (Cloaking region)
phi_c = zeros(cp_c,N); % Neumann data (Cloaking region)

% Get data
for i=1:N
    [U_c, UDX1_c, UDX2_c] = temp(g_c.ctr(:,1),g_c.ctr(:,2),(i-1/2)*dt,0);
    F_c = UDX1_c.*g_c.nv(:,1) +UDX2_c.*g_c.nv(:,2); 
    lambda_c(:,i) = U_c; 
    phi_c(:,i) = F_c; 
end

% Apply SLP and DLP
U1_rec_inc = apply_time_convolution(single_layer_c, phi_c, dt); 
U2_rec_inc = apply_time_convolution(double_layer_c, lambda_c,  dt);
U_cl_all_times = U1_rec_inc - U2_rec_inc; % recreated incident field at all times
U_cl = reshape(U_cl_all_times(:,N),size(X1)); %recreated incident field at final time

% Find exact values of incident field on the boundary of the actual inclusion
U1_rec_inc = apply_time_convolution(single_layer_on_d_omega, phi_c, dt); 
U2_rec_inc = apply_time_convolution(double_layer_on_d_omega, lambda_c,  dt);
U_cl_on_omega = U1_rec_inc - U2_rec_inc; 

% Find scattered field from cloaking field 
% V phi = -1/2 g + K g
% g = -lambda
lambda = -U_cl_on_omega;
U2_bd = apply_time_convolution(double_layer_boundary_a,-lambda,dt); 
% Construct the entire rhs of the calderon system
Y = 1/2*lambda+U2_bd;
nman_data = invert_time_convolution(single_layer_boundary_a,slp_inv_a,Y,dt);

% Apply SLP 
U2_ext = apply_time_convolution(double_layer_a,-lambda,dt);
U1_scat = apply_time_convolution(single_layer_a,nman_data,dt);
Scat_all_times_cl = (U2_ext-U1_scat); % Scattered field from cloaking field at all times
U_scat_cl = reshape(Scat_all_times_cl(:,N),size(X1)); % Scattered field from cloaking field at final time

% Find scattered field from fake inclusion to mimic
lambda = zeros(cp_i,N); % Dirichlet data (fake inclusion)
phi = zeros(cp_i,N); % Neumann data (fake inclusion)
nman_data = zeros(cp_i,N); % Unkown neumann data for scattered field

% Get Data
for i=1:N
    [U, UDX1, UDX2] = temp(g_i.ctr(:,1),g_i.ctr(:,2),(i-1/2)*dt,0);
    F = UDX1.*g_i.nv(:,1) +UDX2.*g_i.nv(:,2); 
    lambda(:,i) = U;  
    phi(:,i) = F; 
end

% Get unknown Neumann data
% V phi = -1/2 g + K g
% g = -lambda
U2_bd = apply_time_convolution(double_layer_boundary_i,-lambda,dt); 
% Construct the entire rhs of the calderon system
Y = 1/2*lambda+U2_bd;

% Invert SLP
nman_data = invert_time_convolution(single_layer_boundary_i,slp_inv_i,Y,dt);

% Apply SLP & DLP (on the boundary of cloaking region)
U2_ext = apply_time_convolution(double_layer_on_d_r,-lambda,dt);
U1_scat = apply_time_convolution(single_layer_on_d_r,nman_data,dt);
mimic_dirichlet_data = (U2_ext-U1_scat); 

% Apply SLP & DLP derivatives (on boundary of cloaking region)
deriv_U2_ext = apply_time_convolution(hyp_sing_on_d_r,-lambda,dt);
deriv_U1_scat = apply_time_convolution(adj_double_layer_on_d_r,nman_data,dt);
mimic_neumann_data = (deriv_U2_ext-deriv_U1_scat);

% Apply SLP & DLP 
m_U1 = apply_time_convolution(single_layer_c, mimic_neumann_data,dt);
m_U2 = apply_time_convolution(double_layer_c, mimic_dirichlet_data,dt);
mimic_all_times = (m_U2-m_U1); % Mimic field at all times
U_mimic = reshape(mimic_all_times(:,N),size(X1)); % Mimic field at final time

% Find scattered field from fake inclusion

% Apply SLP and DLP
U1_scat = apply_time_convolution(single_layer,nman_data,dt);
U2_ext = apply_time_convolution(double_layer,-lambda,dt);
Scat_all_times = (U2_ext-U1_scat); % Scattered field from fake inclusion all times
U_scat_i = reshape(Scat_all_times(:,N),size(X1)); % Scattered field from fake inclusion final time

% Mask for exterior errors
eps = 0;
ext_mask = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 > (radius+eps)^2;

% Calculate errors in scattered fields 
whole_field= U_scat+U_scat_cl+U_mimic;

%Plots

cmin = min(min(U_scat)); 
cmax = -.005;
figure(1); clf();
hold on;
[~,h]=contourf(linspace(0,L1,N1),linspace(0,L2,N2),(U_scat).',[linspace(cmin,cmax,16)]);
set(h,'edgecolor','none');
axis xy;
h= patch(g_a.ctr(:,1),g_a.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
axis tight;
axis xy;
axis off;
axis equal;
xlim([0 L1]);
ylim([0 L2]);
colorbar;
hold off;

cmin = min(min(U_scat_i)); 
cmax = -.005;
figure(2); clf();
hold on;
[~,h]=contourf(linspace(0,L1,N1),linspace(0,L2,N2),(U_scat_i).',[linspace(cmin,cmax,16)]);
set(h,'edgecolor','none');
axis xy;
h= patch(g_i.ctr(:,1),g_i.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
axis tight;
axis xy;
axis off;
axis equal;
xlim([0 L1]);
ylim([0 L2]);
colorbar;
hold off;

cmin = min(min(whole_field)); 
cmax = -.005;
figure(3); clf();
hold on;
[~,h]=contourf(linspace(0,L1,N1),linspace(0,L2,N2),whole_field.',[linspace(cmin,cmax,16)]);
set(h,'edgecolor','none');
plot(g_c.ctr(:,1),g_c.ctr(:,2),'rx');
axis xy;
h= patch(g_a.ctr(:,1),g_a.ctr(:,2),'k'); set(h, 'edgecolor', 'none');
axis tight;
axis xy;
axis off;
axis equal;
xlim([0 L1]);
ylim([0 L2]);
colorbar
hold off;

