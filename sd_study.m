%Study of the sensitivity of numerical field reproductions to different
%numbers of equispaced discretizations of the boundary of a circular
%reproduction region
%
%Input:
%   config.k: thermal diffusivity
%   config.N1: number of grid points in x
%   config.N2: number of gid points in y
%   config.N: number of time steps
%   config.t_eval: final time of evaluation
%   config.Y1: x-coordinate of point source
%   config.Y2: y-coordinate of point source
%Ouput:
%   figure(1): relative errors over time where field is being reproduced
%   figure(2): absolute errors over time where zero field is expected
function[] = sd_study(config)
    L1 = 1; L2 = 1; % dimensions of domain
    
    dt = config.t_eval / config.N; 
    cut = config.N*.1; % cut off time step for plotting errors
    devices = [10,50,100]; % number of points on reproduction regions surface
     
    % Points at which to evaluate
    [X1,X2] = ndgrid(linspace(0,L1,config.N1),linspace(0,L2,config.N2));
    p = [X1(:),X2(:)];
    
    Rcenter = [L1/2,L2/2]; radius = min (L1,L2)/4; % center and radius of reproduction region
    eps= .0125; % buffer size
    
    % Classify source as interior (1) or exterior (0)
    condition = (config.Y1-Rcenter(1))^2+(config.Y2-Rcenter(2))^2<radius^2;
    
    % Temperature function
    temp = @(X1,X2,t,source_time) temp_source(X1,X2,t,source_time, config.Y1, config.Y2, 1, config.k); 
    
    % Masks to give inner and outer region with buffer 
    inner_mask = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 < (radius-eps)^2;
    outer_mask = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 > (radius+eps)^2;

    errors_rec = zeros(config.N-cut+1,length(devices)); % errors where the recreated field is
    errors_zero = zeros(config.N-cut+1,length(devices)); % errors where the field is expected to be zero
    for j = 1:length(devices)
        cp = devices(j); % number of points on boundary of circ region
        geo = circ(cp, Rcenter, radius); %normal vectors, diffs, and centers
        
        % Differences between centers of reproduction and grid points
        [pmc1, pmc2] = point_diff(p,geo.ctr);
        
        % Single and double layer operators
        double_layer = @(t,x ) dlp(x, geo.nv, pmc1, pmc2, 'reg', geo.h, t, config.k)*config.k;
        single_layer = @(t, x) slp(x, pmc1, pmc2, 'reg', geo.h , t,  config.k)*config.k;
        
        % Get data for layer potentials
        lambda = zeros(cp,config.N); % Dirichlet Data
        phi = zeros(cp,config.N); % Neumann Data
        for i=1:config.N
            % evaluate the temperature and derivatives at the centers of
            % the segments on the bdry
            [U, UDX1, UDX2] = temp(geo.ctr(:,1),geo.ctr(:,2),(i-1/2)*dt,0);
            F = UDX1.*geo.nv(:,1) +UDX2.*geo.nv(:,2); 
            lambda(:,i) = U;
            phi(:,i) = F; 
        end
        
        % Evaluate SLP and DLP
        U1 = apply_time_convolution(single_layer, phi, dt);
        U2 = apply_time_convolution(double_layer, lambda, dt);

        if (condition == 1)
            U_rec_all_times = U2-U1; % exterior reproduction
        else
            U_rec_all_times = U1-U2; % interior reproduction
        end
        
        % find the approximate norm of the error/relative error
        for i = cut:config.N
           U_rec = reshape(U_rec_all_times(:,i),size(X1));
           if condition==1 %exterior reproduction
               errors_rec(i-cut+1,j) = nap((temp(X1,X2,i*dt,0)-U_rec).*outer_mask,config.N1,config.N2)/nap(temp(X1,X2,i*dt,0).*outer_mask,config.N1,config.N2);
               errors_zero(i-cut+1,j) = nap(inner_mask.*U_rec,config.N1,config.N2);
           else %interior reproduction
               errors_rec(i-cut+1,j) = nap((temp(X1,X2,i*dt,0)-U_rec).*inner_mask,config.N1,config.N2)/nap(temp(X1,X2,i*dt,0).*inner_mask,config.N1,config.N2);
               errors_zero(i-cut+1,j) = nap(outer_mask.*U_rec,config.N1,config.N2);
           end
        end
    end
    times = dt*cut:dt:config.t_eval;
    
    % plot
    figure(1); clf();
    hold on;
    plot(times,errors_rec(:,1)*100)
    plot(times,errors_rec(:,2)*100)
    plot(times,errors_rec(:,3)*100, 'k--')
    legend('10 points', '50 points', '100 points')
    ylabel('percent L2 error')
    xlabel('time')
    xlim([.1,1])
    hold off;

    figure(2); clf();
    hold on;
    plot(times,errors_zero(:,1))
    plot(times,errors_zero(:,2))
    plot(times,errors_zero(:,3), 'k--')
    legend('10 points', '50 points', '100 points')
    ylabel('L2 error')
    xlim([.1,1])
    xlabel('time')
    hold off;
    
end%function