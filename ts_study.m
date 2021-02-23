%Study of the sensitivity of numerical field reproductions to different
%numbers of time steps
%
%Input:
%   config.k: thermal diffusivity
%   config.N1: number of grid points in x
%   config.N2: number of gid points in y
%   config.t_eval: final time of evaluation
%   config.Y1: x-coordinate of point source
%   config.Y2: y-coordinate of point source
%   config.cp: number of points on the surface of reproduction region
%Ouput:
%   figure(1): relative errors over time where field is being reproduced
%   figure(2): absolute errors over time where zero field is expected
function[] = ts_study(config)
    L1 = 1; L2 = 1; % dimensions of domain

    %Make grid for domain
    [X1,X2] = ndgrid(linspace(0,L1,config.N1),linspace(0,L2,config.N2));
    p = [X1(:),X2(:)];

    %Geometry of region to reproduce field
    Rcenter = [L1/2,L2/2]; radius = min (L1,L2)/4; % center and radius
    geo = circ(config.cp, Rcenter, radius); %normal vectors, diffs, and centers

    %Classify source as interior (1) or exterior (0)
    condition = (config.Y1-Rcenter(1))^2+(config.Y2-Rcenter(2))^2<radius^2;

    %Differences between centers of reproduction region and grid points
    [pmc1, pmc2] = point_diff(p,geo.ctr);

    %Temperature function for a point source
    temp = @(X1,X2,t,source_time) temp_source(X1,X2,t,source_time, config.Y1, config.Y2, 1, config.k);

    %Single and double layer operators
    double_layer = @(t,x ) dlp(x, geo.nv, pmc1, pmc2, 'reg', geo.h, t, config.k)*config.k;
    single_layer = @(t, x) slp(x, pmc1, pmc2, 'reg', geo.h , t,  config.k)*config.k;

    % Masks to give inner and outer region with buffer 
    eps= .0125;
    inner_mask = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 < (radius-eps)^2;
    outer_mask = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 > (radius+eps)^2;

    er = cell(3,1); % errors where the recreated field is
    ez = cell(3,1); % errors where the field is expected to be zero
    times = cell(3,1); % cell containing vectors of times

    timing = [100,500,1000]; % different number of time stpes
    for ii = 1:length(timing)
        N = timing(ii); 
        dt = config.t_eval / N; 
        cut = N*.1; % cut off value for error plots

        % Get data for SLP and DLP
        lambda = zeros(config.cp,N); % Dirichlet
        phi = zeros(config.cp,N); % Neumann
        for i=1:N
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

        % Find the approximate norm of the error at all times
        errors_rec = zeros(N-cut+1,1); % errors where the recreated field is
        errors_zero = zeros(N-cut+1,1); % errors where the field is expected to be zero
        for i = cut:N
           U_rec = reshape(U_rec_all_times(:,i),size(X1));
           if condition==1 %exterior reproduction
               errors_rec(i-cut+1,1) = nap((temp(X1,X2,i*dt,0)-U_rec).*outer_mask,config.N1,config.N2)/nap(temp(X1,X2,i*dt,0).*outer_mask,config.N1,config.N2);
               errors_zero(i-cut+1,1) = nap(inner_mask.*U_rec,config.N1,config.N2);
           else %interior reproduction
               errors_rec(i-cut+1,1) = nap((temp(X1,X2,i*dt,0)-U_rec).*inner_mask,config.N1,config.N2)/nap(temp(X1,X2,i*dt,0).*inner_mask,config.N1,config.N2);
               errors_zero(i-cut+1,1) = nap(outer_mask.*U_rec,config.N1,config.N2);
           end
        end
        er{ii} = errors_rec;
        ez{ii} = errors_zero;
        times{ii}=dt*cut:dt:config.t_eval;
    end

    % Plot
    figure(1); clf();
    hold on;
    plot(times{1},er{1}*100)
    plot(times{2},er{2}*100)
    plot(times{3},er{3}*100)
    legend('100 steps', '500 steps', '1000 steps')
    ylabel('percent L2 error')
    xlabel('time')
    xlim([.1,1])
    hold off;

    figure(2); clf();
    hold on;
    plot(times{1},ez{1})
    plot(times{2},ez{2})
    plot(times{3},ez{3})
    legend('100 steps', '500 steps', '1000 steps')
    ylabel('L2 error')
    xlabel('time')
    xlim([.1,1])
    hold off;
end %function