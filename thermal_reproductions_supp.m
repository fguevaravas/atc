%Reproductions of a field generated by a point source using a
%circular region. Includes the max error over [t_min,t_eval]
%and the absolute difference between this error and the error at t_min
%
%Input:
%   config.k: thermal diffusivity
%   config.N1: number of grid points in x
%   config.N2: number of gid points in y
%   config.N: number of time steps
%   config.t_eval: final time of evaluation
%   config.Y1: x-coordinate of point source
%   config.Y2: y-coordinate of point source
%   config.source_time: time of point source
%   config.cp: number of points on the surface of reproduction region
%   config.t_min: displayed time of evaluation
%   config.errors: switch to calculate with errrors in the data
%   config.perc: percentage of error to make at each time step
%Ouput:
%   figure(1): reproduced field using boundary integral method
%   figure(2): original field from point source
%   figure(3): log plot of errors inside and outside
function[] = thermal_reproductions_supp(config)
    L1 = 1; L2 = 1; % dimensions of domain
    dt = config.t_eval / config.N;
    
    %Make grid for domain
    [X1,X2] = ndgrid(linspace(0,L1,config.N1),linspace(0,L2,config.N2));
    p = [X1(:),X2(:)];
    
    %Geometry of region to reproduce field
    Rcenter = [L1/2,L2/2]; radius = min (L1,L2)/4; % center and radius
    geo = circ(config.cp, Rcenter, radius); %normal vectors, diffs, and centers
    
    %Classify source as interior (1) or exterior (0)
    condition = (config.Y1-Rcenter(1))^2+(config.Y2-Rcenter(2))^2<radius^2; 
    
    %Differences between centers of reproduction and grid points
    [pmc1, pmc2] = point_diff(p,geo.ctr);
    
    %Temperature function for a point source
    temp = @(X1,X2,t,source_time) temp_source(X1,X2,t,source_time, config.Y1, config.Y2, 1, config.k);
    
    %Single and double layer operators
    double_layer = @(t,x ) dlp(x, geo.nv, pmc1, pmc2, 'reg', geo.h, t, config.k)*config.k;
    single_layer = @(t, x) slp(x, pmc1, pmc2, 'reg', geo.h , t,  config.k)*config.k;
    
    %Get data
    lambda = zeros(config.cp,config.N); %Dirichlet data
    phi = zeros(config.cp,config.N); %Neumann data
    for i=1:config.N
        % evaluate the temperature and derivatives at the centers of
        % the segments on the bdry
        [U, UDX1, UDX2] = temp(geo.ctr(:,1),geo.ctr(:,2),(i-1/2)*dt,0);
        F = UDX1.*geo.nv(:,1) +UDX2.*geo.nv(:,2); 
       if config.errors==1
            lambda(:,i) = U+ config.perc*norm(U)*randn(size(U)); 
            phi(:,i) = F+ config.perc*norm(F)*randn(size(F));
       else 
            lambda(:,i) = U;
            phi(:,i) = F;
       end
    end
    
    U1 = apply_time_convolution(single_layer, phi, dt);
    U2 = apply_time_convolution(double_layer, lambda, dt);

    if (condition == 1)
        U_rec_all_times = U2-U1; %Exterior representation formula
    else
        U_rec_all_times = U1-U2; %Interior representation formula
    end
    
    %Recreated field at final time
    U_rec = reshape(U_rec_all_times(:,config.t_min/dt),size(X1));
    
    %Calculate errrors, separate where field is expected to be zero and
    %where field is recreated
    mask_ext = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 > (radius)^2;
    mask_int = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 < (radius)^2;
    if (condition ==1)
        err_1 = log10(abs(U_rec-temp(X1,X2,config.t_min,0))).*mask_ext;
        err_2 = log10(abs(U_rec)).*mask_int;
    else
        err_1 = log10(abs(U_rec)).*mask_ext;
        err_2 = log10(abs(U_rec-temp(X1,X2,config.t_min,0))).*mask_int;
    end
    log_error = err_1 + err_2;
    
    %Calculate max norm of errors on [t_min,t_eval]
    for i = config.t_min/dt:config.t_eval/dt
        t_new = dt*i;
        U_rec = reshape(U_rec_all_times(:,i),size(X1));
        if (condition ==1)
            err_1 = log10(abs(U_rec-temp(X1,X2,t_new,0))).*mask_ext;
            err_2 = log10(abs(U_rec)).*mask_int;
        else
            err_1 = log10(abs(U_rec)).*mask_ext;
            err_2 = log10(abs(U_rec-temp(X1,X2,t_new,0))).*mask_int;
        end
        error_time_step = err_1 + err_2;
        all_errors(:,:, i-config.t_min/dt+1) = error_time_step;
    end
    log_error_sup = max(all_errors, [],3);
    
    %Plot
    figure(1); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),(U_rec));
    plot(geo.ctr(:,1),geo.ctr(:,2),'rx');
    axis xy;
    axis square;
    axis tight;
    axis off;
    hold off;
    xlim([0 L1]);
    ylim([0 L2]);
    colorbar;

    figure(2); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),temp(X1,X2,config.t_eval,0));
    hold off;
    axis xy;
    axis tight;
    axis square;
    axis off;
    xlim([0 L1]);
    ylim([0 L2]);
    colorbar;

    figure(3); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),log_error);
    axis xy;
    axis tight;
    hold off;
    axis square;
    axis off;
    xlim([0 L1]);
    caxis([-10,0])
    ylim([0 L2]);
    colorbar;
    
    figure(4); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),log_error_sup);
    axis xy;
    axis tight;
    hold off;
    axis square;
    axis off;
    xlim([0 L1]);
    caxis([-10,0])
    ylim([0 L2]);
    colorbar;
    
    figure(5); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),abs(log_error_sup-log_error));
    axis xy;
    axis tight;
    hold off;
    axis square;
    axis off;
    xlim([0 L1]);
    caxis([0,2])
    ylim([0 L2]);
    colorbar;
end %function


