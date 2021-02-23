%Reproductions of a field generated by a harmonic IC
%
%Input:
%   config.k: thermal diffusivity
%   config.N1: number of grid points in x
%   config.N2: number of gid points in y
%   config.N: number of time steps
%   config.t_eval: final time of evaluation
%   config.cp: number of points on the surface of reproduction region
%Ouput:
%   figure(1): reproduced field using boundary integral method
%   figure(2): original field from point source
%   figure(3): log plot of errors inside and outside
function[] = thermal_reproductions_supp_IC(config)
    L1 = 1; L2 = 1; % dimensions of domain
    dt = config.t_eval / config.N;
    
    %Make grid for domain
    [X1,X2] = ndgrid(linspace(0,L1,config.N1),linspace(0,L2,config.N2));
    p = [X1(:),X2(:)];
    
    %Geometry of region to reproduce field
    Rcenter = [L1/2,L2/2]; radius = min (L1,L2)/4; % center and radius
    geo = circ(config.cp, Rcenter, radius); %normal vectors, diffs, and centers
    
    %Differences between centers of reproduction and grid points
    [pmc1, pmc2] = point_diff(p,geo.ctr);
    
    %Temperature function for a point source
    temp = @(X1,X2) temp_int(X1,X2);
    
    %Single and double layer operators
    double_layer = @(t,x ) dlp(x, geo.nv, pmc1, pmc2, 'reg', geo.h, t, config.k)*config.k;
    single_layer = @(t, x) slp(x, pmc1, pmc2, 'reg', geo.h , t,  config.k)*config.k;
    
    %volume integral 
    [U, UDX1, UDX2]= temp_int(geo.ctr(:,1),geo.ctr(:,2));
    F = UDX1.*geo.nv(:,1)+UDX2.*geo.nv(:,2);
    
    %Get data
    lambda = zeros(config.cp,config.N); %Dirichlet data
    phi = zeros(config.cp,config.N); %Neumann data
    for i=1:config.N
        % evaluate the temperature and derivatives at the centers of
        % the segments on the bdry
        [U, UDX1, UDX2] = temp(geo.ctr(:,1),geo.ctr(:,2));
        F = UDX1.*geo.nv(:,1) +UDX2.*geo.nv(:,2); 
        lambda(:,i) = U;
        phi(:,i) = F;
    end
    
    U1 = apply_time_convolution(single_layer, phi, dt);
    U2 = apply_time_convolution(double_layer, lambda, dt);
    U_rec_all_times = U1-U2; %Interior representation formula
    
    %Volume integral
    U10 = dlp(U, geo.nv, pmc1, pmc2, 'int', geo.h, config.t_eval, config.k);
    U20 =slp(F, pmc1, pmc2, 'int', geo.h , config.t_eval,  config.k);
    
    initial_pot = reshape(U20 - U10, size(X1));
    %Recreated field at final time
    U_rec = reshape(U_rec_all_times(:,config.t_eval/dt),size(X1));
    %Calculate errrors, separate where field is expected to be zero and
    %where field is recreated
    mask_ext = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 > (radius)^2;
    mask_int = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 < (radius)^2;
    
    err_1 = log10(abs(U_rec+initial_pot)).*mask_ext;
    err_2 = log10(abs(U_rec+initial_pot-temp(X1,X2))).*mask_int;

    log_error = err_1 + err_2;
    
    %%%%Plot
    %t = t_eval
    cmax1 = max(max(U_rec));
    cmin1 = -cmax1;
    figure(1); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),(U_rec));
    axis xy;
    axis square;
    axis tight;
    axis off;
    hold off;
    xlim([0 L1]);
    ylim([0 L2]);
    caxis([cmin1,cmax1])
    colormap(redblue)
    %colorbar('southoutside')
   
    figure(2); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),initial_pot);
    hold off;
    axis xy;
    axis tight;
    axis square;
    axis off;
    xlim([0 L1]);
    ylim([0 L2]);
    caxis([cmin1,cmax1])
    colormap(redblue)
    %colorbar('southoutside')
   
    figure(3); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),U_rec+initial_pot);
     plot(geo.ctr(:,1),geo.ctr(:,2),'rx');
    hold off;
    axis xy;
    axis tight;
    axis square;
    axis off;
    xlim([0 L1]);
    ylim([0 L2]);
    caxis([cmin1,cmax1])
    colormap(redblue)
    %colorbar('southoutside')
    
    figure(4); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),log_error);
    axis xy;
    axis tight;
    hold off;
    axis square;
    axis off;
    xlim([0 L1]);
    caxis([-10,0])
    ylim([0 L2]);
    colorbar('southoutside')
    %%%%t=t_pics(2)
    step = config.t_pics(2)/dt;
    %Volume integral
    U10 = dlp(U, geo.nv, pmc1, pmc2, 'int', geo.h, config.t_pics(2), config.k);
    U20 =slp(F, pmc1, pmc2, 'int', geo.h , config.t_pics(2),  config.k);
    
    initial_pot = reshape(U20 - U10, size(X1));
    %Recreated field at final time
    U_rec = reshape(U_rec_all_times(:,step),size(X1));
    %Calculate errrors, separate where field is expected to be zero and
    %where field is recreated
    mask_ext = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 > (radius)^2;
    mask_int = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 < (radius)^2;
    
    err_1 = log10(abs(U_rec+initial_pot)).*mask_ext;
    err_2 = log10(abs(U_rec+initial_pot-temp(X1,X2))).*mask_int;

    log_error = err_1 + err_2;
    
    figure(5); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),(U_rec));
    axis xy;
    axis square;
    axis tight;
    axis off;
    hold off;
    xlim([0 L1]);
    ylim([0 L2]);
    caxis([cmin1,cmax1])
    colormap(redblue)
    
    figure(6); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),initial_pot);
    hold off;
    axis xy;
    axis tight;
    axis square;
    axis off;
    xlim([0 L1]);
    ylim([0 L2]);
    caxis([cmin1,cmax1])
    colormap(redblue)
    
    figure(7); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),U_rec+initial_pot);
     plot(geo.ctr(:,1),geo.ctr(:,2),'rx');
    hold off;
    axis xy;
    axis tight;
    axis square;
    axis off;
    xlim([0 L1]);
    ylim([0 L2]);
    caxis([cmin1,cmax1])
    colormap(redblue)
    
    figure(8); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),log_error);
    axis xy;
    axis tight;
    hold off;
    axis square;
    axis off;
    xlim([0 L1]);
    caxis([-10,0])
    ylim([0 L2]);
    
    %%%%t=t_pics(1)
    step = config.t_pics(1)/dt;
    %Volume integral
    U10 = dlp(U, geo.nv, pmc1, pmc2, 'int', geo.h, config.t_pics(1), config.k);
    U20 =slp(F, pmc1, pmc2, 'int', geo.h , config.t_pics(1),  config.k);
    
    initial_pot = reshape(U20 - U10, size(X1));
    %Recreated field at final time
    U_rec = reshape(U_rec_all_times(:,step),size(X1));
    %Calculate errrors, separate where field is expected to be zero and
    %where field is recreated
    mask_ext = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 > (radius)^2;
    mask_int = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 < (radius)^2;
    
    err_1 = log10(abs(U_rec+initial_pot)).*mask_ext;
    err_2 = log10(abs(U_rec+initial_pot-temp(X1,X2))).*mask_int;

    log_error = err_1 + err_2;
    
    figure(9); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),(U_rec));
    axis xy;
    axis square;
    axis tight;
    axis off;
    hold off;
    xlim([0 L1]);
    ylim([0 L2]);
    caxis([cmin1,cmax1])
    colormap(redblue)

    figure(10); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),initial_pot);
    hold off;
    axis xy;
    axis tight;
    axis square;
    axis off;
    xlim([0 L1]);
    ylim([0 L2]);
    caxis([cmin1,cmax1])
    colormap(redblue)
    
    figure(11); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),U_rec+initial_pot);
     plot(geo.ctr(:,1),geo.ctr(:,2),'rx');
    hold off;
    axis xy;
    axis tight;
    axis square;
    axis off;
    xlim([0 L1]);
    ylim([0 L2]);
    colormap(redblue)
    caxis([cmin1,cmax1])
    
    figure(12); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),log_error);
    axis xy;
    axis tight;
    hold off;
    axis square;
    axis off;
    xlim([0 L1]);
    caxis([-10,0])
    ylim([0 L2]);
    
     %%%%t=t_pics(3)
    step = config.t_pics(3)/dt;
    %Volume integral
    U10 = dlp(U, geo.nv, pmc1, pmc2, 'int', geo.h, config.t_pics(3), config.k);
    U20 =slp(F, pmc1, pmc2, 'int', geo.h , config.t_pics(3),  config.k);
    
    initial_pot = reshape(U20 - U10, size(X1));
    %Recreated field at final time
    U_rec = reshape(U_rec_all_times(:,step),size(X1));
    %Calculate errrors, separate where field is expected to be zero and
    %where field is recreated
    mask_ext = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 > (radius)^2;
    mask_int = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 < (radius)^2;
    
    err_1 = log10(abs(U_rec+initial_pot)).*mask_ext;
    err_2 = log10(abs(U_rec+initial_pot-temp(X1,X2))).*mask_int;

    log_error = err_1 + err_2;
    
    figure(13); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),(U_rec));
    axis xy;
    axis square;
    axis tight;
    axis off;
    hold off;
    xlim([0 L1]);
    ylim([0 L2]);
    colormap(redblue)
    caxis([cmin1,cmax1])

    figure(14); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),initial_pot);
    hold off;
    axis xy;
    axis tight;
    axis square;
    axis off;
    xlim([0 L1]);
    ylim([0 L2]);
    colormap(redblue)
    caxis([cmin1,cmax1])
    
    figure(15); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),U_rec+initial_pot);
     plot(geo.ctr(:,1),geo.ctr(:,2),'rx');
    hold off;
    axis xy;
    axis tight;
    axis square;
    axis off;
    xlim([0 L1]);
    ylim([0 L2]);
    colormap(redblue)
    caxis([cmin1,cmax1])
    
    figure(16); clf(); hold on;
    imagesc(linspace(0,L1,config.N1),linspace(0,L2,config.N2),log_error);
    axis xy;
    axis tight;
    hold off;
    axis square;
    axis off;
    xlim([0 L1]);
    caxis([-10,0])
    ylim([0 L2]);
    save('harmonic_rec.mat')
end %function



