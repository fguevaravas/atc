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
function[] = therm_rep_1D_Error_study(config)
L1 = 1; 


for ii = 1:length(config.N) 
    dt = config.t_eval/config.N(ii);
    cut = config.N(ii)*.1;
    nv = [1,-1;-1,1]; %interior rep
    
    %Spatial discretization
    
    X1 = linspace(0,L1,config.N1);
    eps = 0.005; %buffer
    %geometry
    geo.cp = config.cp;
    geo.nv = nv;
    
    %operators
    temp = @(X,t) tempsrc1d(config.k,t,X,config.Y1);
    dtemp = @(X,t) dertempsrc1d(config.k,t,X,config.Y1);
    single_layer = @(t, x) slp_1d(x,t,X1,config.k,geo);
    double_layer = @(t, x) dlp_1d(x, t,X1,config.k,geo);
    
    %Get Data for Convolutions
    lambda = zeros(2,config.N(ii));
    phi = zeros(2,config.N(ii));
    for i = 1:config.N(ii)
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
    
    
    %Calculate errors
    for i = cut:config.N(ii)
        U_rec = U_rec_all_times(:,i);
        U_rec = reshape(U_rec,size(X1));
        if config.Y1-config.cp(1) <0||config.Y1-config.cp(2)>0
            ext_err = ((X1<geo.cp(1)-eps)+(X1>geo.cp(2)+eps)).*abs(U_rec);
            int_err = ((X1>geo.cp(1)+eps).*(X1<geo.cp(2)-eps)).*abs(U_rec-temp(X1,i*dt));
            int = ((X1>geo.cp(1)+eps).*(X1<geo.cp(2))-eps).*abs(temp(X1,i*dt));
            errs_abs(i-cut+1,1) = sqrt(sum(ext_err.^2)/config.N1);
            errs_rel(i-cut+1,1) = sqrt(sum(int_err.^2)/config.N1)/sqrt(sum(int.^2)/config.N1);
        else
            ext_err = ((X1<geo.cp(1)-eps)+(X1>geo.cp(2)+eps)).*abs(U_rec-temp(X1,i*dt));
            int_err = ((X1>geo.cp(1)+eps).*(X1<geo.cp(2)-eps)).*abs(U_rec);
            ext = ((X1<geo.cp(1)-eps)+(X1>geo.cp(2)+eps)).*abs(temp(X1,i*dt));
            errs_rel(i-cut+1,1) = sqrt(sum(ext_err.^2)/config.N1)/sqrt(sum(ext.^2)/config.N1);
            errs_abs(i-cut+1,1) = sqrt(sum(int_err.^2)/config.N1);
        end
    end
    err1{ii} =errs_rel;
    err2{ii} =errs_abs;
    times{ii} = dt*cut:dt:config.t_eval;
end
%Plot
if config.Y1-config.cp(1) <0||config.Y1-config.cp(2)>0
    figure(1); clf;
    hold on;
    plot(times{1},err1{1}*100)
    plot(times{2},err1{2}*100)
    xlabel('time')
    ylabel('relerr_-(u;t) (percent)')
    legend('300 steps', '3000 steps')
    xlim([0.1,1])
    ylim([0,max(err1{1})*150])
    hold off;
    
    figure(2); clf;
    hold on;
    plot(times{1},err2{1})
    plot(times{2},err2{2})
    xlabel('time')
    ylabel('err_+(u;t)')
    xlim([0.1,1])
    ylim([0,max(err2{1})*1.5])
    legend('300 steps', '3000 steps')
    hold off;
else
    figure(1); clf;
    hold on
    plot(times{1},err1{1}*100)
    plot(times{2},err1{2}*100)
    xlabel('time')
    ylabel('relerr_+(v;t) (percent)')
    xlim([0.1,1])
    ylim([0,max(err1{1})*150])
    legend('300 steps', '3000 steps')
    hold off
    
    figure(2); clf;
    hold on
    plot(times{1},err2{1})
    plot(times{2},err2{2})
    xlabel('time')
    ylabel('err_-(v;t)')
    ylim([0,max(err2{1})*1.5])
    legend('300 steps', '3000 steps')
    xlim([0.1,1])
    hold off
end


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
