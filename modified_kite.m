% function [x] = modified_kite(t)
%
% Kite scatterer in parametric form, as defined in Fig 3.1 of Colton and Kress,
% Inverse acoustic and electromagnetic scattering.
% 
% Input:
%   t: parameter
% Ouput:
%   geo.h: length of segments
%   geo.ctr: center points of segments
%   geo.nv: normal vectors to center points of segments
function [geo] = modified_kite(t)
  %curve
  x=  1.25*(.08*(cos(t)  +   0.65*cos(2*t)) + 0.465/1.2);
  y =  1.25*(.15*sin(t)+.5/1.25);
  p = [x,y];
  s = [diff(p,1,1); p(1,:) - p(end,:)];
  nv = [ s(:,2) , -s(:,1) ];
  for i=1:length(nv)
    nv(i,:) = nv (i,:)/norm(nv(i,:));
  end
  geo.h = sqrt(s (:,1).^2 + s (:,2).^2);
  geo.ctr = (p+ circshift(p,[-1,0]))/2;
  geo.nv = nv;
end
