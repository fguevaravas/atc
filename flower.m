% function [x] = flower(t)
%
% Flower scatterer in parametric form
% 
% Returns center points, lengths of segmenets, and normal vectors
function [geo] = flower(t)
  %parameters 
  a = 4.8;
  b = -.03;
  c = .5;
  x=  b*(a+cos(8*t)).*cos(t)+c;
  y =  b*(a+cos(8*t)).*sin(t)+c;
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
