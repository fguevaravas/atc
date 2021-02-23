%%%Circle Geometry
%function [geo] = circ(cp,Rcenter,radius)
%returns normal vectors, ctrs, and h for a circle
function[geo] = circ(cp,Rcenter,radius)
circp = zeros(cp,2);
for i=1:cp
    theta = ((i-1)*pi)./(cp*0.5);
    circp(i,1) = (radius*cos(theta)+Rcenter(1)); % xvals on our circle
    circp(i,2) = (radius*sin(theta)+Rcenter(2)); % yvals on our circle
end
% first calculate the normal vectors
% segments
s = [diff(circp,1,1); circp(1,:) - circp(end,:)];

% lengths of the segments
h = sqrt(s (:,1).^2 + s (:,2).^2);

% centers of each segment
ctr = (circp + circshift(circp,[-1,0]))/2;

% normals
nv = [ s(:,2) , -s(:,1) ];
for i=1:length(nv)
    nv(i,:) = nv(i,:)/norm(nv(i,:));
end

geo.nv = nv;
geo.h = h;
geo.ctr = ctr;

end%function