% Figure 1 - plots Peltier elements with Matlab
theta = linspace(0,2*pi,200);
region.r = 2;
plate.a1 = 7; 
plate.a2 = 7;;
kite.x = 6*0.08*(cos(theta) + 0.65*cos(2*theta));
kite.y = 6*(0.15*sin(theta));

figure(); hold on;

% plot plate
x = [-1 1 1 -1]*plate.a1/2;
y = [-1 -1 1 1]*plate.a2/2;
z = [ 0 0 0 0];
h=patch(x,y,z);
set(h,'facecolor','none')

% plot region omega
h=patch(region.r*cos(theta),region.r*sin(theta),0*theta)
set(h,'facecolor','none')

% plot kite
patch(kite.x,kite.y,0*kite.x,'k')

% plot out of plane Peltier elements
dr=0.4; dt=0.4; dz = 0.1; n = 10;
for t=pi/n + (0:n-1)*2*pi/n,
 x = region.r*cos(t); y = region.r*sin(t);
 vr = dr*[cos(t);sin(t);0];
 vt = dt*[-sin(t);cos(t);0];
 vz = [0;0;dz];
 c = [x;y;0] - vr/2 - vt/2;
 h = rect(c,vr,vt,vz);
 set(h,'facecolor','r','edgecolor','none')

 h = rect(c+vz,vr,vt,vz);
 set(h,'facecolor','b','edgecolor','none')
end

% plot in plane Peltier elements
delta_r = region.r*0.05;
dr=0.2; dt=0.4; dz = 0.2; n = 10;
for t=(0:n-1)*2*pi/n,
 x = (region.r + delta_r)*cos(t); y = (region.r + delta_r)*sin(t);
 vr = dr*[cos(t);sin(t);0];
 vt = dt*[-sin(t);cos(t);0];
 vz = [0;0;dz];
 c = [x;y;0] - vr/2 - vt/2 - vz/2;
 h = rect(c,vr,vt,vz);
 set(h,'facecolor','r','edgecolor','none')

 x = (region.r - delta_r)*cos(t); y = (region.r - delta_r)*sin(t);
 c = [x;y;0] - vr/2 - vt/2 - vz/2;
 h = rect(c,vr,vt,vz);
 set(h,'facecolor','b','edgecolor','none')

end


view(3)
axis equal;
axis off;
%light('Position',[1 3 2]);
%light('Position',[-3 -1 3]);
light('position',[0.1,0.1,1])
light('position',[-0.1,-0.1,1])
lighting gouraud 
%text(2,2,0,'\Omega','fontsize',18); % outside
text(0.8,-1,0,'\Omega','fontsize',18); % inside
print('-dpng','peltier.png')
