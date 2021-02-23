% gives a collection of patch objects for a box with corner c
% and along directions u,v,w
function h=plotbox(c,u,v,w)
A = [u,v,w];
verts = [0 0 0
         0 1 0
	 1 1 0
	 1 0 0
	 0 0 1
	 0 1 1
	 1 1 1
	 1 0 1]*A';
verts(:,1) = verts(:,1) + c(1);
verts(:,2) = verts(:,2) + c(2);
verts(:,3) = verts(:,3) + c(3);

faces = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
h = patch('vertices',verts,'faces',faces);
