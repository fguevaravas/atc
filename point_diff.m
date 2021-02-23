%Point differences in x and y between 2 sets of points in 2D
%
%Input:
%   p1: set of 2d points
%   p2: set of 2d points
%Output:
%   x: x-coordinates of p1-p2
%   y: y-coordinates of p1-p2
function[x,y] = point_diff(p1,p2)
    x = p1(:,1)*ones(1,size(p2,1)) - ones(size(p1,1),1)*p2(:,1)';
    y = p1(:,2)*ones(1,size(p2,1)) - ones(size(p1,1),1)*p2(:,2)';
end%function