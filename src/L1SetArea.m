function [area Xg Yg ind] = L1SetArea(f,bounds,dx,dy)
% Calculate area of set approximation obtained with l1 heuristic

% Note this is approximate as it treats each point on a grid as
% having an area rather than checking if the neighbors are also points

if nargin == 0
f = @(x1,x2) 2 - ((x1.^2 + x2.^2)/0.5^2);
bounds.x_lb = [-1;-1];
bounds.x_ub = [+1;+1];

dx = 0.001;
dy = 0.001;
end
xp = [bounds.x_lb(1):dx:bounds.x_ub(1)];
yp = [bounds.x_lb(2):dy:bounds.x_ub(2)];

options.tolerance = 1e-6;
[Xg Yg] = ndgrid(xp,yp);
Vg = NaN*Xg;
ind = [];
for i = 1:numel(Xg)
    if f(Xg(i),Yg(i)) >= 1 - options.tolerance
        ind = [ind i];
        Vg(i) = 1;
    else
        Vg(i) = 0;
    end
end

area = numel(ind)*dx*dy;






