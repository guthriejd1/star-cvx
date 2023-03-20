function [X Y R] = RadialLevelSet2D(F,theta,options)
if nargin == 1
    theta = linspace(0,2*pi,1e3);
end
if nargin == 2
    options.r_min = 0.01;
    options.r_tol = 1e-6;
    options.gamma = 1;
end
if ~iscell(F)
    f = F; clear F;
    F{1} = f;
end

for k = 1:numel(theta)
    r_ub = options.r_min;
    is_outside_set = false;
    while(is_outside_set == false)
        x = r_ub*cos(theta(k));
        y = r_ub*sin(theta(k));
        for nf = 1:numel(F)
            v(nf) = F{nf}(x,y);
        end 
        if(max(v) >= options.gamma)
            is_outside_set = true;
            r_lb = r_ub/2;
        else
            r_ub = r_ub*2;
        end
    end
    while(r_ub - r_lb > options.r_tol)
        r = 0.5*(r_lb+r_ub);
        x = r*cos(theta(k));
        y = r*sin(theta(k));
        for nf = 1:numel(F)
            v(nf) = F{nf}(x,y);
        end
        if(max(v) >= options.gamma)
            r_ub = r;
        else
            r_lb = r;
        end
    end
    X(k) = x;
    Y(k) = y;
    R(k) = r;
end