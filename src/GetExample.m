function [S varargout] = GetExample(mdl)

% Sources:
% [Henrion2012] Henrion, D., & Lasserre, J. B. (2012). Inner approximations for polynomial matrix inequalities and robust stability regions. IEEE Transactions on Automatic Control

% Output:
% Semialgebraic set S
%   S = {x | g{i} <= 1,  i in [n] }

varargout = {};
switch mdl
    case 'Demo'
        g{1} = @(x) x(1)^2 + x(2)^2;
        g{2} = @(x)  x(1)*x(2) - x(1) - x(2);
        g{3} = @(x) -x(1)*x(2) - x(1) + x(2);
        S.g = g;
        S.h = [];
        S.dim = 2;
        
        varargout{1} = [];
    case 'Cerone2012'
        S.volume = 0.9966;
        g{1} = @(x) (x(1)-1)^2 + (x(2)-1)^2;
        g{2} = @(x) x(2) - 0.5*x(1)^2 + 1;
        S.g = g;
        S.h = [];
        S.dim = 2;
        varargout{1} = [];
    case 'Henrion2012_exG'
        % Format: g(x) >= 0
        g{1} = @(x) (1-16*x(1)*x(2))*(1-x(1)^2-x(2)^2) - x(1)^2;
        g{2} = @(x) (1-16*x(1)*x(2));
        g{3} = @(x) (1-x(1)^2-x(2)^2);
        
        % Convert inequality into standard form
        % Original: g(x) >= 0
        % Final: g(x) <= 1
        for i = 1:numel(g)
            g{i} = @(x) -g{i}(x) + 1;
        end
        S.g = g;
        S.h = [];
        S.dim = 2;
        
        varargout{1} = [];
    case 'Henrion2012_ex4p4'
        % [Henrion2012] - Example 4.4
        g{1} = @(x) 1 + 2*x(2);
        g{2} = @(x) 2 - 4*x(1) - 3*x(2);
        g{3} = @(x) 10 - 28*x(1) - 5*x(2) - 24*x(1)*x(2) - 18*x(2)^2;
        g{4} = @(x) 1 - x(2) - 8*x(1)^2 - 2*x(1)*x(2) - x(2)^2 - 8*x(1)^2*x(2) - 6*x(1)*x(2)^2;
        
        % Convert inequality into standard form
        % Original: g(x) >= 0
        % Final: g(x) <= 1
        for i = 1:numel(g)
            g{i} = @(x) -g{i}(x) + 1;
        end
        S.g = g;
        S.h = [];
        S.dim = 2;
        
        varargout{1} = [];
    case 'Henrion2012_ex4C'
        % Vertices of convex hull
        xv1 = [1;3;3];
        xv2 = [-1;-1;1];
        xv3 = [1;-1;-1];
        xv4 = [-1;3;-3];
        Phull = Polyhedron([xv1 xv2 xv3 xv4]');
        Phull.computeHRep;
        Phull.normalize;
        
        for i = 1:numel(Phull.b)
            g{i} = @(x) Phull.b(i) - Phull.A(i,:)*x;
        end
        
        % Hyperbolic paraboloid
        g{end+1} = @(x) 1 - x(1).^2 + x(1).*x(3) - x(2);
        
       % Convert inequality into standard form
       % Original: g(x) >= 0
       % Final: g(x) <= 1
        for i = 1:numel(g)
            g{i} = @(x) -g{i}(x) + 1;
        end
        S.g = g;
        S.h = [];
        S.dim = 3;
        
        varargout{1} = [];
    case 'Henrion2012_ex4C_old'
        % [Henrion2012] - Example IV.C
        % "The boundary of P consists of two triangles and a hyperbolic paraboloid"
        
        % (17) Convex hull of P(x)
        Phull = Polyhedron([-3 3 -1; -1 -1 1; 1 -1 -1; 3 3 1]);
        
        % Note! Henrion2012 uses reverse ordering of coefficients
        P = @(x,y,z) [  1-z^2       x - y*z                 y - x*z;...
                        x - y*z     1 + x^2 - y^2 - z^2     x - y*z;...
                        y - x*z     x - y*z                 1 - y^2];
        varargout{1} = P;
        varargout{2} = Phull;
        
        g{1} = [];
    case 'Henrion2003_ex4A_Sol'
        % [Henrion2003] - Example IV.A
        g{1} = [];
        
        % Ellipsoidal solution
        Pe = [2.3189 0 0.4458; 0 2.0180 0; 0.4458 0 1.7998];
        xc = [0;0.1459;0];
        xv1 = [1;3;3];
        xv2 = [-1;-1;1];
        xv3 = [1;-1;-1];
        xv4 = [-1;3;-3];
        % Triangles 1,2
        T1 = Polyhedron([xv1 xv2 xv3]');
        T2 = Polyhedron([xv2 xv3 xv4]');
        % Paraboloid
        xs = [0;1;0];
        
        close all
        
        plot(T1)
        hold on;
        plot(T2)
        
 

                
        % Solve for hyperbolic paraboloid
        
        H = [1 3^2; 1 1]^-1*(-[2;-2]);
        cx = H(1); cz = H(2);
        
        dx = 0.1; dy = 0.1; dz = 0.1;
        [X,Y,Z] = meshgrid(-4:dx:4,-4:dy:4,-4:dz:4);
        a=1;
        b=1;
        c=1;
        V = cz*(Z).^2 + cx*X.^2 + (Y-1);
        V = 1 - X.^2 + X.*Z - Y;
        p=patch(isosurface(X,Y,Z,V,0));
        set(p,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.5);
        clear p
        daspect([1 1 1])
        view(3);
        camlight
        plot3(0,1,0,'g.','MarkerSize',30)
        plot3(1,3,3,'g.','MarkerSize',30);
       
        xlabel('x0');
        ylabel('x1');
        zlabel('x2');
        
        [X Y Z] = sphere(10);
        Ps = sqrtm(Pe^-1);
        f1 = @(x,y,z) cz*(z).^2 + cx*(x).^2 + (y-1);
        xv = [xv1 xv2 xv3 xv4];
        for i = 1:4
            residual = f1(xv(1,i),xv(2,i),xv(3,i))
        end
        for i = 1:numel(X)
            p(:,i) = Ps*[X(i);Y(i);Z(i)] + xc;
            
            res(i) = f1(p(1,i), p(2,i), p(3,i));
        end
        Pe = Polyhedron(p');
        h = plot(Pe)
        set(h,'FaceAlpha',0.25);

        
        x0 = linspace(-1,1,10);
        x1 = linspace(-1,3,10);
        x2 = linspace(-3,3,10);
        [X0,X1,X2] = meshgrid(x0,x1,x2);
        ps = [];
        for i = 1:numel(X0)
            Gd = tf([1],[1 X2(i) X1(i) X0(i)],1);
            if(isstable(Gd,1))
                ps = [ps [X0(i);X1(i);X2(i)]];
            end
        end
        plot3(ps(1,:),ps(2,:),ps(3,:),'.')
              
                
       
        
        
end
