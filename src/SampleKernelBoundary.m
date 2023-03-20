function [P A b X] = SampleKernelBoundary(S,C,options)

%% Input:
% g_i(x) <= 1
if nargin == 0
    close all
    demo = 2;
    switch demo 
        case 1
            S.dim = 2;
            g{1} = @(x,y) 1 + 2*y;
            g{2} = @(x,y) 2 - 4*x - 3*y;
            g{3} = @(x,y) 10-28*x-5*y - 24*x*y-18*y^2;
            g{4} = @(x,y) 1 - y - 8*x^2 - 2*x*y - y^2 - 8*x^2*y - 6*x*y^2;

            for i = 1:numel(g)
                gi = @(x) -g{i}(x(1),x(2)) + 1;
                S(1).g{i} = gi;
                S(1).dim = 2;
            end

            S.h = [];
            theta = linspace(0,2*pi,2000);
            C = [cos(theta); sin(theta)];

            for i = 1:numel(S.g)
                G{i} = @(x,y) S.g{i}([x;y]);
            end
            [Xp Yp R] = RadialLevelSet2D(G,theta);
            fill(Xp,Yp,0.7*[1 1 1]);
            hold on;
            
            theta = linspace(0,2*pi,100);
            C = [cos(theta); sin(theta)];
            
        case 2
            
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
            %@(x,y,z) 1 - x.^2 + x.*z - y;
            g{end+1} = @(x) 1 - x(1).^2 + x(1).*x(3) - x(2);
            %clear g
            %g{1} = @(x) x(1).^4 + x(2).^4 + x(3).^4 - x(1).^2.*x(2).*x(3);
            S.dim = 3;
            S.g = g;
            S.h = [];
            [Xp Yp Zp] = sphere(floor(sqrt(1000)));
            C = [Xp(:) Yp(:) Zp(:)]';
            
            for i = 1:numel(g)
                gi = @(x) -g{i}(x) + 1;
                S(1).g{i} = gi;
                S(1).dim = 3;
            end

            
            bounds = 3.5*[-1 1 -1 1 -1 1];
            options.ngrid = 50;
            options.color = [0    0.4470    0.7410];
            options.translate = [0;0;0];
            options.color = 'r';
            options.levelset = 1;
            PlotSurface(UnpackVariables(S.g),bounds,options);
            hold on;
            axis(3*[-1 1 -1 1 -1 1])
    end

end

opti = casadi.Opti();
sopts.print_level = 0;
popts = struct('expand',true,'jit',false,'print_time',false,'verbose',false,'compiler','shell');
opti.solver('ipopt',popts,sopts);

t = opti.variable;
c = opti.parameter(S.dim,1);
opti.minimize(-t);
opti.set_initial(t,1);
x = t*c;
for i = 1:numel(S.g)
    opti.subject_to(S.g{i}(x) <= 1);
end
for i = 1:numel(S.h)
    error('Equality not yet supported');
    opti.subject_to(S.h{i}(x) == 1);
end

X = [];
Ind = [];
for i = 1:size(C,2)
    i
    opti.set_value(c,C(:,i));
    try
        opti.solve();
        X = [X opti.value(x)];
        lam = opti.value(opti.lam_g);
        assert(max(lam) > 0);
        [temp ind] = max(lam);
        Ind = [Ind ind];
    catch
        disp('Solver failed');
        disp('y')
    end

end



% Construct Halfplanes
x = casadi.MX.sym('x',S.dim,1);
for i = 1:numel(S.g)
    dg{i} = casadi.Function('dg',{x},{jacobian(S.g{i}(x), x)});
end
for i = 1:size(X,2)
    xi = X(:,i);
    ind = Ind(i);
    A(i,:) = full(dg{ind}(xi));
    b(i,1) = A(i,:)*xi;
end

P = Polyhedron(A,b);

if(nargin == 0)
    switch demo
        case 1
            fill(Xp,Yp,0.7*[1 1 1]);
            hold on;
            plot(X(1,:),X(2,:),'.');
            P = Polyhedron(A,b);
            h = plot(P);
            set(h,'FaceAlpha',0.1);
            hold on
            plot(X(1,:),X(2,:),'.');
        case 2
            plot3(X(1,:),X(2,:),X(3,:),'.');
            hold on;
            P = Polyhedron(A,b);
            h = plot(P);
            set(h,'FaceAlpha',0.1);
            hold on
            save('Henrion_KernelOuterApprox','P');
    end
    
end


end

