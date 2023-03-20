function [X result] = SampleKernel(S,C,options)

%% Input:
% g_i(x) <= 1
if nargin == 0
S.dim = 2;
S.g{1} = @(x) x(1).^2 + x(2).^2;
S.g{2} = @(x) 4*(x(1).^3 + x(2).^2 + x(2).^1);
S.g{3} = @(x) 2*x(1);
S.h = [];
c = [-1;0];
options.s.dmax = 4;
options.s.dmin = 0;
end

%%
for i = 1:size(C,2)
    disp(['Solving ' num2str(i) ' : ' num2str(size(C,2))]);
    yalmip('clear');
    c = C(:,i);
    nx = S.dim;
    x  = sdpvar(nx,1);
    
    impl = 1;
    switch impl
        case 1
            xp = sdpvar(nx,1);
            obj = c'*xp;
            param = xp;
            con = [];
        case 2
            t = sdpvar;
            xp = sdpvar(nx,1);
            obj = -t;
            param = [t;xp];
            con = [t*c == xp];
        case 3
            xp = c;
            obj = 0;
            param = [];
            con = [];
    end
    

    G = [];
    for ng = 1:numel(S.g)
        G = [G; 1-S.g{ng}(x)];
    end

    for k = 1:numel(G)
        grad_f = -jacobian(G(k),x);
        set(k).lhs = grad_f*(x - xp);
        set(k).ce = G(k);
        set(k).ci = G(setdiff(1:numel(G),k));
        set(k).z = [x];
    end
    
    
    n_set = numel(set);
    for k = []
        set(n_set+k).lhs = 1 - S.g{k}(x);
        set(n_set+k).ce = (x-xp)'*(x-xp);
        set(n_set+k).ci = [];
        set(n_set+k).z = [x];
    end
    
    if 0
        xa = [1;3;3];
        xb = [-1;-1;1];
        xc = [1;-1;-1];
        xd = [-1;3;-3];
        Phull = Polyhedron([xa xb xc xd]');
        con = [con; Phull.A*xp <= Phull.b];
    end

    sosprog.set = set;
    sosprog.obj = obj;
    sosprog.con = con;
    sosprog.param = param;

    [result solver_info] = SolveSetContainSOS(sosprog,options);
    solver_info.info;
    xp = value(xp);
    if solver_info.problem == 0
        X(:,i) = value(xp);
    else
        warning('Solver issues')
        solver_info
        X(:,i) = NaN*ones(nx,1);
    end

end
ind = find(~isnan(X(1,:)));
P = Polyhedron(X(:,ind)');

end

