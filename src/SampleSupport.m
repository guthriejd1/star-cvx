function [Alpha P] = SampleSupport(S,C,options)

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
            % Given vector c in Rn
            % min alpha 
            % s.t.
            %   c'x <= alpha forall x in X
            alpha = sdpvar;
            obj = alpha;
            param = alpha;
            con = [];
    end

    G = [];
    for ng = 1:numel(S.g)
        G = [G; 1-S.g{ng}(x)];
    end
    H = [];
    for nh = 1:numel(S.h)
        H = [H; 1-S.h{nh}(x)];
    end

    set(1).lhs = alpha - c'*x;
    set(1).ci = G;
    set(1).ce = H;
    set(1).z = [x];
    
    sosprog.set = set;
    sosprog.obj = obj;
    sosprog.con = con;
    sosprog.param = param;

    [result solver_info] = SolveSetContainSOS(sosprog,options);
    solver_info.info;

    if solver_info.problem == 0
        Alpha(i,1) = value(alpha);
    else
        warning('Solver issues')
        solver_info
        Alpha(i,1) = NaN;
    end

end
P = Polyhedron(C',Alpha);

end

