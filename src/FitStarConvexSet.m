function [result] = FitStarConvexSet(S,dS,zmon,options,Scomp)
yalmip('clear')
nx = dS(1).dim;
x = sdpvar(nx,1);
nm = numel(zmon([x]));

if(~isfield(options, 'P_form'))
    options.P_form = 'Gram';
end

switch options.P_form
    case 'Gram'
        P = sdpvar(nm,nm);
        f_outer = zmon([x])'*P*zmon([x]);
        
        if(isfield(options, 'GramPSD') && options.GramPSD == false)
            con = [];
        else
            con = [ P >= 0];
        end
    case 'coeff'
        con = [];
        
        P = sdpvar(numel(zmon(x)),1);
        f_outer = dot(zmon(x),P);        
    otherwise
        error('Unrecognized');
        
end
param = P(:);

switch options.objective
    case 'geomean_S'
        objective = -geomean(P(2:end,2:end));
        set = [FormConstraint_LessEqualOneOnS(f_outer,x,S);];
    case 'geomean_dS'
        objective = -geomean(P);
        set = [FormConstraint_LessEqualOneOnDS(f_outer,x,dS)];
    case 'traceP_S'
        Pm = P(2:end,2:end);
        objective = -trace(Pm);
        set = [FormConstraint_LessEqualOneOnS(f_outer,x,S);];
    case 'AhmadiTraceGram_S'
        Pm = P(2:end,2:end);
        nv = size(Pm,1);
        
        V = sdpvar(nv,nv,'symmetric');
        con = [con; [V eye(nv); eye(nv) Pm] >= 0];
        param = [param; V(:)];
        objective = trace(V);
        
        set = [FormConstraint_LessEqualOneOnS(f_outer,x,S);];
    case 'geomeanHessian_S'
   
        [constraint H] = FormConstraint_SosCvx(f_outer,x,zmon);
        con = [con constraint];
        param = [param; H(:)];
        objective = -geomean(H);
        
        set = [FormConstraint_LessEqualOneOnS(f_outer,x,S);];
    case {'quasi_cvx_S', 'radial_S'}
        objective = -geomean(P);
        zi = zmon((1+options.kappa)*x);
        f_inner = zi'*P*zi;
        
        epsilon = options.epsilon;
        kappa = options.kappa;
        
        set = [FormConstraint_LessEqualOneOnS(f_outer,x,S);...
               FormConstraint_GreaterEqualOneOnDS(f_inner,x,dS,epsilon)];
    case {'quasi_cvx_dS', 'radial_dS'}
        objective = -geomean(P);
        zi = zmon((1+options.kappa)*x);
        f_inner = zi'*P*zi;
        
        epsilon = options.epsilon;
        kappa = options.kappa;

        set = [FormConstraint_LessEqualOneOnDS(f_outer,x,dS);...
               FormConstraint_GreaterEqualOneOnDS(f_inner,x,dS,epsilon)]; 
           
   case {'radial_S_Scomp_scaleX'}
        
        con = [con; replace(f_outer,x,zeros(nx,1)) <= options.max_value_at_origin];
        
        objective = 0;
        zi = zmon((1+options.kappa)*x);
        switch options.P_form
            case 'Gram'
                f_inner = replace(f_outer,x,(1+options.kappa)*x);
            case 'coeff'
                f_inner = replace(f_outer,x,(1+options.kappa)*x);
            otherwise
                error('Unrecognized');
        end
        
        epsilon = options.epsilon;
        kappa = options.kappa;
        
        set = [FormConstraint_LessEqualOneOnS(f_inner,kappa*x,S);...
           FormConstraint_GreaterEqualOneOnScomp_v2(f_inner,x,S,epsilon)];

        if isfield(options, 'fi_bound')
            set_bound_outer.lhs = f_outer-1;
            set_bound_outer.ci = options.fo_bound(x(1),x(2)) - 1;
            set_bound_outer.ce = [];
            set_bound_outer.z = x;
            
            set_bound_inner.lhs = 1 - f_inner;
            set_bound_inner.ci = 1 - options.fi_bound(x(1),x(2));
            set_bound_inner.ce = [];
            set_bound_inner.z = x;
            
            set = [set; set_bound_outer; set_bound_inner];
        end
    case {'radial_S_Scomp'}
        % Make f_inner be the "main" decision variable
        f_inner = f_outer;
        con = [con; replace(f_inner,x,zeros(nx,1)) <= options.max_value_at_origin];
        
        objective = 0;
        switch options.P_form
            case 'Gram'
                f_outer = replace(f_inner,x,x/options.kappa);
            case 'coeff'
                f_outer = replace(f_inner,x,x/options.kappa);
            otherwise
                error('Unrecognized');
        end
        
        epsilon = options.epsilon;
        kappa = options.kappa;
        
        set = [FormConstraint_LessEqualOneOnS(f_outer,x,S);...
           FormConstraint_GreaterEqualOneOnScomp_v2(f_inner,x,S,epsilon)];

        if isfield(options, 'fi_bound')
            set_bound_outer.lhs = f_outer-1;
            set_bound_outer.ci = options.fo_bound(x(1),x(2)) - 1;
            set_bound_outer.ce = [];
            set_bound_outer.z = x;
            
            set_bound_inner.lhs = 1 - f_inner;
            set_bound_inner.ci = 1 - options.fi_bound(x(1),x(2));
            set_bound_inner.ce = [];
            set_bound_inner.z = x;
            
            set = [set; set_bound_outer; set_bound_inner];
        end
  
    case 'l1_outer'
        x_lb = options.x_lb;
        x_ub = options.x_ub;
        
        objective = int(f_outer,x,x_lb,x_ub);
        
        set = [FormConstraint_PositiveOnB(f_outer,x,x_lb,x_ub);...
               FormConstraint_GreaterEqualOneOnS(f_outer,x,S);];
    case 'l1_inner'
        x_lb = options.x_lb;
        x_ub = options.x_ub;
        
        objective = int(f_outer,x,x_lb,x_ub);
        set = [FormConstraint_PositiveOnB(f_outer,x,x_lb,x_ub);...
            FormConstraint_L1_Inner(f_outer,x,x_lb,x_ub,Scomp,options.epsilon)];
    case 'geomean_S_inner'
        objective = -geomean(P);
        set = [FormConstraint_GreaterEqualOneOnScomp_v3(f_outer,x,Scomp,options.epsilon)];     
end

if options.quasi_cvx && ~strcmp(options.objective,'l1_inner') && ~strcmp(options.objective,'l1_outer')
    set = [set; FormConstraint_QuasiCvx(f_outer,x,zeros(nx,1))];
end

if options.sos_cvx
    [constraint H] = FormConstraint_SosCvx(f_outer,x,zmon);
    con = [con constraint];
    param = [param; H(:)];
end

P_offdiag = P - diag(diag(P));
options.sdpsettings = sdpsettings('solver','mosek','sos.model',1,'verbose',0,'sos.newton',1,'sos.congruence',1,'sos.csp',1,'sos.scale',0);
sosprog.set = set;
if isfield(options,'max_objective')
    % Pose as feasibility problem
    sosprog.obj = 0;
    if(~isnan(options.max_objective))
        sosprog.con = [con; objective <= options.max_objective];
    else
        sosprog.con = con;
    end
else
    sosprog.obj = objective;
    sosprog.con = con;
end
sosprog.param = param;
if(~isfield(options,'s'))
    options.s.dmax = max(degree(f_outer,x)) + 2;
    options.s.dmin = 0;
end

[result solver_info] = SolveSetContainSOS(sosprog,options);
result.obj = value(objective);
result.P = value(P);
solver_info;
P = value(P);
if nx == 2
    switch options.P_form
        case 'Gram'
            result.f0_outer = @(x,y) zmon([x;y])'*P*zmon([x;y]);
        case 'coeff'
            result.f0_inner = @(x,y) P'*zmon([x;y]);
            result.f0_outer = @(x,y) P'*zmon([x/options.kappa; y/options.kappa]);
        otherwise
            error('Unrecognized');
    end
    
else
    error('script only supports 2D currently due to post-processing & visualization needs');
end


% Attempt to clean up YALMIP
yalmip('clear')
end

function [constraint H] = FormConstraint_SosCvx(f,x,zmon)
    HessP = jacobian(jacobian(f,x)',x);
        
    d = degree(zmon(x)) - 1;
    y = sdpvar(numel(x),1);
    xmon = monolist(x,d);
    wmon = unique(kron(y,xmon));

    H = sdpvar(size(wmon,1),size(wmon,1),'symmetric');
    h1 = y'*HessP*y;
    h2 = wmon'*H*wmon;

    [c1 v1] = coefficients(h1,[x;y]);
    [c2 v2] = coefficients(h2,[x;y]);

    constraint = [coefficients(h1 - h2,[x;y]) == 0; H >= 0];         
end

function set = FormConstraint_PositiveOnB(f,x,x_lb,x_ub)
% Note: (Dabbene, 2017) suggests encoding as (x-xlb)(xub-x) >= 0.
% It is probably better to do this as two constraints that are linear: x>=xlb, x<= xub
% However, we follow their implementation.

% p(x) >=0 on B    
set.lhs = f;
set.ce = [];
for n = 1:numel(x)
    set.ci(n,1) = (x(n) - x_lb(n))*(x_ub(n) - x(n));
end
set.z = x;
end

function set = FormConstraint_GreaterEqualOneOnS(f,x,S)
for k = 1:numel(S)
    set(k,1).lhs = f-1;
    set(k,1).ce = []; 
    for ng = 1:numel(S(k).g)
        set(k,1).ci(ng,1) = [1 - S(k).g{ng}(x)];
    end
    if(~isfield(set(k,1),'ci'))
        % Add empty
        set(k,1).ci = [];
    end
    set(k,1).z = [x];
end
end

function set = FormConstraint_GreaterEqualOneOnScomp(f,x,x_lb,x_ub,S)
if(numel(S) > 1)
    error('Currently only support single set')
end

k = 0;
for ng = 1:numel(S.g)
    k = k+1;
    set(k,1).lhs = f - 1;
    set(k,1).ce = [];
    set(k,1).ci(1,1) = -(1-S.g{ng}(x));
    for n = 1:numel(x)
        set(k,1).ci(end+1,1) = (x(n) - x_lb(n))*(x_ub(n) - x(n));
    end
    set(k,1).z = [x];
end

end

function set = FormConstraint_GreaterEqualOneOnScomp_v2(f,x,S,epsilon)
if(numel(S) > 1)
    error('Currently only support single set')
end
    
k = 0;
for ng = 1:numel(S.g)
    k = k+1;
    set(k,1).lhs = f - (1+epsilon);
    set(k,1).ce = [];
    set(k,1).ci(1,1) = -(1-S.g{ng}(x));
    set(k,1).z = [x];
end

end



function set = FormConstraint_GreaterEqualOneOnScomp_v3(f,x,S,epsilon)

k = 0;
for i = 1:numel(S)
    k = k+1;
    set(k,1).lhs = f - (1+epsilon);
    set(k,1).ce = [];
    set(k,1).z = [x];
    for ng = 1:numel(S(i).g)
        set(k,1).ci(ng,1) = (1-S(i).g{ng}(x));     
    end
end
end

function set = FormConstraint_L1_Inner(f,x,x_lb,x_ub,Scomp,epsilon)


k = 0;
for i = 1:numel(Scomp)
    k = k+1;
    set(k,1).lhs = f - (1+epsilon);
    set(k,1).ce = [];
    set(k,1).z = [x];
    for ng = 1:numel(Scomp(i).g)
        set(k,1).ci(ng,1) = (1-Scomp(i).g{ng}(x));     
    end
    for n = 1:numel(x)
        set(k,1).ci(end+1,1) = (x(n) - x_lb(n))*(x_ub(n) - x(n));
    end
end
end


function set = FormConstraint_QuasiCvx(f,x,xc)
% Quasi-Star-Cvx
gf = jacobian(f,x);
set.lhs = gf*(x-xc);
set.ce = [];
set.ci = [];
set.z = [x];
end

function set = FormConstraint_LessEqualOneOnDS(f,x,dS)
% f_outer(x) <= 1 on dS
k = 0;
for i = 1:numel(dS)
    k = k+1;
    set(k,1).lhs = 1-f;
    set(k,1).ce = 1-dS(i).f(x);
    for ng = 1:numel(dS(i).g)
        set(k,1).ci(ng,1) = [1 - dS(i).g{ng}(x)];
    end
    if(~isfield(set(k,1),'ci'))
        % Add empty
        set(k,1).ci = [];
    end
    set(k,1).z = [x];
end
end

function set = FormConstraint_LessEqualOneOnS(f,x,S)
% f_outer(x) <= 1 on S
k = 0;
for i = 1:numel(S)
    k = k+1;
    set(k,1).lhs = 1-f;
    set(k,1).ce = [];
    for ng = 1:numel(S(i).g)
        set(k,1).ci(ng,1) = [1 - S(i).g{ng}(x)];
    end
    if(~isfield(set(k,1),'ci'))
        % Add empty
        set(k,1).ci = [];
    end
    set(k,1).z = [x];
end
end

function set = FormConstraint_GreaterEqualOneOnDS(f,x,dS,epsilon)
% f_inner(x) > 1 on dS
% Encode this as f_innner(x) >= 1+epsilon (epsilon small positive constant)
k = 0;
for i = 1:numel(dS)
    k = k+1;  
    set(k,1).lhs = f-(1+epsilon);
    set(k,1).ce = 1-dS(i).f(x);
    for ng = 1:numel(dS(i).g)
        set(k,1).ci(ng,1) = [1 - dS(i).g{ng}(x)];
    end
    if(~isfield(set(k,1),'ci'))
        % Add empty
        set(k,1).ci = [];
    end
    set(k,1).z = [x];
end
end

