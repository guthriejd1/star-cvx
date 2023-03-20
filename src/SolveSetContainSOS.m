function [result solver_info] = SolveSetContainSOS(sosprog,options)

obj = sosprog.obj;
con = sosprog.con;
param = sosprog.param;
set = sosprog.set;

verbose = false;
max_d = options.s.dmax; min_d = options.s.dmin;
for ns = 1:numel(set)
    if verbose
        disp(['Constructing set constraint ' num2str(ns) ' of ' num2str(numel(set))]);
    end
    
    prog = set(ns); 
    lhs = prog.lhs;
    for ni = 1:numel(prog.ci)
        z = prog.z;
        if isempty(z)
            s = sdpvar;
            c = s;
        else
            [s c m] = polynomial(z,max_d,min_d);  
        end
        con = [con; sos(s)];
        param = [param;c];
        lhs = lhs - s*( prog.ci(ni) );
        
        set(ns).params.ci{ni,1} = c;
    end
    for ne = 1:numel(prog.ce)
        z = prog.z;
       
        if isempty(z)
            p = sdpvar;
            c = p;
        else
            [p c m] = polynomial(z,max_d,min_d);  
        end
       
        param = [param;c];
        lhs = lhs + p*( prog.ce(ne) );

        set(ns).params.ce{ne,1} = c;
    end

    con = [con; sos(lhs)];

end

if ~isfield(options,'sdpsettings')
	options.sdpsettings = sdpsettings('solver','mosek','sos.model',2,'verbose',0,'sos.congruence',1,'sos.scale',1);
    options.sdpsettings.sos.newton = 1;
    options.sdpsettings.sos.csp = 1;
    options.sdpsettings.removeequalities = 0;     
    options.sdpsettings.sos.postprocess = 0;
end

[solver_info,v,Q,res] = solvesos(con,obj,options.sdpsettings,param);

result.param = value(param);
result.Q = Q;
result.obj = value(obj);
result.set = set;
result.solver_info = solver_info;

end

function [s c m] = polynomialCUSTOM(z,max_d,min_d)
c0 = sdpvar;
s = c0;

c = [c0];
m = 1;

for i = 1:numel(z)
    [si ci mi] = polynomial(z(i),max_d,max(min_d,1));
    c = [c;ci];
    s = s+si;
end
end