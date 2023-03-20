% Author: Jay Guthrie
% Description: Replicates Figure 4 of L-CSS paper

clear
clc
close all

addpath('./src')
warning off;

DefineColors;
generate = true;

k = 0; 
for degree = [6]
    for obj_type = {'kappa','det','tracePinv','l1'}
        k = k+1;
        Results{k}.objective = char(obj_type);
        Results{k}.degree = degree;
    end
end

mdl = 'Henrion2012_ex4p4';
save_filename = ['data_' mdl '.mat'];
switch mdl
    case 'Cerone2012'
        scale = 1; xc = zeros(2,1); xc = [1.3877; 0.3525];
        options.x_lb = scale*[-1 -1] + xc';
        options.x_ub = scale*[1 1] + xc';    
        options.x_lb = [-0.9 -0.4];
        options.x_ub = [0.6 1.3];
        options.x_lb = [-0.88 -0.3526];
        options.x_ub = [0.614 1.256];
        [S varargout] = GetExample(mdl);
        axis_bnds = [-2 2 -2 2];
    case 'Henrion2012_exG'
        scale = 1; xc = zeros(2,1); xc = zeros(2,1);
        options.x_lb = scale*[-1 -1];
        options.x_ub = scale*[1 1];    
        
        options.x_lb = [-0.876 -1];
        options.x_ub = [+0.876 +1];
        [S varargout] = GetExample(mdl);
        axis_bnds = [-1.1 1.1 -1.1 1.1];
    case 'Henrion2012_ex4p4'
        scale = 1; xc = zeros(2,1);
        options.x_lb = scale*[-0.8 -0.6];
        options.x_ub = scale*[0.6 1];
        
        axis_bnds = [-0.9 0.7 -0.7 1.1];
        [S varargout] = GetExample(mdl);
end

% Scale and translate
for i = 1:numel(S.g)
    S.g{i} = @(x) S.g{i}(1/scale*x+xc);
end

% Construct dS, Scomp
dS(1).dim = S.dim;
for i = 1:numel(S.g)
    Scomp(i).g{1} = @(x) (2)-S.g{i}(x); % Complement of S
    dS(i).f = Scomp(i).g{1};
    dS(i).g = [];
end

if generate
for k = 1:numel(Results)
    yalmip('clear');
    disp(['Objective: ' Results{k}.objective ' Degree: ' num2str(Results{k}.degree)]);
    options.s.dmin = 0;

    powers = monpowers(S.dim,0.5*Results{k}.degree)';
    xs = sym('xs',[S.dim 1]);
    zmon = matlabFunction( prod(xs.^powers).', 'Vars', {xs});
    options.s.dmax = Results{k}.degree;
    
    options.kappa = NaN;
    options.epsilon = 1e-6;
    options.max_value_at_origin = 1;
    options.quasi_cvx = false;
    options.sos_cvx = false;
    options.GramPSD = true;
    options.P_form = 'Gram';
    switch char(Results{k}.objective)
        case 'det'          
            options.objective = 'geomean_S';
            Results{k}.output = FitStarConvexSet(S,dS,zmon,options,Scomp);
        case 'tracePinv'
            options.objective = 'AhmadiTraceGram_S';
            Results{k}.output = FitStarConvexSet(S,dS,zmon,options,Scomp);
        case 'tracePinv_inner'
            options.objective = 'AhmadiTraceGram_S_Inner';
            Results{k}.output = FitStarConvexSet(S,dS,zmon,options,Scomp);
        case 'traceP'
            options.objective = 'traceP_S';
            Results{k}.output = FitStarConvexSet(S,dS,zmon,options,Scomp);
        case 'l1'
            options.objective = 'l1_outer';
            options.GramPSD = true;
            Results{k}.output = FitStarConvexSet(S,dS,zmon,options,Scomp);
         case 'l1_inner'
            options.objective = 'l1_inner';
            options.GramPSD = true;
            Results{k}.output = FitStarConvexSet(S,dS,zmon,options,Scomp);           
        case {'kappa', 'kappa_psd'}
            if(strcmp(char(Results{k}.objective),'kappa_psd'))
                options.GramPSD = true;   
            else
                options.GramPSD = false;
                options.P_form = 'coeff';
                powers = monpowers(S.dim,Results{k}.degree)';
                xs = sym('xs',[S.dim 1]);
                zmon = matlabFunction( prod(xs.^powers).', 'Vars', {xs});
            end
            

    
            kappa_lb = 1; kappa_ub = 1.5;
            options.objective = 'radial_S_Scomp';
            while(kappa_ub - kappa_lb > 1e-3)
                disp(['Bisection Interval: [' num2str(kappa_lb) ' ' num2str(kappa_ub) ']']);
                options.kappa = 0.5*(kappa_ub + kappa_lb);
                yalmip('clear');
                result = FitStarConvexSet(S,dS,zmon,options,Scomp);
                if(result.solver_info.problem == 0)
                    result_kappa = result;
                    kappa_ub = options.kappa;
                else
                    kappa_lb = options.kappa;
                end
            end
            Results{k}.output = result_kappa;
            Results{k}.output.kappa = kappa_ub;
    end
end

save(save_filename,'Results','S');
end

%% Plot    
load(save_filename)

figure(1)
t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
tile_order = [1:4];
result_order = [1:4];

theta = linspace(0,2*pi,1e3);
for i = 1:numel(S.g)
    G{i} = @(x,y) S.g{i}([x;y]);
end
[Set.x Set.y Set.r] = RadialLevelSet2D(G,theta);

area_X = CalcArea(Set.r, theta);
FontSize = 11;

output_string = 'Percent Error: ';
for k = 1:numel(Results)
    nexttile(tile_order(k))

   
    R = Results{result_order(k)};
    switch R.objective
        case 'det'
            fo = R.output.f0_outer;
            [SetOuter.x SetOuter.y SetOuter.r] = RadialLevelSet2D(fo,theta);
            fill(SetOuter.x, SetOuter.y,MyColors.SetOuter);
            hold on;
            fill(Set.x, Set.y,MyColors.Set);
            lam = eig(R.output.P(2:end,2:end));
            Pinv = inv(R.output.P(2:end,2:end));
            obj = (det(Pinv));
            [area] = CalcArea(SetOuter.r,theta);
            obj = det(R.output.P(2:end,2:end));
            base10_power = floor(log10(obj));
            value = obj/(10^base10_power);
            
            per_error = 100*(area - area_X)/area_X;
            objective_string = ['$\textup{det} P = ' num2str(value,'%10.2f') '\times 10^' num2str(base10_power) '$, Error = '  num2str(per_error,'%.1f') '\%'];
            title(['$\textup{det} P^{-1} = ' num2str(obj,'%.2f') '$' ' Area = ' num2str(area)],'interpreter','latex','FontSize',FontSize);
        case 'traceP'
            fo = R.output.f0_outer;
            [SetOuter.x SetOuter.y SetOuter.r] = RadialLevelSet2D(fo,theta);
            fill(SetOuter.x, SetOuter.y,MyColors.SetOuter);
            hold on;
            fill(Set.x, Set.y,MyColors.Set);
            Pm = R.output.P(2:end,2:end);
            obj = -trace(Pm);
            [area] = CalcArea(SetOuter.r,theta);
            area = area/area_X;
            objective_string = '$-\textup{tr} P$';
            values_string = ['$J = ' num2str(obj,'%.2f') '$, Percent Error = '  num2str(area,'%.2f')];
            title(['-$\textup{tr} P = ' num2str(obj,'%.2f') '$' ' Area = ' num2str(area)],'interpreter','latex','FontSize',FontSize);
        case 'tracePinv'
            fo = R.output.f0_outer;
            [SetOuter.x SetOuter.y SetOuter.r] = RadialLevelSet2D(fo,theta);
            fill(SetOuter.x, SetOuter.y,MyColors.SetOuter);
            hold on;
            fill(Set.x, Set.y,MyColors.Set);
            Pinv = inv(R.output.P(2:end,2:end));
            obj = trace(Pinv);
            [area] = CalcArea(SetOuter.r,theta);
            
            per_error = 100*(area - area_X)/area_X;
            objective_string = ['$\textup{tr} P^{-1}$ = ' num2str(obj,'%.2f')];
            values_string = ['Percent Error = '  num2str(per_error,'%.1f')];
            
            objective_string = ['$\textup{tr} P^{-1}$ = ' num2str(obj,'%.2f') ', Error = '  num2str(per_error,'%.1f') '\%'];
            title(['$\textup{tr} P^{-1} = ' num2str(obj,'%.2f') '$' ' Area = ' num2str(area)],'interpreter','latex','FontSize',FontSize);
          case 'tracePinv_inner'
            fo = R.output.f0_outer;
            [SetOuter.x SetOuter.y SetOuter.r] = RadialLevelSet2D(fo,theta);
            fill(Set.x, Set.y,MyColors.Set);
            hold on;
            fill(SetOuter.x, SetOuter.y,MyColors.SetOuter);
            
            Pinv = inv(R.output.P(2:end,2:end));
            obj = trace(Pinv);
            [area] = CalcArea(SetOuter.r,theta);
            area = area/area_X;
            per_error = 100*(areaOuter - area_X)/area_X;
            objective_string = ['$\textup{tr} P^{-1} = ' num2str(obj,'%.2f') '$'];
            values_string = ['Percent Error = '  num2str(per_error,'%.2f')];
            title(['$\textup{tr} P^{-1} = ' num2str(obj,'%.2f') '$' ' Area = ' num2str(area)],'interpreter','latex','FontSize',FontSize);
        case {'kappa', 'kappa_psd'}
            fo = R.output.f0_outer;
            [SetOuter.x SetOuter.y SetOuter.r] = RadialLevelSet2D(fo,theta);
            fill(SetOuter.x, SetOuter.y,MyColors.SetOuter);
            hold on;
            fill(Set.x, Set.y,MyColors.Set);   
            hold on;
            fi = R.output.f0_inner;
            [SetInner.x SetInner.y SetInner.r] = RadialLevelSet2D(fi,theta);
            % Uncomment to plot inner approximation
            % plot(SetInner.x,SetInner.y);
            % fill(SetInner.x, SetInner.y,MyColors.SetInner);    
            [areaOuter] = CalcArea(SetOuter.r,theta);
            [areaInner] = CalcArea(SetInner.r,theta);
            area = areaOuter;
            areaRatio = areaOuter/areaInner;
            expected_areaRatio = (R.output.kappa)^S.dim;
            assert(abs(areaRatio - expected_areaRatio) < 1e-3);
            per_error = 100*(areaOuter - area_X)/area_X;
            objective_string = ['$s = $ ' num2str(R.output.kappa,'%.2f')];
            values_string = ['Percent Error = '  num2str(per_error,'%.1f')];
            
            objective_string = ['$s = $ ' num2str(R.output.kappa,'%.2f') ', Error = '  num2str(per_error,'%.1f') '\%'];
        case 'l1'
            fo = R.output.f0_outer;
            figure(2)
            h_l1 = fcontour(fo,[-1 1 -1 1]*1.5, 'LevelList',[1],'MeshDensity',500,'LineColor',[1 0 0.5],'LineWidth',2.0);
            hold on;
            dF = h_l1.ContourMatrix(:,3:end);
            close(figure(2))

            figure(1)
            dx = 0.01; dy = 0.01;
            [area Xg Yg ind] = L1SetArea(fo,options,dx,dy);
            plot(Xg(ind),Yg(ind),'.','Color',MyColors.SetOuter,'MarkerSize',4); 
            hold on;
            c = colormap(lines(10));
            plot(dF(1,:),dF(2,:),'.','MarkerSize',3,'Color',c(1,:));
            
            plot([options.x_lb(1) options.x_ub(1) options.x_ub(1) options.x_lb(1) options.x_lb(1)],...
            [options.x_lb(2) options.x_lb(2) options.x_ub(2) options.x_ub(2) options.x_lb(2)], 'Color', 'k','LineWidth',2.0);
        
            fill(Set.x, Set.y,MyColors.Set);  
            
      
            per_error = 100*(area - area_X)/area_X;
            objective_string = ['$l_1$ = ' num2str(R.output.obj,'%.2f')];
            
            objective_string = ['$l_1$ = ' num2str(R.output.obj,'%.2f') ', Error = '  num2str(per_error,'%.1f') '\%'];
            values_string = ['Percent Error = '  num2str(per_error,'%.1f')];
            title(['$l_1 = ' num2str(R.output.obj,'%.2f') '$' ' Area = ' num2str(area)],'interpreter','latex','FontSize',FontSize);
      case 'l1_inner'
            
            figure(2)
            fi = R.output.f0_outer;
            h_l1 = fcontour(fi,[-1 1 -1 1]*1.5, 'LevelList',[1],'MeshDensity',500,'LineColor',[1 0 0.5],'LineWidth',2.0);
            hold on;
            dF = h_l1.ContourMatrix(:,3:end);
            close(figure(2))

            figure(1)
            fill(Set.x, Set.y,MyColors.Set); 
            hold on;
            
            dx = 0.01; dy = 0.01;
            [area Xg Yg ind] = L1SetArea(@(x,y) 2-fi(x,y),options,dx,dy);
            plot(Xg(ind),Yg(ind),'.','Color',MyColors.SetOuter); 
            hold on;
            plot(dF(1,:),dF(2,:),'r.','MarkerSize',2);
            
            plot([options.x_lb(1) options.x_ub(1) options.x_ub(1) options.x_lb(1) options.x_lb(1)],...
            [options.x_lb(2) options.x_lb(2) options.x_ub(2) options.x_ub(2) options.x_lb(2)], 'Color', 'b');
        
             
            area = area/area_X;
            objective_string = '$l_1$';
            values_string = ['$J = ' num2str(R.output.obj,'%.2f') '$, Percent Error = '  num2str(area,'%.2f')];
            title(['$l_1 = ' num2str(R.output.obj,'%.2f') '$' ' Area = ' num2str(area)],'interpreter','latex','FontSize',FontSize);
    
    end
    
    clear titleString
    if(tile_order(k) <= 5)
        titleString{1} = objective_string;
        titleString{2} = values_string;
    else
        titleString{1} = values_string;
    end
    titleString = objective_string;
    title(titleString, 'interpreter','latex','FontSize',FontSize);
    
    axis(axis_bnds);
    pbaspect([axis_bnds(2)-axis_bnds(1) axis_bnds(4)-axis_bnds(3) 1]);
    output_string = [output_string R.objective ' ' num2str(area/area_X, '%.4f') ' '];
end

FontSize = 10;
for i = 1:4
nexttile(i)
xlabel('$x_1$','interpreter','latex','FontSize',FontSize);
ylabel('$x_2$','interpreter','latex','FontSize',FontSize);
end


disp(output_string)
set(figure(1),'Position',[2183 1061 495 547]);
nexttile(2)
legend({'$\mathcal{F}_o$', '$\mathcal{X}$'},'interpreter','latex','FontSize',12,'NumColumns',1);

exportgraphics(figure(1),[mdl '_tiled.eps'],'BackgroundColor','none')

function [area] = CalcArea(R,theta)
    g = @(t) 0.5*interp1(theta,R,t,'linear').^2;
    area = integral(g,0,2*pi);
end