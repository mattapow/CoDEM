function f = muVsI(DP, showPlot)
% Plot the mu(I) rheology

[m, n] = size(DP);

str = {'x', 'v', 'o', 'p', 's', '+', '>', '*', 'd', 'x', 'v', 'o', 'p', 's'};

if showPlot
%% Pre-Formatting
fig=gcf; clf
set(fig, 'WindowStyle', 'normal');
set(fig, 'Units', 'pixels');
set(fig, 'OuterPosition', [100 100 700 600])
set(fig, 'InnerPosition', [100 100 700 600])
set(fig, 'Visible', 'on');
set(fig, 'Color', '[1 1 1]')

ax = axes();
set(ax, 'Box', 'on')
% set(ax, 'XTick', [])
set(ax, 'Units', 'pixels')
set(ax, 'Position', [80 60 600 500]);

ax.FontSize = 20;    
%     ax.YLim = [.1 1];
% ax.YTick = .25:.25:1.75;    
% ax.YMinorTick = 'on';
ax.XLim = [0 .95];
% ax.XScale = 'log';
% ax.XTick = [.005 .01 .05 .1 .2 .3];
% ax.XTick = [0.1,0.5,1,2,3];

hold on
end
mu = []; I = []; 
d=1; rho_g=1;

for i=1:m
    j=1;
%     for j=1:n    
        rheology = Rheology.loadRheology(DP, i, j);
        if isnumeric(rheology)
            continue
        end
        

        if showPlot
%         color = [(j-1)/(n-1) .2 .3];
        color = [0 0 0];
        plot(rheology.shearRate, rheology.stress(3)/-rheology.stress(4), 'Marker', str{j}, 'Color', color, 'MarkerSize', 10)                                           
        end
%         if j==1
            P=-(rheology.stress(4)+rheology.stress(1))/2;
            tau=(rheology.stress(2)+rheology.stress(3))/2;
%             P=-rheology.stress(4);
%             tau=rheology.stress(2);
            mu(i) = tau/P; %#ok<AGROW>
            gm = rheology.shearRate;
            t_i = d*sqrt(rho_g/P);
            I(i) = gm*t_i; %#ok<AGROW>            
%         end
    
%         LH(n-j+1) = plot(nan, nan, 'Color', color', 'Marker', str{j}, 'LineStyle', 'none'); %#ok<AGROW>
%         drawnow
%     end
end
    

%% Fit mu(I) cohesionless
ft = 'a+(b-a)/(c/x+1)';
[f, gof, ~] = fit(I', mu', ft, 'StartPoint', [.3 .9 .3]);
% ft = 'poly1';
% [f, gof, ~] = fit(I', mu', ft);
if showPlot
LH(n+1) = plot(f, '-k');

xlab = xlabel('$$\mathcal{I}$$');
xlab.Interpreter = 'Latex';
xlab.FontSize = 24;

ylab = ylabel('$$\mu$$');
ylab.Interpreter = 'Latex';
ylab.FontSize = 24;
ylab.Rotation = 0;
ylab.HorizontalAlignment = 'right';

% L = {'C=30', 'C=25', 'C=20', 'C=15', 'C=10', 'C=7', 'C=5', 'C=2', 'C=1', 'C=0.5', 'C=0', 'C=0 fit'};
% leg = legend(LH, L);
leg = legend();
leg.Position = [.75 .75 .19 .043];
ttl = title('$$\mu(\mathcal{I}, C=0)$$');
ttl.Interpreter = 'Latex';

leg = legend();
delete(leg);
end



end