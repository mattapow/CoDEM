% function [eta, C, I, rmse] = getViscoGK(plotIt)
function [x,y] = getViscoGK(plotIt)

% 	allSims = 0;
%     DP = simPaths(allSims, '', '');
%     [m, n] = size(DP);

%     dp = '/Users/mattmacaulay/Documents/DEM/DATA/N10000/gm0.1/C0/runFine1';
%     dp = '/Users/mattmacaulay/Documents/DEM/DATA/N10000/gm0.3/C0/run';
    dp = '/Users/mattmacaulay/Documents/DEM/DATA/N10000_1/gm0.1/C0/runExtraFine';
    m=1;
    n=1;

    eta = nan(m,n);
    C = nan(m,n);
    I = nan(m,n);
    rmse = nan(m,n);

    clf
    hold on
	lineStyle = {'-', '--', '-.', ':'};
    L = {};
    jj = 1;

    for i = 1
        for j = jj

%             dp = DP{i,j};
%             if isempty(dp); continue; end
            fp = strcat(dp, '/Rheology.mat');
            fprintf('%s\n', fp)

            try 
                in = load(fp);
                rheology = in.rheology;
            catch
                rheology = Rheology(dp);
                stress = rheology.stress; % includes saving
            end
                        
            rheology = rheology.setViscoGKViscardy();
            save(fp, 'rheology')

            t_i = 1;
            P = -rheology.stress(4);
            
            eta(i,j) = rheology.viscoLimGK / t_i / P;
            rmse(i,j) = rheology.viscoLimRmse / t_i / P;
        
            
            C(i,j)=rheology.cohesion;
            I(i,j)=rheology.shearRate;
%             rheology.viscoGK = rheology.viscoGK .* linspace(1, 0, length(rheology.viscoGK))';
%             rheology.viscoLimGK = trapz(rheology.dtGK, rheology.viscoGK);

            if nargin == 1 && plotIt
                x = rheology.dtGK * rheology.shearRate;
                y = rheology.viscoGK;
                h1 = plot(x, y, '-');
%                 set(h1, 'Color', [0 .3 (i-1)/(m-1)])
                set(h1, 'LineWidth', 1.4)
                set(gca, 'XLim', [0 22])
                set(h1, 'LineStyle', lineStyle{mod(i-1,4)+1})
                drawnow
            end
            

            % if rheology.cohesion == .5
            %     L{end+1} = sprintf('C=%0.1f', rheology.cohesion);  
            % else
            %     L{end+1} = sprintf('C=%0.0f', rheology.cohesion);  
            % end
        end        
        L{end+1} = sprintf('I=%0.3f', rheology.shearRate);

    end   

    if nargin == 1 && plotIt
        leg=legend(L);
        leg.FontSize = 16;
        leg.Location = 'northeast';
        leg.Interpreter = 'latex';
        ylabel('$$\frac{L^x L^y d}{m (\delta v^y)^2}\langle \tau(0) \tau(t) \rangle$$', 'Interpreter', 'latex')
        ylabel('$$\langle \tau(0) \tau(t) \rangle$$', 'Interpreter', 'latex')
        xlabel('$$t \dot \gamma$$', 'Interpreter', 'latex')
        set(gca, 'FontSize', 16)
    end
end


