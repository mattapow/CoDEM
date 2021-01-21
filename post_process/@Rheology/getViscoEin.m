function [eta, C, I, rmse] = getViscoEin(isMu)
    % NB. rmse holds the gof.rsquare value

    DP = simPaths(1.3);
    [m, n] = size(DP);
    str = {'x', 'v', 'o', 'p', 's', '+', '>', '*', 'd', 'v', 'o', 'p', 's'};

    mu = nan(m,n);
    C = nan(m,n);
    I = nan(m,n);
    rmse = nan(m,n);

    clf
    hold on
    L = {};
    jj = 1;

    for i = m:-1:1
        for j = jj

            dp = DP{i,j};
            if isempty(dp); continue; end
            fp = strcat(dp, '/Rheology.mat');
            fprintf('%s\n', fp)

            try 
                in = load(fp);
                rheology = in.rheology;
            catch
                rheology = Rheology(dp);
                stress = rheology.stress; % includes saving
            end

%             rheology = rheology.setViscoGK();
            rheology = rheology.setViscoEin();
            save(fp, 'rheology')

            t_i = 1;
            sr = rheology.shearRate;
            P = -rheology.stress(4);
            if isMu 
                eta(i,j) = rheology.viscoLim * sr / P;
                rmse(i,j) = rheology.viscoLimRmse * sr / P;
            else
                eta(i,j) = rheology.viscoLim / t_i / P;
                rmse(i,j) = rheology.viscoLimRmse / t_i / P;
            end
            
            C(i,j)=rheology.cohesion;
            I(i,j)=rheology.shearRate;
            
            % if rheology.viscoLimRmse > .98
                h1 = plotViscoT(rheology);
%                 % set(h1, 'Color', [(j-1)/(length(jj)-1) .3 .1])
                set(h1, 'Color', [(i-1)/(m-1) .3 .1])
                set(h1, 'LineWidth', 2)
                set(gca, 'XLim', [0 22])
                set(h1, 'Marker', 'none')
            % end
            drawnow

            % if rheology.cohesion == .5
            %     L{end+1} = sprintf('C=%0.1f', rheology.cohesion);  
            % else
            %     L{end+1} = sprintf('C=%0.0f', rheology.cohesion);  
            % end
        end        
        % L{end+1} = sprintf('%0.3f fit', rheology.shearRate);
%         L{end+1} = sprintf('I=%0.3f', rheology.shearRate);
        

    end
    

    % ax=gca;
    % ax.XLim = [0 21];

    leg=legend(L);
    leg.FontSize = 16;
    leg.Location = 'northwest';
    leg.Interpreter = 'latex';

    % set(gca, 'XScale', 'log')
    % set(gca, 'YScale', 'log')
    % set(gca, 'YLim', [1e-10 1e-3])


end



