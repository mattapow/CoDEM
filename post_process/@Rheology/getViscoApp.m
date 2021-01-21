function [eta, C, I] = getViscoApp(isMu)
    % Apparent viscosity in each simulation

    allSims = 0;
    DP = simPaths(allSims, '', '');
    [m, n] = size(DP);

    str = {'x', 'v', 'o', 'p', 's', '+', '>', '*', 'd', 'v', 'o', 'p', 's'};

    mu = nan(m,n);
    C = nan(m,n);
    I = nan(m,n);
    rmse = nan(m,n);
    

    for i = 1:m

        for j = 1:n

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

            % eta(i,j) = rheology.getThisAppVisco;
            t_i = 1;
            P = -rheology.stress(4);
            sr = rheology.shearRate;
            if isMu
                eta(i,j) = rheology.getThisAppVisco / P * sr;
            else
                eta(i,j) = rheology.getThisAppVisco / P / t_i;
            end
            C(i,j)=rheology.cohesion;
            I(i,j)=rheology.shearRate;
            
        end        

    end

end



