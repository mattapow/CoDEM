% [eta, C, rmse] = GETVISCOTAU Get viscosity at each cohesion level
%
% For each cohesion level, fit the shear viscosity according to:
% tau = eta (shearRate) + tau_0
% Output ETA is an array of viscosities with corresponging cohesion array C
% Linear fits have root mean squared errors rmse (standard error)

function [eta, C, I, rmse, isWP] = getViscoTau(rheoOpt, isMu)

    % clf
    % hold on
    
    DP = simPaths(2);
    [m, n] = size(DP);

    eta = nan(m,n);
    C = nan(m,n);
    I = nan(m,n);
    rmse = nan(m,1); 
    isWP = nan(m, n);   

    ii = 1:m;
    jj = 1

    if length(ii)<2 && strcmp(rheoOpt, 'Bingham')
        warning('Need at least 2 inertial numbers to fit to Bingham rheology');
    elseif  length(ii)<3 && strcmp(rheoOpt, 'mu(I)')
        warning('Need at least 3 inertial numbers to fit to mu(I) rheology');
    end
    

    for j = jj

        mu1 = [];
        sr = [];

        for i = ii

            dp = DP{i,j};
            if isempty(dp); continue; end
            fp = strcat(dp, '/Rheology.mat');

            in = load(fp);
            rheology = in.rheology;

            % Values for fitting mu(I) or Netwonian rheology
            mu1=[mu1; rheology.stress(3)/-rheology.stress(4)];
            sr=[sr; rheology.shearRate];

            % Return values
            C(i,j) = rheology.cohesion;
            I(i,j) = rheology.shearRate;

        end

        % Fit shear stress (tau) to shear rate (sr) at this cohesion        
        switch rheoOpt
        case 'mu(I)'
            if length(mu1)<3; continue; end
            ft = 'a+(b-a)/(c/x+1)';
            [f, gof]=fit(sr, mu1, ft, 'StartPoint', [.3 .9 .3]);        
        case 'Bingham'
            if length(mu1)<2; continue; end
            [f, gof]=fit(sr, mu1, 'poly1');
        otherwise
            error('What is your rheology?')
        end


        rmse(j) = gof.rmse;

        % h1=plot(sr, mu1, 'o', 'color', [(j-1)/(length(jj)-1) .3 .1]);
        % set(h1, 'MarkerFaceColor', [(j-1)/(length(jj)-1) .3 .1])
        % if gof.rsquare > .6
        %     h2=plot(f);
        %     set(h2, 'color', [(j-1)/(length(jj)-1) .3 .1])
        %     set(h2, 'LineWidth', 1.5)            
        % else
        %     set(h1, 'LineStyle', '--')
        %     set(h1, 'LineWidth', 1)  
        %     set(h1, 'MarkerFaceColor', 'none')          
        % end

        % legend off

        % xlab=xlabel('$$\mathcal{I}$$');
        % ylab=ylabel('$$\mu = \frac{\tau_{yx}}{\sigma_{yy}}$$');
        % % ylab=ylabel('$$\frac{\sigma_{xy}+\sigma_{yx}}{\sigma_{xx}+\sigma_{yy}}$$');
        % set(xlab, 'Interpreter', 'latex')
        % set(ylab, 'Interpreter', 'latex')
        % set(ylab, 'Rotation', 0)
        % set(ylab, 'HorizontalAlignment', 'right')
        % set(xlab, 'FontSize', 20)
        % set(ylab, 'FontSize', 20)



        % Evaluate values on fitted function
        for i = ii

            dp = DP{i,j};
            if isempty(dp); continue; end
            fp = strcat(dp, '/Rheology.mat');

            in = load(fp);
            rheology = in.rheology;
            sr = rheology.shearRate;
            P = rheology.stress(4);
            t_i = 1;

            switch rheoOpt
            case 'mu(I)'
                % Using formula from \mu(I) rheology paper (Nature 2006) \eta = \mu(I) * P/sr
                if isMu 
%                     eta(i,j) = f(sr);
%                     eta(i,j) = f(sr)*-P/sr;
                    eta(i,j) = f(sr)*-P/sr /-P /t_i; % non-dimensional
                else
                    eta(i,j) = differentiate(f, sr);%*-P;
                end
            case 'Bingham'
                if isMu
                    eta(i,j) = f(sr); % bingham = poly1 fit
                else
                    eta(i,j) = differentiate(f, sr)*-P; % bingham = poly1 fit
                end
            end

            if strcmp(rheoOpt, 'mu(I)')
                % Is linearly well posed?
                isThisWP = Rheology.isWellPosed(f, sr);
                if isThisWP; str = '';
                else; str = 'not '; end
                fprintf('I=%.4f, C=%.1f: Rheology is linearly %swell posed.\n', rheology.shearRate, rheology.cohesion, str)
                isWP(i,j) = isThisWP;
            end


        end

        


    end


end
