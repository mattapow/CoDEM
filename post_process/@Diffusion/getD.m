function obj = getD(obj)
            % obj = GETD(obj)
            % GETD Extract the diffusion coefficient Dy: msd = 2Dt
            % NB: This is not the true diffusion coefficient except in
            % uncorrelated motions?!

            % Need mean square diplacement readings
            if isnan(obj.dt)
                obj = obj.getMsd();
            end
            
            % This isn't quite working yet. Uncomment below to go back to
            % working version with fixed sample period.
%             obj = msdLimFit(obj);

%             if ~(isnan(obj.D))
%                 startPoint = [obj.D 1];
%             else
%                 startPoint = [.1 1];
%             end
%             % Fit msd(t) = 2Dt^v
%             ft = fittype('2*D*x^v');
%             selec = round(numel(obj.dt)/obj.top):numel(obj.dt);
%             [f,gof,~] = fit(obj.dt(selec), obj.msd(selec), ft, 'StartPoint', startPoint);
%             obj.D = f.D; obj.DError = NaN;
%             obj.v = f.v;
%             obj.rSquare = gof.adjrsquare;
%             fprintf('Adjusted R-squared: %f\n', obj.rSquare);

            if ~(isnan(obj.D))
                startPoint = obj.D;
            else
                startPoint = .1;
            end
            
%             % Fit msd(t) = 2Dt
            ft = fittype('2*D*x');
%             selec = round(numel(obj.dt)/obj.top):numel(obj.dt);
            T=numel(obj.dt);
            selec = round(T/obj.top):T;
            [f,gof] = fit(obj.dt(selec), obj.msd(selec), ft, 'StartPoint', startPoint);

            obj.D = f.D;
            % obj.DError = gof.rmse;
            obj.v = NaN;
            obj.rSquare = gof.adjrsquare;
            fprintf('Adjusted R-squared: %f\n', obj.rSquare);

            % Extract power fit
           nPara = 1;
           dfe = numel(obj.dt(selec))-nPara;
           Y1 = 2*f.D*obj.dt(selec);
           Y2 = obj.msd(selec);
           sse = sum(( Y1 - Y2 ).^2);
           rmse = sqrt(sse / dfe);
           obj.DError = rmse;
                      

end
        
function obj = msdLimFit(obj)    
    % Linear fit of for t large

    n = 1;
    N = length(obj.dt);
    
    % 1 Increase lower bound of interval until best fit
    increment=+1;
    removeTop=0;
    [rsquare,n, N] = fitInterval(obj, n, N, increment, removeTop);


    % Make sure the best fit is actually good and fit it
    if rsquare >= .8
        x = obj.dt(n:N);
        y = obj.msd(n:N)/2;
        D = x\y;

%         h=plot(x*obj.shearRate, msd*x, 'g-');
%         h.LineWidth = 4;
    else
        fprintf('Unable (r^2 = %f) to fit linear MSD in %s\n', rsquare, obj.dirPath)
        D = nan;
        rsquare = nan;                
    end

    fprintf('------------------\nt0=%.1f, t1=%.1f, r^2=%f\n------------------\n',n/20, N/20, rsquare)

    obj.DError = rsquare;
    obj.D = D;

end

function [rsquare, n, N] = fitInterval(obj, n, N, increment, removeTop)
    % [rsquare, n, N] = fitInterval(obj, n, N, increment, removeTop)
    % find best interval for a linear fit of
    % obj.dt(n0:N0)  ~ obj.visco(n0:N0)
    % 
    % if removeTop
    %   iterate N+=increment
    %  (in this case make increment negative)
    % if ~removeTop
    %   iterate n+=increment    
    % until r^2 not improving    
    
    % initial rsquare
    y = obj.msd(n:N);
    x = obj.dt(n:N);  
    D2 = x\y;
    yCalc = x*D2;        
    rsquare = 1 - sum((y - yCalc).^2)/sum((y - mean(y)).^2);

    fitImproving = 1;

    while fitImproving

        % If interval too short, no good fit found
        if N-n<2 || N < length(obj.dt) || n < 1
            % rsquare = nan;
            % fitImproving = 0;
            return
        end

        % Linear fit to this sub-sequence
        y = obj.msd(n:N);
        x = obj.dt(n:N);
        D2 = x\y;
        yCalc = x*D2;        
        thisRsquare = 1 - sum((y - yCalc).^2)/sum((y - mean(y)).^2);

        % If this fit isn't improving then stop
        if thisRsquare < rsquare
            fitImproving = 0;
        else
            rsquare=thisRsquare;

            if removeTop
                N=N+increment;
            elseif ~removeTop
                n=n+increment;
            else
                error('invalid removeTop boolean input')
            end
        end

    end
end