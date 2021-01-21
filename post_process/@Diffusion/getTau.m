function [obj, rmse] = getTau(obj)
    % obj = GETTAU(obj) Fit exp(-t/tau) to VACF
    % Only fit to first positive interval
    %
    % tau has units obj.shearRate

    % this should go in get.acf, get.lags
    if any(isnan(obj.acf)) || any(isnan(obj.lags))
        instantaneous = 1;
        [acf, lags] = getVAcf(obj, instantaneous);
        obj.acf = acf;
        obj.lags = lags;
    end

    firstNegIdx = find(obj.acf < 0, 1);
    maxShearLag = obj.lags(firstNegIdx - 1);

    shearLags = obj.lags * obj.shearRate;
    acf = obj.acf/obj.acf(1); % normalise
    acf = acf.*(1-obj.lags/obj.lags(end));

    acfS   = acf(shearLags <= maxShearLag); % Short interval
    shearLagsS = shearLags(shearLags <= maxShearLag);

    ft = fittype('exp(-x/t)');
    startPoint = [.2 .01 1 .001]';

    % ft = fittype('(t*(x-1))^-3');
    % startPoint = [.2 .01 1 .001]';


    % ft = fittype('cos(a*x)*exp(-x/t)');
    % warning('Fit type cos(a*x)*exp(-x/t). Nargout affected.')
    % if obj.cohesion < 8
    %     startPoint = .5*zeros(4,2);
    % else
    %     startPoint = .5*ones(4,2);
    % end
    % startPoint(:,2) = [.2 .01 1 .001]';

    try
        [f,gof] = fit(shearLagsS', acfS', ft, 'StartPoint', startPoint(1,:));
    catch
        try
            [f,gof] = fit(shearLagsS', acfS', ft, 'StartPoint', startPoint(2,:));
        catch
            try
                [f,gof] = fit(shearLagsS', acfS', ft, 'StartPoint', startPoint(3,:));
            catch
                try
                    [f,gof] = fit(shearLagsS', acfS', ft, 'StartPoint', startPoint(4,:));
                catch
                  try
                    [f,gof] = fit(shearLagsS', acfS', ft);
                  catch
                    error('Startpoint guess failed.')
                  end
                end
            end
        end
    end

    rmse = Inf;

    if exist('gof', 'var') && (gof.rsquare >= .8)
        obj.tau = f.t;
        rmse = gof.rmse;
        % a = f.a;
    else
        obj.tau = Inf;
        warning('Bad Fit.');
        % warning('r^2 = %s\n', flag)
    end

end
