function obj = setViscoEin(obj)
% OBJ = SETVISCOEIN(OBJ) Get cell viscosity
%
% Compute the shear viscosity from Einstein-Helfand equation
% 
% n(t) = lim 1/(2kVTt) <[G(t)-G(0)]^2>
% 
%   where G(t) is the Helfand moment
%       G(t) = sum_grains p_x . y
%       where p_x is the grain momentum in the x-direction
%       and y is the position in the y-direction
    
    % Read in grains at each timestep    
    tic
    ts=TimeStep.loadTS(obj.dirPath,obj.nPass+1,obj.tEnd);
    fprintf('Read took %.4f seconds.\n', toc)
    ts = noDrift(ts);

    
    % Compute the viscosity helper array
    % Viscous tensor component. See https://doi.org/10.1080/08927022.2017.1321760
    tnsrCmp = [2,1];
    [visco, mass] = viscoEngine(obj, ts, tnsrCmp);    

    % Prefactors: Volume and energy terms
    E = .5 * mass * obj.dvY^2;
    V = obj.L(1) * obj.L(2) * 1;
    % E = getBeta(ts, mass);

    obj.visco =  visco/(2*E*V);
  


    % Use time difference between first two files to compute all lags
    dt1 = PostP.getDt1(obj.dirPath);
    N = length(obj.visco);
    obj.dt = (1:N)'*dt1;


    % Compute the viscosity as t -> infty
%     obj = viscoLimFit(obj);

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
    y = obj.visco(n:N);
    x = obj.dt(n:N);  
    eta = x\y;
    yCalc = x*eta;        
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
        y = obj.visco(n:N);
        x = obj.dt(n:N);
        eta = x\y;
        yCalc = x*eta;        
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


function obj = viscoLimFit(obj)    
    % Linear fit of for t large

    warning('Fitting lower section')
    n = 1;
    N = 10;
    
    % 1 Decrease upper bound of interval until best fit
    increment=-1;
    removeTop=1;
    [rsquare,n, N] = fitInterval(obj, n, N, increment, removeTop);
    rsquare
    
    % 1a Decrease upper bound of interval until best fit
    % This time start away from the top to avoid getting stuck in the jumpy tail
    N1=5;
    increment=-1;
    removeTop=1;
    [rsquare1, n1, N1] = fitInterval(obj, n, N1, increment, removeTop);
    if rsquare1 > rsquare
        rsquare=rsquare1;        
        n=n1;
        N=N1;
    end
    
    % 1b Decrease upper bound of interval until best fit
    % This time start away from the top to avoid getting stuck in the jumpy tail
    N1=3;
    increment=-1;
    removeTop=1;
    [rsquare1, n1, N1] = fitInterval(obj, n, N1, increment, removeTop);
    rsquare1
    if rsquare1 > rsquare
        rsquare=rsquare1;        
        n=n1;
        N=N1;
    end
    
    
    
    
%     n=20;
%     N=length(obj.dt);    
% 
%     % 1 Increase lower bound of interval until best fit    
%     increment=1;
%     removeTop=0;
%     [rsquare,n, N] = fitInterval(obj, n, N, increment, removeTop);
% 
%     % 2a Decrease upper bound of interval until best fit
%     increment=-1;
%     removeTop=1;
%     [rsquare,n, N] = fitInterval(obj, n, N, increment, removeTop); 
% 
%     % 2b Decrease upper bound of interval until best fit
%     % This time start away from the top to avoid getting stuck in the jumpy tail
%     N1=length(obj.dt) - 20;
%     increment=-1;
%     removeTop=1;
%     [rsquare1, n1, N1] = fitInterval(obj, n, N1, increment, removeTop);    
%     % eta or eta1 may be best, find the best
%     if rsquare1 > rsquare
%         rsquare=rsquare1;        
%         n=n1;
%         N=N1;
%     end
% 
%     % 2c Decrease upper bound of interval until best fit
%     % This time start away from the top to avoid getting stuck in the jumpy tail
%     N1=length(obj.dt) - 50;
%     increment=-1;
%     removeTop=1;
%     [rsquare1, n1, N1] = fitInterval(obj, n, N1, increment, removeTop);    
%     % eta or eta1 may be best, find the best
%     if rsquare1 > rsquare
%         rsquare=rsquare1;        
%         n=n1;
%         N=N1;
%     end
% 
%     % 3a Increase lower bound of interval until best fit
%     n1=20;
%     increment=1;
%     removeTop=0;
%     [rsquare1,n1, N1] = fitInterval(obj, n1, N, increment, removeTop);
%     if rsquare1 > rsquare
%         rsquare=rsquare1;        
%         n=n1;
%         N=N1;
%     end
% 
% 
%     % 3b Increase lower bound of interval until best fit
%     n1=50;
%     increment=1;
%     removeTop=0;
%      [rsquare1,n1, N1] = fitInterval(obj, n1, N, increment, removeTop);
%     if rsquare1 > rsquare
%         rsquare=rsquare1;        
%         n=n1;
%         N=N1;
%     end
% 
%     % 3c Increase lower bound of interval until best fit
%     n1=100;
%     increment=1;
%     removeTop=0;
%      [rsquare1,n1, N1] = fitInterval(obj, n1, N, increment, removeTop);
%     if rsquare1 > rsquare
%         rsquare=rsquare1;        
%         n=n1;
%         N=N1;
%     end
% 
%     % 3d Decrease lower bound of interval until best fit
%     n1=100;
%     increment=-1;
%     removeTop=0;
%      [rsquare1,n1, N1] = fitInterval(obj, n1, N, increment, removeTop);
%     if rsquare1 > rsquare
%         rsquare=rsquare1;        
%         n=n1;
%         N=N1;
%     end


    % Make sure the best fit is actually good and fit it
    if rsquare >= .8
        x = obj.dt(n:N);
        y = obj.visco(n:N);
        eta = x\y;

        h=plot(x*obj.shearRate, eta*x, 'g-');
        h.LineWidth = 4;
    else
        fprintf('Unable (r^2 = %f) to fit linear viscosity in %s\n', rsquare, obj.dirPath)
        eta = nan;
        rsquare = nan;                
    end

    fprintf('------------------\nt0=%.1f, t1=%.1f, r^2=%f\n------------------\n',n/20, N/20, rsquare)

    obj.viscoLim = eta;
    obj.viscoLimRmse = rsquare;

end

function E = getBeta(ts, mass)
    % Determine the energy term in the Einstein-Helfand relation
    % Use E = m/<p^2> as in Appendix A of VISCOSITY AND MICROSCOPIC CHAOS : THE HELFAND-MOMENT APPROACH
    % Sebastien Viscardy

    n = length(ts);
    EX = nan(n,1);
    EY = nan(n,1);

    for i = 1:n
        EX(i) = mass / mean( (ts(i).grain.mass.*ts(i).grain.V(:,1)).^2 );
        EY(i) = mass / mean( (ts(i).grain.mass.*ts(i).grain.V(:,2)).^2 );
    end

    E = mean( [mean(EX), mean(EY)] );

end
