function [visco, mass] = viscoEngine(obj, ts, tnsrCmp)

    % Average over all possible windows provided that there at least minAvg such windows
    minAvg = 40; % 20 files per 1 shear
    N = obj.tEnd - obj.nPass - (minAvg-1) - 1;
    visco = zeros(N,1);

    fprintf('Computing viscosity over %d window sizes...', N)
    fprintf('%6.0f/%6.0f', 0, N)

    % for each window size (difference dt = t2 - t1)
    for fnDt = 1:N
    
        if mod(fnDt, 500) == 0
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%6.0f/%6.0f', fnDt, N)
        end

        fn1 = 1;
        fn2 = fn1 + fnDt;

        if fn2 > N+1; continue; end


        % array for each start time at this window size
        nStarts = N-fnDt+1;
        theseVisco = zeros(nStarts, 1);        

        for k = 1:nStarts            

            % tnsrCmp(1)-component of position
            X2 = ts(fn2).grain.XTrue(:, tnsrCmp(1));
            X1 = ts(fn1).grain.XTrue(:, tnsrCmp(1));

            % tnsrCmp(2)-components of momentum
            % NB: unit density not unit mass
            P2 = ts(fn2).grain.V(:, tnsrCmp(2)) .* ts(fn2).grain.mass;
            P1 = ts(fn1).grain.V(:, tnsrCmp(2)) .* ts(fn1).grain.mass;

%             % tnsrCmp(2)-components of fluctuating component of momentum 
%             if tnsrCmp(2) == 1
%                 U2Mean = ts(fn2).grain.XTrue(:, 2)*obj.shearRate;
%                 U1Mean = ts(fn1).grain.XTrue(:, 2)*obj.shearRate;
%             else
%                 U2Mean = 0;
%                 U1Mean = 0;
%             end
%             P2 = (ts(fn2).grain.V(:, tnsrCmp(2)) - U2Mean) .* ts(fn2).grain.mass;
%             P1 = (ts(fn1).grain.V(:, tnsrCmp(2)) - U1Mean) .* ts(fn1).grain.mass;

            % Einstein-Helfand Moment (Helfand  1960 Eq. 3.14)
            theseVisco(k) = (sum(P2.*X2 - P1.*X1))^2;            

            % shift the window
            fn1 = fn1 + 1; fn2 = fn2 + 1;

        end

        % average over all start times for this interval    
        visco(fnDt) = mean(theseVisco);
    end

    assert(~any(visco==0))
    assert(~any(isnan(visco)))

    % Also return the mean grain mass
    mass = mean(ts(1).grain.mass);

    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bdone\n')

end