function visco = viscoEngine2(obj, tnsrCmp)
    
    xOnly = false;
    ts=TimeStep.loadTS(obj.dirPath,obj.nPass+1,obj.tEnd,xOnly);

    % Average over all possible windows provided that there at least minAvg such windows
    minAvg = 1; % 20 files per 1 shear %   minAvg = 40; % 20 files per 1 shear
%     N = obj.tEnd - obj.nPass - (minAvg-1) - 1;
%     warning('nPass not used')
%     N = obj.tEnd - (minAvg-1) - 1;
    N=numel(ts)-minAvg;
    visco = zeros(N,1);

    fprintf('Computing viscosity over %d window sizes...', N)
    fprintf('%6.0f/%6.0f', 0, N)
    
    % for each window size (difference dt = t2 - t1)
    for tau = 1:N
   
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%6.0f/%6.0f', tau, N)

        fn1 = 1;
        fn2 = fn1 + tau;

        if fn2 > obj.tEnd; continue; end

        % array for each start time at this window size
        nStarts = N-tau+1;
        theseVisco = zeros(nStarts, 1);        

        for t0 = 1:nStarts                        
            
            X1 = ts(fn1).grain.XTrue(:, tnsrCmp(1));
            P1 = ts(fn1).grain.V(:,tnsrCmp(2)) .* ts(fn1).grain.mass;

            X2 = ts(fn2).grain.XTrue(:, tnsrCmp(1));
            P2 = ts(fn2).grain.V(:,tnsrCmp(2)) .* ts(fn1).grain.mass;                        

            theseVisco(t0) = (sum(P2.*X2 - P1.*X1))^2; % Helfand moment     
            
            fn1 = fn1 + 1; fn2 = fn2 + 1; % shift the window

        end

        % average over all start times for this interval    
        visco(tau) = mean(theseVisco);
    end

    assert(~any(visco==0))
    assert(~any(isnan(visco)))

    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bdone\n')

end