function [dX,t] = getLyapunov(dp,tStart,tEnd)
    
    xOnly = true;
    ts=TimeStep.loadTS(dp,tStart,tEnd,xOnly);

    % Average over all possible windows provided that there at least minAvg such windows
    minAvg = 1; % NB: 20,50 files per 1 shear
    N = numel(ts) - minAvg;
    dX = nan(N,2);
           
    fprintf('Computing |dY| over %d window sizes...', N)
    fprintf('%4.0f/%4.0f', 0, N)

    % for each window size (lag)
    for tau = 1:N

        fprintf('\b\b\b\b\b\b\b\b\b%4.0f/%4.0f', tau, N)
        fn0 = 1;
        fn1 = fn0 + tau;

        if fn1 > tEnd; continue; end

        % array for each start time at this window size
        nStarts = N-tau+1;
        theseDY = nan(nStarts, 2); 

        for t0=1:nStarts            
            
            % only get grains in contact at fn1=start
            ID_A=ts(fn0).contact.ID_A+1;
            ID_B=ts(fn0).contact.ID_B+1;
            
            
%             X1 = ts(fn0).grain.XTrue(:,:);
%             X2 = ts(fn1).grain.XTrue(:,:);
%             dX1 = ( (X1(ID_B,1)-X1(ID_A,1)).^2 - ...
%                     (X1(ID_B,2)-X1(ID_A,2)).^2 ).^.5;
%             dX2 = ( (X2(ID_B,1)-X2(ID_A,1)).^2 - ...
%                     (X2(ID_B,2)-X2(ID_A,2)).^2 ).^.5;
%             
%             dy = abs(dX2./dX1);
%             dy(dy==-Inf|dy==Inf)=NaN;
%             theseDY(t0,1) = mean(dy,'omitnan'); % mean over all grains
            
            % X,Y-dir
            for dir=1:2
                X1 = ts(fn0).grain.XTrue(:,dir);
                X2 = ts(fn1).grain.XTrue(:,dir);
                dX1=X1(ID_B)-X1(ID_A);
                dX2=X2(ID_B)-X2(ID_A);
                
                dy = abs(dX2./dX1);
                dy(dy==-Inf|dy==Inf)=NaN;
                theseDY(t0,dir) = mean(dy,'omitnan'); % mean over all grains
            end

            % shift the window
            fn0 = fn0 + 1; fn1 = fn1 + 1;

        end

        % mean over all start times for this interval
        dX(tau,1:2) = mean(theseDY,1); % DIM=1 for when nStarts=1
        
    end

    assert(~any(any(dX==0)))
    assert(~any(any(isnan(dX))))
    fprintf('\b\b\b\b\b\b\b\b\b\n\tdone\n')
    
    % Use time difference between two files to compute all lags
    fprintf('Getting corresponding times\n')
    dt = PostP.getDt1(dp);
    t = (1:N)'*dt;   
    fprintf('\tdone\n')

end