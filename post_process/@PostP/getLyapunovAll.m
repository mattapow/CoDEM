function [dX,t] = getLyapunovAll(dp,tStart,tEnd)
%{
[dX,t] = getLyapunovAll(dp,tStart,tEnd)
Get the Lyapunov exponents for individual grains

The Lyapunov exponent is only defined for each orbit (grain path), so
instead of averaging across all the grains, why not see which grain have
big (well mixed) and small, perhaps even negative (i.e. non-chaotic)
orbits. This may reveal chaotic and non-chaotic regions.

Beacuse size(dX) = [nTimes,direction,nGrains] is large and slow to compute, 
it's saved to file in dp/dXAll.mat for future use. tStart and tEnd a called
when loading the TimeStep object, quite a big object.
%}

fp=strcat(dp,'/dXAll');
try
    in=load(fp);
    dX=in.dX;
    t=in.t;
catch e
    fprintf('Unable to load dX: \n\t%s\n', e.message)
    ts=TimeStep.loadTS(char(dp),tStart,tEnd);
    ts=noDrift(ts);

    % Average over all possible windows provided that there at least minAvg such windows
    minAvg = 1; % NB: 20,50 files per 1 shear
    N = tEnd - tStart - (minAvg-1);
    nG = PostP.nGrains1(dp);
    dX = zeros(N,2,nG);

    fprintf('Computing |dY| over %d window sizes...', N)
    fprintf('%4.0f/%4.0f', 0, N)

    % for each window size (lag)
    for tau = 1:N

        fprintf('\b\b\b\b\b\b\b\b\b%4.0f/%4.0f', tau, N)
        fn1 = 1;
        fn2 = fn1 + tau;

        if fn2 > tEnd; continue; end

        % array for each start time at this window size
        nStarts = N-tau+1;
        theseDY = zeros(nStarts,2,nG);
    
        for t0=1:nStarts            

            % only get grains in contact at fn1=start
            ID_A=ts(fn1).contact.ID_A+1;
            ID_B=ts(fn1).contact.ID_B+1;

            for g=1:nG
                I=find(ID_A==g|ID_B==g);

                X2 = ts(fn2).grain.XTrue(:,:);
                X1 = ts(fn1).grain.XTrue(:,:);
                dX2=abs(X2(ID_B(I),:)-X2(ID_A(I),:));
                dX1=abs(X1(ID_B(I),:)-X1(ID_A(I),:));
                dy=log(dX2)-log(dX1); % dy=log(dX2./dX1);

                dy(dy==-Inf)=NaN; dy(dy==Inf)=NaN;                
                theseDY(t0,:,g) = mean(dy,1,'omitnan'); % mean over initial contacts of a grain
            end

            % shift the window
            fn1 = fn1 + 1; fn2 = fn2 + 1;

        end

        % average over all start times for this interval
        dX(tau,1:2,1:nG) = mean(theseDY,1); % DIM=1 for when nStarts=1

    end
    
    fprintf('\b\b\b\b\b\b\b\b\b\n\tdone\n')

    % Use time difference between two files to compute all lags
    fprintf('Getting lag times\n')
    dt = PostP.getDt1(dp);
    t = (1:N)'*dt;   
    fprintf('\tdone\n')

    fprintf('Saving Output\n')
    save(fp, 'dX', 't')
end

assert(~any(any(any(dX==0))))
if any(any(any(isnan(dX))))
    warning('Some values of dX are NaN')
end

end