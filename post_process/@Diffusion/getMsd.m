function obj = getMsd(obj)

    %% Read in cell and grains at each timestep
    k = obj.tEnd-obj.nPass;

    fprintf('Reading grain true positions in file number...%4.0f/%4.0f', 0, obj.tEnd-obj.nPass)

    for filenum = obj.tEnd:-1:1+obj.nPass

        fprintf('\b\b\b\b\b\b\b\b\b%4.0f/%4.0f', obj.tEnd-k-obj.nPass, obj.tEnd-obj.nPass)

        % initialise Timestep class
        ts(k) = TimeStep(obj.dirPath);

        % Read cell time, length, and shift
        % Required for recomputing XTrue
        fp = strcat(obj.dirPath, '/cell/cell_', string(filenum));
        formatspec = '%f %f %f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f';
        ts(k).cell = ts(k).cell.readCell(fp, formatspec);


         % Read grain XTrue positions
        fp = strcat(obj.dirPath, '/grain/grain_', string(filenum));
        % formatspec = '%d %*f %*f %*f %*f %f %f %*f %*f %*f %*f';
        formatspec = '%d %f %f %*f %*f %f %f %*f %*f %*f %*f'; % add raw X to recompute X True
        ts(k).grain = ts(k).grain.readGrain(fp, obj.nG, formatspec);
        k=k-1;
    end

    fprintf('\b\b\b\b\b\b\b\b\bdone\n')

    fprintf('Recomputing XTrue\n')
    checkV = 0;
    ts = ts.recomputeXTrue(checkV);
    fprintf('Done\n')
    

    %% Average over all possible windows provided that there at least minAvg such windows
    minAvg = 40; % 20 files per 1 shear
    N = obj.tEnd - obj.nPass - (minAvg-1) - 1;
    MSD = zeros(N,1);
    
    fprintf('Computing MSD over %d window sizes...', N)
    fprintf('%4.0f/%4.0f', 0, N)

    % for each window size (difference dt = t2 - t1)
    for fnDt = 1:N

        fprintf('\b\b\b\b\b\b\b\b\b%4.0f/%4.0f', fnDt, N)

        fn1 = 1;
        fn2 = fn1 + fnDt;

        if fn2 > obj.tEnd; continue; end % if fn2 > N+1; continue; end


        % array for each start time at this window size
        nStarts = N-fnDt+1;
        theseMSD = zeros(nStarts, 1);        

        for k = 1:nStarts            

            % y-component of position
            X2 = ts(fn2).grain.XTrue(:, 2);
            X1 = ts(fn1).grain.XTrue(:, 2);
  
            % Position moment
            theseMSD(k) = mean( (X2 - X1).^2 ); % msd
%             theseMSD(k) = (sum(X2 - X1)).^2; % msd

            % shift the window
            fn1 = fn1 + 1; fn2 = fn2 + 1;

        end

        % average over all start times for this interval    
        MSD(fnDt) = mean(theseMSD);
    end

    assert(~any(MSD==0))
    assert(~any(isnan(MSD)))
    obj.msd = MSD;

    fprintf('\b\b\b\b\b\b\b\b\bdone\n')
    
    % Use time difference between two files to compute all lags
    fprintf('Getting corresponding MSD times\n')
    dt = PostP.getDt1(obj.dirPath);
    obj.dt = (1:N)'*dt;


end

