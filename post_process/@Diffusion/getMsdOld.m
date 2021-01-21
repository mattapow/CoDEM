function obj = getMsdOld(obj)
            % obj = GETMSD(obj)
            % Calculate the MSD(t) (mean square vertical displacement)
            % Computes using XTrue, each time interval is averaged over
            % a number of start times
            %
            % NB: fn - file number, fp - file path
            %
            % Skip obj.nPass number of files at the start

            fprintf('Getting mean squared displacement recordings\n')

            % Presize Mean time change and Mean MSD
            fnEnd = obj.tEnd; %#ok<*PROP>
            obj.maxT  = 2; % Careful. This will yield tEnd/2
            fnNStarts = obj.maxT;
            fnNIntervals = obj.maxT;
            obj.dt = zeros(fnNStarts, 1); %#ok<*PROPLC>
            obj.msd = zeros(fnNStarts, 1);


            % Read in cell and grains at each timestep
            k = fnEnd-obj.nPass;

            fprintf('Reading grain true positions in file number...%4.0f/%4.0f', 0, fnEnd-obj.nPass)

            for filenum = fnEnd:-1:1+obj.nPass

                fprintf('\b\b\b\b\b\b\b\b\b%4.0f/%4.0f', fnEnd-k-obj.nPass, fnEnd-obj.nPass)

                % initialise Timestep class
                ts(k) = TimeStep(obj.dirPath);

                % Read cell time and length
                % Only required for recomputing XTrue
                fp = strcat(obj.dirPath, '/cell/cell_', string(filenum));
                formatspec = '%f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
                ts(k).cell = ts(k).cell.readCell(fp, formatspec);


                 % Read grain XTrue positions
                fp = strcat(obj.dirPath, '/grain/grain_', string(filenum));
                % formatspec = '%d %*f %*f %*f %*f %f %f %*f %*f %*f %*f';
                formatspec = '%d %f %f %*f %*f %f %f %*f %*f %*f %*f'; % add raw X to recompute X True
                ts(k).grain = ts(k).grain.readGrain(fp, obj.nG, formatspec);
                k=k-1;
            end

            fprintf('\b\b\b\b\b\b\b\b\bdone\n')

            % X = ts(1).grain.X;
            % X_1 = ts(end).grain.X;
            % XTrue0 = ts(1).grain.XTrue;
            % XTrue0_1 = ts(end).grain.XTrue;
            fprintf('Recomputing XTrue\n')
            ts = ts.recomputeXTrue;
            fprintf('Done\n')
            % XTrue1 = ts(1).grain.XTrue;
            % XTrue1_1 = ts(end).grain.XTrue;
            % disp([X(1:4,:), XTrue0(1:4, :), XTrue0(1:4, :)-X(1:4,:), XTrue1(1:4, :)])
            % X_1(1:4,:), XTrue0_1(1:4, :), 
            % disp([X(1:4, 2), X(1:4, 2)*obj.shearRate XTrue1(1:4, 2)-XTrue0(1:4, 2), XTrue1_1(1:4, 2)-XTrue0_1(1:4, 2)])


            % Compute the square displacements
            fprintf('Computing MSD over %d delays and %d start times...', fnNIntervals, fnNStarts)
            fprintf('%4.0f/%4.0f', 0, fnNIntervals)

            % for each time interval (difference dt = t2 - t1)
            for fnDt = 1:fnNIntervals

                fprintf('\b\b\b\b\b\b\b\b\b%4.0f/%4.0f', fnDt, fnNIntervals)

                fn1 = 1;
                fn2 = fn1 + fnDt;

                % tmp arrays for each start time at this dt
                msd_t = zeros(fnNStarts, 1);

                % For each start time (indexed by k)
                % k = 1;
                % for k = 1+obj.nPass:obj.maxT
                for k = 1:fnNStarts

                    % y-components of XTrue
                    X2 = ts(fn2).grain.XTrue(:, 2);
                    X1 = ts(fn1).grain.XTrue(:, 2);
                    msd_t(k) = mean((X2 - X1).^2);

                    % shift the interval
                    fn1 = fn1 + 1;
                    fn2 = fn2 + 1;
                    % k = k + 1;

                end

                % average over all start times for this interval
                obj.msd(fnDt) = mean(msd_t);

            end
            fprintf('\b\b\b\b\b\b\b\b\bdone\n')

            % Use time difference between two files to compute all lags
            dt = PostP.getDt1(obj.dirPath);
            obj.dt = (1:fnNStarts)'*dt;

        end