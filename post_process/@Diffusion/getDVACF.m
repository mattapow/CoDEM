function obj = getDVACF(obj)
            % GETDVACF Get time dependent diffusion constant computed from
            % velocity autocorrelation function
            % D = int_{0}^{t} (1-s/t)< v(s).v(0)> ds
            % See Utter 2004

            if isnan(obj.nG)
                fprintf('Getting nGrains of PostP object\n')
                obj.nG = PostP.nGrains1(obj.dirPath);
            end

            % Compute the autocorrelation of the instantaneous velocities
            instantaneous = 1;
            [obj.acf, obj.lags] = getVAcf(obj, instantaneous);

            if ~isnan(obj.dvY)
            assert(eps(abs(obj.dvY^2 - obj.acf(1))) < eps(), 'Mismatched autocorrelation at t=0 compared to instantaneous velocity fluctuations')
            end
            
            % % Code to integrate only up to a truncation of lags and acf
            % lags = obj.lags;
            % acf = obj.acf;
            % sr = obj.shearRate;
            % lagsSr = lags*sr;
            % maxShearLag = Inf; % Integrating to less than infty misses data
            % acf = acf(lagsSr <= maxShearLag);
            % lags = lags(lagsSr <= maxShearLag);

            x = obj.lags;
            y = obj.acf.*(1-obj.lags/obj.lags(end));
            obj.DVACF = trapz(x, y);

        end