function [acf, lags] = getVAcf(obj, instantaneous)
            % GETVACF Get the velocity autocorrelation function
            %
            % Time averaged over all grains
            % lags - lag times
            % acf - non-normalised mean (over grains) autocorrelation
            % bool instantaneous set to 1 to use instantaneous velocities, see getDVAll
            %
            % Computation using xcorr from time nPass to tEnd.
            % Updated from using econ:autocorr function which gave bad results
            warning('acf of norm of velocity fluctuations')

            % Number of lags to compute acf
            nFiles = obj.tEnd-obj.nPass-1;
            numLags = nFiles-(1-instantaneous);

            % compute lags from evenly space timesteps dt
            dt = PostP.getDt1(obj.dirPath);
            lags = (0:numLags)*dt;


            if isnan(obj.nG)
                fprintf('Getting nGrains of PostP object\n')
                obj.nG = PostP.nGrains1(obj.dirPath);
            end


            % Get grain velocity fluctuations at all times
            [dU, dV] = getDvAll(obj, instantaneous);            

            fprintf('Using dvY = mean(sqrt(mean(dV.^2, 2)))\n')
            dvY = mean(sqrt(mean(dV.^2, 2)));

            if isnan(obj.dvY)
                obj.dvY = dvY;
            else
                err = dvY - obj.dvY < eps();
                assert(err, 'obj.dvY does not match dV calculated here')
            end

            % compute autocorrelation of velocty fluctuations
            acfs = zeros(obj.nG, numLags+1);

            % for each grain
            for g = 1:obj.nG

                % this grain's fluctuations over time
                gDV = dV(:,g); % gDV = (dU(:,g).^2+dV(:,g).^2).^.5;
                gDV = gDV - mean(gDV);

                % use autocorr=xcorr
                % not autocorr=ifft(fft(X)*conj(fft(X))) which had spectal leakage issues
                C = xcorr(gDV);
                
                if mod(length(C),2) == 0
                    C = C(length(C)/2:end);
                else
                    C = C((length(C)+1)/2:end);
                end
                C = C(1:numLags+1);
                C = C/C(1);
                             
                acfs(g,:) = C;

            end


            acf = mean(acfs, 1); % Time averaged acf. No grain averaged
            acf = acf*obj.dvY^2; % scale by velocity fluctuations
            
            
            

        end