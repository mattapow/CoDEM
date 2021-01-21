classdef ProfileP < PostP
    %ProfileP Post Process profile
    %   Any sort of space profile space-time profile
    %   Also shear related functions (possible to extract this as a class)

    properties

        vXSlice = NaN % V_x of corresponding slice normalised by wall speed
        ySlice = NaN % y of corresponding slice normalised by cell height
        times = NaN % Averaging times used. -1 for all

        % Shear rate profiles
        % These have already been normalised as per vXSlice and ySlice
        dVdy = NaN % Shear rate in slices
        sbWidth = NaN % List of Widths of shear bands
        sbTime = NaN % List of shear band time persistance
        plugWidth = NaN

        inhomFft = NaN % Whether the max power of time frequency is more than 1/2H


    end

    methods

        % Constructor
        function obj = ProfileP(dp)
            obj = obj@PostP(dp);
        end

%         function S = saveobj(obj)
%             % Don't save the super (commented out)
%             % S = saveobj@PostP(obj);
%
%             % Save the sub
%             S.vXSlice = obj.vXSlice;
%             S.ySlice = obj.ySlice;
%             S.times = obj.times;
%             S.dVdy = obj.dVdy;
%             S.sbWidth = obj.sbWidth;
%             S.sbTime = obj.sbTime;
%             S.plugWidth = obj.plugWidth;
%         end
%
%         function obj = reload(obj, S)
%
%             % Super from saved
%             fp = strcat(obj.dirPath, '/PostP.mat');
%             in = load(fp);
%             varObj = in.varObj;
%             obj = reload@PostP(obj,varObj);
%
%             % Sub
%             obj.vXSlice = S.vXSlice;
%             obj.ySlice = S.ySlice;
%             obj.times = S.times;
%             obj.dVdy = S.dVdy;
%             obj.sbWidth = S.sbWidth;
%             obj.sbTime = S.sbTime;
%             obj.plugWidth = S.plugWidth;
%          end


        %% Convergence
        function isConverged(obj, complete, saveIt)
            % if(complete) then also plot space-time diagrams of density
            % and shear rate. Takes a bit longer and must have profiles
            % if(saveIt) then save the plots to the PostP directory path as
            % a .tif file

            clf
            for i = obj.tEnd:-1:1

                % Read cell: phi and stress
                fpCell = strcat(obj.dirPath, '/cell/cell_', string(i));
                formatspec = '%f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %*f %*f';
                cell = Cell();
                cell = cell.readCell(fpCell, formatspec);

                phi(i) = cell.phi;
                sigma(i) = cell.stress(4);
                tau(i) = cell.stress(2);
                time(i) = cell.time;

            end

            nPlot = 2;
            if (nargin > 1) && (complete == 1)
                nPlot = 3;
            end
            disp(obj.inhomX)
            disp((phi(end)-phi(61))/mean(phi(61:end)))
            % Plot of mean density in time
            subplot(1,nPlot,1)
            line(time,phi)

            xlab = xlabel('$$t$$');
            ylab = ylabel('$$\phi$$');
            set(xlab, 'Interpreter','latex');
            set(ylab, 'Interpreter','latex');
            axis tight
            ylim([.55 .85])


            % Plot of normal and shear stresses in time
            subplot(1,nPlot,2)
            line(time ,sigma,'Color','r')

            axis tight
            ax1 = gca;
            ax1.YColor = 'r';
            ax1_pos = ax1.Position;
            xlab = xlabel('$$t$$');
            set(xlab, 'Interpreter','latex');
            ylab = ylabel('$$\sigma_{yy}$$');
            set(ylab, 'Interpreter','latex');
            set(xlab, 'Interpreter','latex');
            ax2 = axes('Position',ax1_pos, 'YAxisLocation','right','Color','none');


            % Plot of and shear stresses in time
            line(time, tau, 'Parent', ax2, 'Color', 'k')
            ylab = ylabel('$$\sigma_{xy} = \tau$$');
            set(ylab, 'Interpreter','latex');
            axis tight


            if (nargin > 1) && (complete == 1)
                % Plot of shear rate in space-time
                subplot(1,nPlot,3)
                obj.plotShearBandsContour();

%                 % Plot of density space-time profile
%                 subplot(1,nPlot,4)
%                 [phiSpaceTime, ySpaceTime] = getDensityProfile(obj);
%                 s=surf(repmat(time'*obj.shearRate, [1,size(phiSpaceTime,2)]), ySpaceTime, phiSpaceTime);
%
%                 s.EdgeColor = 'none';
%                 view([0 90])
%                 title('Density in Space-Time')
%                 xlab = xlabel('$$t \times (\dot{\gamma})$$');
%                 set(xlab, 'Interpreter', 'Latex')
%                 ylabel('Height')
%                 caxis([.55 .8])
%                 axis tight
%                 colorbar
            end

            % Save
            if (nargin > 2) && (saveIt == 1)
                oldFolder = cd(obj.dirPath);
                saveas(gca,'Convergence.tif')
                cd(oldFolder)
            end

            ttl = sprintf('%s', obj.dirPath);
            suptitle(ttl)

        end

        function obj = setInhom(obj)
            % setInhom(obj) set the time averaged, profile inhom.
            % sqrt(2-Norm of velocity fluctuations) computed by
            % Sqrt(Sum over slices (v_x/v_{wall} - y/H)^2)
            % See Flow of wet granular materials by Khamseh, 2015, but
            % using integration limits -H to H, not -H/2 to H/2 
            % Better than the PostP version, which doesn't time average the each height
            
            u = obj.vXSlice - obj.ySlice*obj.shearRate;
            x = obj.ySlice;
            intDu = trapz(x,u.^2);
            obj.inhomX = 3/2 * intDu / obj.L(2)^3 / obj.shearRate^2;

        end

        function obj = getInhomFft(obj)

            if any(isnan(obj.ySlice)) || any(isnan(obj.vXSlice))
                obj = getShearProfile(obj, -1);
            end

            % All non-dimensional
            L = numel(obj.ySlice);      % Number of samples
            T = mean(diff(obj.ySlice)); % Sample Period
            Fs = 1/T;                   % Sample frequency

            if mod(L,2)==1; L = L-1; end

            % Signal: Non-dimensional velocity fluctuations
%             t = obj.ySlice(1:L);
            f = Fs*(0:L/2)/L;
            X = obj.vXSlice(1:L) - obj.ySlice(1:L);

            Y = fft(X);

            P2 = abs(Y/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);

            if max(P1(f>1.1/2)) > max(P1(f<=1.1/2))
                obj.inhomFft = 1;
            else
                obj.inhomFft = 0;
            end

%             if varObj.inhomX > .1
%                 color = [0 1 0];
%             else
%                 color = [0 0 0];
%             end
%
%             subplot(2,1,1)
%             plot(t, X, 'Color', color, 'Marker', '-', 'LineWidth', 1)
%             hold on
%
%             subplot(2,1,2)
%             plot(2*f, P1, 'Color', color, 'Marker', '-', 'LineWidth', 1)
        end

        %% Set Profiles
        function [phiSpaceTime, ySpaceTime] = getDensityProfile(obj)

            for i = obj.tEnd-1:-1:1

                % Read Profile: y, phi
                fpProfile = strcat(obj.dirPath, '/profile/profile_', string(i));
                formatspec = '%f %*f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
                profile = Profile();
                profile = profile.readProfile(fpProfile, formatspec);

                % Read cell: L
                fpCell = strcat(obj.dirPath, '/cell/cell_', string(i));
                formatspec = '%*f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
                cell = Cell();
                cell = cell.readCell(fpCell, formatspec);

                profile = PostP.validateProfile(profile, cell.L(2));

                phiSpaceTime(i,1:length(profile.phi)) = profile.phi;
                ySpaceTime(i,1:length(profile.phi)) = profile.y;
            end
        end

        function obj = getShearProfile(obj, filenums)
            % If filenum = -1 or not given, get the average shear Profile
            % over all timesteps

%             % Don't recalculate if already saved
%             if ~(any(isnan(obj.vXSlice)) || any(isnan(obj.ySlice)))
%                 return;
%             end

            % Which times to calculate
            if nargin == 1 || numel(filenums)==1 && filenums<0
                
                obj.tEnd=500;
%                 obj.nPass=200;
%                 warning('hardcode')                

                % Cannot use end time because DEM post process doesn't compute at last
                % timestep
                obj.times = obj.tEnd-1:-1:obj.nPass;
            else
                obj.times = filenums;
            end

            k=1;
            for i = obj.times
                % Read parameters: P and shear rate
                fpPara = strcat(obj.dirPath, '/para/parameter_', string(i));
                formatspec = '%*f %*f %*f %*f %f %f %*f';
                para = Para();
                para = para.readPara(fpPara, formatspec);

                % Read Profile: y, V
                fpProfile = strcat(obj.dirPath, '/profile/profile_', string(i));
                formatspec = '%f %*f %*f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
                profile = Profile();
                profile = profile.readProfile(fpProfile, formatspec);

                % Read cell: L
                fpCell = strcat(obj.dirPath, '/cell/cell_', string(i));
                formatspec = '%*f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
                cell = Cell();
                cell = cell.readCell(fpCell, formatspec);

                profile = PostP.validateProfile(profile, cell.L(2));

                % Get the horizontal speed profile and normalise by the
                % wall speed
                V_wall = 1; %para.shear_rate * cell.L(2) / 2;
                n = length(profile.V(:, 1));
                V(1:n, k) = profile.V(:, 1) / V_wall; %#ok<AGROW>
                V(n+1:end, k) = nan;


                % Get slice height and normalise by the cell height H
                H = 1; %cell.L(2)/2;
                n = length(profile.y);
                y(1:n, k) = profile.y / H;
                y(n+1:end, k) = nan;
                
                k=k+1;
            end

            % Average each slice in the shear profile over time
            V = mean(V, 2);
            y = mean(y, 2);

            % Trim the ends
            V(isnan(V)) = [];
            y(isnan(y)) = [];
            while (1)
                if y(end) < y(end-1)
                    y(end) = [];
                    V(end) = [];
                else
                    break
                end
            end

            obj.vXSlice = V;
            obj.ySlice = y;

        end

        function obj = getShearRateProfile(obj, smooth)
            % Get the shear rate profile using slices
            % If shear profile exists then it will use that, so be aware
            % what filenum has been used
            % Optional Boolean smooth to smooth data

             % Don't Recalculate
            if ~(any(isnan(obj.dVdy)))
                return;
            end

            % Requires shear Profile aka V_x(y) for all y
            if all(isnan(obj.vXSlice)) || all(isnan(obj.ySlice))
                fprintf('Shear Profile has nan. Must get shear Profile first\n')
                filenum = input('filenums (-1 for all times): ');
                obj = getShearProfile(obj, filenum);
            end

            % Compute derivative
            dVdy = gradient(obj.vXSlice, obj.ySlice); %#ok<PROPLC>

            if (nargin == 2) && (smooth == 1)
                % Smooth this profile
                obj.dVdy = smoothdata(dVdy, 'gaussian', 10); %#ok<PROPLC>
            else
                obj.dVdy = dVdy;%#ok<PROPLC>
            end

        end

        %% Plot Profiles
        function plotShearProfile(obj, color)
            % Plot the shear Profile
            V = obj.vXSlice;
            y = obj.ySlice;

            if nargin == 1
                plot(V, y, '.')
            else
                plot(V, y, 'Color', color, 'LineWidth', 1)
            end

            xlab = xlabel('$$\frac{V_x} {V_x^{Wall}}$$');
            ylab = ylabel('$$\frac{y}{H}$$');
            set(ylab, 'Rotation', 0)
            set(xlab, 'Interpreter','latex');
            set(ylab, 'Interpreter','latex');

            ax=gca;
            ax.XLim = [-1 1];
            ax.YLim = [-1 1];

        end

        function plotShearRateProfile(obj, formatspec)
            % PLOTSHEARRATEPROFILE Plot normalised shear rate profile
            % obj
            % formatspec
            if nargin==2
                plot(obj.dVdy, obj.ySlice, formatspec)
            else
                plot(obj.dVdy, obj.ySlice, '.-')
            end
            xlab = xlabel('$$\frac{dV}{dy} / \frac{V_{wall}}{H}$$');
            ylab = ylabel('$$y / H$$');
            set(xlab, 'Interpreter', 'Latex')
            set(ylab, 'Interpreter', 'Latex')
            set(ylab, 'Rotation', 0)
            ax=gca;
%             ax.XLim = [.5 1.5];
            ax.YLim = [-1 1];
        end

        function plotStressProfile(obj)

            for filenum = obj.tEnd-1:-1:1

                % Read parameters: P and shear rate
                fpPara = strcat(obj.dirPath, '/para/parameter_', string(filenum));
                formatspec = '%*f %*f %*f %*f %f %f %*f';
                para = Para();
                para = para.readPara(fpPara, formatspec);

                % Read Profile: y, stress
                fpProfile = strcat(obj.dirPath, '/profile/profile_', string(filenum));
                formatspec = '%f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f';
                profile = Profile();
                profile = profile.readProfile(fpProfile, formatspec);

                % Read cell: time, L
                fpCell = strcat(obj.dirPath, '/cell/cell_', string(filenum));
                formatspec = '%f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
                cell = Cell();
                cell = cell.readCell(fpCell, formatspec);

%                 profile = PostP.validateProfile(profile, cell.L(2));

                stress_yy(:,filenum) = profile.stress(:,5);
                y(:,filenum) = profile.y;

            end

            plot(mean(stress_yy, 2), mean(y, 2));
            ttl=sprintf('Time-Averaged Stress Profile: C = %.0f\nShear Rate = %.3f after %.1f shears', obj.cohesion, para.shear_rate, para.shear_rate*cell.time);
            title(ttl)
            xlab=xlabel('$$\sigma_{yy}$$');
            ylabel('Height')
            set(xlab, 'Interpreter', 'Latex')
            ax = gca;
            ax.XLim = [0 3];
            ax.YLim = [-cell.L(2)/2 cell.L(2)/2];

        end

        %% Set Shear parameters
%         function obj = getShearInhom(obj)
%             % GETSHEARINHOM A metric of shear profile distance from homogeneous shear
%             % sqrt(2-Norm of velocity fluctuations in x-direction) computed by
%             % Sqrt(Sum over slices (v_x/v_{wall} - y/H)^2)
%             % See Flow of wet granular materials: A numerical study, Khamseh, 2017
%
% %             if ~isnan(obj.shearInhom)
% %                 return
% %             end
%
%             if any(isnan(obj.vXSlice)) || any(isnan(obj.ySlice))
%                 fprintf('Must determine shear profile first\n')
%                 obj = obj.getShearProfile();
%             end
%
%             V = obj.vXSlice(:);
%             y = obj.ySlice;
%
%             diff = V - y;
%             obj.shearInhom = sqrt(sum((diff).^2));
%             assert(obj.shearInhom >= 0, 'Negative Inhomogeneous Shear Metric')
%         end
%
%         function obj = getShearRateInhom(obj)
%             % A measure of how far the shear rate profile is from being
%             % homogeneous where \dot{\gamma} = 1
%             diff = obj.dVdy - 1;
%             diff(isnan(diff)) = [];
%             obj.shearRateInhom = sqrt(sum( (diff).^2) );
%         end

        function obj = getShearBandWidth(obj)

            [gms, dy, ~] = getGms(obj);

            % Find the width of the shear bands shear greater than thresh
            width = [];
            thresh = 1.33;

            for t = obj.tEnd-1:-1:1
                len = 0;
                for i = 1:size(gms, 1)
                    if gms(i, t) > thresh
                        len = len + dy;
                    elseif len > 0
                        if i == size(gms,1)
                            % Connect shear bands across the periodic BC
                            width(1) = width(1) + len;
                        else
                            width(end+1) = len;
                            len = 0;
                        end
                    end
                end
            end

            obj.sbWidth = width;

        end

        function obj = getShearBandTime(obj)
            [gms, ~, dt] = getGms(obj);

            % Find the time persistance of the shear bands greater than
            % thresh
            sbTime = [];
            thresh = 1.33;

            for i = 1:size(gms, 1)
                len = 0;
                for t = 1:obj.tEnd-1
                    if gms(i, t) > thresh
                        len = len + dt;
                    elseif len > 0
                        sbTime(end+1) = len;
                        len = 0;
                    end
                end
            end

            obj.sbTime = sbTime;

        end

        function [gms, dy, dt] = getGms(obj)
            % Get the shear rate profile through time
            % Each column gms(:, t) is a shear rate profile for some time t
            % Also return the profile slice widths dy

            tLast = 0;

            % For each timestep
            for t = obj.tEnd-1:-1:1

                % Get the normalised shear rate profile
                fp = strcat(obj.dirPath, '/profile/profile_', string(t));
                profile = Profile();
                profile = profile.readProfile(fp, '%f %*f %*f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f');

                % Read cell: t, L
                fpCell = strcat(obj.dirPath, '/cell/cell_', string(t));
                formatspec = '%f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
                cell = Cell();
                cell = cell.readCell(fpCell, formatspec);

                profile = PostP.validateProfile(profile, cell.L(2));

                dt = tLast-cell.time;
                tLast = cell.time;

                % Read parameters: P and shear rate
                fpPara = strcat(obj.dirPath, '/para/parameter_', string(t));
                formatspec = '%*f %*f %*f %*f %f %f %*f';
                para = Para();
                para = para.readPara(fpPara, formatspec);

                V = profile.V(:,1);
                y = profile.y;
                n = length(V);

                % Normalise the velocity and height
                V_wall = para.shear_rate * cell.L(1) / 2;
                V(:, 1) = profile.V(:, 1) / V_wall;


                % Get slice height and normalise by the cell height
                % NB: H = L_x / 2 not L_y / 2
                H = cell.L(1)/2;
                y = y / H;

                % Compute derivative and return it as the shear rate
                % profile gms
                gms(:,t) = PostP.deriv(V, y);


%                % Smooth this profile
%                 dVdy = smoothdata(dVdy, 'gaussian', 30);

            end

            dy = mean(y(3:n)-y(1:n-2))/2;

            % Normalised time step
            dt = dt * para.shear_rate;
        end

        function [stress, dvdy]= getGrainStress(obj)
            % Get the stress profile at all timesteps
            for i = obj.tEnd-1:-1:1
                % Read Grain stress
                error('Redo this section')
                fpGrainPp = strcat(obj.dirPath, '/grain_post_process/grain_pp', string(i));
                formatspec = '%f %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f';
                grainPp = GrainPp();
                grainPp = grainPp.readGrainPp(fpGrainPp, formatspec);
                stress(:,i) = grainPp.stress(:, 2);
                dvdy(:,i) = grainPp.gradV(:, 2);
            end

        end

        function obj = getShearplugWidth(obj)

            [gms, dy, ~] = getGms(obj);

            % Find the width of the shear bands shear greater than thresh
            width = [];
            thresh = 2/3;

            for t = obj.tEnd-1:-1:1
                len = 0;
                for i = 1:size(gms, 1)
                    if gms(i, t) < thresh
                        len = len + dy;
                    elseif len > 0
                        width(end+1) = len;
                        len = 0;
                    end
                end
            end

            obj.plugWidth = width;

        end

        function [obj, times, ys, gms, H] = getShearBands(obj)
            % For a single simulation


            % For each timestep
            for t = obj.tEnd-1:-1:1

                % Get the normalised shear rate profile
                fp = strcat(obj.dirPath, '/profile/profile_', string(t));
                profile = Profile();
                profile = profile.readProfile(fp, '%f %*f %*f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f');

                % Read cell: t, L
                fpCell = strcat(obj.dirPath, '/cell/cell_', string(t));
                formatspec = '%f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
                cell = Cell();
                cell = cell.readCell(fpCell, formatspec);

                % Read parameters: P and shear rate
                fpPara = strcat(obj.dirPath, '/para/parameter_', string(t));
                formatspec = '%*f %*f %*f %*f %f %f %*f';
                para = Para();
                para = para.readPara(fpPara, formatspec);

                profile = PostP.validateProfile(profile, cell.L(2));

                V = profile.V(:,1);
                y =  profile.y;

                % Get non-dimensional time
                times(t) = cell.time * para.shear_rate;

                % Normalise the velocity and height
%                 V_wall = para.shear_rate * cell.L(1) / 2;
%                 V(:, 1) = profile.V(:, 1) / V_wall;

                % Get slice height and normalise by the cell height
                H = cell.L(2)/2;
%                 y = y / H;
                ys(1:length(y),t) = y;
                ys(length(y)+1:end,t) = nan;

                % Compute derivative and return it as the shear rate
                % profile gms
                dVdy= gradient(V, ys(1:length(profile.y),t)); %#ok<*PROP>
                gms(1:length(dVdy),t) = dVdy;
                gms(length(dVdy)+1:end,t) = nan;

        %        % Smooth this profile
        %         dVdy = smoothdata(dVdy, 'gaussian', 30);

            end

        end


        %% Plot Shear Parameters
        function plotShearBandWidth(obj)
            histogram(obj.sbWidth, 'Normalization', 'pdf')
            xlab = xlabel('Shear Band Width $$\Delta y / (H)$$');
            set(xlab, 'Interpreter', 'Latex')
            ylabel('pdf')
            ttl = sprintf('Shear Zones Widths. C=%i H=%i', obj.cohesion, obj.L(1));
            title(ttl)
        end

        function plotPlugWidth(obj)
            hist(obj.plugWidth)
            xlab = xlabel('Plug Width $$\Delta y / (H)$$');
            set(xlab, 'Interpreter', 'Latex')
            ylabel('Count')
            ttl = sprintf('Plug Zone Widths. C=%i H=%i', obj.cohesion, obj.L(1));
            title(ttl)
        end

        function plotShearBandTime(obj)
            eps = 1;
            hist(obj.sbTime(obj.sbTime>eps))
            xlab = xlabel('Shear Band Time Persitence $$ t * (\dot{\gamma})$$');
            set(xlab, 'Interpreter', 'Latex')
            ylabel('Count')
            ttl = sprintf('Shear Band Persistance. C=%i H=%i', obj.cohesion, obj.L(1));
            title(ttl)
        end

        function [stress, dvdy] = plotShearVsShearRate(obj)
            % Each point has unique height and time
            [stress, dvdy] = getGrainStress(obj);

            plot(dvdy, stress, '.k')
            xlab = xlabel('$$\dot{\gamma}$$');
            ylab = ylabel('$$\tau$$');
            set(xlab, 'Interpreter', 'Latex')
            set(ylab, 'Interpreter', 'Latex')

        end

        function plotShearBandsContour(obj, times, ys, gms, H)
            if nargin==1
                [obj, times, ys, gms, H] = getShearBands(obj);
            end

            % Plot the filled contours
            s=surf(times,ys,gms);
            s.EdgeColor = 'none';

            colorbar
%             caxis([-.5 2.5])
            gmMax = quantile(gms(:), .9);
            gmMin = quantile(gms(:), .1);
            caxis([gmMin gmMax])
            axis tight
            view(0,90)

            colormap parula

            xlab = xlabel('$$t \dot{\gamma}$$');
            ylabel('y');
            set(xlab, 'Interpreter', 'Latex')
            title('$\dot{\gamma}(t,y)$', 'Interpreter', 'Latex')

        end

    end

    methods (Static)
%         function obj = loadobj(S)
%             if isstruct(S)
%                 obj = ProfileP(S.dp);
%             else
%                 obj = S;
%             end
%             obj = reload(obj, S);
%         end

        function inhomFFT = returnInhomFFT(dp)
            fpProf = strcat(dp, '/ProfileP.mat');
            try
                in = load(fpProf);
                profileP = in.profileP;
            catch
                profileP = ProfileP(dp);
            end

            if isnan(profileP.inhomFft)
                profileP = profileP.refresh();
                profileP = profileP.getInhomFft();
                save(fpProf, 'profileP')
            end

            assert( ~isnan(profileP.inhomFft), 'profileP.inhomFft is NaN\n')
            inhomFFT = profileP.inhomFft;
        end

    end

end
