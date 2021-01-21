classdef PostP
    % POSTP Post process the data of 1 simulation
    %
    % Use as a superclass with only the core information such as
    % directory path, number of files, number of grains,
    % shear rate and level of cohesion.
    %
    % Also contains visualisation tools such as visu, makeMovie and quiverOne

    properties

        dirPath = '~/documents/DEM/DATA/run' % absolute path to simulation directory
        tEnd = NaN % number of files in ./grain folder

        nG = NaN % Number of grains
        nContacts = NaN % Average number of contacts
        dnC = NaN % Standard Error of nContacts


        cohesion = NaN % Cohesion intensity
        shearRate = NaN % Shear Rate
        L = [NaN NaN] % Cell Dimensions (Lx2H) (time averaged)
        dL = NaN % Standard Error of L(2)

        solidFrac = NaN % Solids fraction V_{grains}/V_cell (time averaged)

        dvX = NaN % Velocity fluctuation standard deviation (time-space averaged) in x-direction
        dvY = NaN % Velocity fluctuation standard deviation (time-space averaged) in y-direction
        dvX_sd = NaN % standard deviation of above mean across times
        dvY_sd = NaN % standard deviation of above mean across times

        da = NaN % Velocity fluctuation standard deviation (time-space averaged) in x-direction
        da_sd = NaN % standard deviation of da

        df = NaN % Force fluctuations y-direction (on grains)
        df_sd = NaN

        df_c = [NaN NaN] % Force fluctuations (over contacts)
        df_c_sd = [NaN NaN]

        inhomX = NaN % A measure of the system homogeneity (0=homo., 1=inhom.)
        inhomY = NaN % A measure of the system homogeneity (0=homo., 1=inhom.)

        nPass = 1; % Number of files to skip for all measurements

    end

    methods

        function obj = PostP(dp)
            % POSTP Constructor
            % dp - path to experiment folder

%             warning('nPass=%.0f', obj.nPass)
            if nargin==1

                % If 'PostP.mat' file exists, read from that
                contents = struct2cell(dir(dp));
                if any(strcmp('PostP.mat', contents(1,:)))

                    S = PostP.readPostP(dp);
                    obj = refreshStruct(obj, S);

                else

                    obj.dirPath = dp;
                    obj = nTimesteps(obj);
                    obj = getCohesion(obj);
                    obj = getDimensions(obj, -1);
                    obj = getNG(obj);

                end
            end

        end

        function obj = nTimesteps(obj)

            obj.tEnd = PostP.nFiles(obj.dirPath);

        end

        function obj = setDirPath(obj, dirPath)
            obj.dirPath = dirPath;
        end

        function obj = getCohesion(obj)
            % Get the cohesion and shear rate
            para1Path = strcat(obj.dirPath, '/para/parameter_1');
            para = Para();
            formatspec = '%*f %*f %*f %*f %*f %f %f';
            para = para.readPara(para1Path, formatspec);
            obj.shearRate = para.shear_rate;
            obj.cohesion = para.cohesion;

            assert(obj.cohesion >= 0, 'Negative cohesion value read')
        end

        function obj = getNG(obj)
            % set number of grains
            obj.nG = PostP.nGrains1(obj.dirPath);
        end

        function obj = getDimensions(obj, filenums)
            % GETDIMENSIONS Get the time averaged cell dimensions
            % filenums = array of file numbers to average over, set to -1
            % for all timesteps

            if (double(filenums) == -1)
                filenums = obj.tEnd:-1:obj.nPass;
            elseif ~(size(filenums,1)==1 || size(filenums,2)==1)
                error('filenums must be a vector or -1')
            end
            if isempty(filenums)
                error('filenums between nPass and tEnd must be non-empty')
            end

            n = length(filenums);
            idx = 1;

            for i = n:-1:1
                % Get the cell dimensions
                cellPath = strcat(obj.dirPath, '/cell/cell_', string(filenums(i)));
                cell = Cell();
                formatspec = '%*f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
                cell = cell.readCell(cellPath, formatspec);
                L(idx, :) = cell.L; %#ok<PROPLC>
                idx = idx + 1;
            end
            obj.L = mean(L, 1); %#ok<PROPLC>
            obj.dL = std(L,0,1);
        end

        function shift = getShift(obj, filenum)
            % GETSHIFT Get the cell shift at filenum
            fp = strcat(obj.dirPath, '/cell/cell_', string(filenum));
            cell = Cell();
            formatspec = '%f %f %f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f';
            cell = cell.readCell(fp, formatspec);
            shift = cell.shift;
        end

%         function dvX = get.dvX(obj)
%             if isnan(obj.dvX)
%                 instantaneous = 1;
%                 obj = getDv(obj, instantaneous);
%             end
%             dvX = obj.dvX;
%         end
%
%         function dvY = get.dvY(obj)
%             if isnan(obj.dvY)
%                 instantaneous = 1;
%                 obj = getDv(obj, instantaneous);
%             end
%             dvY = obj.dvY;
%         end

        function obj = getDv(obj, instantaneous)
            % GETDV Return the system time and spatially averaged velocity
            % fluctuiations with units [dv] = L/T. Save all grain velocity
            % fluctuations to file.
            % <dv^2> = <(dv) - <dv>)^2>
            %
            % Instantaneous boolean whether or not to use instantaneous grain
            % velocities (true) or compute grain velocities from
            % differneces in displacement (false)

            [dU, dV] = getDvAll(obj, instantaneous);

            % Take the standard deviation (time averaged) of the velocity (aka
            % the velocity fluctuations)
            % NB: alt way to calculate: mean(sqrt(mean(dU.^2, 2))), mean(sqrt(var(dU, 1, 2)))
%             obj.dvX = mean(sqrt(mean(dU.^2, 2)));
%             obj.dvY = mean(sqrt(mean(dV.^2, 2)));
            obj.dvX = mean(mean(sqrt(dU.^2)));
            obj.dvY = mean(mean(sqrt(dV.^2)));

            obj.dvX_sd = std(sqrt(mean(dU.^2, 2)));
            obj.dvY_sd = std(sqrt(mean(dV.^2, 2)));

        end

        function [dU, dV] = getDvAll(obj, instantaneous)
            % GETDVALL(obj, instantaneous)
            % Get velocity fluctuations for every grain and every timestep
            % Matrix indexed by dv(filenum, grainID)
            %
            % Instantaneous boolean whether or not to use instantaneous grain
            % velocities (true) or compute grain velocities from
            % differneces in displacement (false).

            if instantaneous
                fp = strcat(obj.dirPath, '/dV_Raw');
%                 fprintf('Instantaneous Velocity fluctuations\n')
            elseif ~instantaneous
                fp = strcat(obj.dirPath, '/dV');
                fprintf('Not Instantaneous Velocity fluctuations\n')
            else
                error('Invalid boolean input instantaneous')
            end

            % try loading from matlab file
            try
                in = load(fp);
                dU = in.dU;
                dV = in.dV;
            catch
                nTimesteps = obj.tEnd-obj.nPass;
                fprintf('Reading gain velocities over %.0f times ...', nTimesteps)
                fprintf('%6.0f/%6.0f', 0, nTimesteps)

                formatspec = '%d %f %f %f %f %*f %*f %*f %*f %*f %*f';
                dpGg = strcat(obj.dirPath, '/grain/grain_');

                U = zeros(nTimesteps, obj.nG);
                V = zeros(nTimesteps, obj.nG);
                X = zeros(nTimesteps, obj.nG);
                Y = zeros(nTimesteps, obj.nG);

                for filenum = 1+obj.nPass:obj.tEnd

                    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%6.0f/%6.0f', filenum-obj.nPass, obj.tEnd-obj.nPass)
                    % Read grain position and velocity
                    fpG = strcat(dpGg, string(filenum));
                    grain = Grain();
                    grain = grain.readGrain(fpG, obj.nG, formatspec);

                    U(filenum-obj.nPass, :) = grain.V(:,1);
                    V(filenum-obj.nPass, :) = grain.V(:,2);

                    X(filenum-obj.nPass, :) = grain.X(:,1);
                    Y(filenum-obj.nPass, :) = grain.X(:,2);

                end

                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bdone\n')

                if ~instantaneous
                    % use grain positions to generate velocity

                    dt = PostP.getDt1(obj.dirPath);
                    U = (X(2:end, :) - X(1:end-1, :)) / dt;
                    V = (Y(2:end, :) - Y(1:end-1, :)) / dt;
                    Y = (Y(2:end, :) + Y(1:end-1, :)) / 2;
                end

                fprintf('Correct for drift: ')
                % Correct velocity for homogeneous flow and total drift
                dU = U - obj.shearRate*Y;
                dU = dU - mean(mean(dU));
                fprintf('1/2')
                dV = V - 0;
                dV = dV - mean(mean(V));
                fprintf('\b\b\b done\n')

                % Save to file for reuse speed
                fprintf('Saving velocity fluctuations to file: ')
                save(fp, 'dU', 'dV')
                fprintf('done\n')
            end



        end


        function obj = getInhom(obj)
            % GETINHOM A metric for how inhomogeneous the system is
            % sqrt(2-Norm of velocity fluctuations) computed by
            % Sqrt(Sum over slices (v_x/v_{wall} - y/H)^2)
            % See Flow of wet granular materials: A numerical study, Khamseh, 2017

            error('Depreciated. Use Library/misc/setInhom.m')
            warning('Use ProfileP.setInHom() for a better answer. It has time averaging over each slice, which allows for turbulent flows.')

            if (isnan(obj.dvX))
                fprintf('Getting velocity fluctuations.\n')
                obj = obj.getDv(1);
            end

            obj.inhomX = 12 * obj.dvX^2 / obj.L(2)^2 / obj.shearRate^2;
            obj.inhomY = 12 * obj.dvY^2 / obj.L(2)^2 / obj.shearRate^2;

            if ~((obj.inhomX <= 1) && (obj.inhomX >=0))
                warning('Invalid shear inhomogenious metric in x-dir: %f', obj.inhomX);
            end
            if ~((obj.inhomY <= 1) && (obj.inhomY >=0))
                warning('Invalid shear inhomogenious metric in y-dir: %f', obj.inhomY);
            end
        end

%         function obj = refresh(obj, fp)
%           nb function name clash with rheology
%             % REFRESH loads in PostP file data
%             S = load(fp);
%
%             try
%                 S = S.postP;
%             catch
%                 S = S.varObj;
%             end
%
%             obj.tEnd = S.tEnd;
%             obj.nG = S.nG;
%             obj.cohesion = S.cohesion;
%             obj.shearRate = S.shearRate;
%             obj.L = S.L;
%             obj.dvX = S.dvX;
%             obj.dvY = S.dvY;
%             obj.inhomX = S.inhomX;
%             obj.inhomY = S.inhomY;
%             obj.nPass = S.nPass;
%         end

        function obj = refreshStruct(obj, S)
            % REFRESHSTRUCT Refresh object properties from structure S
            obj.dirPath = S.dirPath;
            obj.tEnd = S.tEnd;
            obj.nG = S.nG;
            obj.cohesion = S.cohesion;
            obj.shearRate = S.shearRate;
            obj.L = S.L;
            obj.dvX = S.dvX;
            obj.dvY = S.dvY;
            obj.inhomX = S.inhomX;
            obj.inhomY = S.inhomY;
            obj.nPass = S.nPass;
        end

        function dsolidFrac = getSolidFracSTD(obj)
            % returns the standard error of solid fraction

           formatspec = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %f';
           grain = Grain();
           fp = strcat(obj.dirPath, '/grain/grain_1');
           grain = readGrain(grain, fp, obj.nG, formatspec);

           % Use total derivative
           areaSum = sum(pi*grain.R.^2);
           dsolidFrac = obj.dL(2) * areaSum / obj.L(1)/obj.L(2)^2;


        end

        function obj = getSolidFrac(obj)
           %  GETVOIDFRAC Returns the average solid fraction

           if isnan(obj.nG)
               obj.nG = getNG(obj);
           end

           if any(isnan(obj.L))
               filenums = -1;
               obj = getDimensions(obj, filenums);
           end

           formatspec = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %f';
           grain = Grain();
           fp = strcat(obj.dirPath, '/grain/grain_1');
           grain = readGrain(grain, fp, obj.nG, formatspec);

           areaSum = sum(pi*grain.R.^2);
           Volume = obj.L(1)*obj.L(2);

           obj.solidFrac = areaSum / Volume;

        end

        function obj = getDf_c(obj)
            T1=obj.nPass+1;
            T2=obj.tEnd;
            T = T2-T1+1;
            contact=Contact();

            fy = nan(1,T); fx=fy;

             for t = 1:T
                fn = T1+t-1;
                fp = strcat(obj.dirPath, '/contact/contact_', string(fn));
                contact = readContact(contact, fp);

                Fy = contact.Force(:,2);
                fy(t) = std(Fy,0,'all','omitnan');
                
                Fx = contact.Force(:,1);
                fx(t) = std(Fx,0,'all','omitnan');
             end

             obj.df_c = [mean(fx,2) mean(fy,2)];
             obj.df_c_sd = [std(fx,0,2) std(fy,0,2)];
        end

        function obj = getDf(obj)
            % Grain averaged(?) force fluctuations. Check it

            T1=obj.nPass+1;
            T2=obj.tEnd;
            T2=T1+10
            T = T2-T1+1;

            formatspec = '%d %*f %*f %*f %*f %*f %*f %*f %f %*f %*f';
            contact = Contact();
            grain = Grain();

            fy = nan(1,T);

             for t = 1:T
                fn = T1+t-1;
                fp = strcat(obj.dirPath, '/contact/contact_', string(fn));
                contact = readContact(contact, fp);

                Fy = full(sparse(double(contact.ID_A)+1,double(contact.ID_B)+1,contact.Force(:,2), obj.nG, obj.nG));
                Fy=Fy-Fy';
                

                fp = strcat(obj.dirPath, '/grain/grain_', string(fn));
                grain = readGrain(grain, fp, obj.nG, formatspec);
                fy(t) = mean( sum(Fy.^2,2)).^.5;
             end

             obj.df = mean(fy,2);
             obj.df_sd = std(fy,0,2);
        end

        function [Fx,Fy] = getMeanForce(obj)

            T1=obj.nPass+1;
            T2=obj.tEnd;
            T2=obj.nPass+3;
            T = T2-T1+1;

            contact = Contact();

            fy = nan(1,T);
            fx = nan(1,T);

             for t = 1:T
                fn = T1+t-1;
                fp = strcat(obj.dirPath, '/contact/contact_', string(fn));
                contact = readContact(contact, fp);

                Fy = full(sparse(double(contact.ID_A)+1,double(contact.ID_B)+1,contact.Force(:,2), obj.nG, obj.nG));
                Fy=Fy-Fy';
                Fx = full(sparse(double(contact.ID_A)+1,double(contact.ID_B)+1,contact.Force(:,1), obj.nG, obj.nG));
                Fx=Fx-Fx';

                fy(t) = mean( sum(Fy,2));
                fx(t) = mean( sum(Fx,2));
             end

             Fy = mean(fy,2);
             Fx = mean(fx,2);

        end

        function obj = getDa(obj)
            % Get the acceleration fluctuations parameter da in x-direction
            T1=obj.nPass+1;
            T2=obj.tEnd;
            N = obj.nG;

            a1=[];
            a = nan(T2-T1+1,1);
            formatspec = '%d %*f %*f %*f %*f %*f %*f %*f %f %*f %*f';
            contact = Contact();
            grain = Grain();

            i=1;
             for t = T1:T2
                fp = strcat(obj.dirPath, '/contact/contact_', string(t));
                contact = readContact(contact, fp);
                Fy = sparse(double(contact.ID_A)+1,double(contact.ID_B)+1,contact.Force(:,2), N, N);
                Fy=Fy-Fy';
                Fx = sparse(double(contact.ID_A)+1,double(contact.ID_B)+1,contact.Force(:,1), N, N);
                Fx=Fx-Fx';

                fp = strcat(obj.dirPath, '/grain/grain_', string(t));
                grain = readGrain(grain, fp, N, formatspec);
                m = grain.mass;

                % (sum(Fx,2)./m) has size nGrain x 1, so it's std over the grains
                at = (sum(Fx,2).^2 + sum(Fy,2).^2).^.5./m;
                at(abs(at)<eps())=[];
                a(i) = std(at);
                i=i+1;
                
                a1=[a1; at];
                
             end
%              obj.da = mean(a);
%              obj.da_sd = std(a);

            fprintf('\nAcceleration computations: %f = %f\n', full(mean(a)), full(std(a1)))
             
            obj.da=std(a1);
            warning('PostP.da_sd is NaN.')
            
        end

        function ts = readTimeSteps(obj, isStress, isXTrue)
            % Read in parameters at all timesteps (except last) into cell array of timestep objects
            %
            % If every grain at every timestep is read into the working memory (as in this function),
            % it can be slow and could cause memory issues. (Worked for 2440 timesteps on uni computer).
            % Do you need all this information at once?
            %
            % optional input isStress to read grain stresses only

            T1 = obj.nPass+1;
            T2 = obj.tEnd-1;
            N = T2-T1+1;
            k = N;
            ts = TimeStep(obj.dirPath, N);

            if nargin == 2 && isStress

                formatspecGP = '%f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
                fprintf('Reading grain stresses in file number %2.0f + %6.0f/%6.0f', obj.nPass, 0, T2)
                for fn = T1:T2
%                     fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%6.0f/%6.0f', fn, T2)

                    % Read grainPp
                    fp = strcat(obj.dirPath, '/grain_post_process/grain_pp', string(fn));
                    ts(k).grainPp = ts(k).grainPp.readGrainPp(fp, obj.nG, formatspecGP);
                    k=k-1;
                end
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bdone\n')

                assert(~any(any(isnan(ts(1).grainPp.stress))))
                assert(~any(any(isnan(ts(end-1).grainPp.stress))))

            else

                formatspecG = '%d %f %f %f %f %*f %*f %*f %f %*f %*f'; % read grain ID X V
                formatspecC = '%*f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f'; % read cell


                fprintf('Reading grain positions and momentum in file number %2.0f + %6.0f/%6.0f', obj.nPass, 0, T2)
                for fn = T1:T2
%                     fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%6.0f/%6.0f', fn, T2)

                    % Read grains
                    fp = strcat(obj.dirPath, '/grain/grain_', string(fn));
                    ts(k).grain = ts(k).grain.readGrain(fp, obj.nG, formatspecG);

                    % Read cell (required for XTrue)
                    fp = strcat(obj.dirPath, '/cell/cell_', string(fn));
                    ts(k).cell = ts(k).cell.readCell(fp, formatspecC);

                    % Read cell
                    fp = strcat(obj.dirPath, '/cell/cell_', string(fn));
                    ts(k).cell = ts(k).cell.readCell(fp);

                    k=k-1;
                end
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bdone\n')


                if isXTrue
                    % Recopmute XTrue and VTrue
                    checkV = 1;
                    ts = recomputeXTrue(ts, checkV);

%                     Print to screen to see if it's working
%                     for i = 1:N
%                         fprintf('%10.5f, %10.5f; ', ts(i).grain.XTrue(2, :))
%                         fprintf('%10.5f, %10.5f; ', ts(i).grain.XTrue(3, :))
%                         fprintf('%10.5f, %10.5f; ', ts(i).grain.XTrue(4, :))
%                         fprintf('%10.5f, %10.5f\n', ts(i).grain.XTrue(5, :))
%                     end
                    assert(~any(any(isnan(ts(1).grain.XTrue))))
                    assert(~any(any(isnan(ts(end).grain.XTrue))))
                else
                    assert(~any(any(isnan(ts(1).grain.X))))
                    assert(~any(any(isnan(ts(end).grain.X))))
                end


                assert(~any(any(isnan(ts(1).grain.V))))
                assert(~any(any(isnan(ts(1).grain.mass))))
                assert(~any(any(isnan(ts(end).grain.V))))
                assert(~any(any(isnan(ts(end).grain.mass))))



            end



            % Parfor is slower!
            % parfor k = 1:N
            %     fn = k + obj.nPass;
            %     fp = strcat(obj.dirPath, '/grain/grain_', string(fn));
            %     ts(k).grain = ts(k).grain.readGrain(fp, obj.nG, formatspec);
            % end


        end

        function ts = readTimeStepsContacts(obj)
            % Read in parameters at all timesteps (except last) into cell array of timestep objects

            T1 = obj.nPass+1;
            T2 = obj.tEnd-1;
            N = T2-T1+1;
            k = N;
            ts = TimeStep(obj.dirPath, N);

            fprintf('Reading contact ID and forces in file number %2.0f + %6.0f/%6.0f', obj.nPass, 0, T2)
            for fn = T1:T2
%                 fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%6.0f/%6.0f', fn, T2)

                % Read contacts
                fp = strcat(obj.dirPath, '/contact/contact_', string(fn));
                ts(k).contact = ts(k).contact.readContact(fp);

                k=k-1;
            end
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bdone\n')

        end

        function obj = setNContacts(obj)
            [obj.nContacts, obj.dnC] = getNContacts(obj);
        end

        function [Nc, dNc] = getNContacts(obj)
            % Average number of contacts in the simulation

            T1 = obj.nPass+1;
            T2 = obj.tEnd;
%             warning('Sampling over 40 times')
%             T2 = T1+40;

            NNc = nan(T2-T1+1, 1);

            k=0;
            for fn = T1:T2
                k=k+1;
                fp = strcat(obj.dirPath, '/contact/contact_', string(fn));
                NNc(k) = PostP.nLines1(fp);
            end

            Nc = mean(NNc);
            dNc = std(NNc);

        end

        function Z = getCoordingation(obj)
            % Z = getCoordingation(obj) get mean coordination number
            % Z = 2 Nc / Ng

            if isnan(obj.nContacts)
                obj.nContacts = getNContacts(obj);
            end

            if isnan(obj.nG)
                obj.nG = getNG(obj);
            end

            Z = 2 * obj.nContacts / obj.nG;
        end

        tc = getTc(obj)

        % Visu
        function trackID = getTrackID(obj, tStart)
            % Get ID of particles in the middle half at time fileNum=tStart

            % Read in grains at tStart: ID, position
            fpGrain = strcat(obj.dirPath, '/grain/grain_', string(tStart));
            formatspec = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %*f';
            grain = Grain();
            grain = grain.readGrain(fpGrain, obj.nG, formatspec);

            % Read in the cell at tStart: dimensions
            fpCell = strcat(obj.dirPath, '/cell/cell_', string(tStart));
            cell = Cell();
            formatspec = '%*f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
            cell = cell.readCell(fpCell, formatspec);

            % If the grains are in the middle half, append grain.ID to trackID
            track = abs(grain.X(:,2)) < cell.L(2)/4;
            trackID=grain.ID(track);

        end

        function visu(obj, fileNum, contacts, trackID)
            % VISU(obj, fileNum, contacts, trackID)
            % bool contacts - whether or not to draw contacts
            % trackID - an array with the indexes of grains to colour
            % Uses speed of patch to plot grains, contacts are still slow

            assert(~isnan(obj.nG))
            ax=gca;
            PostP.visuOne(ax, obj.dirPath, fileNum, trackID, contacts, obj.nG)

        end

        function movieStress(obj)
           % Make movie of shear stress in grains
           f=figure(1);
           f.Visible = 'off';
           ax=gca;

           nFrames = obj.tEnd - obj.nPass + 1;
           nGrains=obj.nG;

           fp = strcat(obj.dirPath, '/stressVid');
           vid = VideoWriter(fp);
           vid.FrameRate = 10;
           open(vid)

            for  fileNum = 1:nFrames
                cla(ax)
                PostP.visuOneStress(ax, obj.dirPath, fileNum, nGrains)
                frame = getframe(f);
                writeVideo(vid, frame);
                fprintf('%d%% complete\n', round(fileNum/nFrames*100))
            end

            close(f)
            close(vid)
        end

        function makeMovie(obj)
            % MAKEMOVIE Visulise each timestep and save as avi file
            % Animate in background
%             f=figure(1);
            f=figure(1);
            f.Visible = 'off';
            ax=gca;

            postP = PostP(obj.dirPath);
            nFrames = postP.tEnd;

            % Videowriter object
            fp =  strcat(obj.dirPath, '/vid');
            vid = VideoWriter(fp);
            vid.FrameRate = 10;  % Matlab default is 30

            % Don't ask to plot contacts
            contacts = 0;
            trackID = [];

            open(vid)

            for fileNum = 1:nFrames
                cla(ax)
                postP.visuOne(ax, obj.dirPath, fileNum, trackID, contacts, obj.nG)

                frame = getframe(f);
                writeVideo(vid, frame);
                fprintf('%d%% complete\n', round(fileNum/nFrames*100))
            end

            close(f)
            close(vid)
        end

        function makeMovieTrack(obj)
            % MAKEMOVIE Visulise each timestep and save as avi file
            % Animate in background and color middle grains

            f=figure;
            set(f, 'WindowStyle', 'normal');
            set(f, 'OuterPosition', [0 0 1000 1100])
            set(f, 'InnerPosition', [0 0 1000 1000])
            set(f, 'Visible', 'off');
            set(f, 'Color', '[1 1 1]')

            ax=gca;
            set(ax, 'box', 'on')
            set(ax, 'XTick', [])
            set(ax, 'YTick', [])
            set(ax, 'Units', 'pixels')
            set(ax, 'Position', get(f, 'InnerPosition'));


%             postP = PostP(obj.dirPath); // This was a silly line of code
            nFrames = obj.tEnd-obj.nPass;
%             nFrames = 101;

            trackID = getTrackID(obj, obj.nPass+1);

            % Videowriter object
            fp =  strcat(obj.dirPath, '/vidTrackLFull');

%             vid = VideoWriter(fp, 'MPEG-4');
%             vid.Quality = 75;
            vid = VideoWriter(fp, 'Motion JPEG AVI');
            vid.Quality = 100;

            vid.FrameRate = 10;  % Matlab default is 30

            % Don't ask to plot contacts
            contacts = 0;

            nGrains = PostP.nGrains1(obj.dirPath);

            open(vid)

            for fileNum = [1:nFrames]+obj.nPass
                cla(ax)
                PostP.visuOne(ax, obj.dirPath, fileNum, trackID, contacts, nGrains)

                frame = getframe(f);
                writeVideo(vid, frame);
                fprintf('%d%% complete\n', round((fileNum-obj.nPass)/nFrames*100))
            end

            close(f)
            close(vid)
        end


    end

    methods (Static)

        function visuSingleG()
            ax=gca;
            dirPath =  '/Users/mattmacaulay/Documents/DEM/DATA/N10000/gm0.01/C40/run';
            trackID = 1;
            contacts = 0;
            nGrains = 1e4;

            for fileNum = 60:100
                cla
                PostP.visuOne(ax, dirPath, fileNum, trackID, contacts, nGrains)
                drawnow
            end
        end

        visuOne(ax, dirPath, fileNum, trackID, cPlot, nGrains)

        visuOne3D(ax, dirPath, fileNum, nGrains)

        visuOneFxLy(ax, dirPath, fileNum, nGrains)

        visuSweep()

        function visuOneStress(ax, dirPath, fileNum, nGrains)
            % VISUONE(dirPath, fileNum, trackID, contacts)
            % Plot grains showing shear stresses
            %
            % ax - axes to plot in
            % dirOne - directory of experiment
            % fileNum - array of fileNumbers to render
            %
            % Uses speed of patch function to plot grains

%             cla(ax)

            % Read grains
            fpGrain = strcat(dirPath, '/grain/grain_', string(fileNum));
            formatspec = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %f';
            grain = Grain();
            grain = grain.readGrain(fpGrain, nGrains, formatspec);

            % Read grains Post Processing (stress)
            fpGrainPp = strcat(dirPath, '/grain_post_process/grain_pp', string(fileNum));
            formatspec = '%f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
            grainPp = GrainPp();
            grainPp = grainPp.readGrainPp(fpGrainPp, nGrains, formatspec);

            % Read in Cell time, length and shift
            cell = Cell();
            fpCell = strcat(dirPath, '/cell/cell_', string(fileNum));
            cell = cell.readCell(fpCell, '%f %f %f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f');

            % Set limits
            set(ax, 'XLim', [0-1 cell.L(1)+1]);
            set(ax, 'YLim', [-cell.L(2)/2-1 +cell.L(2)/2+1]);
            % axis equal


%             % Read in the Pressure, shear Rate, cohesion
%             para = Para();
%             fpPara = strcat(dirPath, '/para/parameter_', string(fileNum));
%             para = para.readPara(fpPara, '%*f %*f %*f %*f %f %f %f');
%
%             inertial = 'I';
%             time = '$t \dot{\gamma}$';
%             ttl = sprintf('%s = %.3f, C = %.0f\n%s = %.1f', inertial, para.shear_rate, para.cohesion, time, cell.time*para.shear_rate-3);
%             ttl = title(ax, ttl, 'Interpreter','latex');
%             set(ttl, 'FontSize', 20)

            % Points of a unit circle
            nPoints = 30;
            circX = cos(2*pi*(1:nPoints)/nPoints)';
            circY = sin(2*pi*(1:nPoints)/nPoints)';

            % Get all circle coordinates
            X = grain.X(:,1)' + circX .* grain.R(:)';
            Y = grain.X(:,2)' + circY .* grain.R(:)';

            % Set all cicle colours
            colormap winter
%             cmap = colormap;

            q = grainPp.stress(:,3)'+grainPp.stress(:,2)';
            p = -grainPp.stress(:,4)'+grainPp.stress(:,1)';
            mu=q./p;
            mu(isnan(mu)) = 0;
            mu(mu>0) = mu(mu>0)/max(mu(mu>0));
            mu(mu<0) = mu(mu<0)/min(mu(mu<0));
%             mu = 64*(mu+1)/2;

            C = zeros(nPoints, length(grain.R), 3);
            C(:,:,1) = 1-repmat(mu,  nPoints, 1);
            C(:,:,2) = ones(nPoints, length(grain.R), 1);
            C(:,:,3) = ones(nPoints, length(grain.R), 1);

            % Plot grains
            patch(ax,X,Y,C, 'EdgeColor', 'none')

            % Plot the contacts
            contact = Contact();
            fpContact = strcat(dirPath, '/contact/contact_', string(fileNum));
            contact = contact.readContact(fpContact);
            contact.plotContact(grain, cell);

        end

        function quiverOne(ax, dirPath, fileNum)
            % QUIVERONE(dirPath, fileNum, trackID, contacts)
            % Quiver plot of fluctuating grain velocities
            %
            % dirOne - directory of experiment
            % fileNum - array of fileNumbers to render
            % ax - axes to plot in

            % Read grains
            fpGrain = strcat(dirPath, '/grain/grain_', string(fileNum));
            formatspec = '%d %f %f %f %f %*f %*f %*f %*f %*f %*f';
            grain = Grain();
            nG = PostP.nGrains(fpGrain);
            grain = grain.readGrain(fpGrain, nG, formatspec);


            % Read in Cell time, length and shift
            cell = Cell();
            fpCell = strcat(dirPath, '/cell/cell_', string(fileNum));
            cell = cell.readCell(fpCell, '%f %f %f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f');

            % Read in the Pressure, shear Rate, cohesion
            para = Para();
            fpPara = strcat(dirPath, '/para/parameter_', string(fileNum));
            para = para.readPara(fpPara, '%*f %*f %*f %*f %f %f %f');

%             inertial = 'I';
%             time = '$t \dot{\gamma}$';
%             ttl = sprintf('%s = %.3f, C = %.0f\n%s = %.1f', inertial, para.shear_rate, para.cohesion, time, cell.time*para.shear_rate-3);
%             ttl = title(ax, ttl, 'Interpreter','latex');
%             set(ttl, 'FontSize', 20)


            U = grain.V(:,1);
            V = grain.V(:,2);

            X = grain.X(:,1);
            Y = grain.X(:,2);

            dU = U - para.shear_rate*Y - mean(mean(U - para.shear_rate*Y));
            dV = V - 0                 - mean(mean(V));

%             cla(ax)

%             try
%                 addpath('/Users/mattmacaulay/Documents/MATLAB/FileExchange/quiverc')
%             catch
%                 error('Need quiverc (from File Exchange) in PATH')
%             end
%             quiverc(X,Y,dU,dV); % slow but arrows coloured by size
            h = quiver(X,Y,dU,dV);
            set(h, 'LineWidth', .8)
            set(h, 'MarkerSize', 10)
            set(h, 'Color', 'k')

% %             set(ax, 'XLim', [0-1 cell.L(1)+1]);
%             set(ax, 'YLim', [-cell.L(2)/2-1 +cell.L(2)/2+1]);
%             set(ax, 'box', 'on')
%             set(ax, 'XTick', [])
%             set(ax, 'YTick', [])
%             set(ax, 'Units', 'pixels')
%
%             axis equal

        end

        function I = getInertial(dirPath, fileNum)
            % I = getInertial(fileNum)
            %  Assumes: m=d=1 so that the density of a grain = 4/pi

            % FilePath for parameters
            fpPara = strcat(dirPath, '/para/parameter_', string(fileNum));

            % Read in Parameters
            formatspec = '%*f %*f %*f %*f %f %f %*f';
            para = Para();
            para = para.readPara(fpPara, formatspec);

            rho = 4/pi;
            shearRate = para.shear_rate;
            P = para.P;
            d = 1;

            I = shearRate * d * sqrt( rho / P);
        end

        function dydx = deriv(y, x)
            % Derivative in periodic signal
            % for regular derivative use gradient

            if length(x) > 1
                dx = mean(diff(x));
            elseif length(x) == 1
                dx = x;
            else
                error('Could not compute derivative\n')
            end

            dydx = gradient(y, dx);

            % Periodic side cases of derivative - still using central difference
            dydx(1) = (y(2) + y(end)) / ((x(2) + x(end)));
            dydx(end) = (y(1) + y(end-1)) / (x(1) + x(end-1));

        end

        function profile = validateProfile(profile, H)

            % remove invalid slices - those with
            valid = abs(profile.y) < H/2 - 1;

            if ~any(isnan(profile.y));      profile.y = profile.y(valid,:); end
            if ~any(isnan(profile.y_H));    profile.y_H = profile.y_H(valid,:); end
            if ~any(isnan(profile.phi));    profile.phi = profile.phi(valid,:); end
            if ~any(isnan(profile.V));      profile.V = profile.V(valid,:); end
            if ~any(isnan(profile.gradV));  profile.gradV = profile.gradV(valid,:); end
            if ~any(isnan(profile.stress)); profile.stress = profile.stress(valid,:); end
            if ~any(isnan(profile.V_H));    profile.V_H = profile.V_H(valid); end
            if ~any(isnan(profile.shear_rate_normalised)); profile.shear_rate_normalised = profile.shear_rate_normalised(valid); end
            if ~any(isnan(profile.ome));    profile.ome = profile.ome(valid,:); end
            if ~any(isnan(profile.dV2));    profile.dV2 = profile.dV2(valid,:); end

        end

        function n = nFiles(dirPath)
            % how many saved timesteps are there
            grainPath = strcat(dirPath, '/grain');
            contents = dir(grainPath);
            contents = contents(~ismember({contents.name},{'.','..','.DS_Store'}));
            contents = struct2cell(contents);
            contents = cell2mat(contents(5,:));
            contents = contents==0;
            n = sum(contents);
        end

        function n = nGrains1(dirPath)
            % NGRAINS Return the number of grains found in the first grain file

            fp = strcat(dirPath, '/grain/grain_1');
            n = PostP.nGrains(fp);

        end

        function n = nGrains(fp)
            % NGRAINS Return the number of grains found in this file

%             formatspec = '%d %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
%             grain = Grain();
%             grain = grain.readGrain(fp, formatspec);
%
%             n = numel(grain.ID);

            fid = fopen(fp);

            n = 0;
            tline = fgetl(fid);
            while ischar(tline)
              tline = fgetl(fid);
              n = n+1;
            end
            fclose(fid);
        end

        function time = getTime(dirPath)
            % Return array of times in each cell file
            nFiles = PostP.nFiles(dirPath);
            for t = nFiles:-1:1
                fp = strcat(dirPath, '/cell/cell_', string(t));
                cell = Cell();
                cell = cell.readCell(fp, '%f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f');
                time(t) = cell.time;
            end
        end

        function dt = getDt1(dirPath)
            % Return timestep between first and second files
            % Only read 2 files but remove truncation error. Using PostP.getTime leads to truncation error

            for t = 2:-1:1
                fp = strcat(dirPath, '/cell/cell_', string(t));
                cell = Cell();
                cell = cell.readCell(fp, '%f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f');
                time(t) = cell.time;
            end

            dt = mean(diff(time));

            % Account for one particluar case of rounding error used in these simulations
            if dt == .167
               dt = 1/6;
            elseif mod(dt, .001) ~= 0
               warning('possible truncation error\n')
           end

        end

        function S = readPostP(dp)

            fp = strcat(dp, '/PostP.mat');
            in = load(fp);
            postP = in.postP;

            S = PostP();

            S.dirPath = postP.dirPath;
            S.tEnd = postP.tEnd;
            S.nG = postP.nG;
            S.cohesion = postP.cohesion;
            S.shearRate = postP.shearRate;
            S.L = postP.L;
            S.dvX = postP.dvX;
            S.dvY = postP.dvY;
            S.inhomX = postP.inhomX;
            S.inhomY = postP.inhomY;
            S.nPass = postP.nPass;
        end

        function makeFourMovie()

            fp = '/Users/mattmacaulay/Dropbox (Sydney Uni)/Papers/2020_JFM_Cohesive_Rheology/figure/movie/out/track';

            dp = {
                '/Users/mattmacaulay/Documents/DEM/DATA/N10000/gm0.005/C0/run',...
                '/Users/mattmacaulay/Documents/DEM/DATA/N10000/gm0.005/C40/run2',...
                '/Users/mattmacaulay/Documents/DEM/DATA/N10000/gm0.3/C0/run',...
                '/Users/mattmacaulay/Documents/DEM/DATA/N10000/gm0.3/C40/run4',...
                };

            % Use number of frames in first axes for all plots
%             tEnd = obj.nPass+101;
%             fileNums = obj.nPass+1:tEnd;
            fileNums = 200:500;
            nFrames = numel(fileNums);

            f=gcf;
            clf(f)
            set(f, 'WindowStyle', 'normal');
            set(f, 'OuterPosition', [0 0 1000 1000])
            set(f, 'InnerPosition', [0 0 1000 1000])
            set(f, 'Visible', 'off'); % Faster for getFrame
            set(f, 'Color', '[1 1 1]')

            % Get ID of blue grains in centre
            for i = 4:-1:1
                postP = PostP(dp{i});
                thisTrackID = getTrackID(postP, 201);
                nMiddle(i) = numel(thisTrackID);
                trackID(i,1:nMiddle(i)) = thisTrackID;

                ax{i}=subplot(2,2,i);
                set(ax{i}, 'box', 'on')
                set(ax{i}, 'XTick', [])
                set(ax{i}, 'YTick', [])
            end


            width = 0.5;
            height = 0.425;
            set(ax{1}, 'Position', [.0 .5  width height]);
            set(ax{2}, 'Position', [.5 .5  width height]);
            set(ax{3}, 'Position', [.0 0   width height]);
            set(ax{4}, 'Position', [.5 0   width height]);


            contacts = 0; % Ask to not plot contacts


            % Videowriter object
            vid = VideoWriter(fp, 'Motion JPEG AVI');
            set(vid, 'Quality', 75);
            set(vid, 'FrameRate', 10);  % Matlab default is 30
            fprintf('VideoWriter object created.\n')

            open(vid)

            for fileNum = fileNums
                tic
                for i = 1:4
                    PostP.visuOne(ax{i}, dp{i}, fileNum, trackID(i,1:nMiddle(i)), contacts)
                end

                frame = getframe(f);
                writeVideo(vid, frame);
                fprintf('%3.0f%% complete. Filenum: %d. Frame time: %5.3f\n', (fileNum-nPass)/nFrames*100, fileNum, toc)

            end

            close(f)
            close(vid)
        end

        function [dvXAll, dvYAll] = smoothDv(dvXAll, dvYAll)
            % Get velocity fluctuations for every grain and every timestep
            % Forward average: not instantaneous fluctuation but
            % the average fluctuation between two files
            %
            % this was a silly way to do it. Rather than V = [v(t+1)+v(t)]/2
            % use V = [x(t+1) - x(t)] / dt as in getDvAll

            dvXAll = ( dvXAll(2:end, :) + dvXAll(1:end-1, :) ) / 2;
            dvYAll = ( dvYAll(2:end, :) + dvYAll(1:end-1, :) ) / 2;

        end

        function shearRate = getShearRate(dp)
            % Static method to get the shear rate in the first parameter file of the provided dir
            para1Path = strcat(dp, '/para/parameter_1');
            para = Para();
            formatspec = '%*f %*f %*f %*f %*f %f %f';
            para = para.readPara(para1Path, formatspec);
            shearRate = para.shear_rate;
        end

        N = nLines1(fp)

        [dX,t] = getLyapunovAll(dp,tStart,tEnd)
        [dX,t] = getLyapunov(dp,tStart,tEnd)

        blue = tag(dp,nGrains,t0,t1)

    end

end
