classdef TimeStep
% ts = TimeStep() Timestep class
% ts = TimeStep(dirPath, _)
% ts = TimeStep(length) 
% 
% Use length to create array of timestep objects
    
    properties
        grain = Grain()
        grainPp = GrainPp()
        contact = Contact()
        cell = Cell()
        
        dirPath
        tStart
        tEnd
    end
    
    methods

        function obj = TimeStep(dirPath, N)
            % Constructor
            % obj = TimeStep()
            % obj = TimeStep(dirPath)
            % obj = TimeStep(dirPath, N)

            if  nargin == 1
                obj.dirPath = dirPath;
            end

            if nargin == 2

                obj(N) = TimeStep();
                for i = 2:N-1
                    obj(i) = TimeStep();
                end
                obj(1) = TimeStep(dirPath);

            end

        end
        
        function ts = read(obj, tStart, tEnd,xOnly)
            
            % if only need timestep for position data, don't read in
            % velocity or mass data
            if xOnly
                fsGr = '%d %f %f %*f %*f %f %f %*f %*f %*f %*f'; 
            else
                fsGr = '%d %f %f %f %f %*f %*f %*f %f %*f %*f';
            end
            
            % Read in cell and grains at each timestep
            k = tEnd-tStart+1;
            nG = PostP.nGrains1(obj(1).dirPath);
            fsCe = '%f %f %f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f';

            fprintf('\tReading data. Progress: %4.0f/%4.0f', tEnd-k-tStart+1, tEnd-tStart+1)
            for filenum=tEnd:-1:tStart
                fprintf('\b\b\b\b\b\b\b\b\b%4.0f/%4.0f', tEnd-k-tStart+1, tEnd-tStart+1)
                ts(k) = TimeStep(obj(1).dirPath);

                % reads time, L shift
                fp = strcat(obj(1).dirPath, '/cell/cell_', string(filenum));
                ts(k).cell = ts(k).cell.readCell(fp, fsCe);

                % Read grain ID, X, XTrue %X,V,mass
                fp = strcat(obj.dirPath, '/grain/grain_', string(filenum));
                ts(k).grain = ts(k).grain.readGrain(fp, nG, fsGr);

                % Read contacts
                fp = strcat(obj(1).dirPath, '/contact/contact_', string(filenum)); 
                ts(k).contact = ts(k).contact.readContact(fp);

                k=k-1;
            end
            fprintf('\b\b\b\b\b\b\b\b\bdone\n')

            fprintf('Computing XTrue (and VTrue if it exists):\n')
            if xOnly; checkV=0; else; checkV = 1; end
            ts = ts.recomputeXTrue(checkV);
            fprintf('\tDone\n') 
            
        end
        
        function obj = recomputeXTrue(obj, checkV)
        % obj = recomputeXTrue(obj)
        % Recompute XTrue with by checking for particle jumps at each each timestep.
        % If obj.grains have velocity read the recompute V:=VTrue                        
            if checkV
                fprintf('\tRecomputing XTrue and VTrue with a non-periodic metric\n')
            else
                fprintf('\tRecomputing XTrue (and not V) with a non-periodic metric\n')
            end

            nT = numel(obj);

            for i = 1:nT   
                assert(~any(any(isnan(obj(i).grain.X))), 'Need to have all grain.X')
                assert(~any(any(isnan(obj(i).cell.L))), 'Need to have all cell.L')
                assert(~any(any(isnan(obj(i).cell.shift))), 'Need to have all cell.Shift')
                if checkV 
                    assert(~any(any(isnan(obj(i).grain.V))), 'Need to have all grain.V if recomputing V')
                end
            end            

            fp = strcat(obj(1).dirPath, '/grain/grain_1');
            nG = PostP.nGrains(fp);
            shearRate = PostP.getShearRate(obj(1).dirPath);
            
            % set boundary count
            XCount=zeros(nG,1);
            YCount=zeros(nG,1);
            
            for i = 1:nT
                
                if i==1
                    obj(i).grain.XTrue(1:nG,1) = obj(i).grain.X(:,1);
                    obj(i).grain.XTrue(1:nG,2) = obj(i).grain.X(:,2);
                    continue;
                end
                
                L = obj(i).cell.L;
                shift = obj(i).cell.shift;
                
                thisX = obj(i).grain.X(:,1);
                thisY = obj(i).grain.X(:,2);                
                Xlast = obj(i-1).grain.X(:,1);
                Ylast = obj(i-1).grain.X(:,2);
                
                right=(thisX-Xlast) < -L(1)/2;
                left=(thisX-Xlast) > L(1)/2;
                XCount=XCount+right-left;
                
                up=(thisY-Ylast) < -L(2)/2;
                down=(thisY-Ylast) > L(2)/2;
                YCount=YCount+up-down;
                
                obj(i).grain.XTrue(1:nG,1) = thisX+XCount*L(1)+YCount*shift;
                obj(i).grain.XTrue(1:nG,2) = thisY+YCount*L(2);
                                
                if checkV
                    obj(i).grain.V(:,1) = obj(i).grain.V(:,1) + YCount*L(2)*shearRate;
                end
                
            end
            

        end

        function obj = noDrift(obj)
            % OBJ = OBJ.NODRIFT() removes total vertical drift at each
            % timestep. y -> y - mean(y) for each t
            % grain.X(:,2) -> grain.X(:,2) - mean(grain.X(:,2))
            
            n = length(obj);
            if n == 0 || isnan(n) 
                error('No timesteps in Timestep object.')
            end
            
            for t=1:n
                y = obj(t).grain.X(:,2);
                obj(t).grain.X(:,2) = y - mean(y);
            end            
            
        end
    end
    
    
    
    
    methods (Static)
        
       function timeStep=loadTS(dirPath,tStart,tEnd,xOnly)
            
            if xOnly
                fp=strcat(dirPath, '/timeStep.mat');
            else
                fp=strcat(dirPath, '/timeStep2.mat');
            end
            
            
            try 
                error('Force to read timestep from scratch')
                in=load(fp);
                timeStep=in.timeStep;
                fprintf('Loading TimeStep with tStart: %d\t tEnd: %d\n',timeStep(1).tStart, timeStep(1).tEnd);
                assert(timeStep(1).tStart==tStart)
                assert(timeStep(1).tEnd==tEnd)
    %                 assert(any(any(~isnan(timeStep(1).XTrue))))
            catch e
                fprintf('Unable to load timeStep data file with requested tStart=%d, tEnd=%d.\n Error: %s\n', tStart,tEnd,e.message)
                N=tEnd-tStart+1;
                timeStep=TimeStep(dirPath, N);
                timeStep= read(timeStep, tStart, tEnd,xOnly);

                timeStep(1).tStart=tStart;
                timeStep(1).tEnd=tEnd;

                fprintf('Saving to timeStep.mat\n')
                save(fp,'timeStep');
                fprintf('\tDone.\n')
            end
        end 
    end
end

