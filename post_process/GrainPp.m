classdef GrainPp
    %GrainPp 
    % Grain Post process
    
    properties
        
        stress = [NaN NaN NaN NaN]
        gradV = [NaN NaN NaN NaN]
        V = [NaN NaN NaN]
        
        voro_area = NaN
        voro_area_n = NaN
        solid_area_n = NaN
        voro_neighbour_ID = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]
    end
    
    methods

        function obj = readGrainPp(obj, fp, nGrain, formatspec)
             % Read the stress
             % '%f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f'
             % Read the gradV
             % '%*f %*f %*f %*f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f'
             % Read the gradV, V
             % '%*f %*f %*f %*f %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f'
             % Read the voro_neighbour_ID
             % '%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %f %f %f'            
            
            if nargin == 3
                formatspec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
            end
            
            
            fileID = fopen(fp,'r');
            g = textscan(fileID, formatspec, nGrain, 'ReturnOnError', 0);
            fclose(fileID);
            

            if isempty(g)
                obj = Grain();
                return;
            end

            
            if nargin == 3           

                error('Do you really need all this data?\n If so actually write this code\n')

            elseif nargin==4 

                if strcmp(formatspec, '%f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    obj.stress = cell2mat(g(1:4));
                elseif strcmp(formatspec, '%*f %*f %*f %*f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    obj.gradV  = cell2mat(g(1:4));
                elseif strcmp(formatspec, '%*f %*f %*f %*f %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    obj.gradV  = cell2mat(g(1:4));
                    obj.V  = cell2mat(g(5:7));
                elseif strcmp(formatspec, '%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %f %f %f')
                    % NB: indexing in cpp starts at 0, in matlab from 1
                    obj.voro_neighbour_ID = cell2mat(g(1:12)) + 1;
                else
                    error('Invalid formatspec')
                end

            end

        end



         function obj = readGrainPpSlow(obj, fp, formatspec)
             % Read the stress
             % '%f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f'
             % Read the gradV
             % '%*f %*f %*f %*f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f'
             % Read the gradV, V
             % '%*f %*f %*f %*f %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f'
             % Read the voro_neighbour_ID
             % '%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %f %f %f'
             
            fileID = fopen(fp,'r');             
            
            if nargin == 3

                if strcmp(formatspec, '%f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    sizeG = [4 Inf];
                elseif strcmp(formatspec, '%*f %*f %*f %*f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                   sizeG = [4 Inf];
                elseif strcmp(formatspec, '%*f %*f %*f %*f %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                   sizeG = [7 Inf];
                elseif strcmp(formatspec, '%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %f %f %f')
                   sizeG = [12 Inf];
                else                    
                    error('Invalid formatspec to read grains')
                end
                g=fscanf(fileID, formatspec, sizeG);
                
            elseif nargin == 2
                
                sizeG = [11  Inf];
                g=fscanf(fileID, '%f', sizeG);  
                
            else     
                error('Too many or few input arguments\n')     
            end
            
            fclose(fileID);
            g=g';

            if isempty(g)
                obj = Grain();
                return; 
            end
            
            if nargin == 2                
                fprintf('Do you really need all this data?\n If so actually write this code\n')
            elseif nargin==3 
                if strcmp(formatspec, '%f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    obj.stress = g(:,1:4);
                elseif strcmp(formatspec, '%*f %*f %*f %*f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    obj.gradV  = g(:, 1:4);
                elseif strcmp(formatspec, '%*f %*f %*f %*f %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    obj.gradV  = g(:, 1:4);
                    obj.V  = g(:, 5:7);
                elseif strcmp(formatspec, '%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %f %f %f')
                    % NB: indexing in cpp starts at 0, in matlab from 1
                    obj.voro_neighbour_ID = g(:, 1:12) + 1;
                end
            end

         end
    end
    
end

