classdef Para
    % Class for all the parameters at a single time
    % Read in from ./para/parameter_[i]
    
    properties
        E = NaN          % Young's modulus
        stiff = NaN      % Stiffness (E/unit_length)
        e = NaN          % Coefficient of resitution
        visco = NaN      % Viscosity
        P = NaN          % Imposed Pressure
        shear_rate = NaN % Imposed Shear Rate
        cohesion = NaN   % Tensile force between contacts
    end
    
    methods        
        function obj = Para(E, stiff, e, visco, P, shear_rate, cohesion)    
            % Constructor
            if nargin > 0
                obj.E = E;
                obj.stiff = stiff;
                obj.e = e;
                obj.visco = visco;
                obj.P = P;
                obj.shear_rate = shear_rate;
                obj.cohesion = cohesion;
            end
        end
        
        function obj = readPara(obj, fp, formatspec)  
            % Read in grain from file
            % formatspec to read E, stiff, e, visco, P, shear_rate, cohesion
            %
            % '%*f %*f %*f %*f %f %f %*f'
            %  reads P and shear rate
            %
            % '%*f %*f %*f %*f %f %f %f'
            %  reads P, shear rate and cohesion
            %
            % '%*f %*f %*f %*f %*f %f %f'
            %  reads shear rate and cohesion
            
            fileID = fopen(fp,'r');
            if nargin == 3
                
                if strcmp(formatspec, '%*f %*f %*f %*f %f %f %*f')
                    sizeP = [2 Inf];  
                elseif strcmp(formatspec, '%*f %*f %*f %*f %f %f %f')
                    sizeP = [3 Inf];  
                elseif strcmp(formatspec, '%*f %*f %*f %*f %*f %f %f')
                    sizeP = [2 Inf];  
                end
                try
                    p=fscanf(fileID, formatspec, sizeP);
                catch
                    error('Unable to read file: %s', fp)
                end
                
            elseif nargin == 2
                
                sizeP = [7  Inf];
                try
                    p=fscanf(fileID, '%f', sizeP);
                catch
                    error('Unable to read file: %s', fp)
                end
                
            else
                error('Too many or few input arguments\n')
            end
            
            fclose(fileID);
            p=p';

            if isempty(p)
                obj = para();
                return; 
            end
            
            if nargin == 2
                obj.E          = p(:, 1);
                obj.stiff      = p(:, 2);
                obj.e          = p(:, 3);
                obj.visco      = p(:, 4);
                obj.P          = p(:, 5);
                obj.shear_rate = p(:, 6);
                obj.cohesion   = p(:, 7);
            elseif nargin==3 
                if strcmp(formatspec,'%*f %*f %*f %*f %f %f %*f')
                    obj.P          = p(:, 1);                    
                    obj.shear_rate = p(:, 2);
                elseif strcmp(formatspec, '%*f %*f %*f %*f %f %f %f')
                    obj.P          = p(:, 1);                    
                    obj.shear_rate = p(:, 2);
                    obj.cohesion   = p(:, 3);
                elseif strcmp(formatspec, '%*f %*f %*f %*f %*f %f %f')
                    obj.shear_rate = p(:, 1);
                    obj.cohesion   = p(:, 2);
                end
            end
        end
        
    end
    
end

