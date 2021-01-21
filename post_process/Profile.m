classdef Profile
    %Profile
    % Class for all the profile information at a single time
    % Read in from ./profile/profile_[i]
    
    properties
        y = NaN
        y_H = NaN
        phi = NaN
        V = [NaN NaN NaN]
        gradV = [NaN NaN NaN]
        stress = [NaN NaN NaN NaN NaN NaN NaN NaN NaN]
        V_H = [NaN NaN NaN]
        shear_rate_normalised = NaN
        ome = NaN
        dV2 = NaN
    end
    
    methods
        function obj = Profile(y, y_H, phi, V, gradV, stress, V_H, shear_rate_normalised, ome, dV2)
            % Constructor
            if nargin > 0
                obj.y = y;
                obj.y_H = y_H;
                obj.phi = phi;
                obj.V = V;
                obj.gradV = gradV;
                obj.stress = stress;
                obj.V_H = V_H;
                obj.shear_rate_normalised = shear_rate_normalised;
                obj.ome = ome;
                obj.dV2 = dV2;
            end
        end
        
        function obj = readProfile(obj, fp, formatspec)
            % Read in profile from file
            %
            % '%f %*f %*f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f'
            % Reads y, V
            %
            % '%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %*f %*f %*f %f %*f %*f'
            % Reads stress and shear_rate_normalised            
            %
            % '%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f'
            % Reads stress
            %
            % '%f %*f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f'
            % Reads y, phi
            %
            % '%f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f'
            % Reads y, stress
            
            %% Read saved variables to a matrix     
            fileID = fopen(fp,'r');
            
            if nargin == 3

                if strcmp(formatspec, '%f %*f %*f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    sizeP = [4  Inf];   
                elseif strcmp(formatspec, '%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %*f %*f %*f %f %*f %*f')
                    sizeP = [10  Inf];   
                elseif strcmp(formatspec, '%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f')
                    sizeP = [9  Inf];   
                elseif strcmp(formatspec, '%f %*f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    sizeP = [2  Inf];   
                elseif strcmp(formatspec, '%f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f')
                    sizeP = [10  Inf];   
                else
                    error('Invalid formatspec to read profile')
                end                
                p=fscanf(fileID, formatspec, sizeP);  

      
            elseif nargin == 2                   
                sizeP = [30  1];
                p=fscanf(fileID, '%f', sizeP); 
            else               
                error('Too many or few input arguments\n')
            end
            
            fclose(fileID);
            p=p';
            if isempty(p)
                obj = Profile();
                return; 
            end
            
            %% Assign read variables to object            
             if nargin == 2             
                
                obj.y                       = p(:,1);
                obj.y_H                     = p(:,2);
                obj.phi                     = p(:,3);
                obj.V                       = p(:,4:6);
                obj.gradV                   = p(:,7:15);
                obj.stress                  = p(:,16:24);
                obj.V_H                     = p(:,25:27);
                obj.shear_rate_normalised   = p(:,28);
                obj.ome                     = p(:,29);
                obj.dV2                     = p(:,30);                  
                
            elseif nargin==3 
                
                if strcmp(formatspec, '%f %*f %*f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    obj.y = p(:, 1);
                    obj.V = p(:, 2:4);                
                elseif strcmp(formatspec, '%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %*f %*f %*f %f %*f %*f')
                    obj.stress                  = p(:, 1:9);
                    obj.shear_rate_normalised   = p(:, 10);                                   
                elseif strcmp(formatspec, '%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f')
                    obj.stress                  = p(:, 1:9);
                elseif strcmp(formatspec, '%f %*f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    obj.y   = p(:, 1);
                    obj.phi = p(:, 2);
                elseif strcmp(formatspec, '%f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f')
                    obj.y      = p(:, 1);
                    obj.stress = p(:, 2:10);
                end
                
            end
            
        end
        
    end
    
end

