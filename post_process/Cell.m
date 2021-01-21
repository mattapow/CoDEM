classdef Cell
    % Class for all the Cell information at a single time
    % Read in from ./cell/cell_[i]
    
    properties
        time = NaN
        L = [NaN NaN]
        V = [NaN NaN]
        mass = NaN
        shift = NaN
        phi = NaN
        stress = [NaN NaN NaN NaN]
        boundary_normal_stress = NaN
        boundary_shear_stress = NaN
    end
    
    methods
        
        function obj = Cell(time, L, V, mass, shift, phi, stress, boundary_normal_stress, boundary_shear_stress)   
            % Constructor
            if nargin > 0
                obj.time = time;
                obj.L = L;
                obj.V = V;
                obj.mass = mass;
                obj.shift = shift;
                obj.phi = phi;
                obj.stress = stress;
                obj.boundary_normal_stress = boundary_normal_stress;
                obj.boundary_shear_stress = boundary_shear_stress;
            end
        end
        
        function obj = readCell(obj, fp, formatspec)
            % obj = readCell(obj, fp, formatspec)
            % Read in cell from file 
            % optional formatspec to read =time, L, V, mass, shift, phi,
            % stress, boundary_normal_stress, boundary_shear_stress
            % Only certain combinations coded: 
            %   '%f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f'
            %   reads time, L
            %
            %   '%f %f %f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f'
            %   reads time, L shift
            %
            %   '%f %*f %*f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f'
            %   reads time and shift
            %
            %   '%*f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f'
            %    reads L
            %
            %   '%*f %f %f %*f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f'
            %    reads L and phi
            %            
            %   '%*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %*f %*f'
            %    reads stress
            %
            %   '%f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %*f %*f'
            %    reads t phi and stress
            %
            
                        
            %% Read saved variables to a matrix
            fileID = fopen(fp,'r'); 
            
            if nargin == 3
                
                if strcmp(formatspec, '%f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    sizeC = [3  Inf];   
                elseif strcmp(formatspec, '%f %f %f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f')
                    sizeC = [4 Inf];
                elseif strcmp(formatspec, '%f %*f %*f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f')
                    sizeC = [2 Ing];
                elseif strcmp(formatspec, '%*f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    sizeC = [2 Inf];
                elseif strcmp(formatspec, '%*f %f %f %*f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f')
                    sizeC = [3 Inf];
                elseif strcmp(formatspec, '%*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %*f %*f')
                    sizeC = [4 Inf];
                elseif strcmp(formatspec, '%f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %*f %*f')
                    sizeC = [6 Inf];
                else
                    error('Invalid formatspec to read grains')
                end                
                c=fscanf(fileID, formatspec, sizeC);           
            elseif nargin == 2   
                
                sizeC = [14  1];
                c=fscanf(fileID, '%f', sizeC);              
            else               
                error('Too many or few input arguments\n')
            end
            
            fclose(fileID);
            c=c';
            if isempty(c)
                obj = Cell();
                return; 
            end
            
            %% Assign read variables to object
            if nargin == 2             
                
                obj.time                    = c(:, 1);
                obj.L                       = c(:, 2:3);
                obj.V                       = c(:, 4:5);
                obj.mass                    = c(:, 6);
                obj.shift                   = c(:, 7);
                obj.phi                     = c(:, 8);
                obj.stress                  = c(:, 9:12);
                obj.boundary_normal_stress  = c(:, 13);
                obj.boundary_shear_stress   = c(:, 14);
                
            elseif nargin==3 
                
                if strcmp(formatspec, '%f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    obj.time = c(:, 1);
                    obj.L    = c(:, 2:3);
                elseif strcmp(formatspec, '%f %f %f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f')
                    obj.time = c(:, 1);
                    obj.L    = c(:, 2:3);
                    obj.shift  = c(:, 4);
                elseif strcmp(formatspec, '%f %*f %*f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f')
                    obj.time = c(:, 1);
                    obj.shift = c(:,2);
                elseif strcmp(formatspec, '%*f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    obj.L    = c(:, 1:2);
                elseif strcmp(formatspec, '%*f %f %f %*f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f')
                    obj.L    = c(:, 1:2);
                    obj.phi  = c(:, 3);
                elseif strcmp(formatspec, '%*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %*f %*f')
                    obj.stress = c(:, 1:4);
                 elseif strcmp(formatspec, '%f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %*f %*f')
                     obj.time = c(:, 1);
                     obj.phi  = c(:, 2);
                     obj.stress = c(:, 3:6);
                end
                
            end

        end
        
        function plotCell(obj)
            % Plot cell boundaries
            xMax = obj.L(1);
            
            yMin = -obj.L(2)/2;
            yMax = +obj.L(2)/2;
            
            line([0 xMax], [yMin yMin])
            line([0 xMax], [yMax yMax])
            line([0 0], [yMin yMax])
            line([xMax xMax], [yMin yMax])
        end
        
        function [v1, v2] = CLP(obj, v1, v2) 
            % [v1, v2] = CLP(obj, v1, v2) 
            % Adjust vector V = [v1, v2] for periodic Lees-Edwards boundary
            % v1, v2 can be matrixes with each pair {v1(i,j), v2{i,j)} a
            % vector
            
            L = obj.L(1); %#ok<*PROPLC>
            H = obj.L(2);
            shift = obj.shift;            
            
            v1 = v1 - (v1>L/2).*L + (v1<-L/2).*L - (v2>=H./2).*shift + (v2<-H./2).*shift;
            v2 = v2 - (v2>=H./2).*H + (v2<-H./2).*H;
            v1 = v1 - (v1>L/2).*L + (v1<-L/2).*L;            

        end        
        
        function L = getL(obj)
            L = obj.L;
        end
    end 
    
end

