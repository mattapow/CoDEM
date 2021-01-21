classdef Grain
    % Class for all the grains at a single time
    % Read in from ./grain/grain_[i]

    properties
        ID = NaN
        X = [NaN, NaN]
        V = [NaN, NaN]
        XTrue = [NaN, NaN]
        ome = NaN
        mass = NaN
        Imass = NaN
        R = NaN
    end

    methods

        function obj = Grain(ID, X, V, XTrue, ome, mass, Imass, R)
            % Constructor

            if nargin > 0
                obj.ID = ID;
                obj.X = X;
                obj.V = V;
                obj.XTrue = XTrue;
                obj.ome = ome;
                obj.mass = mass;
                obj.Imass = Imass;
                obj.R = R;
            end
        end        

        function obj = readGrain(obj, fp, nGrain, formatspec)
            %
            % Read in grain from file
            % formatspec to read ID, X, X, V, V, XTrue, XTrue, ome, mass,
            % Imass, R. Only certain combinations allowed:
            %   '%d %f %f %*f %*f %*f %*f %*f %*f %*f %f' reads ID, X, R
            %   '%d %f %f %f %f %*f %*f %*f %*f %*f %f' reads ID, X, V, R
            %   '%d %f %f %*f %*f %*f %*f %*f %*f %*f %*f' reads ID X
            %   '%d %*f %*f %*f %*f %*f %*f %*f %f %*f %*f' reads ID mass
            %   '%d %f %f %*f %*f %*f %*f %*f %f %*f %*f' reads ID X mass
            %   '%d %*f %*f %*f %*f %f %f %*f %*f %*f %*f' reads ID XTrue
            %   '%d %f %f %f %f %*f %*f %*f %*f %*f %*f' reads ID X V
            %   '%d %f %f %f %f %*f %*f %*f %f %*f %*f' reads ID X V mass
            %   '%d %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f' reads ID
            %   '%d %*f %*f %f %f %f %f %*f %f %*f %*f' reads ID, V, XTrue, mass
            %   '%d %f %f %*f %*f %f %f %*f %*f %*f %*f' reads ID, X, XTrue

            if nargin == 3
                formatspec = '%d %f %f %f %f %f %f %f %f %f %f';
            end
            
            fileID = fopen(fp,'r');
            g = textscan(fileID, formatspec, nGrain, 'ReturnOnError', 0);
            fclose(fileID);

            if isempty(g)
                obj = Grain();
                return;
            end

            if nargin == 3
                obj.ID          = cell2mat(g(1));
                obj.X           = cell2mat([g(2), g(3)]);
                obj.V           = cell2mat([g(4), g(5)]);
                obj.XTrue       = cell2mat([g(6), g(7)]);
                obj.ome         = cell2mat(g(8));
                obj.mass        = cell2mat(g(9));
                obj.Imass       = cell2mat(g(10));
                obj.R           = cell2mat(g(11));
            elseif nargin==4
                if strcmp(formatspec,     '%d %f %f %*f %*f %*f %*f %*f %*f %*f %f')
                    obj.ID          = cell2mat(g(1));
                    obj.X           = cell2mat([g(2), g(3)]);
                    obj.R           = cell2mat(g(4));
                elseif strcmp(formatspec, '%d %f %f %f %f %*f %*f %*f %*f %*f %f')
                    obj.ID          = cell2mat(g(1));
                    obj.X           = cell2mat([g(2), g(3)]);
                    obj.V           = cell2mat([g(4), g(5)]);
                    obj.R           = cell2mat(g(6));
                elseif strcmp(formatspec, '%d %f %f %*f %*f %*f %*f %*f %*f %*f %*f')
                    obj.ID          = cell2mat(g(1));
                    obj.X           = cell2mat([g(2), g(3)]);
                elseif strcmp(formatspec, '%d %f %f %f %f %*f %*f %*f %f %*f %*f')
                    obj.ID          = cell2mat(g(1));
                    obj.X           = cell2mat([g(2), g(3)]);
                    obj.V           = cell2mat([g(4), g(5)]);
                    obj.mass        = cell2mat(g(6));
                elseif strcmp(formatspec, '%d %*f %*f %*f %*f %*f %*f %*f %f %*f %*f')
                    obj.ID          = cell2mat(g(1));
                    obj.mass        = cell2mat(g(2));
                elseif strcmp(formatspec, '%d %*f %*f %*f %*f %f %f %*f %*f %*f %*f')
                    obj.ID          = cell2mat(g(1));
                    obj.XTrue       = cell2mat([g(2), g(3)]);
                elseif strcmp(formatspec, '%d %f %f %f %f %*f %*f %*f %*f %*f %*f')
                    obj.ID          = cell2mat(g(1));
                    obj.X           = cell2mat([g(2), g(3)]);
                    obj.V           = cell2mat([g(4), g(5)]);
                elseif strcmp(formatspec, '%d %f %f %f %f %*f %*f %*f %f %*f %*f')
                    obj.ID          = cell2mat(g(1));
                    obj.X           = cell2mat([g(2), g(3)]);
                    obj.V           = cell2mat([g(4), g(5)]);
                    obj.mass        = cell2mat([g(6)]);                    
                elseif strcmp(formatspec, '%d %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f')
                    obj.ID          = cell2mat(g(1));
                elseif strcmp(formatspec, '%d %*f %*f %f %f %f %f %*f %f %*f %*f')
                    obj.ID          = cell2mat(g(1));
                    obj.V           = cell2mat([g(2), g(3)]);
                    obj.XTrue       = cell2mat([g(4), g(5)]);
                    obj.mass        = cell2mat(g(6));
                elseif(strcmp(formatspec, '%d %f %f %f %f %*f %*f %*f %f %*f %*f'))
                    obj.ID          = cell2mat(g(1));
                    obj.X           = cell2mat([g(2), g(3)]);
                    obj.V           = cell2mat([g(4), g(5)]);
                    obj.mass        = cell2mat(g(6));
                elseif strcmp(formatspec, '%d %f %f %*f %*f %f %f %*f %*f %*f %*f')
                    obj.ID          = cell2mat(g(1));
                    obj.X           = cell2mat([g(2), g(3)]);
                    obj.XTrue       = cell2mat([g(4), g(5)]);
                else
                    error('Invalid formatspec')
                end
            end

        end
        
        function obj = readGrainUnfold(obj, fp, nGrain)
            % obj = READGRAINUNFOLD(obj, fp, nGrain, formatspec) Read and unfold XTrue
            % XTrue = X + XTrue*L
            
            obj = readGrain(obj, fp, nGrain);            
            assert(all(mod(obj.XTrue,1)==0, 'all'))
            
            dp = split(fp, 'grain/grain_');
            fp = strcat(dp{1}, '/cell/cell_', dp{2});
            cell = Cell();
            formatspec = '%*f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
            cell = readCell(cell, fp, formatspec);
            
            dp = split(fp, '/cell/cell_');
            fp = strcat(dp{1}, '/para/parameter_', dp{2});
            para = Para();
            formatspec='%*f %*f %*f %*f %f %f %*f';
            para = readPara(para, fp, formatspec);
                        
            obj.V(:,1) = obj.V(:,1) + obj.XTrue(:,2) * cell.L(2) * para.shear_rate;
            
            obj.XTrue(:,1) = obj.X(:,1) + obj.XTrue(:,1) * cell.L(1)+obj.XTrue(:,2)*cell.shift;
            obj.XTrue(:,2) = obj.X(:,2) + obj.XTrue(:,2) * cell.L(2);

            
            % remove affine
%             obj.XTrue(:,1) = obj.XTrue(:,1) - obj.XTrue(:,2) * cell.L(2) * para.shear_rate;
            
%             % normalise by cell size and shift up to unit grid
%             obj.XTrue(:,1) = obj.XTrue(:,1) / cell.L(1);
%             obj.XTrue(:,2) = obj.XTrue(:,2) / cell.L(2) + .5;
%             Leff = (cell.L(1)+cell.L(2))/2;
%             obj.R = obj.R / Leff;
            
        end
        
        function plotGrain3d(obj,ax,fp,nPoints)
            % Plot each grain in this object
            % Looks nicer plotting spheres but is slower
            nG = length(obj.ID);
            [x,y,z] = sphere(nPoints);
%             cmap=colormap(brewermap([],'YlOrRd'))
%             cmap=colormap('copper');
%             cmap=colormap(brewermap([],'Reds'));            
%             cmap=colormap('parula');
%             cmap = colormap(flipud(customCols()));
            cmap = colormap(plasma);
%             cmap=colormap('jet');            
            nCmap = size(cmap, 1);
            
%             nCmap = 64;
%             cmap = [[0:nCmap-1]' zeros(nCmap,2)]/(nCmap-1);
%             cmap = colormap(cmap);
            
            
            contact = Contact();    
            contact = readContact(contact, fp);

            Fy = sparse(double(contact.ID_A)+1,double(contact.ID_B)+1,contact.Force(:,2), nG, nG);
            Fy = Fy-Fy';
            Fx = sparse(double(contact.ID_A)+1,double(contact.ID_B)+1,contact.Force(:,1), nG, nG);
            Fx = Fx-Fx';
            da = (sum(Fx).^2 + sum(Fy).^2).^.5;
            CMax=150; %132 75.0013;
            if max(da)>CMax; da(da>CMax) = CMax; end
%             warning('normalisation hardcoded. da max = %f', CMax)
            CMax = log(CMax*1+1);
            C = log(da*1+1);
            
            colIndx=ceil(eps()+C/CMax * nCmap);
            col = cmap(colIndx,:);
            
            hold(ax,'on')
%             fprintf('%4.0f%%\n',0/nG*100)
            
            epsilon=4; % grains outisde frame
            for i = 1:nG
%                 fprintf('\b\b\b\b\b%4.0f%%',i/nG*100)
                r = obj.R(i);
                X_x = obj.X(i,1);
                X_y = obj.X(i,2);
                if X_x < ax.XLim(1)-epsilon || X_x > ax.XLim(2)+epsilon ||...
                   X_y < ax.YLim(1)-epsilon || X_y > ax.YLim(2)+epsilon
                    continue; 
                end
            
                color=col(i)*ones(size(z)); % color=C(i)/max(C)*nCmap*ones(size(z));
                surf(ax,(X_x-r*x), (X_y-r*y), r*z-1, color,...
                    'EdgeColor', 'none',...
                    'FaceColor','flat',...
                    'Facelighting','gouraud',...
                    'AmbientStrength',.5)
%                     'BackFaceLighting','lit',...
%                     'DiffuseStrength',1,...
%                     'SpecularStrength',1,...
%                     'SpecularExponent',1,...
%                     'SpecularColorReflectance',1)
            end
%             fprintf('\n')
            hold(ax,'off')


        end

        function plotGrain(obj, tag, outside)
            % Plot each grain in this object as circles
            % outside=1 -> use XTrue
            
            
            n = length(obj.ID);
            if(nargin==1)
                tag=[];
            end
            
            if(outside==1 && any(any(isnan(obj.XTrue))))
                error('XTrue is nan')
            end

            hold on
            for i = 1:n

                r = obj.R(i);
                if outside
                    X_x = obj.XTrue(i,1);
                    X_y = obj.XTrue(i,2);
                else
                    X_x = obj.X(i,1);
                    X_y = obj.X(i,2);
                end
                pos = [X_x-r X_y-r 2*r 2*r];

                % Color the grain if it's ID is in trackID
                if ismember(obj.ID(i), tag)
                    facecolor = [204/256 0/256 0/256];
                else
                    facecolor = [142/256 162/256 144/256];
%                     facecolor='none';
                end

                % Draw a circle
                rectangle('Position',pos,'Curvature',[1 1], ...
                    'FaceColor', facecolor, 'Edgecolor', 'black', 'LineWidth', .5)

            end
            hold off

        end
        
        
        

        function plotGrainClust(obj, clustID)
            % Plot each grain in this object as circles
            % Show non-singleton clusters in different colors
            n = length(obj.ID);
            nClust = max(clustID);
%             nClust = 1000;

            hold on
            for i = 1:n

                r = obj.R(i);
                X_x = obj.X(i,1);
                X_y = obj.X(i,2);
                pos = [X_x-r X_y-r 2*r 2*r];

                colormap colorcube
                cm = colormap; % returns the current color map

                if (sum(clustID == clustID(i))==1)
                    % White singleton clusters
                    facecolor = [1 1 1];
%                 elseif  clustID(i) == 1
%                     % Grey background cluster. NB: it may not be 1
%                     facecolor = [.5 .5 .5];
                else
                    colorID = max(1, sum(clustID(i)/nClust > [0:1/length(cm(:,1)):1]));
                    facecolor = cm(colorID, :); % returns your color
                end

                % Draw a circle
                rectangle('Position',pos,'Curvature',[1 1], ...
                    'FaceColor', facecolor, 'Edgecolor', 'black', 'LineWidth', .5)

            end
            hold off

        end

        function plotVonoroi(obj)
            % Plot the vonoroi network of grains
            voronoi(obj.X(:,1),obj.X(:,2), 'k')
        end
        
    end
    
    methods (Static)
       function plotTagGrain(ax,dp, NG, fn, tag,C)
            % Plot each grain in this object as circles
            % outside=1 -> use XTrue
            % omitflag=1 -> don't plot particles that aren't tagged
            % grainID tag indexed from 1
            
            % Points of a unit circle
            nPoints = 30;
            circX = cos(2*pi*(1:nPoints)/nPoints)';
            circY = sin(2*pi*(1:nPoints)/nPoints)';
            
%             % leave colored trail of particle position
%             if fn>1
%                 grain0=Grain();
%                 fp = strcat(dp, '/grain/grain_', string(fn-1));
%                 formatspec = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %f';
%                 grain0 = readGrain(grain0, fp, NG,formatspec);
%                 
%                 grain0.ID=grain0.ID(tag);
%                 grain0.X=grain0.X(tag,:);
%                 grain0.R=grain0.R(tag);
%                 
%                 XX0 = grain0.X(:,1)' + circX .* grain0.R(:)';
%                 YY0 = grain0.X(:,2)' + circY .* grain0.R(:)';
%                 
%                 patch(ax,XX0,YY0,C,...
%                 'EdgeColor', 'none',...
%                 'FaceAlpha', 1)    
%             end
            
            
            
            grain=Grain();
            fp = strcat(dp, '/grain/grain_', string(fn));
            formatspec = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %f';
            grain = readGrain(grain, fp, NG,formatspec);
            
            grain.ID=grain.ID(tag);
            grain.X=grain.X(tag,:);
            grain.R=grain.R(tag);

            
            % Set axis limits
            cell = Cell();
            fpCell = strcat(dp, '/cell/cell_', string(fn));
            cell = cell.readCell(fpCell, '%f %f %f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f');
            set(ax, 'XLim', [0-1 cell.L(1)+1]);
            set(ax, 'YLim', [-cell.L(2)/2-1 +cell.L(2)/2+1]);
            
            

            % Get all circle coordinates
            XX = grain.X(:,1)' + circX .* grain.R(:)';
            YY = grain.X(:,2)' + circY .* grain.R(:)';
            
            
             patch(ax,XX,YY,C,...
                'EdgeColor', 'none',...
                'FaceAlpha', 1)    
       end 
        
       function tagged = tagEven(dp, nGrains, t0,eps)
           % tagged = tagEven(dp, nGrains, t0,eps) Evely tag grains
           %
           % Tag grains in the centre vertical strip with even spacing
           % eps = size of square to search for grains. e.g. eps=.5
    
            % Read grains
            fpGrain = strcat(dp, '/grain/grain_', string(t0));
            formatspec = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %f';
            grain = Grain();
            grain = grain.readGrain(fpGrain, nGrains, formatspec);

             % Read in Cell time, length and shift
            cell = Cell();
            fpCell = strcat(dp, '/cell/cell_', string(t0));
            cell = cell.readCell(fpCell, '%f %f %f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f');

            % dx strip in middle
            X=grain.X(:,1);
            tag1 = grain.ID(abs(X-50)<eps)+1;

            % y strips at 10d spacings
            dh=10;
            h = 0:dh:cell.L(2);
            h = [-fliplr(h) h]';
            Y=grain.X(:,2);
            tag2 = [];
            for i=1:numel(h)
                tag2 = [tag2; grain.ID(abs(Y-h(i))<eps)+1];
            end
            tagged = intersect(tag1,tag2);
            
            [~,I]=sort(Y(tagged));
            tagged=tagged(I);
            
            dY=diff(Y(tagged));
            tagged(abs(dY)<1)=[];

        end

    end


end
