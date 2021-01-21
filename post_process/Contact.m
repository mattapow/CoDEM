classdef Contact
    % Class for all the contact at a single time
    % Read in from ./contact/contact_[i]
    
    properties
        ID_A = NaN 
        ID_B = NaN
        Force = [NaN NaN]
        fp = '~/documents/MATLAB'
    end
    
    methods        
        function obj= Contact(ID_A, ID_B, Force)
            % Constructor
            if nargin > 0
                obj.ID_A = ID_A;
                obj.ID_B = ID_B;
                obj.Force = Force;
            end
        end
                
        function obj = readContact(obj, fp, formatspec)
            % Read in contacts from file           
            % Don't bother with formatspec, haven't written

            if nargin == 2
                formatspec = '%d %d %f %f';
            end
                                    
            nContact = PostP.nLines1(fp);
            
            fileID = fopen(fp,'r'); 
            c = textscan(fileID, formatspec, nContact, 'ReturnOnError', 0);
            fclose(fileID);

            if isempty(c)
                obj = Contact();
                return;
            end

            if nargin == 2
                obj.ID_A      = cell2mat(c(1));
                obj.ID_B      = cell2mat(c(2));
                obj.Force     = cell2mat([c(3), c(4)]);
            else
                error('No code written for formatspec input')
            end
            
        end
        
        function plotContact(obj, grain, cell)
            % Plot the contacts in 2D
            % Grain has ID and X only
            % Cell has time, L and shift
            
            % if no contacts, do nothing
            if isnan(obj.ID_A)
                return;
            end
            
            % Normalise contact size by the maximum contact force
            scale=3; %1.7 works well for cohesionless
            scale2=1;%50 works well for all cohesion
            Fmax=30; %1000; %5318;
            Fmax = scale*log(Fmax*scale2+1);
%             Fmax = max(hypot(obj.Force(:,1), obj.Force(:,2)));

%             cmap=colormap(brewermap([],'YlOrRd'));
%             cmap = colormap('hot');            
%             cmap = colormap(cmap(1:48,:));
            cmap=colormap(plasma);
            nCmap = size(cmap, 1);
            
            % sort so that bigger are on top
            [~,fIdx]=sort(obj.Force(:,1).^2+obj.Force(:,2).^2);
            
            
            nC = length(obj.ID_A);
            for i = 1:nC
                c=fIdx(i);                
                x1 = grain.X(obj.ID_A(c)+1, 1);
                x2 = grain.X(obj.ID_B(c)+1, 1);
                y1 = grain.X(obj.ID_A(c)+1, 2);
                y2 = grain.X(obj.ID_B(c)+1, 2);
                
                dX = [x2 - x1, y2 - y1];
                
                v1 = dX(:,1); v2 = dX(:,2);
                [v1, v2] = cell.CLP(v1, v2);
                dX = [v1 v2];
                
                x2 = x1 + dX(1);
                y2 = y1 + dX(2);                
                
                nF = hypot(obj.Force(c, 1), obj.Force(c, 2));                
                nF = scale*log(nF*scale2+1);
                colIndx=min(ceil(nF/Fmax * nCmap),nCmap);
                col = cmap(colIndx,:);

                        
                % Show blue if net attractive force
                % Use dot product of force direction and dX direction
%                 ndX = norm(dX); 
%                 dirX = dX ./ ndX;
%                 F = obj.Force(i, :);                                
%                 dirF = F ./ nF;                
%                 if (dot(dirX,dirF) > 0)
%                     col='blue';
%                 else                    
%                     col='red';
%                 end
                
                line([x1 x2], [y1 y2],...
                    'color', col,...
                    'LineWidth', nF,...
                    'Marker','o',...
                    'MarkerSize',nF,...
                    'MarkerFaceColor',col,...
                    'MarkerEdgeColor','none')
                
            end            
            
        end
        
        function plotContactLight(obj, grain, cell)
            % Plot the contacts as cylinders
            % Grain has ID and X only
            % Cell has time, L and shift
            
            % if no contacts, do nothing
            if isnan(obj.ID_A); return; end
            
            % Normalise contact size by the maximum contact force
            Fmax = max(hypot(obj.Force(:,1), obj.Force(:,2)));            
            
            hold on
%             lgt = camlight(0, 90);
%             lighting phong

            points = 20;
            t = 0:pi/10:2*pi;
            R=1;
            [X,Y,Z] = cylinder(R,points);
            
            
            n = length(obj.ID_A);
            
            for i = 1:n
                % Master grain: 1, slave grain: 2
                x1 = grain.X(obj.ID_A(i)+1, 1);
                x2 = grain.X(obj.ID_B(i)+1, 1);
                y1 = grain.X(obj.ID_A(i)+1, 2);
                y2 = grain.X(obj.ID_B(i)+1, 2);
                
                dX = [x2 - x1, y2 - y1];
                
                v1 = dX(:,1); v2 = dX(:,2);
                [v1, v2] = cell.CLP(v1, v2);
                dX = [v1 v2];                
                                
                F = obj.Force(i, :);                                
                nF = norm(F); %/Fmean;                
                
                cmap = colormap('parula');
                nCmap = size(cmap, 1);
                colIndx=ceil(nF/Fmax * nCmap);
                col = cmap(colIndx,:);               

%                 R = abs(nF)/Fmax*.5;
                R = max(log(abs(nF/Fmax)), .1)*.6;                                
                d = (dX(1)^2+dX(2)^2)^.5;                
                h=mesh(X,Y,Z);
                rotate(h,[1 0 0],270)
                h.YData = h.YData + 1/2;   
                h.XData = h.XData*R + x1;
                h.YData = h.YData*d + y1;
                h.ZData = h.ZData*R;                

                rotate(h,[0 0 -1],atand(dX(1)/dX(2)),[x1,y1,0])
                surf(h.XData, h.YData, h.ZData, ...
                    'FaceColor', col, 'FaceLighting','gouraud', ...
                    'LineStyle', 'none')
                
                                             
                if mod(i,100)==0                    
                    drawnow                  
                end
                
            end

            
            
        end
    end
    
    
    
end

