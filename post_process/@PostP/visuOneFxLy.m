function visuOneFxLy(ax, dirPath, fileNum, nGrains)
% VISUONEFXLY(ax, dirPath, fileNum, nGrains)
% Plot grains and colour them by the sum over their contacts of Fx*Ly

% Read grains
fpGrain = strcat(dirPath, '/grain/grain_', string(fileNum));
formatspec = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %f';
grain = Grain();
grain = grain.readGrain(fpGrain, nGrains, formatspec);

% Read contacts
fpContact = strcat(dirPath, '/contact/contact_', string(fileNum));
contact = Contact();
contact = contact.readContact(fpContact);

% Read in Cell time, length and shift
fpCell = strcat(dirPath, '/cell/cell_', string(fileNum));
cell = Cell();
cell = cell.readCell(fpCell, '%f %f %f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f');

% Set limits
set(ax, 'XLim', [0-1 cell.L(1)+1]);
set(ax, 'YLim', [-cell.L(2)/2-1 +cell.L(2)/2+1]);

% Points of a unit circle
nPoints = 30;
circX = cos(2*pi*(1:nPoints)/nPoints)';
circY = sin(2*pi*(1:nPoints)/nPoints)';

% Get all circle coordinates
X = grain.X(:,1)' + circX .* grain.R(:)';
Y = grain.X(:,2)' + circY .* grain.R(:)';

% Set weight for colours Fx*Ly
weight = setWeight(contact, grain, nGrains, cell);

% Set all cicle colours
C = .95*ones(nPoints, nGrains, 3);          

% red negative
C(:,weight<0,1) = .95+.05*ones(nPoints, sum(weight<0)).*weight(weight<0)/min(weight);
C(:,weight<0,2) = .95-.95*ones(nPoints, sum(weight<0)).*weight(weight<0)/min(weight);
C(:,weight<0,3) = .95-.95*ones(nPoints, sum(weight<0)).*weight(weight<0)/min(weight);

% blue positive
C(:,weight>0,1) = .95-.95*ones(nPoints, sum(weight>0)).*weight(weight>0)/max(weight);
C(:,weight>0,2) = .95-.95*ones(nPoints, sum(weight>0)).*weight(weight>0)/max(weight);
C(:,weight>0,3) = .95+.05*ones(nPoints, sum(weight>0)).*weight(weight>0)/max(weight);

% C(:,weight==0,:) = 1;

% weight = abs(weight);
% C = weight;
% colormap(flipud(pink))

% Plot grains
patch(ax,X,Y,C, 'EdgeColor', 'none')  
            
end

function weight = setWeight(contact, grain, nGrains, cell)
% Set grain weight to sum over its contacts of (Fx*Ly) 
% 
% You could insert Fx = Fx - Fx' to get the full force matrix
% then sum the upper triangular part to get the same result

    ID_A = double(contact.ID_A(:))+1;
    ID_B = double(contact.ID_B(:))+1;

    Fx = sparse(ID_A,ID_B,contact.Force(:,1),nGrains,nGrains);

    x = grain.X(:,1);
    y = grain.X(:,2);

    thisLx = x(ID_B)-x(ID_A);
    thisLy = y(ID_B)-y(ID_A);       
    [~, thisLy] = CLP(cell, thisLx, thisLy);
    Ly = sparse(ID_A,ID_B,thisLy,nGrains,nGrains);
    
    warning('Unsure about the sign of Ly')
    FXLY = Fx .* Ly;
    weight = sum(FXLY);
    
%     warning('corr not sum')
%     weight=nan(nGrains);
%     for i = 1:nGrains
%         if sum(Fx(:,i)~=0) > 1
%             weight(i) = corr(Fx(:,i), Ly(:,i));
%         end
    end
end