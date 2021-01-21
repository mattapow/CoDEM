function visuOne(ax, dirPath, fileNum, trackID, cPlot, nGrains)
    % VISUONE(ax, dirPath, fileNum, trackID, cPlot, nGrains)
    %
    % ax - axis object
    % dirPath - directory of experiment
    % fileNum - array of fileNumbers to render
    % cPlot - 0: don't plot contacts, 1 plot all contacts, 2 plot inter-species contacts only
    % trackID - an array with the indexes (indexed from 0) of grains to colour
    % nGrains - number of grains
    %
    % Uses speed of patch function to plot grains, plotting contacts is still slow

    % Read grains
    fpGrain = strcat(dirPath, '/grain/grain_', string(fileNum));
    formatspec = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %f';
    grain = Grain();
    grain = grain.readGrain(fpGrain, nGrains, formatspec);


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
    if isempty(trackID)
        % uniform colour
%                 C = [207 207 207]/255; % Light grey
        C = [127 127 127]/255; % Medium grey
        
%     elseif cPlot==2
%         % plot the grains according to their proportion of interparticle
%         % contacts
%         
%         contact = Contact();
%         fpContact = strcat(dirPath, '/contact/contact_', string(fileNum));
%         contact = contact.readContact(fpContact);
%         Z = sum(sparse(double(contact.ID_A)+1, double(contact.ID_B)+1, 1,nGrains,nGrains));
%         
%         contact = interContacts(contact,trackID);
%         contact = removeBoundary(contact, grain,cell);
%         
%         ZInter = sum(sparse(double(contact.ID_A)+1, double(contact.ID_B)+1, 1,nGrains,nGrains));        
%         
%         colormap(brewermap(64,'OrRd'))
%         C = ZInter./Z;
%         C(Z==0)=0;
%         
    else
        % two colours by trackID
        
        C = zeros(nPoints, length(grain.R), 3);
        C(:,:,2) = ones(nPoints, length(grain.R), 1)/2;
   

%         C(:,trackID+1,2) = zeros(nPoints, numel(trackID), 1);
        C(:,trackID+1,3) = ones(nPoints, numel(trackID), 1);
%         C(:,trackID+1,2) = ones(nPoints, numel(trackID), 1);
    end    

    % Plot grains
    patch(ax,X,Y,C,...
        'EdgeColor', 'none',...
        'FaceAlpha', .6)

    % Plot the contacts
    if (cPlot==1 || cPlot==2 )
        
        contact = Contact();
        fpContact = strcat(dirPath, '/contact/contact_', string(fileNum));
        contact = contact.readContact(fpContact);
        
        if cPlot==2
            contact = interContacts(contact,trackID);
            contact = removeBoundary(contact, grain,cell);
        end
        
        contact.plotContact(ax,grain, cell);
        
    end

end

function contact = interContacts(contact,trackID)

    % Only keep inter-particle contacts, by trackID
    n=length(contact.ID_A);
    for c=n:-1:1
        if (ismember(contact.ID_A(c), trackID) && ismember(contact.ID_B(c), trackID)) ||...
            (~ismember(contact.ID_A(c), trackID) && ~ismember(contact.ID_B(c), trackID))
            contact.ID_A(c)=[];
            contact.ID_B(c)=[];
            contact.Force(c,:)=[];
        end
    end


end
    
function contact = removeBoundary(contact, grain,cell)
 % remove contacts accross boundary
 
 n=length(contact.ID_A);
    for c=n:-1:1        
        dY = grain.X(contact.ID_A(c)+1,2) - grain.X(contact.ID_B(c)+1,2);
        if abs(dY) > cell.L(2)/2
            contact.ID_A(c)=[];
            contact.ID_B(c)=[];
            contact.Force(c,:)=[];
        end
    end

end