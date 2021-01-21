function visuOne3D(ax, dirPath, fileNum, nGrains)
    % VISUONE3S(dirPath, fileNum, trackID, nGrains)
    %
    % dirOne - directory of experiment
    % fileNum - array of fileNumbers to render
    % bool contacts - whether or not to draw contacts
    % trackID - an array with the indexes of grains to colour
    % ax - axes to plot in
    %
    % Uses speed of patch function to plot grains, contacts are still slow

%             cla(ax)

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
    axis equal


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

    % Points of a unit sphere
    n = 20;
    phi=2*pi/n*[1:n];
    theta = pi/n*[1:n];
    circX = sin(theta).*cos(phi);
    circY = sin(theta).*sin(phi);
    circZ = cos(theta);

    % Get all circle coordinates    
    X = grain.X(:,1)' + circX' .* grain.R(:)';
    Y = grain.X(:,2)' + circY' .* grain.R(:)';
    Z = circZ' .* grain.R(:)';

    C = [.8 .8 1];

    % Plot grains
    patch(ax,X,Y,Z,C, 'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'BackFaceLighting', 'lit')


end