function obj = covFL(obj, isCorr)
% average (over N grains and T times) covariance between contact x-force and y-length.
% optional flag isCorr to find correlation instead of covariance


    T1 = obj.nPass+1;
    T2 = obj.tEnd;

    c = nan(T2-T1+1, 1);

    dpCon = strcat(obj.dirPath, '/contact/contact_');
    contact = Contact();
    
    dpGg = strcat(obj.dirPath, '/grain/grain_');
    grain = Grain();
    formatspecG = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %*f';
    
    dpCell = strcat(obj.dirPath, '/cell/cell_');
    cell = Cell();
    formatspecCell = '%f %f %f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f';
    
    warning('Changed to Fy, lx')
    
    for fn = T1:T2
        
        % Read contacts
        fpC = strcat(dpCon, string(fn));    
        contact = contact.readContact(fpC);
        
        % Read grain positions
        fpG = strcat(dpGg, string(fn));        
        grain = grain.readGrain(fpG, obj.nG, formatspecG);
        
        % Read cell dimensions at this time
        fpCell = strcat(dpCell, string(fn));
        cell = cell.readCell(fpCell, formatspecCell);
        
        
        x = grain.X(:,1);
        y = grain.X(:,2);
        ID_A = contact.ID_A(:)+1;
        ID_B = contact.ID_B(:)+1;
        
        Lx = x(ID_B)-x(ID_A);
        Ly = y(ID_B)-y(ID_A);        
        [Lx, Ly] = CLP(cell, Lx, Ly);
                
%         Fx = contact.Force(:,1);                
%         if nargin==2 && isCorr
%             c(fn-T1+1) = corr(Fx, Ly);
%         else
%             covMat = cov(Fx, Ly);
%             c(fn-T1+1) = covMat(1,2);
%         end
        
        Fy = contact.Force(:,2);
        if nargin==2 && isCorr
            c(fn-T1+1) = corr(Fy, Lx);
        else
            covMat = cov(Fy, Lx);
            c(fn-T1+1) = covMat(1,2);
        end
        
    end

    if nargin==2 && isCorr
        obj.corFL = mean(c);
        obj.dcorFL = std(c)/sqrt(length(c));
    else
        obj.covFL = mean(c);
        obj.dcovFL = std(c)/sqrt(length(c));
    end
    
   
end