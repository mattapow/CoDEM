function tc = getTc(postP)

    T1=postP.nPass;
    T2=postP.tEnd;
    T2=T1+10;
    
    warning('T2=%d',T2)
    T = T2-T1;
    counts = nan(T,1);
    
    for t = 1
        fnStart=T1+t;
        meanCount = getDurations(postP, fnStart);
        counts(t)=meanCount;
    end
    
    tc=mean(counts, 'omitnan');
    

end


function meanCounts = getDurations(postP, fnStart)
    % meanCounts = getDurations(postP, fnStart)
    % Find the mean contact duration for contacts starting at fnStart>1
    % meanCounts is the number of files
    if fnStart==1; error('fnStart<=1'); end

    fp = strcat(postP.dirPath, '/contact/contact_', string(fnStart-1));
    c0 = Contact();
    c0 = readContact(c0, fp);
    
    fp = strcat(postP.dirPath, '/contact/contact_', string(fnStart));
    c1 = Contact();
    c1 = readContact(c1, fp);
    
    newContacts = setdiff([c1.ID_A,c1.ID_B], [c0.ID_A,c0.ID_B], 'rows');

    tcs = zeros(postP.tEnd-1,1);
    
    for t = fnStart:postP.tEnd-1
        fp = strcat(postP.dirPath, '/contact/contact_', string(t));
        c0 = readContact(c0, fp);
        fp = strcat(postP.dirPath, '/contact/contact_', string(t+1));
        c1 = readContact(c1, fp);
        
        coi = intersect(newContacts, [c0.ID_A,c0.ID_B], 'rows');% contacts of interest
        if isempty(coi); break; end
        
        % number of dead contacts of interest at this time
        tcs(t-fnStart+1) = numel(setdiff(coi, [c1.ID_A,c1.ID_B], 'rows'))/2; 
        
    end

    N=length(tcs);
    dt=PostP.getDt1(postP.dirPath);
    lags=(0:N-1)*dt;
    gm=postP.shearRate;
    plot(lags*gm,tcs)
    meanCounts=trapz(lags,tcs'.*lags);
end
        