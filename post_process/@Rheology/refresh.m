function obj = refresh(obj)
% obj = refresh(obj)
% Update instance values of rheology class
% Do manually below

fprintf('Refreshing rheology object\n')

sv=0;

if any(isnan(obj.dStress))
    obj = setStress(obj);
    sv=1;
end

if isnan(obj.kappa)
    f=figure();
    df_da_fit;
    close(f)
end

% nPass=201;


% if obj.nPass~=nPass
%     obj.nPass = nPass;
    
%     %1
%     obj = getDa(obj);
%     %2
%     obj = getDf_c(obj);
%     %3
%     obj = getSolidFrac(obj);
%     %4
%     obj = setNContacts(obj);
%     %5
%     obj = setStress(obj);
%     %6
%     isCorr=1;
%     obj = getCovFL(obj, isCorr);    
% end

if sv
    fp = strcat(obj.dirPath, '/Rheology.mat');
    rheology=obj;
    save(fp, 'rheology')
end

end