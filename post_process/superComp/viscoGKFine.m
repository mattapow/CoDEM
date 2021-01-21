function viscoGKFine()

%     dp = '/Users/mattmacaulay/Documents/DEM/DATA/N10000/gm0.1/C0/run';
%     dp = '/Users/mattmacaulay/Documents/DEM/DATA/N10000/gm0.1/C0/runFine1';
    dp = '/project/RDS-FEI-CDEM-RW/mmac6772';   

    fp = strcat(dp, '/Rheology.mat');
    fprintf('%s\n', fp)

    try 
        in = load(fp);
        rheology = in.rheology;
    catch
        rheology = Rheology(dp);
        rheology.nPass=60;    
        stress = rheology.stress; % includes saving
    end
    rheology.dirPath = dp;

    
    tic
    rheology = rheology.setViscoGKViscardy();
    fprintf('Crunch time took: %f.0 senconds\nSaving: ', toc)    
    save(fp, 'rheology')                        
    fprintf('done.\n')

end


