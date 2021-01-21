function viscoEinFine()

    dp = '/project/CDEM/DATA/FINE/C0/gm01';   

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
    rheology = rheology.setViscoEin();
    fprintf('Crunch time took: %f.0 senconds\nSaving: ', toc)    
    save(fp, 'rheology')                        
    fprintf('done.\n')

end


