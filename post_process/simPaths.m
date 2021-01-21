function DP = simPaths(ID, doRefresh)
    % Generate simulation Paths
    % DP = SIMPATHS(ID) 
    % DP = SIMPATHS(ID, doRefresh) 
    %
    % Output DP, a nxm cell array containing all the simulation paths.
    % Store thise paths as a matlab binary file in MATLAB/DEM PP/Library
    %
    % allSims: 
    %     0    - only converged simulations
    %     0.1  - all N10000/C*
    %     0.2  - all N10000/C*, using run2 for I=.3, C>=15, ...
    %     0.21 - Converged subset of 0.2
    %     0.22 - Set 0.21 but with max(C)=10. Only a few data points above
    %     
    %     1.0  - all N10000_1/gm*/C0 (no cohesion)
    %     1.1  - N10000_1/gm*/C0 (no cohesion) for .002 <= I <= .9    
    %     1.2  - N10000_1/gm*/C0/runFineExtra for .001<=I<=.9
    %     1.3  - N10000_1/gm*/C0/runFineExtra for .001<=I<=.9 with run2 for I>=.2     
    %     1.4  - N10000_1/gm*/CO/run3   I=.001, .02, .5
    %
    %     2    - N10000/gm*/C10/cJump5
    %     
    %   update the file holding simulation paths  
    %
    % NB: Some entries of DP may be empty, e.g. bad simulations using allSims=0.

    validID=[0 .1 .2 .21 .22 1 1.1 1.2 1.3 1.4 2];  
    if ~isempty(intersect(ID,validID))
        fp = strcat('/Users/mattmacaulay/Documents/MATLAB/DEM PP/Library/simPaths', char(string(ID)),'.mat');
    else
        error('Invalid input allSims')
    end
       
    if nargin==1 || (nargin==2 && ~doRefresh)        
        try
            in = load(fp);
            DP = in.DP;
        catch
            DP = getSimPaths(ID);
            save(fp, 'DP')       
            fprintf('Paths written to %s\n', fp)
        end
    elseif nargin==2 && doRefresh
        DP = getSimPaths(ID);
        save(fp, 'DP')       
        fprintf('Paths written to %s\n', fp)        
    else
        error('Incorrect number of inputs')
    end
    
end



function DP = getSimPaths(ID)
    rootPath = '/Users/mattmacaulay/Documents/DEM/DATA/N10000';
    
    if floor(ID)==1
        if ID==1.0      
            gm = { '0.00005', '0.0001', '0.0005', '0.001', '0.002', '0.005', '0.01', '0.02', '0.05', '0.1' , '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9',  '1', '1.2', '1.5', '2', '3'};
        elseif ID==1.1
            gm = {'0.002', '0.005', '0.01', '0.02', '0.05', '0.1' , '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'};
        elseif ID==1.2 || ID==1.3
            gm = {'0.001', '0.002', '0.005', '0.01', '0.02', '0.05', '0.1' , '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'};
        elseif ID==1.4
            gm = {'0.001', '0.02', '0.5'};
        end
        
        m=numel(gm);
        DP = cell(m,1);
        
        for i = 1:m
            if ID==1.0 || ID==1.1
                tail='/C0/run';
            elseif ID==1.2
                tail='/C0/runExtraFine';        
            elseif ID==1.3                
                if any((intersect(double(string(gm(i))), [.2 .3 .4 .5 .6 .7 .8 .9])))
                    tail='/C0/run2';
                else
                    tail='/C0/runExtraFine';        
                end
            elseif ID==1.4
                tail='/C0/run3';
            end
            dp =  strcat(rootPath, '_1/gm', gm{i}, tail);
            DP{i} = dp;
        end
        
    elseif floor(ID)==2.0
        gm = {'0.005', '0.01', '0.05', '0.1' , '0.2', '0.3'};
        coh='10';
        
        m=numel(gm);
        DP = cell(m,1);
        for i = 1:m
            dp =  strcat(rootPath, '/gm', gm{i}, '/C', coh, '/cJump5');
            DP{i,1} = dp;
        end
        
        
    elseif floor(ID)==0
        
        gm = {'0.005', '0.01', '0.05', '0.1' , '0.2', '0.3'};
        coh = {'0', '0.5', '1', '2', '5', '7', '10', '15', '20', '25', '30', '35', '40'};

        % NB. n,m are switched from normal
        n = numel(gm);
        m = numel(coh);
        DP = cell(n,m);
        
        if ID==.21||ID==.22; warning('Ensure /Library/misc/setInhom.m has been run.'); end

        for i = n:-1:1
        for j = m:-1:1

            if ID==.1
                dp =  strcat(rootPath, '/gm', gm{i}, '/C', coh{j}, '/run');
                DP{i,j} = dp;
                continue
            end
            
            if ID==.2 || ID==.21 || ID==.22
                if i==n && j>=8; tail='4';
                elseif i==2 && j==m; tail='2';  
                elseif i<=2 && j>=9; tail='2';
                elseif i==5 && (j==m-2 || j==m-3); tail='2';
                else; tail=''; 
                end 
                
                dp =  strcat(rootPath, '/gm', gm{i}, '/C', coh{j}, '/run',tail);
                DP{i,j} = dp;

                % filter converged
                if ID==.21 || ID==.22

                    r = Rheology.loadRheology(DP,i,j);

                    P1=-r.stress(1);
                    P2=-r.stress(4);
                    tau1=r.stress(2);
                    tau2=r.stress(3);

                    if r.inhomX^.5 > .05; DP{i,j} = ''; continue; end % velocity profile
                    if abs(P1-1)/P1 > .1; DP{i,j}=''; continue; end % mean pressure in x-dir
                    if abs(P2-1)/P2 > .1; DP{i,j}=''; continue; end % mean pressure in y-dir
                    if abs(tau2-tau1)/(tau2+tau1) > .05; DP{i,j}=''; continue; end % equal shear stresses
                    
%                     dStress = max(r.dStress./r.stress);
%                     if dStress > 1; DP{i,j}=''; continue; end % stress fluctuation
%                     if r.getSolidFracSTD() > .005; DP{i,j} = ''; continue; end % solid fraction standard deviation

                end
                continue
            end
            

            % ID==0
            try

                dp =  strcat(rootPath, '/gm', gm{i}, '/C', coh{j}, '/run');
                fpDif = strcat(dp, '/DiffusionLL.mat');
                in = load(fpDif);
                diffusion = in.diffusion;            

            catch

                fprintf('Initialising Diffusion object in %s\n', fpDif)
                diffusion = Diffusion(dp);

                fprintf('Simpaths fetching number of grains:...')
                diffusion = diffusion.getNG();
                fprintf('done\n')

                fprintf('Simpaths fetching DMSD:...\n')
                diffusion = diffusion.getD();

                fprintf('Simpaths fetching inhomX:...\n')
                diffusion = diffusion.getInhom;

                fprintf('Simpaths fetching DACF:...\n')
                diffusion = diffusion.getDVACF();

                save(fpDif, 'diffusion')

            end

            % Remove simulations with bad convergence
            if diffusion.inhomX < .25 && diffusion.rSquare > .98
                DP{i,j} = dp;
            end


        end        
        end
        
        if ID==.22
            DP(:,8:end) = [];
        end
        
    else
        error('Allsims not caught')        
    end
    fprintf('Paths acquired.\n')
end
