classdef Rheology < PostP
    %Rheology Post Processing

    properties

        % time averaged stress
        stress = [NaN NaN NaN NaN]
        dStress = [NaN NaN NaN NaN]

        % Einstein-Helfand viscosity (time dependent)
        % Helper mean square array, not actual viscosity
        visco = NaN
        dt = NaN        

        % Einstein-Helfand viscosity (lim t -> infty)
        viscoLim = NaN
        viscoLimRmse = NaN

        
        % Green-Kubo viscosity (time dependent)
        % Helper mean square array, not actual viscosity
        viscoGK = NaN
        dtGK = NaN

        viscoLimGK = NaN
        viscoLimGKRmse = NaN
        
        % Force-length correlations - mean
        covFL = NaN
        corFL = NaN
        
        % Force-length correlations - standard error
        dcovFL = NaN
        dcorFL = NaN
        
        % line between accel-force fluctuations
        % (aquired from fitting data from multiple class instances)
        kappa % slope
        f0 % intercept
        dkf0 % rmse error

    end

    methods

        function obj = Rheology(dirPath)
            obj = obj@PostP(dirPath);
        end

        function visco = getThisAppVisco(obj)
            % Get the apparent viscosity \tau / \dot \gamma
            visco = obj.stress(3) / obj.shearRate;
        end        

        function stress = get.stress(obj)

%             if any(isnan(obj.stress))
%                 obj = setStress(obj);
%                 rheology = obj;
%                 fp = strcat(obj.dirPath, '/Rheology.mat');
%                 save(fp, 'rheology')
%                 fprintf('Saved rheology stresses in %s\n', fp)
%             end

            stress = obj.stress;

        end


        function obj = setStress(obj)

            % Read in cell stress
            for filenum = obj.tEnd:-1:obj.nPass
              fp = strcat(obj.dirPath, '/cell/cell_', string(filenum));
              cell = Cell();
              formatspec = '%*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %*f %*f';
              cell = readCell(cell, fp, formatspec);
              obj.stress(filenum-obj.nPass+1, :) = cell.stress;
              
            end

            assert(~any(any(obj.stress==0)), 'Preallocation failed. Zeros stresses detected')
            obj.dStress = std(obj.stress, 0, 1);    
            obj.stress = mean(obj.stress, 1);
            

        end

        obj = setViscoEin(obj);
        obj = refresh(obj);


    end

    methods (Static)
        [eta, C, I, rmse, isWP] = getViscoTau(rheoOpt, isMu);
        [eta, C, I, rmse] = getViscoEin(isMu);
        [eta, C, I, rmse] = getViscoGK(plotIt)
        [eta, C, I] = getViscoApp(isMu);
        isWP = isWellPosed(f, I);
        rheology = loadRheology(DP, i, j);
        f = muVsI(DP,showPlot);
    end

end
