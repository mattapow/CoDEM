classdef Mixing
    %MIXING mixing properties
    % Everything we need for only mixing
    
    properties
    end
    
    methods
        function obj = Mixing(); end
    end
    
    methods (Static)
        function M=getIntensity(dp,t0,t1)
            % IC=getIntensity(dp,t0,t1) get the intensity of mixing
            % See Lacey 1954 (and Dunckwerts 1952)

            nGrains = PostP.nGrains1(dp);
            blue = PostP.tag(dp,nGrains,t0,t1);
            blue=blue+1;

            fpGrain = strcat(dp, '/grain/grain_', string(t1));
            formatspec = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %f';
            grain = Grain();
            grain = grain.readGrain(fpGrain, nGrains, formatspec);
            
            fpCell = strcat(dp,'/cell/cell_', string(t1));
            formatspec='%*f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
            cell=Cell();
            cell=cell.readCell(fpCell,formatspec);

            % normalised y position
            y_blue = grain.X(blue,2)/cell.L(2)*2;

            % get concentration of blue
            f=figure('Visible','off');
            h=histogram(y_blue,50,'Normalization','pdf','BinLimits',[-1 1],'Visible','off');
            C=h.Values;
            close(f)

            av=mean(C);
            n=numel(C);
            
            S2=sum((C - av).^2) / (n-1);
            S02=av*(1-av); % value of S for unmixed. Polydispersion may affect value.
            SR2=S02/n;
            
%             M=1-S/S0; % Dunckwerts
%             M=(SR2/S2)^.5; %Lacey's first proposal
            M=(S02-S2)/(S02-SR2); % Lacey mixing index
                
            assert(av-.5<.00001, 'Mean value is %f', av)
%             assert(IC<=1,'IC value is greater than 1: %f', IC)

            end
        
    end
end

