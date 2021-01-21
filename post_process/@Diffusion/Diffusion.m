classdef Diffusion < PostP
    % DIFFUSION - Particle Self-Diffusion in one experiment
    %
    % Measure diffusion related properties of the simulation

    properties

        % Mean squared Displacement
        dt = NaN % MSD Mean time change dt
        maxT = NaN % MSD Maximum timescale to look at for computing msd
        msd = NaN % MSD Mean msddiffu(dt)

        D = NaN % diffusion coefficient fitted to msd=2Dt^v
        DError = NaN % diffusion coefficient error in fit

        v = NaN % anomalous diffusion power, if relevant
        vError = NaN % anomalous diffusion power error, if relevant

        rSquare = NaN % r-Squared value of linear fit for long term Diffusion
        top = 2 % take the top: [n/top:n] points to fit the normal diffusion coefficient

        % Velocity autocorrelation
        lags = NaN % ACF: Velocity autocorrelation delays
        acf = NaN % ACF: Non-normalised velocity autocorrelation values at delays
        tau = NaN % characteristic decay time of acf

        DVACF = NaN % diffusion coeficient from D = int(<v(t)v(0)>)dt

        Lyapunov = [NaN NaN] % Lyapunov exponent x,y direction
        dLyapunov = [NaN NaN] % Lyapunov exponents error (rmse)

        mad = [NaN NaN] % mean log absolute displacement (MAD) [x y]
        madT = NaN % MAD lags

    end

    properties (Dependent)
         % Diffusion length^2 scale =D/gm
        l = NaN
    end

    methods

        % Constructor
        function obj = Diffusion(dirPath)
            obj = obj@PostP(dirPath);
            obj.maxT = floor((obj.tEnd-obj.nPass)/2);
        end

        % Set functions
        function obj = set.maxT(obj, divide)
            % fprintf('Careful, setting maxT is not intuitive\n')
            obj.maxT = floor((obj.tEnd-obj.nPass)/divide);
        end

        function obj = set.top(obj, top)
            obj.top = top;
        end

        function obj = set.l(obj, ~)
            obj.l = obj.D / obj.shearRate;
        end

    end

    methods (Static)

    diffusion = loadDiffusion(DP, i, j)

    [diffusion, t_i]=loadLyapunov(dp,reCompute)

%         function obj = loadobj(s)
%             if isstruct(s)
%                 newObj = Diffusion();
%                 newObj.dt = s.dt;
%                 newObj.msd = s.msd;
%                 newObj.maxT = s.maxT;
%                 newObj.D = s.D;
%                 newObj.DError = s.DError;
%                 newObj.v = s.v;
%                 newObj.vError = s.vError;
%                 newObj.rSquare = s.rSquare;
%                 newObj.top = s.top;
%                 newObj.lags = s.lags;
%                 newObj.acf = s.acf;
%                 newObj.DVACF = s.DVACF;
%                 newObj.tau = s.tau;
%                 newObj.l = s.l;
%                 obj = newObj;
%             else
%                 obj = s;
%             end
%
% %             try
% %                 fp = strcat(obj.dirPath, '/PostPL.mat');
% %                 obj = obj.refresh(fp);
% %             catch
% %                 fp = strcat(obj.dirPath, '/PostP.mat');
% %                 obj = obj.refresh(fp);
% %             end
% %             fprintf('Refreshed from %s\n', fp)
%         end

    end


end
