function [psi, dPsi] = getLifeTime(obj)
dPsi=nan;
%% Method 1 - integral of acf
% Not really clean data, two trends
    psi =  (obj.D/(obj.dvY^2));        % Use total derivative and assume d(DVACF)=0
    dPsi = 2*obj.DVACF/obj.dvY^3 * obj.dvY_sd; % assuming d(DVACF)=0

%% Method 2 - mean of acf
% worked nicely except for I<.1 where the simulation save timestep was
% picked up
% acfN = obj.acf / trapz(obj.lags, obj.acf); % normalised acf
% t = obj.lags;
% 
% psi = trapz(t, t.*acfN);

%% Method 3 - time at acf=eps
% % worked nicely except for I<.1 where the simulation save timestep was
% % picked up

% eps = .4;
% acfN = obj.acf / obj.acf(1); % normalised acf
% t = obj.lags;
% % clf; plot(t,acfN,'k-'); drawnow; input('')
% 
% t0 = t(acfN<eps);
% psi = t0(1);

%% Method 4 - mean of first part of acf
% kinda noisy and roduces negative lifetimes

% % clf; plot(obj.lags*obj.shearRate, obj.acf); drawnow; error('')
% acfN = obj.acf / obj.acf(1); % acf as pdf
% t = obj.lags;
% 
% % remove all t(i),acfN(i)
% % for all i st |acfN(i)|<eps for all i>N
% eps=.02;
% V = abs(acfN)<eps;
% D = diff(diff(V)==1);
% B = find([true,D>0]);
% E = find([D<0,true])+1;
% [~,idx] = max(E-B);
% 
% acfN = acfN(1:B(idx)-1);
% t = t(1:B(idx)-1);
% if isempty(t); error('No lifetime found'); end
% % t(end)*obj.shearRate
% psi = trapz(t, t.*acfN); % since acf is continuous not disrete
% % psi = t(end);

%% Method 5 - subtract eps towards zero
% acf = obj.acf/obj.acf(1);
% t = obj.lags;
% 
% eps=.01;
% acf(acf>0) = max(acf(acf>0)-eps, 0);
% acf(acf<0) = min(acf(acf<0)+eps, 0);
% 
% psi = trapz(t, t.*acf); % since acf is continuous not disrete

%% Method 6 - mean over certain number of shear deformations

% acf = obj.acf;% / trapz(obj.lags, obj.acf); % normalised acf
% t = obj.lags;
% 
% eps = 2;
% tgm = t * obj.shearRate;
% t = t(tgm<eps); acf = acf(tgm<eps);
% acf = acf / trapz(t, acf); % normalised acf
% warning('Check acf normalisation (normalise full or truncated part?)')
% 
% psi = trapz(t, t.*acf); % since acf is continuous not disrete
% 
% if psi<0; warning('Negative lifetime computed'); end

%% Method 7 - integral of non-negative acf

% acfN = obj.acf / trapz(obj.lags, obj.acf); % normalised acf
% eps=-inf;
% acfN(acfN<eps) = 0;
% t = obj.lags;
% 
% psi = trapz(t, acfN);

%% Method 8 - use dv_x and dv_y to compute acf
% instantaneous=1;
% [acf, lags] = getVAcf(obj, instantaneous);
% acf=acf/acf(1);
% psi=trapz(lags,acf);


end