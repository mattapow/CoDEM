function [acf, dt] = gAutoCorr(obj)
	% Autocorrelation of G = y v_x

	% Read in grains at each timestep    	
	tic
	isStress=1;
	ts = readTimeSteps(obj, isStress); 
	% ts = readTimeSteps(obj); 
	fprintf('Read took %.4f seconds.\n', toc)    

	nG = PostP.nGrains1(obj.dirPath);	
	nT=length(ts);
	if isStress; nT = nT - 1; end
	momentums = zeros(nT, nG);

	dt1 = PostP.getDt1(obj.dirPath);	
	dt = (0:nT-1)'*dt1;

	% fprintf('Computing mean grain momementum autocorrelation\n')
	fprintf('Computing mean grain stress autocorrelation\n')
	% fprintf('Computing mean grain y dv_x \n')

	for g = 1
		% fprintf('%d\n', g)
		thisG = zeros(nT,1);
		% y = zeros(nT, 1);
		% for t = 1:nT
		% 	y(t) = ts(t).grain.XTrue(g,2);
		% end

		for t = 1:nT
			% G(t) = ts(t).grain.XTrue(g, 2) * ts(t).grain.XTrue(g, 2);
			% G(t) = ts(t).grain.XTrue(g, 2) * (ts(t).grain.V(g, 1) - ts(t).grain.XTrue(g, 2)*obj.shearRate); % y dv            
			thisG(t) = ts(t).grainPp.stress(g, 3)/-ts(t).grainPp.stress(g, 4);
            if isnan(thisG(t))
                thisG(t) = 0;
            end
			
			% xx = dt(1:t)';
			% yy = y(1:t)';
			% G(t) = ts(t).grain.V(g, 2) * obj.shearRate * trapz(xx, yy);
		end

		% XTrue = [ts.grain.XTrue(g, 2)]
		% Grain = ts{:}.grain
		% XTrue = {ts.grain.XTrue(g, 2)}
		% XTrue = cell2mat( cellfun(@(x) x.grain.XTrue(g, 2), ts, 'UniformOutput', false));
		% V = cell2mat( cellfun( @(x) x.grain.V(g, 1), ts, 'UniformOutput', false) );
		% G = XTrue .* (V - XTrue*obj.shearRate);


		%% this grain's acf        
		thisacf = xcorr(thisG);
		thisacf = thisacf( nT:end );
% 	 	thisacf = thisacf/thisacf(1);
	 	momentums(:,g) = thisacf;
	 
	 	% Raw G
	 	momentums(:,g) = G;


	end

	acf = mean(momentums, 2);


	


end
