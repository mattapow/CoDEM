function obj = setViscoGK(obj)
%% Read in grains at each timestep    
tic
isStress=1;
tsS = readTimeSteps(obj, isStress);         % stress
fprintf('Read took %.4f seconds.\n', toc)

% isStress=0;
% tsX = readTimeSteps(obj, isStress);         % position
% fprintf('Read took %.4f seconds.\n', toc)

% nG = PostP.nGrains1(obj);
nG = 1e4;
nT=length(tsS);

% P contains the autocorrelations
P = zeros(nT, nG);

for g = nG:-1:1
    
    for t = 1:nT
        tau(t) = tsS(t).grainPp.stress(g,2);
    end
    
    tau = tau-mean(tau);
    thisAcf = xcorr(tau);
    thisAcf = thisAcf(nT:end)';
    P(:, g) = thisAcf;
    
end

dt1 = PostP.getDt1(obj.dirPath);
obj.dtGK = (1:nT)'*dt1;

acf = mean(P, 2);
linear = linspace(1,0,nT)';
acf = acf .* linear;


T = obj.dvY^2;
V = obj.L(1) * obj.L(2) * 1;
acf = V / T * acf;

obj.viscoGK = acf;
obj.viscoLimGK = trapz(obj.dtGK, obj.viscoGK);

end