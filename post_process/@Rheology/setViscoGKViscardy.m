function obj = setViscoGKViscardy(obj)
% Follows Viscardy 2007 equation (11), assuming sum of contact 
%% Read in grains at each timestep    
tic
tsC = readTimeStepsContacts(obj);           % contacts
fprintf('Read took %.1f seconds.\n', toc)

tic
isStress=1; isXTrue=0;
tsS = readTimeSteps(obj, isStress, isXTrue);         % stress
fprintf('Read took %.1f seconds.\n', toc)

tic
isStress=0; isXTrue=0;
tsX = readTimeSteps(obj, isStress, isXTrue);         % position
fprintf('Read took %.1f seconds.\n', toc)


warning('viscoConv = 0;')

% nG = PostP.nGrains1(obj);
nT=length(tsS);

fprintf('Computing momentum flux at time: %6.0f / %6.0f', nT, nT)
% clf
% hold on
% FX = nan(nT, 180000);
% LY = nan(nT, 180000);

for t = 1:nT
% parfor t = 1:nT

    %% a) Momentum flux: grain motions (convection)
%     m = tsX(t).grain.mass(:);
%     Vx = tsX(t).grain.V(:,1);
%     Vy = tsX(t).grain.V(:,2);
%     
%     % remove y drift
%     Vx = Vx - mean(Vx);
%     Vy = Vy - mean(Vy);

%     viscoConv = sum(m.*Vx.*Vy);   
% %     viscoConv = (sum(m.*Vx.^2.*Vy.^2))^.5;   
    viscoConv = 0;   

    %% b) Momentim flux: force (conduction)
    Fy = tsC(t).contact.Force(:,2);
    Fx = tsC(t).contact.Force(:,1);

    ID_A = tsC(t).contact.ID_A(:)+1;
    ID_B = tsC(t).contact.ID_B(:)+1;
    Xx = tsX(t).grain.X(ID_B,1) - tsX(t).grain.X(ID_A,1);
    Xy = tsX(t).grain.X(ID_B,2) - tsX(t).grain.X(ID_A,2);
    
    % Validate contact lengths
    [Xx,Xy] = CLP(tsX(t).cell,Xx,Xy);
%     assert(all(all(abs(Xy) <= 1.2)), 'Max |Xy| = %f', max(max(abs(Xx))));
    
%     warning('Check directions. Previously Fx, Xy')
    viscoCond = (sum(Fy.*Xx) + sum(Fx.*Xy))/2;
%     viscoCond = (sum(Fx.^2.*Xy.^2))^.5;

%     scatter(Fx, Xy, '.k')  
%     FX(t, 1:length(Fx)) = Fx;
%     LY(t, 1:length(Xy)) = Xy;    
    
    %%
    visco(t) = viscoConv + viscoCond;
%     fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%6.0f / %6.0f', t, nT)
    
end

% LY = LY(:);
% FX = FX(:);
% LY(isnan(LY)) = [];
% FX(isnan(FX)) = [];
% 
% hist3([FX LY], 'Nbins', [200 200], 'CDataMode','auto')
% var(LY)
% var(FX)
% corr(LY, FX)
% set(gca, 'FontSize', 16)
% xlabel('$$F^x / Pd^2$$', 'Interpreter', 'Latex')
% ylabel('$$l^y / d$$', 'Interpreter', 'Latex')
% drawnow
% error('end here plz')
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bdone\nComputing correlation\n')


% minAvg = 40; % 20 files per 1 shear
% maxLag = obj.tEnd - obj.nPass - (minAvg-1) - 1;
% acf = xcorr(visco, maxLag);
visco=visco-mean(visco);
plot(visco)
error('Stop here plz')
acf = xcorr(visco);
acf = acf(nT:end);
% acf = acf .* linspace(1, 0, length(acf));


dt1 = PostP.getDt1(obj.dirPath);
obj.dtGK = (1:nT)'*dt1;

if isnan(obj.dvY)
    instantaneous = 1;
    obj = getDv(obj, instantaneous);
end

if isnan(obj.L(1) * obj.L(2))
    filenums = (obj.nPass+1):obj.tEnd;
    obj = getDimensions(obj, filenums);
end

m = 1;
T = obj.dvY^2;
V = obj.L(1) * obj.L(2) * 1;
acf = acf / V / T / m;

obj.viscoGK = acf';
obj.viscoLimGK = trapz(obj.dtGK, obj.viscoGK);

fprintf('Max acf is: %f\n', max(obj.viscoGK))
fprintf('Timestep: %f\n', obj.dtGK(1))
fprintf('Viscosity Integral: %f\n', obj.viscoLimGK)

end