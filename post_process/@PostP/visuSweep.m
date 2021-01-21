clf
hold on
ax=gca;
trackID = 1;
contacts=0;
dp = '~/Documents/DEM/DATA/N10000/gm0.1/C0/run';
nGrains=1e4;

% plot all the grains at mid-time
whiteGrains=0;
PostP.visuOne(ax, dp, 100, [], contacts, nGrains, whiteGrains)

% plot one grain over time
whiteGrains=1;
for t = 1:201
    PostP.visuOne(ax, dp, t, trackID, contacts, nGrains, whiteGrains)
end

drawnow
