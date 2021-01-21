function rheology = loadRheology(DP, i,j)
% rheology = loadRheology(DP, i)
% rheology = loadRheology(DP, i, j)
% Given cell array of directory paths in DP, returns a rheology object
% located in the value of DP{i,j}
% 
% Could extend this to accomodate loading on different computers with
% different root paths

if nargin==2 % && any(size(DP)==1)
    dp = DP{i,1};
elseif nargin==3
    dp = DP{i,j};
elseif nargin==1
    dp = char(DP);
else
    error('Incorrect number of inputs %d', nargin)
end

if isempty(dp)
    fprintf('Empty entry in DP\n')
    rheology = -1;
    return
end

fp = strcat(char(dp), '/Rheology.mat');

try
    in = load(fp);
    rheology = in.rheology;
catch
    rheology = Rheology(dp);
    rheology = rheology.setStress;
    save(fp, 'rheology')
end

% warning('remove this')
% n = PostP.nFiles(rheology.dirPath);
% rheology.tEnd = n;
% save(fp, 'rheology')

end