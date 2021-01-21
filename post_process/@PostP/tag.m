function blue = tag(dp,nGrains,t0,t1)
% returns grain ID of bottom half grains
%
% Accounts for grains passing through the boundary
% ID indexed from 0.

assert(t0==1, 'To use grain.XTrue, the inital time must be 1')

% Read grains at t0
fpGrain = strcat(dp, '/grain/grain_', string(t0));
grain = Grain();
grain = grain.readGrain(fpGrain, nGrains);

% tag bottom half blue
blue = grain.ID(grain.X(:,2)<0);
green =  grain.ID(grain.X(:,2)>=0);


% Read grains at t1
fpGrain = strcat(dp, '/grain/grain_', string(t1));
grain = Grain();
grain = grain.readGrain(fpGrain, nGrains);

% turn green to blue if passed up through boundary
up=rem(grain.XTrue(:,2),2)==1;
gAdd = intersect(grain.ID(up), green);
blue = [blue; gAdd];

% turn blue grains red if passed down through boundary
down=rem(grain.XTrue(:,2),2)==-1;
gSub = intersect(grain.ID(down), blue);
for i=length(blue):-1:1
    if ismember(blue(i), gSub)
       blue(i)=[]; 
    end
end


end