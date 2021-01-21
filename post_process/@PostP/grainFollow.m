function grainFollow(obj,gm,dt,ID)

T0=obj.nPass;
T1=obj.tEnd;

formatspec = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %f';
grain = Grain();
nGrains=1e4;

X=zeros(T1-T0+1,2);

for t=T0:T1

	fpGrain = strcat(obj.dirPath, '/grain/grain_', string(t));
	grain = grain.readGrain(fpGrain, nGrains, formatspec);

	X(t-T0+1,:)=grain.X(ID+1,:);
    
    XAff = X(ID+1,2)*gm*dt;
    X(t-T0+1,1)=X(t-T0+1,1)-XAff;

end

% move intital point to origin
X(:,2)=X(:,2)-X(1,2);


line(X(:,1),X(:,2))

end