function obj = getMsdQ(obj)
    % Quickly calculate the MSD(t) (mean square vertical displacement)
    % Uses XTrue

    % Not finished,  quick because only a selection of random pairs
    % of snapshots in time
    obj = obj;
    error('Code not done\n')

    %% Presize Mean time change and Mean MSD
    % Sample 100 random intervals
    n = 100;
    obj.dt = zeros(n, 1); %#ok<*PROPLC>
    obj.msd = zeros(n, 1);

    for i = 1:n
        t1 = round((obj.tEnd-1)*rand()+1);
        t2 = round((obj.tEnd-1)*rand()+1);
        if t2 < t1; t11=t2; t2=t1; t1=t11; end

        %% More work to be done
    end

end