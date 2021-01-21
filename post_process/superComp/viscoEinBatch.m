%viscoEinBatch.m
%A script to call viscoEinFine.m via the batch command

%%REMOTE MULTI CORE NAIVE
fprintf('Running viscoEinFine remotely...\n')

tic %Let's time how long the entire process takes

job1 = batch(...
    @viscoEinFine,...       % Function Name to send to HPC
    0,...                   % Number of variables returned by fcn
    {},...                  % Inputs to use in function
    'Profile','Artemis',... % Name of the Cluster Profile to use
    'Pool',2,...            % Number of cpus to use in pool (max of 11) NB: if mem/pool > 20GB => high memory job
    'CurrentFolder','/project/CDEM/DATA/FINE/C0/gm01');

%This prevents us from using the terminal while we wait for 
%results from the batch job. 
wait(job1)

           
%Print the time it took
fprintf('Total job submission and queue-wait = %f s.\n\n',toc) %Finish timing of the entire process.
% delete(job1)

