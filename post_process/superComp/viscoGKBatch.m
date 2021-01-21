%viscoGKBatch.m
%A script to call viscoGKFine.m via the batch command

%%REMOTE MULTI CORE NAIVE
fprintf('Running Remotely Multi Core...\n')

tic %Let's time how long the entire process takes

job1 = batch(...
    @viscoGKFine,...        % Function Name to send to HPC
    0,...                   % Number of variables returned by fcn
    {},...                  % Inputs to use in function
    'Profile','Artemis',... % Name of the Cluster Profile to use
    'Pool',4,...            % Number of cpus to use in pool (max of 11)
    'CurrentFolder','/project/RDS-FEI-CDEM-RW/mmac6772');   
%     'CurrentFolder','/home/mmac6772');   % Change into Artemis home
    %  'CurrentFolder','.');   % Don't attempt to change into cwd

%This prevents us from using the terminal while we wait for 
%results from the batch job. 
wait(job1)  
           
%Print the time it took
fprintf('Total job submission and queue-wait = %f s.\n\n',toc) %Finish timing of the entire process.


