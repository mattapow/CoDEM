function [diffusion, t_i]=loadLyapunov(dp,reCompute)
    
    diffusion = Diffusion.loadDiffusion(dp);
    rheology=Rheology.loadRheology(dp);
    
    d=1;
    rho_g=1;
    tStart=diffusion.nPass+1;
    tEnd=diffusion.tEnd;
    P=(rheology.stress(4)+rheology.stress(1))/2;
    t_i=d*sqrt(-rho_g/P);

    if reCompute
        fprintf('Recomputing Lyapunov exponents\n')
        [dX,t] = diffusion.getLyapunov(dp,tStart,tEnd);
        diffusion.mad=dX;
        diffusion.madT=t;

        fp = strcat(dp, '/diffusionLL.mat');
        fprintf('Saving output to %s\n', fp)
        save(fp, 'diffusion')
        fprintf('\tdone\n')
    end
end
