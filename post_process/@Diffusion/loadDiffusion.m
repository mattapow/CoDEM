function diffusion = loadDiffusion(DP, i, j)

if nargin==1
    dp = char(DP);
elseif nargin==2 % && size(DP,2)==1    
    dp = DP{i,1};
elseif nargin==3   
    dp = DP{i,j};    
end

if isempty(dp)
    warning('Empty entry in DP')
    diffusion = -1;
    return
end
fp = strcat(dp, '/diffusionLL.mat');

try
    in = load(fp);
    diffusion = in.diffusion;
    
    if diffusion.nPass==1
        diffusion.nPass=100;
        save(fp, 'diffusion')
    end
catch
    diffusion = Diffusion(dp);
%     diffusion = diffusion.getD();
%     diffusion = diffusion.getDVACF();
    save(fp, 'diffusion')
end

end