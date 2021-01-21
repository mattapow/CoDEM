function obj = setViscoEin2(obj)
% OBJ = SETVISCOEIN2(OBJ) Get cell viscosity
% Doesn't use timestep object (saves overhead)
%
% Compute the shear viscosity from Einstein-Helfand equation
% 
% n(t) = lim 1/(2kVTt) <[G(t)-G(0)]^2>
% 
%   where G(t) is the Helfand moment
%       G(t) = sum_grains p_x . y
%       where p_x is the unfolded grain momentum in the x-direction
%       and y is the unfolded position in the y-direction

    if isnan(obj.dvY)
        instantaneous=1;
        obj = getDv(obj, instantaneous);
    end
        
    
    % Compute the viscosity helper array
    % Viscous tensor component. See https://doi.org/10.1080/08927022.2017.1321760
    tnsrCmp = [2,1];
    visco = viscoEngine2(obj, tnsrCmp);    

    % Prefactors: Volume and energy terms
    mass = 1;
    E = .5 * mass * obj.dvY^2;
    V = obj.L(1) * obj.L(2) * 1;

    obj.visco =  visco/(2*E*V);
  
    % Use time difference between first two files to compute all lags
    dt1 = PostP.getDt1(obj.dirPath);
    N = length(obj.visco);
    obj.dt = (1:N)'*dt1;

end
