% obj = GETMU(obj) Get the time averaged mu
%
% Find the ratio of shear to vertical stresses: sigma_xy / sigma yy
%


function mu = getMu(obj)

    mu = obj.stress(3) / - obj.stress(4);

end
