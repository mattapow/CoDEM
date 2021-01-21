function plotConstitutive(obj)
    % Aiming to plot log(shear Stress) Vs log(shear Rate)

    obj = getShearRateProfile(obj, 0);

    % Read in the time averaged profile
    profile = Profile();
    fp = strcat(obj.dirPath, '/profile/profile_average');
    formatspec = '%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %*f %*f %*f %f %*f %*f';
    profile = profile.readProfile(fp, formatspec);


    shearRate = obj.dVdy(~isnan(obj.dVdy));
    n=length(shearRate);
    mu = profile.stress(1:n,2) ./ profile.stress(1:n,5);

    edge=6;
    scatter(shearRate([1:edge end-edge+1:end]), mu([1:edge end-edge+1:end]), 'or')
    hold on
    scatter(shearRate(edge+1:end-edge), mu(edge+1:end-edge), 'ob')
    hold off
    xlab = xlabel('$$\dot{\gamma}$$');
    ylab = ylabel('$$\mu = \frac{\sigma_{xy}}{\sigma_{yy}}$$');
    set(xlab, 'Interpreter', 'Latex')
    set(ylab, 'Interpreter', 'Latex')
    leg = legend('Profiles near the boundary', 'Profiles in the Middle');
    leg.Location = 'northwest';
    title('Rheology based on time average stress-strain profile')

%           % this was for not time averaged rheology
%             for i = obj.tEnd-1:-1:1
%
%                 obj = getShearProfile(obj, i);
%                 dVdy = gradient(obj.vXSlice);
%
%                 profile = Profile();
%                 fp = strcat(obj.dirPath, '/profile/profile_', string(i));
%                 formatspec = '%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %*f %*f %*f %f %*f %*f';
%                 profile = profile.readProfile(fp, formatspec);
%
%                 n=length(dVdy);
%                 shearRate(1:n, i) = dVdy(~isnan(dVdy));
%                 shearStress(1:n, i) = profile.stress(1:n,2);
%
%             end
%             scatter(shearRate(:), shearStress(:), '.')

end
