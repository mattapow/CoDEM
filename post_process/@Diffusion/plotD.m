function plotD(obj)
    % Plot a line of the fitted msd vs t
    tf = ishold;

    if isnan(obj.D)
        obj = getD.obj;
    end

    if ~isnan(obj.v)
        error('Cannot handle anomalous diffusion\n')
    end

%     selec = round(numel(obj.dt)/obj.top):numel(obj.dt);
    selec = obj.maxT:numel(obj.dt);
    time = obj.dt(selec);
    msd = 2*obj.D*time;
    hold on
    line(time*obj.shearRate, msd, 'Color', 'b', 'LineWidth', 3)

    if ~tf
        hold off
    end

end
