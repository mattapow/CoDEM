function plotMsd(obj, color, marker)
    %% Plot time Vs the Msd log-log
    if nargin==1
        color = 'k';
        marker = 'v';
    end

    hold on

    plot(obj.dt*obj.shearRate, obj.msd, marker, 'Color', color, 'MarkerSize', 4)

    xlab = xlabel('$$ t \dot{\gamma}$$');
    ylab = ylabel('$$\frac{\langle \Delta y^2 \rangle}{d^2}$$');
    set(xlab, 'Interpreter','latex');
    set(ylab, 'Interpreter','latex');
    set(ylab, 'Rotation', 0);
    set(ylab, 'HorizontalAlignment', 'right');

    ax = gca;
    ax.XScale = 'log';
    ax.YScale = 'log';

    hold off

end