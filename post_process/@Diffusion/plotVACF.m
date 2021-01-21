function plotVACF(obj, color)
    % PLOTVACF Plot the velocity autocorrelation function
    %
    % averaged over all grains of this simulation with optional
    % linespec input color

    assert(~(any(isnan(obj.lags)) || any(isnan(obj.acf)) || any(isnan(obj.dvY))), 'obj has nans\n')

    x = obj.lags;
    y = obj.acf/obj.dvY^2;

    if nargin == 1
        plot(x, y, 'k-')
    elseif nargin == 2
        plot(x, y, 'Color', color)
    end

    xlabel('t');
    ylabel('v acf')

end