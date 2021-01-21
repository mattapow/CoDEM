function argout = plotViscoT(obj) 

     % For non-dimensional plotting only
    x = obj.dt*obj.shearRate;
    y = obj.visco;
    
    % y = y/y(10);
    % fprintf('y axis normalised by initial value\n')
   
    h1=plot(x, y, '.-');

    % if rem(length(varargin),2)==1
    %     error('Linespec must come as Name, Value pairs')
    % end
    % while length(varargin)>0
    %     if strcmp(varargin(1), 'Color')
    %         set(h1, 'Color', varargin{2});
    %         varargin(1:2) = [];
    %     elseif strcmp(varargin(1), 'Marker')
    %         set(h1, 'Marker', varargin{2});
    %         varargin(1:2) = [];
    %     elseif strcmp(varargin(1), 'MarkerFaceColor')
    %         set(h1, 'MarkerFaceColor', varargin{2})
    %         varargin(1:2) = [];
    %     % elseif strcmp(varargin(1), 'MarkerEdgeColor')
    %     %     set(h1, 'MarkerEdgeColor', varargin{2})
    %     %     varargin(1:2) = [];
    %     else 
    %         error('Linespec: only Color and Marker')
    %     end
    % end


    % Fit y=p1*x+p2
    % [f, gof] = fit(x, y, 'poly1');
    % viscoLim = f.p1;
    % h2=plot(f);    
    % set(h2, 'Color', h1.Color)        

    xlab=xlabel('$$t \dot{\gamma}$$');
    ylab=ylabel('$$\frac{1}{mv_0^2V} \langle [G(t) - G(0)]^2 \rangle (\frac{d}{m})$$');
%     ylab=ylabel('$$\frac{1}{mv_0^2V} \langle  \rangle (\frac{d}{m})$$');
    set(xlab, 'Interpreter', 'Latex')
    set(ylab, 'Interpreter', 'Latex')
    set(ylab, 'Rotation', 90)


    % uicontrol('Style', 'text',...
    %         'String', 'Slope is viscosity',...
    %         'Position', [20 10 100 20])    

    legend off
    if nargout==1
          argout = h1;
    end

end