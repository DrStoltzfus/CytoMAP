function Invert_Plot_Axes(app)
    k = 0;
    v = 0;
    Constants.plt_menu_y0
    p = get(gcf,'position')
    reverse_x = uicontrol(gcf,'Style', 'push', 'String', ...
        'Invert X Axis','Position', [220 p(4)-30-Constants.plt_menu_y0 80 30],'CallBack', @Rotation_x);
    reverse_y = uicontrol(gcf,'Style', 'push', 'String', ...
        'Invert Y Axis','Position', [300 p(4)-30-Constants.plt_menu_y0 80 30],'CallBack', @Rotation_y);
    set(gca, 'CameraViewAngleMode', 'manual', 'CameraTargetMode', 'manual', ...
                 'CameraPositionMode', 'manual');   
    function Rotation_y(source,event)
        k = mod(k + 1, 2);
        Rotation();
    end
    
	function Rotation_x(source,event)
        v = mod(v + 1, 2);
        Rotation();        
    end

    function Rotation()
        if k==0 && v==0
            set(gca, 'xdir', 'normal');
            set(gca, 'ydir', 'normal');
            set(gca, 'CameraUpVector', [sin(0), cos(0), 0]);
        elseif k == 0 && v==1
            set(gca, 'xdir', 'reverse');
            set(gca, 'ydir', 'normal');
            set(gca, 'CameraUpVector', [sin(0), cos(0), 0]);
        elseif k == 1 && v==1
            set(gca, 'xdir', 'reverse');
            set(gca, 'ydir', 'reverse');
            set(gca, 'CameraUpVector', [sin(pi), cos(pi), 0]);
        elseif k == 1 && v==0
            set(gca, 'xdir', 'normal');
            set(gca, 'ydir', 'reverse');
            set(gca, 'CameraUpVector', [sin(pi), cos(pi), 0]);
        end
    end
end