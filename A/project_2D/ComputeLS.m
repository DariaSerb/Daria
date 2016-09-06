function [Nodes] = ComputeLS(Nodes,Type_LS,Data_LS,Graphic_display)

% LS = InitializeLS(Nodes,Elements,Type_LS,Data_LS)
% Nodal values of the level set function

if strcmp(Type_LS, 'Circle')
    
    x_pos_center = Data_LS(1);
    y_pos_center = Data_LS(2);
    Radius       = Data_LS(3);
    tolerance    = Data_LS(4);
    
    r = sqrt((Nodes(:,2)-x_pos_center).^2 + (Nodes(:,3)-y_pos_center).^2);
    LS = r - Radius;
    LS(abs(LS) < tolerance) = 0;
    
    Nodes(:,4) = LS;
end

xmin = min(Nodes(:,2));
xmax = max(Nodes(:,2));
ymin = min(Nodes(:,3));
ymax = max(Nodes(:,3));

if strcmp(upper(Graphic_display),'YES')
  
    figure(1);
    hold on;
    if strcmp(Type_LS, 'Circle')
        Angle = [0:2*pi/100:2*pi];
        x     = Radius*cos(Angle) + x_pos_center;
        y     = Radius*sin(Angle) + y_pos_center;
        plot(x,y,'k.-');
        axis([xmin xmax ymin ymax]) 
    end
    hold off;
end



