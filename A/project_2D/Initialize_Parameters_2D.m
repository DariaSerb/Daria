function P = Initialize_Parameters_2D()
% system parameters value

    P.rho     = 7800;      % (kg/m3)
    P.e       = 1.1e-3;    % (m)
    P.E       = 2.1e11;    % (Pa)
    P.nu      = 0.3;       % (-)
    P.alpha   = 0;         % = 0 Plane Stress = 1 Plane Strain    
    
    % Number of frequencies and modes of interest
    P.ModeCnt = 20;
    % Number of estimated frequencies and modes
    P.ModeEst = 3;
        
    % Basic characteristics of the LS-function
    pos_x_center  = 0.00;
    pos_y_center  = 0.00;
    radius        = 0.04;
    delta         = 0.08;
    radiusout     = radius + delta;
    
    P.pos_x_center  = pos_x_center;
    P.pos_y_center  = pos_y_center;
    P.radius        = radius;
    P.delta         = delta;
    P.radiusout     = radiusout;
    
    % the identification of plate's sizes
    flag = 0;
    
    if flag == 0
    % domain min(x) max(x) min(y) max(y) 
    a = 0.3;
    b = 0.2;
    domain = [-a, a, -b, b];
    else
    a = 0.4;
    b = 0.2;
    domain = [ 0, a,  0, b];    
    end
    lx = domain(2) - domain(1);
    ly = domain(4) - domain(3);
    
    P.domain = domain;
    P.lx     = lx;
    P.ly     = ly;
    
    tolerance = 1e-7;
    P.tolerance = tolerance;
    P.Data_LS   = [pos_x_center pos_y_center radius tolerance];
end
