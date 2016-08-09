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
    pos_x_center  = 0;
    pos_y_center  = 0;
    radius        = 0.04;
    
    P.pos_x_center  = pos_x_center;
    P.pos_y_center  = pos_y_center;
    P.radius        = radius;
    
    %domain min(x) max(x) min(y) max(y) 
    a = 0.3;
    b = 0.2;
    tolerance = 1e-7;
    
    P.domain          = [-a, a, -b, b];
    P.tolerance       = tolerance;
    P.Data_LS         = [pos_x_center pos_y_center radius tolerance];
end

