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
    P.pos_x_center  = 0.00;
    P.pos_y_center  = 0.00;
    P.radius        = 0.04;
    P.delta         = 0.08;
    P.radiusout     = P.radius + P.delta;
           
    % the identification of plate's sizes
    flag = 1;
    
    if flag == 0
    % domain [min(x) max(x) min(y) max(y)]
    a = 0.3;
    b = 0.2;
    P.domain = [-a, a, -b, b];
    else
    a = 0.4;
    b = 0.2;
    P.domain = [ 0, a,  0, b];    
    end
    P.lx = P.domain(2) - P.domain(1);
    P.ly = P.domain(4) - P.domain(3);
           
    P.tolerance = 1e-7;
    P.Data_LS   = [P.pos_x_center P.pos_y_center P.radius P.tolerance];
 end
