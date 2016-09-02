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
    P.ModeEst = 5;
    
    % dr is a scalar parameter which allows to follow the evolution of the structure     
    P.dr   = [0.000; 0.001; 0.005; 0.006; 0.007; 0.010; 0.015; 0.030; 0.040; 0.060];
    
    % dTau is a scalar parameter which allows to follow the evolution of the structure     
    P.dTau = [0.000; 0.001; 0.005; 0.006; 0.007; 0.010; 0.015; 0.030; 0.040; 0.060];
    
    % Basic characteristics of the LS-function and hole
    
    P.radius        = 0.04;
    P.delta         = 0.08;
    P.radiusout     = P.radius + P.delta;
    
    % the identification of plate's sizes and position of the hole
    flag = 2;
    
    if flag == 0
    P.pos_x_center  = 0.00;
    P.pos_y_center  = 0.00;    
    % domain [min(x) max(x) min(y) max(y)]
    a = 0.3;
    b = 0.2;
    P.domain = [-a, a, -b, b];
    end
    if flag == 1
    P.pos_x_center  = 0.20;
    P.pos_y_center  = 0.10;       
    a = 0.4;
    b = 0.2;
    P.domain = [ 0, a,  0, b];   
    end
    if flag == 2
    P.pos_x_center  = 0.30;
    P.pos_y_center  = 0.20;       
    a = 0.6;
    b = 0.4;
    P.domain = [ 0, a,  0, b];    
    end
    P.lx = P.domain(2) - P.domain(1);
    P.ly = P.domain(4) - P.domain(3);
           
    P.tolerance = 1e-7;
    P.Data_LS   = [P.pos_x_center P.pos_y_center P.radius P.tolerance];
 end
