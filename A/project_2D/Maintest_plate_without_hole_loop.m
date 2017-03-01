function Maintest_plate_without_hole_loop()
    % 27/02/2017    
        
% Read coordinates (in m) of Nodes end Elements 
 Graphic_display  = 'yes';
 Name_GMSH         = 'C:\Users\Dasha\Documents\MATLAB\myTest_2D\new\meshes_others\2Dmesh_0_4_0_2_N_258.msh';

%domain min(x) max(x) min(y) max(y) 
domain = [0.0, 0.400, 0.0, 0.2];
lx     = domain(2) - domain(1);
ly     = domain(4) - domain(3);

% Read of GMSH information
[Nodes,Elements]  = ReadGMSH(Name_GMSH,Graphic_display);
 
global rho;
global e;
global E;
global nu;
global alpha;

% Read the material file
rho     = 7800;      % (kg/m3)
e       = 1.1e-3;    % (m)
E       = 2.1e11;    % (Pa)
nu      = 0.3;       % (-)
alpha   = 0;         % = 0 Plane Stress = 1 Plane Strain

dTau  = linspace(0,0.050,50)';
N_tau = length(dTau);

global ModeCnt;
global ModeEst;

% Number of frequencies and modes of interest
ModeCnt = 20;
% Number of estimated frequencies and modes
ModeEst = 5;

LambdaEst = zeros(N_tau,ModeEst);

[lambda,V,DerivLambda,DerivShapes] = calc_DD_eigensolution(Nodes,Elements,domain);

for nt = 1:N_tau
valid_F = {};

stiff_det_F      = zeros(2*length(Nodes),2*length(Nodes));
stiff_det_F1     = zeros(2*length(Nodes),2*length(Nodes));
mass_det_F       = zeros(2*length(Nodes),2*length(Nodes));
mass_det_F_lump  = zeros(2*length(Nodes),2*length(Nodes));

stiffCurrent     = zeros(2*size(Nodes,1),2*size(Nodes,1));
massCurrent      = zeros(2*size(Nodes,1),2*size(Nodes,1));
massCurrent_lump = zeros(2*size(Nodes,1),2*size(Nodes,1));

stiffintegCurrent     = zeros(2*size(Nodes,1),2*size(Nodes,1));
massintegCurrent      = zeros(2*size(Nodes,1),2*size(Nodes,1));
massintegCurrent_lump = zeros(2*size(Nodes,1),2*size(Nodes,1));

%--------------------------------------------------------------------------
% Integration weights and shape functions
% Number of Gauss points
npg         = 7;
% type of element / triangle /
elem        = 1;
[zgp, wgp]  = Quadrature(elem, npg);
[N, Ns, Nt] = ShapeFunction_T3(zgp); 

% Assemble de stiffness and mass matrices

for in = 1:size(Elements,1)

    Te = Elements(in,2:4);
    nedofT = 6;
    Te_dof = reshape([2*Te-1; 2*Te],1,nedofT);
    Xe     = Nodes(Te,2:3);
         
    % numerical integration using Gauss-points
    [Ke_T3_det_F,Ke_T3_det_F1,Me_T3_det_F,Me_T3_det_F_lump,valid_F_inv] = calc_det_F(Xe,N,Ns,Nt,rho,alpha,nu,E,nedofT,npg,wgp,e,dTau(nt));
    % Calculation the [K_det_F] and [M_det_F] matrices using numerical integration   
    stiff_det_F(Te_dof,Te_dof) = stiff_det_F(Te_dof,Te_dof) + Ke_T3_det_F;
    stiff_det_F1(Te_dof,Te_dof) = stiff_det_F1(Te_dof,Te_dof) + Ke_T3_det_F1;
    mass_det_F(Te_dof,Te_dof)  = mass_det_F(Te_dof,Te_dof)  + Me_T3_det_F;   
    mass_det_F_lump(Te_dof,Te_dof)  = mass_det_F_lump(Te_dof,Te_dof) + Me_T3_det_F_lump;
    
    valid_F{in} =  valid_F_inv;  
end

NodesCurrent = change_configuration(Nodes,dTau(nt));

for in=1:size(Elements,1)
    Te = Elements(in,2:4);
    nedofT = 6;
    Te_dof = reshape([2*Te-1; 2*Te],1,nedofT);
    Xe     = NodesCurrent(Te,2:3);
    
    [Ke_T3,~] = Ke_T3_analytique(NodesCurrent(Elements(in,2),2),NodesCurrent(Elements(in,3),2),...
                                    NodesCurrent(Elements(in,4),2),NodesCurrent(Elements(in,2),3),...
                                    NodesCurrent(Elements(in,3),3),NodesCurrent(Elements(in,4),3),...
                                    E,nu,alpha,e);    

    [Me_T3,Me_T3_lump] = Me_T3_analytique(NodesCurrent(Elements(in,2),2),NodesCurrent(Elements(in,3),2),...
                                    NodesCurrent(Elements(in,4),2),NodesCurrent(Elements(in,2),3),...
                                    NodesCurrent(Elements(in,3),3),NodesCurrent(Elements(in,4),3),...
                                    rho,e);
    % Assemble de stiffness and mass matrices
    stiffCurrent(Te_dof,Te_dof)     = stiffCurrent(Te_dof,Te_dof) + Ke_T3;
    massCurrent(Te_dof,Te_dof)      = massCurrent(Te_dof,Te_dof)  + Me_T3;     
    massCurrent_lump(Te_dof,Te_dof) = massCurrent_lump(Te_dof,Te_dof)  + Me_T3_lump;  
    % numerical integration using Gauss-points
    KetestCurrent_T3 = Elementary_matrix_Ke(Xe,N,Ns,Nt,E,nedofT,npg,wgp,nu,alpha,e);
    [MetestCurrent_T3,MetestCurrent_T3_lump] = Elementary_matrix_Me(Xe,N,Ns,Nt,rho,nedofT,npg,wgp,e);
    % Assemble de stiffness and mass matrices
    % Calculation the stiffness and mass matrices using numerical integration       
    stiffintegCurrent(Te_dof,Te_dof)     = stiffintegCurrent(Te_dof,Te_dof) + KetestCurrent_T3;
    massintegCurrent(Te_dof,Te_dof)      = massintegCurrent(Te_dof,Te_dof) + MetestCurrent_T3;
    massintegCurrent_lump(Te_dof,Te_dof) = massintegCurrent_lump(Te_dof,Te_dof) + MetestCurrent_T3_lump;     
end

free_dofs      = [];
iter           =  1;
for in = 1:size(stiffCurrent,1)
    
    if all(stiffCurrent(in,:)  == 0)
        free_dofs(iter)  = in;
        iter             = iter + 1;
    end
end

dof_out          = setdiff(1:2*length(Nodes),free_dofs);
stiff_det_F      = stiff_det_F(dof_out,dof_out);
stiff_det_F1     = stiff_det_F1(dof_out,dof_out);
mass_det_F       = mass_det_F(dof_out,dof_out);
mass_det_F_lump  = mass_det_F_lump(dof_out,dof_out);

stiffCurrent     = stiffCurrent(dof_out,dof_out);
massCurrent      = massCurrent(dof_out,dof_out);
massCurrent_lump = massCurrent_lump(dof_out,dof_out);

stiffintegCurrent     = stiffintegCurrent(dof_out,dof_out);
massintegCurrent      = massintegCurrent(dof_out,dof_out);
massintegCurrent_lump = massintegCurrent_lump(dof_out,dof_out);

% Apply Boundary conditions
C    = DirBCmodif(Nodes,domain,1e-9);
C    = C';

% suppression of the free dofs of the matrices
increment = 0;
for j = 1:size(C,2)
    position = C(1,j)-increment;
     
    stiff_det_F(position,:)  = [];
    stiff_det_F(:,position)  = [];
    stiff_det_F1(position,:) = [];
    stiff_det_F1(:,position) = [];
    mass_det_F(position,:)   = [];
    mass_det_F(:,position)   = [];
    mass_det_F_lump(position,:) = [];
    mass_det_F_lump(:,position) = [];
    
    stiffCurrent(position,:) = [];
    stiffCurrent(:,position) = [];
    massCurrent(position,:)  = [];
    massCurrent(:,position)  = [];
    massCurrent_lump(position,:) = [];
    massCurrent_lump(:,position) = [];
    
    stiffintegCurrent(position,:) = [];
    stiffintegCurrent(:,position) = [];
    massintegCurrent(position,:) = [];
    massintegCurrent(:,position) = [];
    massintegCurrent_lump(position,:) = [];
    massintegCurrent_lump(:,position) = [];
       
    increment = increment + 1;
end

% to predict the estimated eigensolution
% to evaluate of satisfaction of eigensolution estimated on initial configuration

     for n = 1:ModeEst
        LambdaEst(nt,n) = dTau(nt)*DerivLambda(1,n) + lambda(1,n);   
        uesttemp(:,n) = V(:,n) + dTau(nt)*DerivShapes(:,n);

        value  = stiff_det_F1 - LambdaEst(nt,n) * mass_det_F; 
        ratio(:,n) = value * uesttemp(:,n);
        val1 = stiff_det_F1 * uesttemp(:,n); 
        val2 = LambdaEst(nt,n) * mass_det_F * uesttemp(:,n);
        val  = val1 - val2;
% to calculate the norm of the matrices [Kdet_F1]*Uest    
        [Knorm(nt,:),icK(nt,:)] = max(abs(val1));
        [Mnorm(nt,:),icM(nt,:)] = max(abs(val2));
% to calculate the norm of [Kdet_F]*Uest-Lambda*[Mdet_F]*Uest 
        [val_norm(nt,:),icval(nt,:)] = max(abs(val));
     end
     freqEst(nt,:) = sqrt(LambdaEst(nt,:))/(2*pi);
     uest(:,:,nt)  = uesttemp;
     max_val_ratio(nt,:) = val_norm(nt,:)/Knorm(nt,:);
end
end
