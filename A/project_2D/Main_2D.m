clear all;
close all;
clc;

% Calculation eigenvalue problem in 2D using XFEM for plate with hole
% Init parameters
Graphic_display = 'yes';
Name_Mat        = 'material_properties.txt';

% Read the material file
% [Properties] = ReadMaterial(Name_Mat);
rho     = 7.8e-9;    %(t/mm3)
e       = 1.1;       % (mm)
E       = 210000;    % (MPa)
nu      = 0.3;       % (-)
alpha   = 0;         % = 0 Plane Stress = 1 Plane Strain

% Number of frequencies and modes of interest
ModeCnt = 10;
Elements(:,1) = {1:length(el),1);
Elements(:,2:4) = el;
Elements(:,5) = eps;

% comparison the q - function
% [q1, qm1, qs1, div_q1]  = q_calc_funcTest(Nodes(:,2:3), Data_LS);
% [q2, qm2, qs2, div_q2]  = q_calc_func(Nodes(:,2:3), Data_LS);
% [q3, qm3, qs3, div_q3]  = q_calc_funcPolynom(Nodes(:,2:3), Data_LS);
[q4, qm4, qs4, div_q4]    = q_calc_funcPolynom3(vert);

stiff                = zeros(2*length(vert),2*length(vert));
stiffinteg           = zeros(2*length(vert),2*length(vert));
%mass                = zeros(2*length(vert),2*length(vert));
mass_lump            = zeros(2*length(vert),2*length(vert));
massinteg            = zeros(2*length(vert),2*length(vert));
massinteg_lump       = zeros(2*length(vert),2*length(vert));
%--------------------------------------------------------------------------
% Integration weights and shape functions
% Number of Gauss points
npg         = 7;
% type of element / triangle /
elem        = 1;
[zgp, wgp]  = Quadrature(elem, npg);
[N, Ns, Nt] = ShapeFunction_T3(zgp); 
%--------------------------------------------------------------------------

% Assemble de stiffness and mass matrices
% Calculation the stiffness and mass matrices using analytical integration
for in=1:size(Elements,1)
   
    Te     = Elements(in,2:4);
    nedofT = 6;
    Te_dof = reshape([2*Te-1; 2*Te],1,nedofT);
    Xe     = vert(Te,:);
            
    % check FEM elements
    if Elements(in,5) == 0
        [Ke_T3,detJ] = Ke_T3_analytique(vert(Elements(in,2),2),...
                       vert(Elements(in,3),2),vert(Elements(in,4),2),...
                       vert(Elements(in,2),3),vert(Elements(in,3),3),...
                       vert(Elements(in,4),3),E,nu,alpha,e);
        
        [Me_T3,Me_T3_lump] = Me_T3_analytique(vert(Elements(in,2),2),...
                       vert(Elements(in,3),2),vert(Elements(in,4),2),...
                       vert(Elements(in,2),3),vert(Elements(in,3),3),...
                       vert(Elements(in,4),3),rho,e);
%--------------------------------------------------------------------------            
% numerical integration using Gauss-points
    Ketest_T3                   = Elementary_matrix_Ke(Xe, N, Ns, Nt, E, nedofT, npg, wgp, nu, alpha, e);
    [Metest_T3, Metest_T3_lump] = Elementary_matrix_Me(Xe, N, Ns, Nt, rho, nedofT, npg, wgp, e);
%--------------------------------------------------------------------------           
% Calculation the stiffness and mass matrices using analytical integration      
      stiff(Te_dof,Te_dof)           = stiff(Te_dof,Te_dof) + Ke_T3;
      %mass(Te_dof,Te_dof)           = mass(Te_dof,Te_dof)  + Me_T3;
      mass_lump(Te_dof,Te_dof)       = mass_lump(Te_dof,Te_dof) + Me_T3_lump;
% Calculation the stiffness and mass matrices using numerical integration       
      stiffinteg(Te_dof,Te_dof)      = stiffinteg(Te_dof,Te_dof) + Ketest_T3;
      massinteg(Te_dof,Te_dof)       = massinteg(Te_dof,Te_dof) + Metest_T3;
      massinteg_lump(Te_dof,Te_dof)  = massinteg_lump(Te_dof,Te_dof) + Metest_T3_lump;
   end
    
    if Elements(in,5) == 1
        
        [Ke_XT3,detJ] = Ke_T3_analytique(vert(Elements(in,2),2),...
                       vert(Elements(in,3),2),vert(Elements(in,4),2),...
                       vert(Elements(in,2),3),vert(Elements(in,3),3),...
                       vert(Elements(in,4),3),E,nu,alpha,e);
        
        Ke_XT3    = Elements(in,6) * Ke_XT3;
                
        [Me_XT3,Me_XT3_lump] = Me_T3_analytique(vert(Elements(in,2),2),...
                       vert(Elements(in,3),2),vert(Elements(in,4),2),...
                       vert(Elements(in,2),3),vert(Elements(in,3),3),...
                       vert(Elements(in,4),3),rho,e);
        
        %Me_XT3 = Elements(in,6) * Me_XT3;
        Me_XT3_lump = Elements(in,6) * Me_XT3_lump;
        
%--------------------------------------------------------------------------            
% numerical integration using Gauss-points
    Ketest_T3                   = Elementary_matrix_Ke(Xe, N, Ns, Nt, E, nedofT, npg, wgp, nu, alpha, e);
    [Metest_T3, Metest_T3_lump] = Elementary_matrix_Me(Xe, N, Ns, Nt, rho, nedofT, npg, wgp, e);
    Ketest_XT3                  = Elements(in,6) * Ketest_T3;
    Metest_XT3                  = Elements(in,6) * Metest_T3;
    Metest_XT3_lump             = Elements(in,6) * Metest_T3_lump;
%--------------------------------------------------------------------------
% Calculation the stiffness and mass matrices using analytical integration        
      stiff(Te_dof,Te_dof)           = stiff(Te_dof,Te_dof) + Ke_XT3;
      %mass(Te_dof,Te_dof)           = mass(Te_dof,Te_dof)  + Me_XT3;
      mass_lump(Te_dof,Te_dof)       = mass_lump(Te_dof,Te_dof) + Me_XT3_lump;
% Calculation the stiffness and mass matrices using numerical integration        
      stiffinteg(Te_dof,Te_dof)      = stiffinteg(Te_dof,Te_dof) + Ketest_XT3;
      massinteg(Te_dof,Te_dof)       = massinteg(Te_dof,Te_dof) + Metest_XT3;
      massinteg_lump(Te_dof,Te_dof)  = massinteg_lump(Te_dof,Te_dof) + Metest_XT3_lump;
    end
end

save 'stiffmassFullMatrix.mat' stiff  mass_lump  stiffinteg  massinteg_lump;
  
%--------------------------------------------------------------------------
% q-function that determines the transition between configurations

figure(1)
trimesh(el,vert(:,2),vert(:,3),q(:,2))
title('The determination q(X)-function');
