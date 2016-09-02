clear all;
close all;
clc;

% Calculation eigenvalue problem in 2D using XFEM for plate with hole
% Init parameters
Graphic_display = 'yes';
Name_Mat        = 'material_properties.txt';

% Init parameters in meters
P = Initialize_Parameters_2D();

% Read the material file
ModeCnt = P.ModeCnt;
ModeEst = P.ModeEst;
domain  = P.domain;
Data_LS = P.Data_LS; 
Elements(:,1) = {1:length(el),1);
Elements(:,2:4) = el;
Elements(:,5) = eps;

% Validation of initial parameters
out = Validation(P);
if out == 0
 return
end

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
