clear all;
close all;
clc;
    disp('******************************************************'); 
    disp('Integration of [K] and [M] by gauss points') 
    disp('using the modification of gauss points for cut element') 
    disp('******************************************************'); 
% Integration of [K] and [M] by gauss points using the modification of gauss points for cut element
% Introduction of Reference Element
% Begin Timer
    tic
    
% Init parameters in meters
Graphic_display  = 'yes';
Type_LS          = 'line';
P                = Initialize_Parameters_2D();
ReferenceElement = SetReferenceElement();

% Read of GMSH information
Nodes         = [1, 0.0214, 0.0202; 2, 0.0413, 0.0357; 3, 0.0574, 0.0101];
Elements      = [1, 1, 2, 3];
Nodes(:,4)    = [-0.0106; 0.0146; 0.0183]; 
Elements(:,5) = 1;
Elements(:,6) = 0.7878;

figure(1);
    hold on;
    for j = 1:size(Elements,1)
        tab(1,1:2) = [Nodes(Elements(j,2),2) Nodes(Elements(j,2),3)];
        tab(2,1:2) = [Nodes(Elements(j,3),2) Nodes(Elements(j,3),3)];
        tab(3,1:2) = [Nodes(Elements(j,4),2) Nodes(Elements(j,4),3)];
        tab(4,1:2) = [Nodes(Elements(j,2),2) Nodes(Elements(j,2),3)];
        plot(tab(:,1),tab(:,2),'r.-', 'LineWidth', 3) 
    end

% Read the material file
rho     = P.rho;      
e       = P.e;        
E       = P.E;        
nu      = P.nu;     
alpha   = P.alpha;   

stiff                = zeros(2*size(Nodes,1),2*size(Nodes,1));
mass                 = zeros(2*size(Nodes,1),2*size(Nodes,1));
mass_lump            = zeros(2*size(Nodes,1),2*size(Nodes,1));

stiffinteg           = zeros(2*size(Nodes,1),2*size(Nodes,1));
massinteg            = zeros(2*size(Nodes,1),2*size(Nodes,1));
massinteg_lump       = zeros(2*size(Nodes,1),2*size(Nodes,1));
%--------------------------------------------------------------------------
% Assemble de stiffness and mass matrices
% Calculation the stiffness and mass matrices using analytical integration
for in=1:size(Elements,1)
   % Integration weights and shape functions
     
      zgp = ReferenceElement.zgp; 
      wgp = ReferenceElement.wgp; 
    
    npg = ReferenceElement.ngaus;
    nen = ReferenceElement.nen; 
    Xe_ref = ReferenceElement.Xe_ref;
    
    Te     = Elements(in,2:4);
    nedofT = 2 * nen;
    Te_dof = reshape([2*Te-1; 2*Te],1,nedofT);
    Xe     = Nodes(Te,2:3);
    LSe    = Nodes(Te, 4);
  
   figure(2);
    hold on;
    plot(Xe_ref(:,1), Xe_ref(:,2), 'm.-','LineWidth', 3) 
    plot(zgp(:,1), zgp(:,2), 'k.', 'LineWidth', 4)
    hold off; 
 
   [zgp,wgp]   = ModifyQuadrature(LSe,ReferenceElement,zgp,wgp);
   [N, Ns, Nt] = ShapeFunction_T3(zgp);       
     
        
        [Ke_XT3,detJ]       = Ke_T3_analytique(Nodes(Elements(in,2),2),...
                       Nodes(Elements(in,3),2),Nodes(Elements(in,4),2),...
                       Nodes(Elements(in,2),3),Nodes(Elements(in,3),3),...
                       Nodes(Elements(in,4),3),E,nu,alpha,e);
        
        Ke_XT3    = Elements(in,6) * Ke_XT3;
                
        [Me_XT3,Me_XT3_lump] = Me_T3_analytique(Nodes(Elements(in,2),2),...
                       Nodes(Elements(in,3),2),Nodes(Elements(in,4),2),...
                       Nodes(Elements(in,2),3),Nodes(Elements(in,3),3),...
                       Nodes(Elements(in,4),3),rho,e);
        
        Me_XT3 = Elements(in,6) * Me_XT3;
        Me_XT3_lump = Elements(in,6) * Me_XT3_lump;
%--------------------------------------------------------------------------            
% numerical integration using Gauss-points
%    Ketest_XT3                    = Elementary_matrix_Ke(Xe, N, Ns, Nt, E, nedofT, length(wgp), wgp, nu, alpha, e);
%    [Metest_XT3, Metest_XT3_lump] = Elementary_matrix_Me(Xe, N, Ns, Nt, rho, nedofT, length(wgp), wgp, e);
   Ketest1_XT3                     = Elementary_matrix_Ke(Xe, N(1:7,:), Ns(1:7,:), Nt(1:7,:), E, nedofT, length(wgp), wgp(1,:), nu, alpha, e);
   [Metest1_XT3, Metest1_XT3_lump] = Elementary_matrix_Me(Xe, N(1:7,:), Ns(1:7,:), Nt(1:7,:), rho, nedofT, length(wgp), wgp(1,:), e);
   
   Ketest2_XT3                     = Elementary_matrix_Ke(Xe, N(8:14,:), Ns(8:14,:), Nt(8:14,:), E, nedofT, length(wgp), wgp(2,:), nu, alpha, e);
   [Metest2_XT3, Metest2_XT3_lump] = Elementary_matrix_Me(Xe, N(8:14,:), Ns(8:14,:), Nt(8:14,:), rho, nedofT, length(wgp), wgp(2,:), e);
   
   Ketest3_XT3                     = Elementary_matrix_Ke(Xe, N(15:21,:), Ns(15:21,:), Nt(15:21,:), E, nedofT, length(wgp), wgp(3,:), nu, alpha, e);
   [Metest3_XT3, Metest3_XT3_lump] = Elementary_matrix_Me(Xe, N(15:21,:), Ns(15:21,:), Nt(15:21,:), rho, nedofT, length(wgp), wgp(3,:), e);
   
%    Ketest_XT3 = Ketest1_XT3 + Ketest2_XT3 + Ketest3_XT3;
%    Metest_XT3 = Metest1_XT3 + Metest2_XT3 + Metest3_XT3;

   Ketest_XT3 = Ketest1_XT3 + Ketest2_XT3;
   Metest_XT3 = Metest1_XT3 + Metest2_XT3;
end

% End Timer
computation_time = toc;

disp(['Total computation time: ',num2str(computation_time),' seconds']); 
disp('**************************************************'); 
