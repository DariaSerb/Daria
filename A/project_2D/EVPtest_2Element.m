clear all;
close all;
clc;

% 20/09/2016    
% Begin Timer
tic

% Read coordinates (in m) of Nodes end Elements 
Graphic_display  = 'yes';
Name_INP         = 'plate_400_200.inp';

% Nodes    = [1, 0.0, 0.0; 2, 0.400, 0.0; 3, 0.4, 0.2; 4, 0.0, 0.2];
% Elements = [1, 1, 2, 4; 2, 2, 3, 4];

%domain min(x) max(x) min(y) max(y) 
% domain = [0.0, 0.400, 0.0, 0.2];
domain   = [0, 400, 0, 200];

% Read of GMSH information
[Nodes,Elements]   = ReadINP(Name_INP,Graphic_display);
% [Nodes,Elements] = ReadGMSH(Name_GMSH,Graphic_display);

% figure(1);
%     hold on;
%     for j = 1:size(Elements,1)
%         tab(1,1:2) = [Nodes(Elements(j,2),2) Nodes(Elements(j,2),3)];
%         tab(2,1:2) = [Nodes(Elements(j,3),2) Nodes(Elements(j,3),3)];
%         tab(3,1:2) = [Nodes(Elements(j,4),2) Nodes(Elements(j,4),3)];
%         tab(4,1:2) = [Nodes(Elements(j,2),2) Nodes(Elements(j,2),3)];
%         plot(tab(:,1),tab(:,2),'r.-', 'LineWidth', 3) 
%     end
 
% % Read the material file in m.
% rho     = 7800;      % (kg/m3)
% e       = 1.1e-3;    % (m)
% E       = 2.1e11;    % (Pa)
% nu      = 0.3;       % (-)
% alpha   = 0;         % = 0 Plane Stress = 1 Plane Strain

% Read the material file in mm.
rho     = 7.8e-9;    % (t/mm3)
e       = 1.1;       % (mm)
E       = 210000;    % (MPa)
nu      = 0.3;       % (-)
alpha   = 0;         % = 0 Plane Stress = 1 Plane Strain

% Number of frequencies and modes of interest
ModeCnt = 10;

stiff     = zeros(2*size(Nodes,1),2*size(Nodes,1));
mass      = zeros(2*size(Nodes,1),2*size(Nodes,1));
mass_lump = zeros(2*size(Nodes,1),2*size(Nodes,1));
%--------------------------------------------------------------------------
for in=1:size(Elements,1)
    
    Te     = Elements(in,2:4);
    nedofT = 6;
    Te_dof = reshape([2*Te-1; 2*Te],1,nedofT);
    
    % check FEM elements
    
     [Ke_T3,detJ] = Ke_T3_analytique(Nodes(Elements(in,2),2),Nodes(Elements(in,3),2),...
                                    Nodes(Elements(in,4),2),Nodes(Elements(in,2),3),...
                                    Nodes(Elements(in,3),3),Nodes(Elements(in,4),3),...
                                    E,nu,alpha,e);    

     [Me_T3,Me_T3_lump] = Me_T3_analytique(Nodes(Elements(in,2),2),Nodes(Elements(in,3),2),...
                                    Nodes(Elements(in,4),2),Nodes(Elements(in,2),3),...
                                    Nodes(Elements(in,3),3),Nodes(Elements(in,4),3),...
                                    rho,e);
     
     stiff(Te_dof,Te_dof)     = stiff(Te_dof,Te_dof) + Ke_T3;
     mass(Te_dof,Te_dof)      = mass(Te_dof,Te_dof)  + Me_T3;     
     mass_lump(Te_dof,Te_dof) = mass_lump(Te_dof,Te_dof)  + Me_T3_lump;     
end

free_dofs      = [];
iter           =  1;
for in = 1:size(stiff,1)
    if all(stiff(in,:)  == 0)
        free_dofs(iter)  = in;
        iter             = iter + 1;
    end
end

dof_out        = setdiff(1:2*length(Nodes),free_dofs);
stiff          = stiff(dof_out,dof_out);
mass           = mass(dof_out,dof_out);
mass_lump      = mass_lump(dof_out,dof_out);

% Apply Boundary conditions
C    = DirBCmodif(Nodes,domain,1e-9);
C    = C';
 
%suppression of the free dofs of the matrices
increment = 0;
for j = 1:size(C,2)
    position = C(1,j)-increment;
    stiff(position,:) = [];
    stiff(:,position) = [];
    mass(position,:)  = [];
    mass(:,position)  = [];
    mass_lump(position,:) = [];
    mass_lump(:,position) = [];
                   
    increment = increment + 1;
end

% compute EigenValue problem
% [V,D] = eig(A,B) returns diagonal matrix D of generalized eigenvalues 
% and full matrix V whose columns are the corresponding right eigenvectors, 
% so that A*V = B*V*D.
% sorted natural angular frequencies [rad/s] 
[V,D]   = eig(stiff, mass_lump); 

lambda    = zeros(length(D),1);
freq      = zeros(length(D),1);

  for in = 1:length(D)
    lambda(in) = abs(D(in,in)); 
    freq(in)   = sqrt(abs(D(in,in)))/(2*pi);    
  end
  
lambda_sorted = sort(lambda);  
freq = sort(freq);
freq = freq(1:length(D));
lambda_sorted  = lambda_sorted(1:length(D));

mui  = diag(V' * mass * V);

eigvec_FEM = V';
% MaxMode = zeros(size(eigvec_FEM,1),1);  
%   for nn = 1:size(eigvec_FEM,1)
%    MaxMode(nn)  = max(eigvec_FEM(:,nn));
%    eigvec_FEM_norm(:,nn) = eigvec_FEM(:,nn)/MaxMode(nn);
%   end

% Restoring full matrix of degrees of freedom  
u               = zeros(2*length(Nodes),size(V,1));  
dof_out_BC_hole = setdiff(1:2*length(Nodes),[free_dofs C]);

% u(dof_out_BC_hole,:)   = V';  
u(dof_out_BC_hole,:)   = eigvec_FEM;  
% the values in nodes are equal to 0 (BC)
u(C,:) = u(C,:) * 0;
%-------------------------------------------------------------------------- 
utemp  = u(:,1);
utempx = utemp(1:2:end);
utempy = utemp(2:2:end);

 coef = zeros(length(Nodes),2);
 for in=1:length(Nodes)
    coef(in,1) = Nodes(in,2); 
    coef(in,2) = (ly/lx) * Nodes(in,2); 
    NodesNew(in,2) = Nodes(in,2) + coef(in,1) * utempx(in,1);
    NodesNew(in,3) = Nodes(in,3) + coef(in,2) * utempy(in,1);
 end

% coef = scalefactor(utempx,utempy);

% NodesNew(:,2) = Nodes(:,2) + coef(1) * utempx(:,1);
% NodesNew(:,3) = Nodes(:,3) + coef(2) * utempy(:,1);

figure;
    hold on;
    for j = 1:size(Elements,1)
        tab(1,1:2) = [NodesNew(Elements(j,2),2) NodesNew(Elements(j,2),3)];
        tab(2,1:2) = [NodesNew(Elements(j,3),2) NodesNew(Elements(j,3),3)];
        tab(3,1:2) = [NodesNew(Elements(j,4),2) NodesNew(Elements(j,4),3)];
        tab(4,1:2) = [NodesNew(Elements(j,2),2) NodesNew(Elements(j,2),3)];
        plot(tab(:,1),tab(:,2),'b.-') 
        axis([-50 450 -50 250]) 
    end


% End Timer
 computation_time = toc;

disp(['Total computation time: ',num2str(computation_time),' seconds']); 
disp('**************************************************'); 
