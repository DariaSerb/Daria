clear all;
close all;
clc;
% Calculation eigenvalue problem in 2D using XFEM for plate without hole
% Init parameters
Graphic_display = 'yes';
Name_GMSH       = '2Dmesh.msh';
Name_Mat        = 'material_properties.txt';
Type_LS         = 'None';

pos_x_center    = 0;
pos_y_center    = 0;
radius          = 50;
domain          = [-100, 100, -100, 100];
tolerance       = 1e-4;
Data_LS         = [pos_x_center pos_y_center radius tolerance];

% Read of GMSH information
[Nodes,Elements] = ReadGMSH(Name_GMSH,Graphic_display);

% Read the material file
% [Properties] = ReadMaterial(Name_Mat);
rho     = 7.8e-9;    %(t/mm3)
e       = 1.1;       % (mm)
E       = 210000;    % (MPa)
nu      = 0.3;       % (-)
alpha   = 0;         % = 0 Plane Stress = 1 Plane Strain

ModeCnt = 10;

% Initialization and computing Level Set
[Nodes] = ComputeLS(Nodes,Type_LS,Data_LS,Graphic_display);

% Separation the elements in three sets
if strcmp(Type_LS, 'Circle')
    [Elements] = SeparateElements(Elements,Nodes,Graphic_display);

    % Separation the elements inside the void
    ind_In  = 1;
    for ii = 1:size(Elements,1)
        if Elements(ii,5) == -1;
            Elems_In(ind_In) = ii;
            ind_In = ind_In + 1;
        end
    end
else
    Elements(:,5) = zeros(size(Elements,1),1);
end

stiff      = zeros(2*length(Nodes),2*length(Nodes));
mass       = zeros(2*length(Nodes),2*length(Nodes));

% Assemble de stiffness and mass matrices

for in=1:size(Elements,1)
    Te     = Elements(in,2:4);
    nedofT = 6;
    Te_dof = reshape([2*Te-1; 2*Te],1,nedofT);
         
    % check FEM elements
    if Elements(in,5) == 0
        [Ke_T3,detJ] = Ke_T3_analytique(Nodes(Elements(in,2),2),Nodes(Elements(in,3),2),...
            Nodes(Elements(in,4),2),Nodes(Elements(in,2),3),...
            Nodes(Elements(in,3),3),Nodes(Elements(in,4),3),...
            E,nu,alpha,e);
        
        [Me_T3,Me_T3_lump] = Me_T3_analytique(Nodes(Elements(in,2),2),Nodes(Elements(in,3),2),...
            Nodes(Elements(in,4),2),Nodes(Elements(in,2),3),...
            Nodes(Elements(in,3),3),Nodes(Elements(in,4),3),...
            rho,e);
        
        stiff(Te_dof,Te_dof)      = stiff(Te_dof,Te_dof) + Ke_T3;
        mass(Te_dof,Te_dof)       = mass(Te_dof,Te_dof)  + Me_T3;
    end
    
    if Elements(in,5) == 1
        [Ke_XT3,detJ] = Ke_T3_analytique(Nodes(Elements(in,2),2),Nodes(Elements(in,3),2),...
            Nodes(Elements(in,4),2),Nodes(Elements(in,2),3),...
            Nodes(Elements(in,3),3),Nodes(Elements(in,4),3),...
            E,nu,alpha,e);
        
        Ke_XT3 = Elements(in,6) * Ke_XT3;
        
        [Me_XT3,Me_XT3_lump] = Me_T3_analytique(Nodes(Elements(in,2),2),Nodes(Elements(in,3),2),...
            Nodes(Elements(in,4),2),Nodes(Elements(in,2),3),...
            Nodes(Elements(in,3),3),Nodes(Elements(in,4),3),...
            rho,e);
        
        Me_XT3 = Elements(in,6) * Me_XT3;
        Me_XT3_lump = Elements(in,6) * Me_XT3_lump;
        
        stiff(Te_dof,Te_dof)      = stiff(Te_dof,Te_dof) + Ke_XT3;
        mass(Te_dof,Te_dof)       = mass(Te_dof,Te_dof)  + Me_XT3;
    end
end

save 'stiffmassFullMatrix.mat' stiff mass;

free_dofs      = [];
iter           =  1;
for in = 1:size(stiff,1)
    
    if all(stiff(in,:)  == 0)
        free_dofs(iter)  = in;
        iter             = iter + 1;
    end
end

% convectivities table of Nodes
if ~isempty(free_dofs)
    Free_nodes(:,1) = (free_dofs(2:2:end))/2;
dofs = setdiff(Nodes(:,1),Free_nodes(:,1));
Connectivity_table(:,1) = dofs;
Connectivity_table(:,2) = (1:length(dofs));
else
    Connectivity_table(:,1) = Nodes(:,1);
    Connectivity_table(:,2) = Nodes(:,2);
end

% suppression of the free dofs of the matrices
increment = 0;
for j = 1:size(free_dofs,2)
    position = free_dofs(1,j) - increment;
    stiff(position,:) = [];
    stiff(:,position) = [];
    mass(position,:)  = [];
    mass(:,position)  = [];
    increment = increment + 1;
end

% Apply Boundary conditions
C = DirBC(Nodes,domain,1e-6);
C = C';

% suppression of the free dofs of the matrices
increment = 0;
for j = 1:size(C,2)
    position = C(1,j)-increment;
    stiff(position,:) = [];
    stiff(:,position) = [];
    mass(position,:)  = [];
    mass(:,position)  = [];
    increment = increment + 1;
end

% compute EigenValue problem
% [V,D] = eig(A,B) returns diagonal matrix D of generalized eigenvalues 
% and full matrix V whose columns are the corresponding right eigenvectors, so that A*V = B*V*D.
% sorted natural angular frequencies [rad/s] 
[V,D] = eig(stiff,mass); 
freq    = zeros(length(D),1);
eigvec  = zeros(length(V),ModeCnt);
MaxMode = zeros(length(D),1);
% condnumber: condition number of V (matrix of eigenvectors) 
condnumber = cond(V);

  for in = 1:length(D)
    freq(in) = sqrt(D(in,in))/(2*pi);    
  end
% sorted natural angular frequencies [Hz]  
[freq_sorted, sorted_index] = sort(freq);
freq_sorted  = freq_sorted(1:ModeCnt);
sorted_index = sorted_index(1:ModeCnt);
  for in = 1:ModeCnt
    eigvec(:,in) = V(:,sorted_index(in));
  end
  eigvec_FEM = eigvec';
  
  for nn = 1:size(eigvec_FEM,2)
   MaxMode(nn)            = max(eigvec_FEM(:,nn));
   eigvec_norm_FEM(:,nn) = eigvec_FEM(:,nn)/MaxMode(nn);
  end
