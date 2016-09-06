function [freq_sorted,Vanal,muianal] = eig_val(Nodes,Elements)

P                = Initialize_Parameters_2D();
ReferenceElement = SetReferenceElement();

rho     = P.rho;      
e       = P.e;        
E       = P.E;        
nu      = P.nu;     
alpha   = P.alpha;   
ModeCnt = P.ModeCnt;
domain  = P.domain;

stiff                = zeros(2*length(Nodes),2*length(Nodes));
mass                 = zeros(2*length(Nodes),2*length(Nodes));
mass_lump            = zeros(2*length(Nodes),2*length(Nodes));

stiffanal            = zeros(2*length(Nodes),2*length(Nodes));
massanal             = zeros(2*length(Nodes),2*length(Nodes));
massanal_lump        = zeros(2*length(Nodes),2*length(Nodes));
%--------------------------------------------------------------------------
% Assemble de stiffness and mass matrices
% Calculation the stiffness and mass matrices using analytical integration
for in=1:size(Elements,1)
   % Integration weights and shape functions
   npg = ReferenceElement.ngaus; 
   zgp = ReferenceElement.zgp; 
   wgp = ReferenceElement.wgp; 
   nen = ReferenceElement.nen; 
   N   = ReferenceElement.N; 
   Ns  = ReferenceElement.Ns; 
   Nt  = ReferenceElement.Nt; 
        
    Te     = Elements(in,2:4);
    nedofT = 2 * nen;
    Te_dof = reshape([2*Te-1; 2*Te],1,nedofT);
    Xe     = Nodes(Te,2:3);
    LSe    = Nodes(Te, 4);
            
% check FEM elements
   if Elements(in,5) == 0
% analytical integration      
     [Ke_T3,detJ]       = Ke_T3_analytique(Nodes(Elements(in,2),2),Nodes(Elements(in,3),2),...
                                    Nodes(Elements(in,4),2),Nodes(Elements(in,2),3),...
                                    Nodes(Elements(in,3),3),Nodes(Elements(in,4),3),...
                                    E,nu,alpha,e);    

     [Me_T3,Me_T3_lump] = Me_T3_analytique(Nodes(Elements(in,2),2),Nodes(Elements(in,3),2),...
                                    Nodes(Elements(in,4),2),Nodes(Elements(in,2),3),...
                                    Nodes(Elements(in,3),3),Nodes(Elements(in,4),3),...
                                    rho,e);      
% numerical integration using Gauss-points
 Ketest_T3 = Elementary_matrix_Ke(Xe, N, Ns, Nt, E, nedofT, npg, wgp, nu, alpha, e);
 [Metest_T3, Metest_T3_lump] = Elementary_matrix_Me(Xe, N, Ns, Nt, rho, nedofT, npg, wgp, e);
%--------------------------------------------------------------------------           
% Calculation the stiffness and mass matrices using analytical integration  
     stiffanal(Te_dof,Te_dof)     = stiffanal(Te_dof,Te_dof) + Ke_T3;
     massanal(Te_dof,Te_dof)      = massanal(Te_dof,Te_dof)  + Me_T3;     
     massanal_lump(Te_dof,Te_dof) = massanal_lump(Te_dof,Te_dof)  + Me_T3_lump;  
% Calculation the stiffness and mass matrices using numerical integration 
      stiff(Te_dof,Te_dof)      = stiff(Te_dof,Te_dof) + Ketest_T3;
      mass(Te_dof,Te_dof)       = mass(Te_dof,Te_dof) + Metest_T3;
      mass_lump(Te_dof,Te_dof)  = mass_lump(Te_dof,Te_dof) + Metest_T3_lump;
  end
    
    if Elements(in,5) == 1
     [Ke_XT3,detJ] = Ke_T3_analytique(Nodes(Elements(in,2),2),...
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
       
% the modification of gauss points for cut elements           
   [zgp,wgp]   = ModifyQuadrature(LSe,ReferenceElement,zgp,wgp);
   [N, Ns, Nt] = ShapeFunction_T3(zgp);            
%--------------------------------------------------------------------------   
% numerical integration using Gauss-points
  Ketest_XT3 = Elementary_matrix_Ke(Xe, N, Ns, Nt, E, nedofT, length(wgp), wgp, nu, alpha, e);
  [Metest_XT3, Metest_XT3_lump] = Elementary_matrix_Me(Xe, N, Ns, Nt, rho, nedofT, length(wgp), wgp, e);
% --------------------------------------------------------------------------
% Calculation the stiffness and mass matrices using analytical integration  
     stiffanal(Te_dof,Te_dof)     = stiffanal(Te_dof,Te_dof) + Ke_XT3;
     massanal(Te_dof,Te_dof)      = massanal(Te_dof,Te_dof)  + Me_XT3;     
     massanal_lump(Te_dof,Te_dof) = massanal_lump(Te_dof,Te_dof)  + Me_XT3_lump; 
% Calculation the stiffness and mass matrices using numerical integration 
      stiff(Te_dof,Te_dof)      = stiff(Te_dof,Te_dof) + Ketest_XT3;
      mass(Te_dof,Te_dof)       = mass(Te_dof,Te_dof) + Metest_XT3;
      mass_lump(Te_dof,Te_dof)  = mass_lump(Te_dof,Te_dof) + Metest_XT3_lump;
   end
end
%--------------------------------------------------------------------------
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
    mass_lump(position,:) = [];
    mass_lump(:,position) = [];
    
    stiffanal(position,:) = [];
    stiffanal(:,position) = [];
    massanal(position,:)  = [];
    massanal(:,position)  = [];
    massanal_lump(position,:) = [];
    massanal_lump(:,position) = [];
    
    increment = increment + 1;
end

% Apply Boundary conditions
C    = DirBCmodif(Nodes,domain,1e-9);
C    = C';

% suppression of the free dofs of the matrices
increment = 0;
for j = 1:size(C,2)
    position = C(1,j)-increment;
        
    stiff(position,:) = [];
    stiff(:,position) = [];
    mass(position,:)  = [];
    mass(:,position)  = [];
    mass_lump(position,:) = [];
    mass_lump(:,position) = [];
    
    stiffanal(position,:) = [];
    stiffanal(:,position) = [];
    massanal(position,:)  = [];
    massanal(:,position)  = [];
    massanal_lump(position,:) = [];
    massanal_lump(:,position) = [];
                
    increment = increment + 1;
end
%--------------------------------------------------------------------------
% compute EigenValue problem
% [V,D] = eig(A,B) returns diagonal matrix D of generalized eigenvalues 
% and full matrix V whose columns are the corresponding right eigenvectors, 
% so that A*V = B*V*D.
% sorted natural angular frequencies [rad/s] 
% [V,D]         = eig(stiff,mass_lump); 
  [V,D]         = eig(stiff,mass); 
  [Vanal,Danal] = eig(stiffanal,massanal); 

lambda  = zeros(length(D),1);
freq    = zeros(length(D),1);
eigvec  = zeros(length(V),ModeCnt);

% condnumber: condition number of V (matrix of eigenvectors) 
condnumber = cond(V);

  for in = 1:length(D)
    lambda(in)      = D(in,in); 
    freq(in)        = sqrt(D(in,in))/(2*pi); 
  end
% sorted natural angular frequencies [Hz]  
[freq_sorted, sorted_index]          = sort(freq);
[lambda_sorted, sorted_lambda_index] = sort(lambda);

freq_sorted                          = freq_sorted(1:ModeCnt);
lambda_sorted                        = lambda_sorted(1:ModeCnt);

sorted_index        = sorted_index(1:ModeCnt);
sorted_lambda_index = sorted_lambda_index(1:ModeCnt);
  for in = 1:ModeCnt
    eigvec(:,in)  = V(:,sorted_index(in));
  end
eigvec_XFEM = eigvec';

mui           = diag(V' * mass * V);
muianal       = diag(Vanal' * massanal * Vanal);

% eigvector normalization 
MaxMode = zeros(size(eigvec_XFEM,1),1);  
  for nn = 1:size(eigvec_XFEM,1)
   MaxMode(nn)            = max(eigvec_XFEM(nn,:));
   eigvec_norm_XFEM(nn,:) = eigvec_XFEM(nn,:)/MaxMode(nn);
  end

% Restoring full matrix of degrees of freedom  
u               = zeros(2*length(Nodes),size(V,1));  
dof_out_BC_hole = setdiff(1:2*length(Nodes),[free_dofs C]);

u(dof_out_BC_hole,:) = Vanal';  
% the values in nodes are equal to 0 (BC)
u(C,:) = u(C,:)*0;
%-------------------------------------------------------------------------- 
utemp = u(:,1);

Graphic_display  = 'NO';

   if strcmp(upper(Graphic_display),'YES')
    % plot eigen vector
    figure(2)
    trimesh(Elements(:,2:4), Nodes(:,2), Nodes(:,3), utemp(1:2:end),'FaceColor','Interp','EdgeColor','none')
    axis([-inf inf -inf inf])
    title('The determination eigen vector');
   end
end
