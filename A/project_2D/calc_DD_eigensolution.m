function [lambda_sorted,V,DerivLambda,Du] = calc_DD_eigensolution(Nodes,Elements,domain)
% 27/02/2017
global rho;
global e;
global E;
global nu;
global alpha;

global ModeCnt;
global ModeEst;

    stiff     = zeros(2*size(Nodes,1),2*size(Nodes,1));
    mass      = zeros(2*size(Nodes,1),2*size(Nodes,1));
    mass_lump = zeros(2*size(Nodes,1),2*size(Nodes,1));

    stiffinteg     = zeros(2*size(Nodes,1),2*size(Nodes,1));
    massinteg      = zeros(2*size(Nodes,1),2*size(Nodes,1));
    massinteg_lump = zeros(2*size(Nodes,1),2*size(Nodes,1));

    stiffdq        = zeros(2*size(Nodes,1),2*size(Nodes,1));
    stiffdivq      = zeros(2*size(Nodes,1),2*size(Nodes,1));
    massdq         = zeros(2*size(Nodes,1),2*size(Nodes,1));
    massdivq       = zeros(2*size(Nodes,1),2*size(Nodes,1));
   
% Integration weights and shape functions /initial configuration/
% Number of Gauss points
    npg         = 7;
% type of element / triangle /
    elem        = 1;
    [zgp, wgp]  = Quadrature(elem, npg);
    [N, Ns, Nt] = ShapeFunction_T3(zgp); 
   
   % to calculate [K],[M],[Kdq],[Mdq] and [Mdivq] for DD 
for in=1:size(Elements,1)
    Te = Elements(in,2:4);
    nedofT = 6;
    Te_dof = reshape([2*Te-1; 2*Te],1,nedofT);
    Xe     = Nodes(Te,2:3);
    
    [Ke_T3,~] = Ke_T3_analytique(Nodes(Elements(in,2),2),Nodes(Elements(in,3),2),...
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
     
    % numerical integration using Gauss-points
    [Ketest_T3, Ketest_T3_dq, Ketest_T3_divq] = Elementary_matrix_Ke_Kdq(Xe,N,Ns,Nt,E,nedofT,npg,wgp,nu,alpha,e);
    [Metest_T3, Metest_T3_lump, Metest_T3_dq, Metest_T3_divq] = Elementary_matrix_Me_Mdq(Xe,N,Ns,Nt,rho,nedofT,npg,wgp,e);
    % Assemble de stiffness and mass matrices
     stiffinteg(Te_dof,Te_dof)     = stiffinteg(Te_dof,Te_dof) + Ketest_T3;
     massinteg(Te_dof,Te_dof)      = massinteg(Te_dof,Te_dof) + Metest_T3;
     massinteg_lump(Te_dof,Te_dof) = massinteg_lump(Te_dof,Te_dof) + Metest_T3_lump;     
     % the definition of stiffdq, stiffdivq, massdq and massdivq
     stiffdq(Te_dof,Te_dof)         = stiffdq(Te_dof,Te_dof) + Ketest_T3_dq;
     stiffdivq(Te_dof,Te_dof)       = stiffdivq(Te_dof,Te_dof) + Ketest_T3_divq;
     massdq(Te_dof,Te_dof)          = massdq(Te_dof,Te_dof) + Metest_T3_dq;
     massdivq(Te_dof,Te_dof)        = massdivq(Te_dof,Te_dof) + Metest_T3_divq;
 end
%--------------------------------------------------------------------------
free_dofs      = [];
iter           =  1;
for in = 1:size(stiffinteg,1)
    
    if all(stiffinteg(in,:)  == 0)
        free_dofs(iter)  = in;
        iter             = iter + 1;
    end
end
dof_out        = setdiff(1:2*length(Nodes),free_dofs);

stiff      = stiff(dof_out,dof_out);
mass       = mass(dof_out,dof_out);
mass_lump  = mass_lump(dof_out,dof_out);

stiffinteg     = stiffinteg(dof_out,dof_out);
massinteg      = massinteg(dof_out,dof_out);
massinteg_lump = massinteg_lump(dof_out,dof_out);

stiffdq        = stiffdq(dof_out,dof_out);
stiffdivq      = stiffdivq(dof_out,dof_out);
massdq         = massdq(dof_out,dof_out);
massdivq       = massdivq(dof_out,dof_out);

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
    
    stiffinteg(position,:) = [];
    stiffinteg(:,position) = [];
    massinteg(position,:)  = [];
    massinteg(:,position)  = [];
    massinteg_lump(position,:) = [];
    massinteg_lump(:,position) = [];
    
    stiffdq(position,:) = [];
    stiffdq(:,position) = [];
    stiffdivq(position,:) = [];
    stiffdivq(:,position) = [];
    massdivq(position,:)  = [];
    massdivq(:,position)  = [];
              
    increment = increment + 1;
end

xkdq = 2*stiffdq - stiffdivq;
xmdq = massdivq;  

% compute EigenValue problem
% [V,D] = eig(A,B) returns diagonal matrix D of generalized eigenvalues 
% and full matrix V whose columns are the corresponding right eigenvectors, 
% so that A*V = B*V*D.
% sorted natural angular frequencies [rad/s] 
[V,D] = eig(stiff,mass); 
[Vinteg,Dinteg] = eig(stiffinteg,massinteg); 

lambda    = zeros(length(D),1);
freq      = zeros(length(D),1);

lambdainteg = zeros(length(Dinteg),1);
freqinteg   = zeros(length(Dinteg),1);

  for in = 1:length(Dinteg)
    lambda(in) = abs(D(in,in)); 
    freq(in)   = sqrt(abs(D(in,in)))/(2*pi);    
      
    lambdainteg(in) = abs(Dinteg(in,in)); 
    freqinteg(in)   = sqrt(abs(Dinteg(in,in)))/(2*pi);  
  end
lambda_sorted = sort(lambda);  
lambda_sorted  = lambda_sorted(1:length(D));
lambda_sorted = lambda_sorted';  
  
lambdainteg_sorted = sort(lambdainteg);  
lambdainteg_sorted  = lambdainteg_sorted(1:length(Dinteg));
lambdainteg_sorted = lambdainteg_sorted';

freq = sort(freq);
freq = freq(1:length(D));
freqinteg = sort(freqinteg);
freqinteg = freqinteg(1:length(Dinteg));
mui  = diag(V' * mass * V);
mui_integ  = diag(Vinteg' * massinteg * Vinteg);
alp_num     = zeros(ModeCnt,ModeEst);
alp_numtemp = zeros(ModeCnt,ModeEst);
 
for n = 1:ModeEst
    var = (xkdq + lambda_sorted(n) * xmdq) * V(:,n);
    numer = (V(:,n))' * var;
    DerivLambda(1,n) = - numer/mui(n);
    for m = 1:ModeCnt
        if n == m
            continue;
        end
        numer = (V(:,m))' * var; 
        alp_num(m,n) = numer/(mui(m) * (lambda_sorted(m) - lambda_sorted(n)));
    end
end
Du = V(:,1:ModeCnt)*alp_num;

for n = 1:ModeEst
    vartemp = (xkdq + lambdainteg_sorted(n) * xmdq) * Vinteg(:,n);
    numertemp = (Vinteg(:,n))' * vartemp;
    DerivLambdatemp(1,n) = - numertemp/mui_integ(n);
    for m = 1:ModeCnt
        if n == m
            continue;
        end
        numertemp = (Vinteg(:,m))' * vartemp; 
        alp_numtemp(m,n) = numertemp/(mui_integ(m) * (lambdainteg_sorted(m) - lambdainteg_sorted(n)));
    end
end
Dutemp = V(:,1:ModeCnt)*alp_numtemp;

 end 