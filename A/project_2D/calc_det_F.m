function [Ketest_T3_det_F,Ketest_T3_det_F1,Metest_T3_det_F,Metest_T3_det_F_lump,valid_F_inv] = calc_det_F(Xe,N,Ns,Nt,rho,alpha,nu,E,nedof,ngaus,wgp,e,dTau)
% 06/02/2017
% determination of the determinant of the gradient of transformation

M     = 2;
q     = zeros(ngaus, 1);
qm    = zeros(ngaus, M);
qs    = zeros(ngaus, 2 * M);
div_q = zeros(ngaus, 1);

B               = zeros(3, nedof); 
Ketest_T3_det_F = zeros(nedof);
Metest_T3_det_F = zeros(nedof);

Bnew1            = zeros(3, nedof); 
Ketest_T3_det_F1 = zeros(nedof);

xe    = Xe(:,2);
nen   = length(xe); 
valid = zeros(length(ngaus));

% det_F     = 0 + 1;
% I         = eye(2);
% F         = 0 + I;
% F_inverse = 0 + I;

  for ig = 1:ngaus
  
   det_F     = 1;
   I         = eye(2);
   F_inverse = I;   
      
  % plane stress (alpha = 0) plane strain (alpha = 1)
   C = [1 - alpha * nu, nu, 0; nu, 1 - alpha * nu, 0; 0, 0, (1 - nu - alpha * nu)/2];
   C = C * E/(1 - nu - alpha * nu)/(1 + nu);   
      
   N_ig    = N(ig,:);
   Ns_ig   = Ns(ig,:);
   Nt_ig   = Nt(ig,:);
   
   Xst(1,:) = Xe(1,:) + (Xe(2,:) - Xe(1,:)) * N_ig(:,2) + (Xe(3,:) - Xe(1,:)) * N_ig(:,3);

   % q - function selection 
   [qm(ig,:), qs(ig,:), div_q(ig,:)] = q_calc_funcPolynom4(Xst);   
   
   qsm = (reshape(qs(ig,:),2,2))'; 
   det_q = det(qsm);
      
   % the calculation of det_F  and F^(-1) 
   det_F     = det_F + (dTau * div_q(ig,:) + dTau * dTau * det_q);  
   F         = I + dTau * qsm;
   F_inv     = inv(F);
   F_inverse = F_inverse + (-dTau * qsm + dTau * dTau * qsm * qsm);
   valid_F_inv(ig) = norm(F_inv - F_inverse);
   
   Jacob = [Ns_ig(1:nen); Nt_ig(1:nen)] * Xe;
   J     = det(Jacob); 
   dvolu = wgp(ig) * abs(J); 
   
   res = Jacob\[Ns_ig; Nt_ig];
   nx  = res(1,:); 
   ny  = res(2,:); 
    
   B(1,1:2:end) = nx;
   B(2,2:2:end) = ny;
   B(3,1:2:end) = ny;
   B(3,2:2:end) = nx;
    
   nxnew  = F_inverse(1,1) * nx + F_inverse(2,1) * ny; 
   nynew  = F_inverse(1,2) * nx + F_inverse(2,2) * ny; 
   
   nxnew1  = F_inv(1,1) * nx + F_inv(2,1) * ny; 
   nynew1  = F_inv(1,2) * nx + F_inv(2,2) * ny; 
        
   Bnew(1,1:2:end) = nxnew;
   Bnew(2,2:2:end) = nynew;
   Bnew(3,1:2:end) = nynew;
   Bnew(3,2:2:end) = nxnew;
   
   Bnew1(1,1:2:end) = nxnew1;
   Bnew1(2,2:2:end) = nynew1;
   Bnew1(3,1:2:end) = nynew1;
   Bnew1(3,2:2:end) = nxnew1;
   
   N_T_ig    = [reshape([1;0] * N_ig,1,nedof); reshape([0;1] * N_ig,1,nedof)];
   
   Ketest_T3_det_F = Ketest_T3_det_F + e * Bnew' * C * Bnew * dvolu * det_F; 
   Metest_T3_det_F = Metest_T3_det_F + e * N_T_ig' * rho * N_T_ig * dvolu * det_F; 
   
   Ketest_T3_det_F1 = Ketest_T3_det_F1 + e * Bnew1' * C * Bnew1 * dvolu * det_F; 
  end
  
  Metest_T3_det_F_lump = zeros(nedof);
  % the definition of lump mass matrix  
    for in=1:nedof
      Metest_T3_det_F_lump(in,in) = sum(Metest_T3_det_F(in,:));
    end
end
