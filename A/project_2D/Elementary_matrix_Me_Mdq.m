function  [Metest_T3,Metest_T3_lump,Metest_T3_dq,Metest_T3_divq] = Elementary_matrix_Me_Mdq(Xe,N,Ns,Nt,rho,nedof,ngaus,wgp,e)
% Elementary matrix for a 2D elasticity problem / plane stress /
% Xe:           nodal coordinates
% N, Ns, Nt:    shape functions and their derivatives
% nedof:        number of degrees of freedom in the element
% ngaus, wgp:   number of integration points and integration weights
% nu, e, alpha: material property

Metest_T3      = zeros(nedof); 
Metest_T3_lump = zeros(nedof);
Metest_T3_dq   = zeros(nedof); 
Metest_T3_divq = zeros(nedof);

M       = 2;
q       = zeros(ngaus, 1);
qm      = zeros(ngaus, M);
qs      = zeros(ngaus, 2 * M);
%qsm    = zeros(M, M);
div_q   = zeros(ngaus, 1);

% Nm = zeros(3, nedof); 
xe        = Xe(:,2);
nen       = length(xe); 

 for ig = 1:ngaus
    % plane stress (alpha = 0) plane strain (alpha = 1)
    N_ig    = N(ig,:);
    Ns_ig   = Ns(ig,:); 
    Nt_ig   = Nt(ig,:);
    
    Xst(1,:) = Xe(1,:) + (Xe(2,:) - Xe(1,:)) * N_ig(:,2) + (Xe(3,:) - Xe(1,:)) * N_ig(:,3);

   % q - function selection 
   [q(ig,:), qm(ig,:), qs(ig,:), div_q(ig,:)] = q_calc_funcMain(Xst);  
%  [qm(ig,:), qs(ig,:), div_q(ig,:)]          = q_calc_funcPolynom4(Xst);   

    qsm = (reshape(qs(ig,:),2, 2))'; 
    

    Jacob = [Ns_ig(1:nen); Nt_ig(1:nen)] * Xe;
    J     = det(Jacob); 
    dvolu = wgp(ig) * abs(J); 
          
    N_T_ig         = [reshape([1;0] * N_ig,1,nedof); reshape([0;1] * N_ig,1,nedof)];
    
    Metest_T3      = Metest_T3      + e * N_T_ig' * rho * N_T_ig * dvolu; 
    Metest_T3_dq   = Metest_T3_dq   + e * N_T_ig' * rho * qsm * N_T_ig * dvolu; 
    Metest_T3_divq = Metest_T3_divq + e * N_T_ig' * rho * div_q(ig,:)* N_T_ig * dvolu; 
 end
 
% the definition of lump mass matrix  
    for in=1:nedof
      Metest_T3_lump(in,in)=sum(Metest_T3(in,:));
    end
end



