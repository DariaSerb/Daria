function  [Ketest_T3,Ketest_T3_dq,Ketest_T3_divq] = Elementary_matrix_Ke_Kdq(Xe,N,Ns,Nt,E,nedof,ngaus,wgp,nu,alpha,e)
% Elementary matrix for a 2D elasticity problem / plane stress /
% Xe:           nodal coordinates
% N, Ns, Nt:    shape functions and their derivatives
% nedof:        number of degrees of freedom in the element
% ngaus, wgp:   number of integration points and integration weights
% nu, e, alpha: material property

Ketest_T3      = zeros(nedof); 
Ketest_T3_dq   = zeros(nedof); 
Ketest_T3_divq = zeros(nedof);

B            = zeros(3, nedof); 
Bnew         = zeros(3, nedof); 
xe           = Xe(:, 2);
nen          = length(xe); 
M            = 2;
q            = zeros(ngaus, 1);
qm           = zeros(ngaus, M);
qs           = zeros(ngaus, 2 * M);
div_q        = zeros(ngaus, 1);

 for ig = 1:ngaus
    % plane stress (alpha = 0) plane strain (alpha = 1)
    C = [1 - alpha * nu, nu, 0; nu, 1 - alpha * nu, 0; 0, 0, (1 - nu - alpha * nu)/2];
    C = C * E/(1 - nu - alpha * nu)/(1 + nu);
    
    % shape functions in gauss points
    N_ig    = N(ig,:);
    Ns_ig   = Ns(ig,:);
    Nt_ig   = Nt(ig,:);
    
    % the transformation of gauss points's coordinates from local to reference system 
    Xst(1,:) = Xe(1,:) + (Xe(2,:) - Xe(1,:)) * N_ig(:,2) + (Xe(3,:) - Xe(1,:)) * N_ig(:,3);
   
    % q - function selection    
    [q(ig,:), qm(ig,:), qs(ig,:), div_q(ig,:)] = q_calc_funcMain(Xst);   
%   [qm(ig,:), qs(ig,:), div_q(ig,:)]          = q_calc_funcPolynom4(Xst);   
    
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
    
    nxnew  = qs(ig, 1) * nx + qs(ig, 3) * ny; 
    nynew  = qs(ig, 2) * nx + qs(ig, 4) * ny; 
        
    Bnew(1,1:2:end) = nxnew;
    Bnew(2,2:2:end) = nynew;
    Bnew(3,1:2:end) = nynew;
    Bnew(3,2:2:end) = nxnew;
    
    Ketest_T3       = Ketest_T3 + e * B' * C * B * dvolu; 
    Ketest_T3_dq    = Ketest_T3_dq + e * B' * C * Bnew * dvolu; 
    Ketest_T3_divq  = Ketest_T3_divq + e * B' * C * B * dvolu * div_q(ig,:); 
    
%   Ketest_T3_dq    = Ketest_T3_dq + 0.5*(e * Bnew' * C * Bnew * dvolu + (e * Bnew' * C * Bnew * dvolu)');     
%   Ketest_T3_divq  = Ketest_T3_divq + 0.5*(e * Bnew' * C * Bnew * dvolu * div_q(ig,:) + (e * Bnew' * C * Bnew * dvolu * div_q(ig,:))'); 
 end
end



