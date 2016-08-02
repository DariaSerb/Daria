function  Ketest_T3 = Elementary_matrix_Ke(Xe, N, Ns, Nt, E, nedof, ngaus, wgp, nu, alpha, e)
% [Ke_T3,Me_T3] = Elementary_matrices(Xe, N, Ns, Nt, E, nedof, ngaus, wgp, nu, alpha, e)
% Elementary matrix for a 2D elasticity problem / plane stress /
% Xe:           nodal coordinates
% N, Ns, Nt:    shape functions and their derivatives
% nedof:        number of degrees of freedom in the element
% ngaus, wgp:   number of integration points and integration weights
% nu, e, alpha: material property

Ketest_T3 = zeros(nedof); 
B         = zeros(3,nedof); 
xe        = Xe(:,2);
nen       = length(xe); 

 for ig = 1:ngaus
    % plane stress (alpha = 0) plane strain (alpha = 1)
    C = [1 - alpha * nu, nu, 0; nu, 1 - alpha * nu, 0; 0, 0, (1 - nu - alpha * nu)/2];
    C = C * E/(1 - nu - alpha * nu)/(1 + nu);
   
    N_ig    = N(ig,:);
    Ns_ig   = Ns(ig,:);
    Nt_ig   = Nt(ig,:);

    Jacob = [Ns_ig(1:nen); Nt_ig(1:nen)] * Xe;
    J     = det(Jacob); 
    dvolu = wgp(ig) * J; 
    
    res = Jacob\[Ns_ig; Nt_ig];
    nx  = res(1,:); 
    ny  = res(2,:); 
    
    B(1,1:2:end) = nx;
    B(2,2:2:end) = ny;
    B(3,1:2:end) = ny;
    B(3,2:2:end) = nx;
    Ketest_T3    = Ketest_T3 + e * B' * C * B * abs(dvolu); 
%   Ketest_T3    = Ketest_T3 + e * B' * C * B * dvolu; 
 end
end
