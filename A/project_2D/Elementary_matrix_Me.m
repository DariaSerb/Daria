function  [Metest_T3, Metest_T3_lump] = Elementary_matrix_Me(Xe, N, Ns, Nt, rho, nedof, ngaus, wgp, e)
% [Metest_T3, Metest_T3_lump] = Elementary_matrix_Me(Xe, N, rho, nedof, ngaus, wgp, e);
% Elementary matrix for a 2D elasticity problem / plane stress /
% Xe:           nodal coordinates
% N, Ns, Nt:    shape functions and their derivatives
% nedof:        number of degrees of freedom in the element
% ngaus, wgp:   number of integration points and integration weights
% nu, e, alpha: material property

Metest_T3 = zeros(nedof); 

% Nm   = zeros(3, nedof); 
xe        = Xe(:,2);
nen       = length(xe); 

 for ig = 1:ngaus
    % plane stress (alpha = 0) plane strain (alpha = 1)
    N_ig    = N(ig,:);
    Ns_ig   = Ns(ig,:);
    Nt_ig   = Nt(ig,:);

    Jacob = [Ns_ig(1:nen); Nt_ig(1:nen)] * Xe;
    J     = det(Jacob); 
    dvolu = wgp(ig) * J; 
    
    N_T_ig    = [reshape([1;0] * N_ig,1,nedof); reshape([0;1] * N_ig,1,nedof)];
    Metest_T3 = Metest_T3 + e * N_T_ig' * rho * N_T_ig * abs(dvolu); 
%   Metest_T3 = Metest_T3 + e * N_T_ig' * rho * N_T_ig * dvolu; 
 end

   Metest_T3_lump = zeros(nedof);
 for in=1:nedof
   Metest_T3_lump(in,in)=sum(Metest_T3(in,:));
 end
end
