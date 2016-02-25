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
    
    N_T_ig    = [reshape([1;0]*N_ig,1,nedof); reshape([0;1]*N_ig,1,nedof)];
    Metest_T3 = Metest_T3 + e * N_T_ig' * rho * N_T_ig * dvolu; 
 end

 Metest_T3_lump = zeros(nedof);
 for in=1:nedof
   Metest_T3_lump(in,in)=sum(Metest_T3(in,:));
 end
%  Metest_T3_lump(1,1) = Metest_T3(1,1)+Metest_T3(1,2)+Metest_T3(1,3)+Metest_T3(1,4)+Metest_T3(1,5)+Metest_T3(1,6);
%  Metest_T3_lump(2,2) = Metest_T3(2,1)+Metest_T3(2,2)+Metest_T3(2,3)+Metest_T3(2,4)+Metest_T3(2,5)+Metest_T3(2,6);
%  Metest_T3_lump(3,3) = Metest_T3(3,1)+Metest_T3(3,2)+Metest_T3(3,3)+Metest_T3(3,4)+Metest_T3(3,5)+Metest_T3(3,6);
%  Metest_T3_lump(4,4) = Metest_T3(4,1)+Metest_T3(4,2)+Metest_T3(4,3)+Metest_T3(4,4)+Metest_T3(4,5)+Metest_T3(4,6);
%  Metest_T3_lump(5,5) = Metest_T3(5,1)+Metest_T3(5,2)+Metest_T3(5,3)+Metest_T3(5,4)+Metest_T3(5,5)+Metest_T3(5,6);
%  Metest_T3_lump(6,6) = Metest_T3(6,1)+Metest_T3(6,2)+Metest_T3(6,3)+Metest_T3(6,4)+Metest_T3(6,5)+Metest_T3(6,6);
end



