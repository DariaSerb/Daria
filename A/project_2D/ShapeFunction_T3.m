function [N, Ns, Nt] = ShapeFunction_T3(zgp) 

% [N,Ns,Nt] = ShapeFunc(zgp) 
% Input:    
% zgp:  coordinates of Gauss points in the reference element
% Output:
% N, Ns, Nt: matrices storing the values of the shape functions on the Gauss points
% of the reference element. Each row concerns to a Gauss point

s = zgp(:,1); 
t = zgp(:,2); 

vect0 = zeros(size(s));
vect1 = ones(size(s));

N    = [1 - s - t, s, t]; 
Ns   = [-vect1, vect1, vect0]; 
Nt   = [-vect1, vect0, vect1]; 

end