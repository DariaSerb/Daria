function coef = scalefactor(v,w)
%SCALEFACTOR is the coefficient for multiplication in eigenshapes
% v is the eigenvector (x - direction)
% w is the eigenvector (y - direction)

 coef(1) = 1/sum(v(:,1).^2);
 coef(2) = 1/sum(w(:,1).^2);
end

