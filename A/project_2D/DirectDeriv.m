function [lambda,V,DLambda,Du,dof_out_BC_hole,C] = DirectDeriv(Nodes,Elements)

% the calculation of DD of the first order
P = Initialize_Parameters_2D();

ModeCnt = P.ModeCnt;
ModeEst = P.ModeEst;
iorder  = 1;

[stiff,mass,stiffdq,stiffdivq,massdivq,lambda,~,V,mui,dof_out_BC_hole,C] = eig_val_r0(Nodes,Elements);

alp_num = zeros(ModeCnt,ModeEst);
 
for n = 1:ModeEst
    var = (stiffdq + stiffdq - stiffdivq + lambda(n) * massdivq) * V(:,n);
    numer = (V(:,n))' * var;
    DLambda(1,n) = - numer/mui(n);
   
    % the calculation of DD of the second order
      if iorder == 2;
      term1 = 2 * DLambda(1,n) * massdivq;    
      sum   = (1/2) * (- stiff + lambda(n) * mass) + 2 * stiff;
      term2 = 2 * sum ;
      term3 = stiffdq + term4;
      var2 = (term1 + term2 - term3) * V(:,n); 
      numer2 = (V(:,n))' * var2;
      D2Lambda(1,n) = - numer2/mui(n);
      end
    
    for m = 1:ModeCnt
        if n == m
            continue;
        end
        numer = (V(:,m))' * var; 
        alp_num(m,n) = numer/(mui(m) * (lambda(m) - lambda(n)));
    end
end

Du = V(:,1:ModeCnt)*alp_num;
end