function error_indicator(Nodes,Elements,Type_LS)
%This function is to plot the Modal Assurance Criterion (MAC) matrix between identified mode shapes

P = Initialize_Parameters_2D();

% the calculation of Mnorm for each time step
for n = 1:P.N_tau
  [~,~,~,eigvecnorm,~,Knorm,Mnorm,rationorm] = EigenFunctionEst_XFEM_DD(Nodes,Elements,Type_LS,P.dTau(n));
  eigvecnormTotal(:,:,n) = eigvecnorm;
  KnormTotal(:,:,n) = Knorm;
  MnormTotal(:,:,n) = Mnorm;
  rationormTotal(:,:,n) = rationorm;
end
% save('eigvecnormTotal.mat','eigvecnormTotal') 
% save('KnormTotal.mat','KnormTotal') 
% save('MnormTotal.mat','MnormTotal') 
% save('rationormTotal.mat','rationormTotal') 

save('norm.mat','KnormTotal','MnormTotal','rationormTotal')
end

