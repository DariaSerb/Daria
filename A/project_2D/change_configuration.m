function NodCur = change_configuration(Nodes,dTau)
% the definition of displacement of the cut element

qm          = zeros(size(Nodes,1),2);
NodCur      = zeros(size(Nodes,1),3);
NodCur(:,1) = (1:length(Nodes));

for in = 1:size(Nodes,1);
  [~,qm(in,:),~,~] = q_calc_funcMain(Nodes(in,2:3));   
% [qm(in,:),~, ~]  = q_calc_funcPolynom4(Nodes(in,2:3));   
 
  NodCur(in,2) = Nodes(in,2) + qm(in,1)*dTau;
  NodCur(in,3) = Nodes(in,3) + qm(in,2)*dTau;
end
end
