function [Elements] = SeparateElements(Elements,Nodes)

LS = Nodes(:,4);

ind_In  = 1;
ind_Out = 1;
ind_Cut = 1;
for i = 1:size(Elements,1)
    Te = Elements(i,2:end); 
    LSe = LS(Te); 
    if all(LSe <= 0) 
        Elems_In(ind_In) = i; 
        ind_In = ind_In + 1; 
    elseif all(LSe > 0)
        Elems_Out(ind_Out) = i; 
        ind_Out = ind_Out + 1; 
    else
        Elems_Cut(ind_Cut) = i; 
        ind_Cut = ind_Cut + 1; 
    end
end

Elements(:,5) = zeros(size(Elements,1),1);
Elements(:,6) = ones(size(Elements,1),1);

for i=1:size(Elems_In,2)
    Elements(Elems_In(i),5) = - 1;
    Epsilon_matter = 0;
    Elements(Elems_In(i),6) = Epsilon_matter;
end

for i=1:size(Elems_Cut,2)
    Elements(Elems_Cut(i),5) = 1;
    Epsilon_matter = AreaRatio(Elements,Elems_Cut(i),Nodes,LS);
    Elements(Elems_Cut(i),6) = Epsilon_matter;
end


