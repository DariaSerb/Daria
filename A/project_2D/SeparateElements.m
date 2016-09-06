function [Elements] = SeparateElements(Elements,Nodes,Graphic_display)

% Separates the elements in three sets, depending on whether they are
% completely inside the levelset (where the level-set function is < 0), or
% outside the levelset (where the level-set function is > 0), or cut by the interface
% Elements.In, Elements.Out and Elements.Cut
% Elements: connectivities matrix
% LS: values of the level-set functions at the nodes

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

% Add 0 to end of element table in order to identify FEM element
% Add 1 to end + 1 of element table in order to specify epsilon (full matter)
Elements(:,5) = zeros(size(Elements,1),1);
Elements(:,6) = ones(size(Elements,1),1);

% Add -1 to end of element table in order to identify a "no matter element"
% Add 0 to end+1 of element table in order to specify no matter
for i=1:size(Elems_In,2)
    Elements(Elems_In(i),5) = - 1;
    Epsilon_matter = 0;
    Elements(Elems_In(i),6) = Epsilon_matter;
end

% Add 1 to end of element table in order to identify a "XFEM element"
for i=1:size(Elems_Cut,2)
    Elements(Elems_Cut(i),5) = 1;
    Epsilon_matter = AreaRatio(Elements,Elems_Cut(i),Nodes,LS);
    Elements(Elems_Cut(i),6) = Epsilon_matter;
end

if strcmp(upper(Graphic_display),'YES')
    
    figure(1);
    hold on;
    for j = 1:size(Elements,1)
        % Plot the FEM elements
        if Elements(j,5) == 0
            tab(1,1:2) = [Nodes(Elements(j,2),2) Nodes(Elements(j,2),3)];
            tab(2,1:2) = [Nodes(Elements(j,3),2) Nodes(Elements(j,3),3)];
            tab(3,1:2) = [Nodes(Elements(j,4),2) Nodes(Elements(j,4),3)];
            tab(4,1:2) = [Nodes(Elements(j,2),2) Nodes(Elements(j,2),3)];
            plot(tab(:,1),tab(:,2),'b.-')
        end
        % Plot the XFEM elements
         if Elements(j,5) == 1
            tab(1,1:2) = [Nodes(Elements(j,2),2) Nodes(Elements(j,2),3)];
            tab(2,1:2) = [Nodes(Elements(j,3),2) Nodes(Elements(j,3),3)];
            tab(3,1:2) = [Nodes(Elements(j,4),2) Nodes(Elements(j,4),3)];
            tab(4,1:2) = [Nodes(Elements(j,2),2) Nodes(Elements(j,2),3)];
            plot(tab(:,1),tab(:,2),'r.-')
        end
        % Plot the "non matter elements"
         if Elements(j,5) == - 1
            tab(1,1:2) = [Nodes(Elements(j,2),2) Nodes(Elements(j,2),3)];
            tab(2,1:2) = [Nodes(Elements(j,3),2) Nodes(Elements(j,3),3)];
            tab(3,1:2) = [Nodes(Elements(j,4),2) Nodes(Elements(j,4),3)];
            tab(4,1:2) = [Nodes(Elements(j,2),2) Nodes(Elements(j,2),3)];
            plot(tab(:,1),tab(:,2),'w.-')
        end
    end
    hold off;
end

