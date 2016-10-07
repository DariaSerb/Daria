clear all;
close all;
clc;
  
% Begin Timer
tic

Graphic_display  = 'yes';
Name_GMSH     = '2Dmesh_N_321_W_0.6_H_0.4.msh';
% Name_GMSH   = '2Dmesh_N_1211_W_0.6_H_0.4.msh'; 	

Type_LS          = 'Circle';
P                = Initialize_Parameters_2D();

% Read of GMSH information
[Nodes,Elements]  = ReadGMSH(Name_GMSH,Graphic_display);
% Initialization and computing Level Set
[Nodes] = ComputeLS(Nodes,Type_LS,P.Data_LS,Graphic_display);

% Separation the elements in three sets
if strcmp(Type_LS, 'Circle')
    [Elements] = SeparateElements(Elements,Nodes,Graphic_display);

    % Separation the elements inside the void
    ind_In  = 1;
    for ii = 1:size(Elements,1)
        if Elements(ii,5) == -1;
            Elems_In(ind_In) = ii;
            ind_In = ind_In + 1;
        end
    end
else
    Elements(:,5) = zeros(size(Elements,1),1);
end

% the definition of displacement of the cut element
for in=1:size(Elements,1)
 if Elements(in,5) == 1
 
    Te     = Elements(in,2:4);
    Xe     = Nodes(Te,2:3);
    LSe    = Nodes(Te,4); 
    [p1 p2] = FindPointLS0test(Xe,LSe);
    node0 = find(LSe < 0);  
    if length(node0) == 1
            node1 = node0+1; if(node1) == 4, node1 = 1; end
            node2 = node0-1; if(node2) == 0, node2 = 3; end
            P0tr = Te(node0);  
            P1tr = Te(node1);  
            P2tr = Te(node2); 
              if p1new(1,:)> p2new(1,:) && node0 == 2   
                 tab(1,1:2) = [p1(1,:) p1(2,:)];
                 tab(2,1:2) = [p2(1,:) p2(2,:)];
                 tab(3,1:2) = [Nodes(P1tr,2) Nodes(P1tr,3)];
                 tab(4,1:2) = [Nodes(P2tr,2) Nodes(P2tr,3)];
                 tab(5,1:2) = [p1(1,:) p1(2,:)];
                 plot(tab(1:5,1),tab(1:5,2),'c.-') 
              else 
                 tab(1,1:2) = [p1(1,:) p1(2,:)];
                 tab(2,1:2) = [Nodes(P1tr,2) Nodes(P1tr,3)];
                 tab(3,1:2) = [Nodes(P2tr,2) Nodes(P2tr,3)];
                 tab(4,1:2) = [p2(1,:) p2(2,:)];
                 tab(5,1:2) = [p1(1,:) p1(2,:)];
                 plot(tab(1:5,1),tab(1:5,2),'c.-') 
              end
      elseif length(node0) == 2 
            node1 = setdiff([1 2 3],node0);
            P0tr = Te(node1); 
            P1tr = Te(node0(1)); 
            P2tr = Te(node0(2));
            if p1(1,:)> p2(1,:)            
                tab(1,1:2) = [p2(1,:) p2(2,:)];
                tab(2,1:2) = [Nodes(P0tr,2) Nodes(P0tr,3)];
                tab(3,1:2) = [p1(1,:) p1(2,:)];
                tab(4,1:2) = [p2(1,:) p2(2,:)];
                plot(tab(1:4,1),tab(1:4,2),'c.-') 
            else
                tab(1,1:2) = [p1(1,:) p1(2,:)];
                tab(2,1:2) = [Nodes(P0tr,2) Nodes(P0tr,3)];
                tab(3,1:2) = [p2(1,:) p2(2,:)];
                tab(4,1:2) = [p1(1,:) p1(2,:)];
                plot(tab(1:4,1),tab(1:4,2),'c.-') 
            end
    end
 end
end
